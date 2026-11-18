"""R1 共用的樣本、分群與效果量。

從 S10、S17 接上量表。K-means 前先固定 subject 列序，再做 median imputation
與標準化；這個順序要和原分析一致，否則標籤與 bootstrap 不一定能重現。"""

from __future__ import annotations
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import adjusted_rand_score

import sqlite3

from _dynamic_n import DB_PATH, ISI_ITEMS, SAMPLE, get_labels

BASE = Path(__file__).resolve().parent
AN = BASE / "outputs" / "補充資料" / "analysis"
F_S10 = AN / "S10_fooof_full" / "03_subject_brain_avg.csv"
F_S10_CH = AN / "S10_fooof_full" / "01_specparam_features.csv"
F_S17 = AN / "S17_relative_power" / "01_relpower_subject.csv"

OUTCOME = "beta_pw_BRAIN"  # 主要指標：SpecParam periodic beta peak power（全頭皮平均）
SEED = 42
COV = "AGE + C(SEX) + EDU"
COV_COLS = ["AGE", "SEX", "EDU"]
MIDLINE = ["FZ", "CZ", "PZ"]
# 各樣本定義下預期的人數（防呆：資料或排除名單被動到時會報錯）
EXPECTED = {
    "submitted": dict(N=275, K2=[104, 61], K4=[67, 35, 34, 29]),
    "audited": dict(N=271, K2=[102, 60], K4=[61, 40, 31, 30]),
}


CLUSTER_INPUTS: dict = {}  # load_all() 最近一次計算的分群輸入：{sample: {'X_ins': 標準化後的七題矩陣, 'km': {K: 標籤}}}


def load_data(sample: str | None = None) -> pd.DataFrame:
    """EEG∩ISI 樣本（人數依 sample 而定，見 EXPECTED）：subject、原投稿的三群標籤、ISI、人口學、量表、periodic beta。"""
    df = get_labels(sample).merge(
        pd.read_csv(F_S10)[["subject", OUTCOME]], on="subject", how="left"
    )
    # PSQI 總分只作表 1 的描述，不進任何模型（2026-11-16 加）
    with sqlite3.connect(DB_PATH) as con:
        psqi = pd.read_sql(
            "SELECT NEW_NUMBER AS subject, PSQI_T1 FROM isi_raw_data WHERE NEW_NUMBER IS NOT NULL",
            con,
        )
    df = df.merge(psqi, on="subject", how="left")
    for c in [OUTCOME, "AGE", "ISI_T1", "BDI_T1", "BAI_T1", "PSQI_T1"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_all(sample: str | None = None) -> pd.DataFrame:
    """依 subject 升冪（人數依 sample 而定）；含 K=3 標籤、現算 K=2／K=4 標籤、所有 EEG 指標。"""
    sample = sample or SAMPLE
    exp = EXPECTED[sample]
    df = load_data(sample).copy()
    df["group"] = df["group"].astype(str)
    s10 = pd.read_csv(F_S10)
    s10 = s10[[c for c in s10.columns if c == "subject" or c not in df.columns]]
    s17 = pd.read_csv(F_S17)
    s17 = s17[[c for c in s17.columns if c == "subject" or c not in df.columns]]
    ch = pd.read_csv(F_S10_CH)
    mid = (
        ch[ch["channel"].isin(MIDLINE)]
        .groupby("subject")["beta_pw"]
        .apply(lambda s: s.dropna().mean())
        .rename("beta_pw_MIDLINE")
        .reset_index()
    )
    df = (
        df.merge(s10, on="subject", how="left")
        .merge(s17, on="subject", how="left")
        .merge(mid, on="subject", how="left")
        .sort_values("subject", kind="mergesort")
        .reset_index(drop=True)
    )
    assert df["subject"].is_unique and (exp["N"] is None or len(df) == exp["N"]), f"N = {len(df)}"

    ins_mask = df["group"] != "HC"
    ins = df.loc[ins_mask]
    X = StandardScaler().fit_transform(
        SimpleImputer(strategy="median").fit_transform(ins[ISI_ITEMS])
    )
    km = {k: KMeans(n_clusters=k, random_state=SEED, n_init=100).fit_predict(X) for k in (2, 3, 4)}
    assert adjusted_rand_score(ins["group"], km[3]) > 0.9999, "K=3 無法重現主分析標籤"

    for k in (2, 4):
        raw = pd.Series(km[k], index=ins.index)
        order = (
            ins.groupby(raw)["ISI_T1"].mean().sort_values().index.tolist()
        )  # 依 ISI 總分由低到高編號
        name = {c: f"K{k}_C{i + 1}" for i, c in enumerate(order)}
        df[f"k{k}"] = "HC"
        df.loc[ins_mask, f"k{k}"] = raw.map(name).values
    k2n = sorted(df.loc[ins_mask, "k2"].value_counts().tolist(), reverse=True)
    assert exp["K2"] is None or k2n == exp["K2"], k2n
    k4n = sorted(df.loc[ins_mask, "k4"].value_counts().tolist(), reverse=True)
    assert exp["K4"] is None or k4n == exp["K4"], k4n
    df["k2"] = df["k2"].replace({"K2_C1": "Lower", "K2_C2": "Higher"})  # K=2：低分群／高分群
    # 分群用的矩陣與標籤不掛在 df.attrs（attrs 內有陣列時，pandas 在合併兩個衍生資料表時會比較失敗）
    CLUSTER_INPUTS[sample] = dict(X_ins=X, km=km)
    return df
