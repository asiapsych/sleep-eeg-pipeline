"""EEG 與 ISI 的共同樣本、排除名單及初投三群標籤。

R1 預設用 audited；submitted 留給回覆信的初投樣本對照。
群名沿用 MEW／MDI／SCA，新的兩群分析在 _r1_common.load_all。"""

from __future__ import annotations
import json
import os
import warnings
from pathlib import Path
import sqlite3
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

_PAPER_ROOT = Path(__file__).parent
S17_RELPOWER_CSV = (
    _PAPER_ROOT
    / "outputs"
    / "補充資料"
    / "analysis"
    / "S17_relative_power"
    / "01_relpower_subject.csv"
)
# 受試者編號與本機路徑放在 local_settings.json（只留在本機，不公開；格式見 local_settings.example.json）。
_LOCAL = _PAPER_ROOT / "local_settings.json"
_CFG = json.loads(_LOCAL.read_text(encoding="utf-8")) if _LOCAL.exists() else {}
if not _CFG:
    warnings.warn("找不到 local_settings.json：沒有排除任何受試者，資料路徑要用環境變數指定。")


def local_path(key: str, env: str, default: str = "") -> str:
    """資料位置：環境變數優先，其次是 local_settings.json。"""
    return os.environ.get(env) or _CFG.get("paths", {}).get(key) or default


DB_PATH = local_path("db", "PSYML_DB", "data/psy_ml_isi.db")
ISI_ITEMS = [f"ISI_T1_Q{i}" for i in range(1, 8)]
EXCLUDED = list(_CFG.get("excluded", []))

# R1 資料核對（2026-11-16）：逐筆核對收案紀錄後移除，不進任何分析
# 符合原訂排除條件（n = 4）
EXCLUDED_ELIGIBILITY = list(_CFG.get("excluded_eligibility", []))
# 樣本：預設為核對後的樣本；設 PSYML_SAMPLE=submitted 可重現初投版的 N = 275（回覆信對照用）
SAMPLE = os.environ.get("PSYML_SAMPLE", "audited")


def excluded_ids(sample: str | None = None) -> list:
    """回傳該樣本定義下要排除的受試者編號。"""
    sample = sample or SAMPLE
    assert sample in ("audited", "submitted"), sample
    return EXCLUDED + (EXCLUDED_ELIGIBILITY if sample == "audited" else [])


SEED = 42
DISPLAY_GROUPS = ["HC", "MEW", "MDI", "SCA"]

_cache = {}


def _build_labels(sample: str | None = None) -> pd.DataFrame:
    """Build subject + group + ISI dataframe on the v2.2 EEG∩ISI joint sample."""
    if not S17_RELPOWER_CSV.exists():
        raise FileNotFoundError(
            f"_dynamic_n: 找不到 {S17_RELPOWER_CSV}\n"
            "請先跑 `python3 task_f_recompute_relpower.py`。"
        )
    rel = pd.read_csv(S17_RELPOWER_CSV)[["subject"]].drop_duplicates()

    cols = ["NEW_NUMBER", "AGE", "SEX", "ISI_T1", "BDI_T1", "BAI_T1", "EDU"] + ISI_ITEMS
    conn = sqlite3.connect(DB_PATH)
    db = pd.read_sql(
        f"SELECT {','.join(cols)} FROM isi_raw_data "
        f"WHERE ISI_T1 IS NOT NULL AND ISI_T1_Q1 IS NOT NULL",
        conn,
    )
    conn.close()
    for c in db.columns:
        if c != "NEW_NUMBER":
            db[c] = pd.to_numeric(db[c], errors="coerce")
    db = db[~db["NEW_NUMBER"].isin(excluded_ids(sample))].copy()

    merged = rel.merge(db.rename(columns={"NEW_NUMBER": "subject"}), on="subject", how="inner")
    merged = merged.dropna(subset=ISI_ITEMS).copy()

    ins = merged[merged["ISI_T1"] >= 9].copy()
    ctrl = merged[merged["ISI_T1"] < 9].copy()

    X = StandardScaler().fit_transform(
        SimpleImputer(strategy="median").fit_transform(ins[ISI_ITEMS])
    )
    km = KMeans(n_clusters=3, random_state=SEED, n_init=100)
    ins["cluster_raw"] = km.fit_predict(X)

    order = ins.groupby("cluster_raw")["ISI_T1"].mean().sort_values().index.tolist()
    q23 = {
        ci: ins.loc[ins["cluster_raw"] == ci, ["ISI_T1_Q2", "ISI_T1_Q3"]].mean().sum()
        for ci in order
    }
    maint = max(q23, key=q23.get)
    severe = order[2]
    if maint == severe:
        q57 = {
            ci: ins.loc[ins["cluster_raw"] == ci, ["ISI_T1_Q5", "ISI_T1_Q7"]].mean().sum()
            for ci in order
        }
        severe = max(q57, key=q57.get)
    light = [c for c in order if c not in (maint, severe)][0]

    ins["group"] = ins["cluster_raw"].map({maint: "MEW", light: "MDI", severe: "SCA"})
    ctrl["group"] = "HC"
    labels = pd.concat([ctrl, ins], ignore_index=True)[
        ["subject", "group", "ISI_T1", "AGE", "SEX", "BDI_T1", "BAI_T1", "EDU"] + ISI_ITEMS
    ]
    return labels


def get_labels(sample: str | None = None) -> pd.DataFrame:
    """Cached subject + group + covariates dataframe（v2.2 EEG∩ISI 樣本；sample 見 excluded_ids）。"""
    sample = sample or SAMPLE
    if sample not in _cache:
        _cache[sample] = _build_labels(sample)
    return _cache[sample]


def load_group_n() -> dict:
    """{group_display: n} 字典，固定四群順序 HC/MEW/MDI/SCA。"""
    df = get_labels()
    vc = df["group"].value_counts().to_dict()
    return {g: int(vc.get(g, 0)) for g in DISPLAY_GROUPS}


def load_n_string(sep: str = ", ") -> str:
    """產生目前 v2.2 群組樣本數字串，例如 'HC = 110, MEW = 47, MDI = 68, SCA = 50'。"""
    gn = load_group_n()
    return sep.join(f"{g} = {gn[g]}" for g in DISPLAY_GROUPS)


def total_n() -> int:
    return sum(load_group_n().values())


def insomnia_n() -> int:
    gn = load_group_n()
    return gn["MEW"] + gn["MDI"] + gn["SCA"]


def restrict_to_v22(
    df: pd.DataFrame, id_col: str = "NEW_NUMBER", quiet: bool = False
) -> pd.DataFrame:
    """限制 DB DataFrame 到 v2.2 EEG 有效樣本（n=275）。

    用法：
        from _dynamic_n import restrict_to_v22
        df = restrict_to_v22(df)  # 在 K-means 之前
    """
    if not S17_RELPOWER_CSV.exists():
        if not quiet:
            print("  [WARN] S17 CSV not found, restrict_to_v22 skipped")
        return df
    valid_sids = set(pd.read_csv(S17_RELPOWER_CSV)["subject"].drop_duplicates().tolist())
    before = len(df)
    out = df[df[id_col].isin(valid_sids)].copy()
    if not quiet:
        print(f"  [v2.2 new methodology] restrict_to_v22: {before} → {len(out)}")
    return out
