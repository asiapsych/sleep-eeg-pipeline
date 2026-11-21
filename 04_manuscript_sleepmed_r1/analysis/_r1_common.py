"""Merge the EEG and ISI data and split participants with insomnia symptoms into two severity groups."""

from __future__ import annotations

import sqlite3
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import pandas as pd
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

from _dynamic_n import DB_PATH, ISI_ITEMS, SAMPLE, SEED, get_labels

BASE = Path(__file__).resolve().parent
AN = BASE / "outputs" / "analysis"
F_S10 = AN / "S10_fooof_full" / "03_subject_brain_avg.csv"
F_S10_CH = AN / "S10_fooof_full" / "01_specparam_features.csv"
F_S17 = AN / "S17_relative_power" / "01_relpower_subject.csv"

OUTCOME = "beta_pw_BRAIN"  # primary outcome: periodic beta peak power (SpecParam), whole-scalp mean
MIDLINE = ["FZ", "CZ", "PZ"]
# expected sizes, to catch a wrong exclusion list
EXPECTED = {
    "submitted": dict(N=275, K2=[104, 61]),
    "audited": dict(N=271, K2=[102, 60]),
}


CLUSTER_INPUTS: dict = {}


def prepare_isi_items(items):
    """Median-impute the ISI items. Returns the filled scores and their z scores."""
    filled = SimpleImputer(strategy="median").fit_transform(items)
    return filled, StandardScaler().fit_transform(filled)


def load_data(sample: str | None = None) -> pd.DataFrame:
    """Joint EEG and ISI sample with group (HC/INS), ISI, demographics, questionnaires and periodic beta."""
    df = get_labels(sample).merge(
        pd.read_csv(F_S10)[["subject", OUTCOME]], on="subject", how="left"
    )
    # PSQI total is only used to describe the sample.
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
    """Sorted by subject, with the K = 2 labels (Lower/Higher) and all EEG measures."""
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
    assert df["subject"].is_unique and len(df) == exp["N"], f"N = {len(df)}"

    ins_mask = df["group"] != "HC"
    ins = df.loc[ins_mask]
    # K-means initialisation depends on row order; keep the original subject order.
    _, items = prepare_isi_items(ins[ISI_ITEMS])
    labels = KMeans(n_clusters=2, random_state=SEED, n_init=100).fit_predict(items)
    raw = pd.Series(labels, index=ins.index)
    order = ins.groupby(raw)["ISI_T1"].mean().sort_values().index
    names = dict(zip(order, ("Lower", "Higher")))
    df["k2"] = "HC"
    df.loc[ins_mask, "k2"] = raw.map(names).values
    k2n = sorted(df.loc[ins_mask, "k2"].value_counts().tolist(), reverse=True)
    assert k2n == exp["K2"], k2n
    # not in df.attrs: arrays there break attrs comparison when pandas merges tables
    CLUSTER_INPUTS[sample] = {"X_ins": items, "km": {2: labels}}
    return df
