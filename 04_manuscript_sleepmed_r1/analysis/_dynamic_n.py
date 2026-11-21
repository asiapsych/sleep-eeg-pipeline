"""Joint EEG and ISI sample and the exclusion lists.

Participants with ISI >= 9 are labelled INS and the rest HC. The two severity
groups are assigned in _r1_common.load_all."""

from __future__ import annotations

import json
import os
import sqlite3
import warnings
from pathlib import Path

import pandas as pd

_PAPER_ROOT = Path(__file__).parent
S17_RELPOWER_CSV = (
    _PAPER_ROOT
    / "outputs"
    / "analysis"
    / "S17_relative_power"
    / "01_relpower_subject.csv"
)
# Participant IDs and local paths are kept in local_settings.json, which is not
# shared (format: local_settings.example.json).
_LOCAL = _PAPER_ROOT / "local_settings.json"
_CFG = json.loads(_LOCAL.read_text(encoding="utf-8")) if _LOCAL.exists() else {}
if not _CFG:
    warnings.warn("local_settings.json not found: nobody is excluded and data paths must come from environment variables.")


def local_path(key: str, env: str, default: str = "") -> str:
    """Environment variable first, then local_settings.json."""
    return os.environ.get(env) or _CFG.get("paths", {}).get(key) or default


DB_PATH = local_path("db", "PSYML_DB", "data/psy_ml_isi.db")
ISI_ITEMS = [f"ISI_T1_Q{i}" for i in range(1, 8)]
EXCLUDED = list(_CFG.get("excluded", []))

# exclusions added in the revision
EXCLUDED_ELIGIBILITY = list(_CFG.get("excluded_eligibility", []))
# Revised sample by default; PSYML_SAMPLE=submitted gives the submitted sample (N = 275).
SAMPLE = os.environ.get("PSYML_SAMPLE", "audited")


def excluded_ids(sample: str | None = None) -> list:
    """IDs excluded under the given sample definition."""
    sample = sample or SAMPLE
    assert sample in ("audited", "submitted"), sample
    return EXCLUDED + (EXCLUDED_ELIGIBILITY if sample == "audited" else [])


SEED = 42

_cache = {}


def _build_labels(sample: str | None = None) -> pd.DataFrame:
    if not S17_RELPOWER_CSV.exists():
        raise FileNotFoundError(
            f"_dynamic_n: {S17_RELPOWER_CSV} not found.\nPut the relative band-power table there first."
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

    ins["group"] = "INS"
    ctrl["group"] = "HC"
    labels = pd.concat([ctrl, ins], ignore_index=True)[
        ["subject", "group", "ISI_T1", "AGE", "SEX", "BDI_T1", "BAI_T1", "EDU"] + ISI_ITEMS
    ]
    return labels


def get_labels(sample: str | None = None) -> pd.DataFrame:
    """Read the sample once per exclusion set."""
    sample = sample or SAMPLE
    if sample not in _cache:
        _cache[sample] = _build_labels(sample)
    return _cache[sample]
