"""
compute_features.py — PSD → EEG features
========================================

Read the per-subject PSD ``.mat`` files exported by stage 02 and compute,
from a *single* PSD with a *single* band definition:

  1. Traditional band power
       - absolute & relative power (Delta/Theta/Alpha/Beta/Gamma + sub-bands)
       - common band ratios (TBR, ABR, ...)
  2. Periodic / aperiodic decomposition (FOOOF / specparam)
       - aperiodic exponent & offset
       - periodic peak (center freq / power / bandwidth) per band

This script contains **only feature extraction**. Group assignment, clinical
scales and statistics belong to each study's private analysis and are not
included here. No database, no subject grouping, no hypothesis tests.

Input
-----
A folder of ``psd_<ID>.mat`` files, each containing:
    freqs     : [n_freq]            frequency axis (Hz)
    psd       : [n_ch x n_freq]     power spectral density
    ch_names  : [n_ch]              channel labels
    n_ep      : scalar (optional)   number of epochs averaged

Output (CSV, written to --out-dir)
----------------------------------
    traditional_power_subject.csv   subject x band power (BRAIN + regions + ratios)
    specparam_channel.csv           subject x channel x aperiodic/periodic params
    specparam_subject.csv           subject-level BRAIN/region averages
    fit_quality.csv                 FOOOF fit R^2 / error per subject

Usage
-----
    python compute_features.py --psd-dir /path/to/psd_mat --out-dir /path/to/features
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.io import loadmat

try:
    from fooof import FOOOFGroup
except ImportError:                                              # pragma: no cover
    FOOOFGroup = None


# ── Band definitions (literature standard) ───────────────────────────────
# PRIMARY bands are continuous & non-overlapping; their sum (1-45 Hz) is the
# denominator of relative power. DERIVED sub-bands are reported as a share of
# that same total but do NOT enter the denominator.
PRIMARY_BANDS = {
    "Delta": (1.0, 4.0),
    "Theta": (4.0, 8.0),
    "Alpha": (8.0, 13.0),
    "Beta":  (13.0, 30.0),
    "Gamma": (30.0, 45.0),
}
DERIVED_BANDS = {
    "Alpha1":   (8.0, 10.0),
    "Alpha2":   (10.0, 13.0),
    "Beta1":    (13.0, 18.0),
    "Beta2":    (18.0, 24.0),
    "Beta3":    (24.0, 30.0),
    "HighBeta": (30.0, 40.0),
}
ALL_BANDS = {**PRIMARY_BANDS, **DERIVED_BANDS}

# Bands within which to report a FOOOF periodic peak.
PEAK_BANDS = {
    "theta": (4.0, 8.0),
    "alpha": (8.0, 13.0),
    "beta":  (13.0, 30.0),
}

# Normalize any legacy 10-20 temporal labels to the modern convention that
# stage 02 produces (so region maps work regardless of input convention).
CH_NORMALIZE = {"T3": "T7", "T4": "T8", "T5": "P7", "T6": "P8"}

# 19-channel regions in modern 10-20 naming.
BRAIN_REGIONS = {
    "FRONTAL":   ["FP1", "FP2", "F3", "F4", "F7", "F8", "FZ"],
    "CENTRAL":   ["C3", "C4", "CZ"],
    "PARIETAL":  ["P3", "P4", "PZ"],
    "TEMPORAL":  ["T7", "T8", "P7", "P8"],
    "OCCIPITAL": ["O1", "O2"],
}
ALL_19 = sorted({c for chs in BRAIN_REGIONS.values() for c in chs})

# FOOOF settings (fixed aperiodic mode).
FOOOF_KW = dict(
    peak_width_limits=(1.0, 12.0),
    max_n_peaks=8,
    min_peak_height=0.1,
    peak_threshold=2.0,
    aperiodic_mode="fixed",
    verbose=False,
)


# ── PSD loading ──────────────────────────────────────────────────────────

def load_psd(path):
    """Return (freqs, psd[n_ch x n_freq], ch_names) from a PSD .mat."""
    m = loadmat(str(path), squeeze_me=True)
    freqs = np.asarray(m["freqs"], dtype=float).ravel()
    psd = np.asarray(m["psd"], dtype=float)
    ch_raw = m["ch_names"]
    if isinstance(ch_raw, np.ndarray):
        ch = [str(c).strip() for c in ch_raw]
    else:
        ch = [str(ch_raw).strip()]
    ch = [CH_NORMALIZE.get(c.upper(), c.upper()) for c in ch]
    n_ep = int(np.atleast_1d(m.get("n_ep", np.nan))[0]) if "n_ep" in m else np.nan
    return freqs, psd, ch, n_ep


# ── Traditional band power ───────────────────────────────────────────────

def band_abs_power(freqs, psd_1d, lo, hi):
    """Trapezoidal integration of PSD over [lo, hi]."""
    mask = (freqs >= lo) & (freqs <= hi)
    if mask.sum() < 2:
        return np.nan
    return float(np.trapz(psd_1d[mask], freqs[mask]))


def compute_traditional(sid, freqs, psd, ch_names):
    """Subject-level absolute/relative band power + ratios."""
    ch_abs = {}                                   # ch -> {band: abs, _total: ...}
    for ci, ch in enumerate(ch_names):
        bp = {b: band_abs_power(freqs, psd[ci], lo, hi)
              for b, (lo, hi) in ALL_BANDS.items()}
        total = float(np.nansum([bp[b] for b in PRIMARY_BANDS]))
        ch_abs[ch] = {"_total": total, **bp}

    def region_rel(chs, band):
        """Pooled relative power: mean(band_abs) / mean(total) * 100 over chs."""
        absvals, totvals = [], []
        for ch in chs:
            d = ch_abs.get(ch)
            if d and not np.isnan(d[band]) and d["_total"] > 0:
                absvals.append(d[band])
                totvals.append(d["_total"])
        if not absvals:
            return np.nan
        return float(np.mean(absvals) / np.mean(totvals) * 100)

    def region_abs(chs, band):
        vals = [ch_abs[ch][band] for ch in chs
                if ch in ch_abs and not np.isnan(ch_abs[ch][band])]
        return float(np.mean(vals)) if vals else np.nan

    row = {"subject": sid}
    for band in ALL_BANDS:
        row[f"{band}_rel_BRAIN"] = region_rel(ALL_19, band)
        row[f"{band}_abs_BRAIN"] = region_abs(ALL_19, band)
    for region, chs in BRAIN_REGIONS.items():
        for band in PRIMARY_BANDS:
            row[f"{band}_rel_{region}"] = region_rel(chs, band)

    a = {b: region_abs(ALL_19, b) for b in ALL_BANDS}

    def safe(n, d):
        return float(n / d) if (d and not np.isnan(n) and not np.isnan(d) and d != 0) else np.nan

    row["TBR_BRAIN"] = safe(a["Theta"], a["Beta"])      # theta/beta
    row["ABR_BRAIN"] = safe(a["Alpha"], a["Beta"])      # alpha/beta
    row["TAR_BRAIN"] = safe(a["Theta"], a["Alpha"])     # theta/alpha
    row["BAR_BRAIN"] = safe(a["Beta"], a["Alpha"])      # beta/alpha
    row["BDR_BRAIN"] = safe(a["Beta"], a["Delta"])      # beta/delta
    row["DBR_BRAIN"] = safe(a["Delta"], a["Beta"])      # delta/beta
    row["DTABR_BRAIN"] = safe(a["Delta"] + a["Theta"], a["Alpha"] + a["Beta"])
    return row


# ── FOOOF / specparam ────────────────────────────────────────────────────

def _peak_in_band(peaks, band):
    """Strongest (highest-power) peak whose center freq is within band."""
    if peaks is None or len(peaks) == 0:
        return (np.nan, np.nan, np.nan)
    in_band = [(cf, pw, bw) for cf, pw, bw in peaks if band[0] <= cf <= band[1]]
    if not in_band:
        return (np.nan, np.nan, np.nan)
    return max(in_band, key=lambda x: x[1])


def fit_fooof(sid, freqs, psd, ch_names, n_ep, freq_range):
    """Fit FOOOF per channel; return list of per-channel rows."""
    fg = FOOOFGroup(**FOOOF_KW)
    fg.fit(freqs, psd, freq_range)
    rows = []
    for ci, ch in enumerate(ch_names):
        try:
            fm = fg.get_fooof(ci)
            offset, exponent = fm.aperiodic_params_[0], fm.aperiodic_params_[1]
            peaks, r2, err = fm.peak_params_, fm.r_squared_, fm.error_
        except Exception:                                       # noqa: BLE001
            offset = exponent = r2 = err = np.nan
            peaks = None
        row = {"subject": sid, "channel": ch, "n_epochs": n_ep,
               "exponent": exponent, "offset": offset, "r2": r2, "fit_error": err}
        for band, rng in PEAK_BANDS.items():
            cf, pw, bw = _peak_in_band(peaks, rng)
            row[f"{band}_cf"], row[f"{band}_pw"], row[f"{band}_bw"] = cf, pw, bw
        rows.append(row)
    return rows


def aggregate_specparam(df_long):
    """Brain- and region-average the per-channel specparam features."""
    feat_cols = [c for c in df_long.columns
                 if c not in ("subject", "channel", "n_epochs", "r2", "fit_error")]
    rows = []
    for sid, g in df_long.groupby("subject"):
        row = {"subject": sid, "n_epochs": g["n_epochs"].iloc[0]}
        for f in feat_cols:
            row[f"{f}_BRAIN"] = g[f].dropna().mean()
        for region, chs in BRAIN_REGIONS.items():
            gr = g[g["channel"].isin(chs)]
            for f in feat_cols:
                row[f"{f}_{region}"] = gr[f].dropna().mean()
        row["r2_BRAIN_mean"] = g["r2"].mean()
        row["r2_BRAIN_min"] = g["r2"].min()
        row["fit_error_BRAIN_mean"] = g["fit_error"].mean()
        rows.append(row)
    return pd.DataFrame(rows)


# ── Main ─────────────────────────────────────────────────────────────────

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute traditional band power + FOOOF features from PSD .mat.")
    parser.add_argument("--psd-dir", required=True,
                        help="Folder of psd_<ID>.mat files")
    parser.add_argument("--out-dir", required=True,
                        help="Output folder for feature CSVs")
    parser.add_argument("--pattern", default="psd_*.mat",
                        help="Glob for PSD files (default: psd_*.mat)")
    parser.add_argument("--exclude", nargs="*", default=[],
                        metavar="ID", help="Subject IDs to skip")
    parser.add_argument("--freq-min", type=float, default=1.0,
                        help="Lower bound for FOOOF fit & power total (Hz)")
    parser.add_argument("--freq-max", type=float, default=45.0,
                        help="Upper bound for FOOOF fit & power total (Hz)")
    args = parser.parse_args(argv)

    if FOOOFGroup is None:
        parser.error("fooof is not installed. Run: pip install fooof==1.1.0")

    psd_dir = Path(args.psd_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    freq_range = (args.freq_min, args.freq_max)

    def sid_of(p):
        stem = p.stem
        return stem[4:] if stem.startswith("psd_") else stem

    psd_files = sorted(p for p in psd_dir.glob(args.pattern)
                       if sid_of(p) not in set(args.exclude))
    if not psd_files:
        print(f"No PSD files matching {args.pattern} in {psd_dir}")
        return 1
    print(f"Found {len(psd_files)} PSD file(s).")

    trad_rows, specparam_long, failed = [], [], []
    for i, p in enumerate(psd_files, 1):
        sid = sid_of(p)
        try:
            freqs, psd, ch, n_ep = load_psd(p)
            trad_rows.append(compute_traditional(sid, freqs, psd, ch))
            specparam_long.extend(fit_fooof(sid, freqs, psd, ch, n_ep, freq_range))
            if i % 25 == 0 or i == len(psd_files):
                print(f"  [{i:3d}/{len(psd_files)}] {sid} ok")
        except Exception as e:                                  # noqa: BLE001
            failed.append((sid, str(e)))
            print(f"  [{i:3d}/{len(psd_files)}] {sid} FAILED: {e}")

    # Traditional power
    trad = pd.DataFrame(trad_rows)
    trad.to_csv(out_dir / "traditional_power_subject.csv", index=False)
    print(f"\nTraditional power: {trad.shape} -> traditional_power_subject.csv")

    # Specparam (channel-level + subject-level + fit quality)
    df_long = pd.DataFrame(specparam_long)
    df_long.to_csv(out_dir / "specparam_channel.csv", index=False)
    df_subj = aggregate_specparam(df_long)
    df_subj.to_csv(out_dir / "specparam_subject.csv", index=False)
    print(f"Specparam: channel {df_long.shape} -> specparam_channel.csv | "
          f"subject {df_subj.shape} -> specparam_subject.csv")

    fq = df_long.groupby("subject").agg(
        r2_mean=("r2", "mean"), r2_min=("r2", "min"),
        fit_error_mean=("fit_error", "mean"), n_epochs=("n_epochs", "first"),
    ).reset_index()
    fq.to_csv(out_dir / "fit_quality.csv", index=False)
    print(f"Fit quality: mean R^2 = {fq['r2_mean'].mean():.4f} | "
          f"min R^2 = {fq['r2_min'].min():.4f} -> fit_quality.csv")

    if failed:
        print(f"\nFailed: {len(failed)}")
        for sid, e in failed[:10]:
            print(f"  {sid}: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
