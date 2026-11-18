"""逐電極擬合 FOOOF，再彙整腦區及組間比較。

讀 S10_fooof_psd_input 的 MATLAB 頻譜；擬合設定沿用 Bai 2025。
此處的診斷用 Hedges' g 正值代表 b 高於 a，方向與 R1 共用函式相反。"""

from __future__ import annotations

import sqlite3
import warnings
from itertools import combinations
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from fooof import FOOOFGroup
from scipy import stats
from scipy.io import loadmat
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

BASE_DIR = Path(__file__).parent
PSD_DIR = BASE_DIR / "outputs" / "補充資料" / "analysis" / "S10_fooof_psd_input"
OUT_DIR = BASE_DIR / "outputs" / "補充資料" / "analysis" / "S10_fooof_full"
OUT_DIR.mkdir(parents=True, exist_ok=True)

from _dynamic_n import DB_PATH as _DB, EXCLUDED

DB_PATH = Path(_DB)
SEED = 42

ISI_ITEMS = [f"ISI_T1_Q{i}" for i in range(1, 8)]
GROUP_ORDER = ["HC", "MEW", "MD", "SD"]
GROUP_ORDER_CN = ["非失眠對照", "睡眠維持/早醒型", "輕-日常型", "嚴重-高困擾型"]
CODE_MAP = dict(zip(GROUP_ORDER_CN, GROUP_ORDER))

# Bai 2025
FREQ_RANGE = (1.0, 45.0)
PEAK_WIDTH_LIMITS = (1.0, 12.0)
MAX_N_PEAKS = 8
MIN_PEAK_HEIGHT = 0.1
PEAK_THRESHOLD = 2.0
APERIODIC_MODE = "fixed"

# 電極名稱沿用資料庫的 10–20 命名，T7／T8 等別名先轉換。
CH_RENAME = {"T7": "T3", "T8": "T4", "P7": "T5", "P8": "T6"}

BRAIN_REGIONS = {
    "FRONTAL": ["FP1", "FP2", "F3", "F4", "F7", "F8", "FZ"],
    "CENTRAL": ["C3", "C4", "CZ"],
    "PARIETAL": ["P3", "P4", "PZ"],
    "TEMPORAL": ["T3", "T4", "T5", "T6"],
    "OCCIPITAL": ["O1", "O2"],
}

PEAK_BANDS = {
    # 2026/5/13 對齊文獻標準頻段定義（與 task_f 一致）；aperiodic 不受影響
    "theta": (4.0, 8.0),
    "alpha": (8.0, 13.0),
    "beta": (13.0, 30.0),
}


def load_psd_mat(path):
    m = loadmat(str(path), squeeze_me=True)
    freqs = np.asarray(m["freqs"], dtype=float).ravel()
    psd = np.asarray(m["psd"], dtype=float)  # [n_ch × n_freqs]
    ch_raw = m["ch_names"]
    if isinstance(ch_raw, np.ndarray):
        ch_names = [str(c).strip() for c in ch_raw]
    else:
        ch_names = [str(ch_raw).strip()]
    ch_names = [CH_RENAME.get(c.upper(), c.upper()) for c in ch_names]
    n_ep = int(np.atleast_1d(m.get("n_ep", 60))[0])
    return freqs, psd, ch_names, n_ep


def extract_peak_in_band(peaks, band):
    if peaks is None or len(peaks) == 0:
        return (np.nan, np.nan, np.nan)
    in_band = [(cf, pw, bw) for cf, pw, bw in peaks if band[0] <= cf <= band[1]]
    if not in_band:
        return (np.nan, np.nan, np.nan)
    return max(in_band, key=lambda x: x[1])


def fit_subject(subject_id, freqs, psd, ch_names, n_ep):
    fg = FOOOFGroup(
        peak_width_limits=PEAK_WIDTH_LIMITS,
        max_n_peaks=MAX_N_PEAKS,
        min_peak_height=MIN_PEAK_HEIGHT,
        peak_threshold=PEAK_THRESHOLD,
        aperiodic_mode=APERIODIC_MODE,
        verbose=False,
    )
    fg.fit(freqs, psd, FREQ_RANGE)
    rows = []
    for ci, ch in enumerate(ch_names):
        try:
            fm = fg.get_fooof(ci)
            offset, exponent = fm.aperiodic_params_[0], fm.aperiodic_params_[1]
            peaks = fm.peak_params_
            r2 = fm.r_squared_
            err = fm.error_
        except Exception:
            offset = exponent = r2 = err = np.nan
            peaks = None
        row = {
            "subject": subject_id,
            "channel": ch,
            "n_epochs": n_ep,
            "exponent": exponent,
            "offset": offset,
            "r2": r2,
            "fit_error": err,
        }
        for band, rng in PEAK_BANDS.items():
            cf, pw, bw = extract_peak_in_band(peaks, rng)
            row[f"{band}_cf"] = cf
            row[f"{band}_pw"] = pw
            row[f"{band}_bw"] = bw
        rows.append(row)
    return rows


def aggregate_features(df_long):
    feat_cols = [
        c for c in df_long.columns if c not in ("subject", "channel", "n_epochs", "r2", "fit_error")
    ]
    rows = []
    for sid, g in df_long.groupby("subject"):
        row = {"subject": sid, "n_epochs": int(g["n_epochs"].iloc[0])}
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


def load_groups():
    needed = ["NEW_NUMBER", "AGE", "SEX", "ISI_T1", "BDI_T1", "BAI_T1"] + ISI_ITEMS
    quoted = ", ".join(f'"{c}"' for c in needed)
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql(
        f"SELECT {quoted} FROM isi_raw_data WHERE ISI_T1 IS NOT NULL AND ISI_T1_Q1 IS NOT NULL",
        conn,
    )
    conn.close()
    for c in df.columns:
        if c != "NEW_NUMBER":
            df[c] = pd.to_numeric(df[c], errors="coerce")

    from _dynamic_n import restrict_to_v22

    df = restrict_to_v22(df)
    df_ins = df[df["ISI_T1"] >= 9].dropna(subset=ISI_ITEMS).copy()
    df_ctrl = df[df["ISI_T1"] < 9].dropna(subset=["ISI_T1"]).copy()

    imp = SimpleImputer(strategy="median")
    scl = StandardScaler()
    X = scl.fit_transform(imp.fit_transform(df_ins[ISI_ITEMS]))
    km = KMeans(n_clusters=3, random_state=SEED, n_init=100)
    df_ins["cluster_raw"] = km.fit_predict(X)

    order = df_ins.groupby("cluster_raw")["ISI_T1"].mean().sort_values().index.tolist()
    q23 = {
        ci: df_ins.loc[df_ins["cluster_raw"] == ci, ["ISI_T1_Q2", "ISI_T1_Q3"]].mean().sum()
        for ci in order
    }
    maint = max(q23, key=q23.get)
    sev = order[2]
    if maint == sev:
        q57 = {
            ci: df_ins.loc[df_ins["cluster_raw"] == ci, ["ISI_T1_Q5", "ISI_T1_Q7"]].mean().sum()
            for ci in order
        }
        sev = max(q57, key=q57.get)
    light = [ci for ci in order if ci not in (maint, sev)][0]

    cn_map = {maint: "睡眠維持/早醒型", light: "輕-日常型", sev: "嚴重-高困擾型"}
    df_ins["group_cn"] = df_ins["cluster_raw"].map(cn_map)
    df_ctrl["group_cn"] = "非失眠對照"
    out = pd.concat([df_ctrl, df_ins], ignore_index=True)
    out["group"] = out["group_cn"].map(CODE_MAP)
    return out[["NEW_NUMBER", "group", "AGE", "SEX", "BDI_T1", "BAI_T1"]]


def fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    if len(pvals) == 0:
        return np.array([])
    mask = ~np.isnan(pvals)
    out = np.full_like(pvals, np.nan)
    if mask.sum() == 0:
        return out
    out[mask] = multipletests(pvals[mask], method="fdr_bh")[1]
    return out


def hedges_g(a, b):
    """Hedges' g (weighted pooled SD + J correction). Sign: positive means b > a.

    Note (v2.2): 與 make_v22_stats_package.py 同公式 — paper / 補充資料正式
    效果量來源是 S30 master，此處輸出僅作為診斷用途，請勿直接引用此處數字。
    """
    a, b = pd.Series(a).dropna().astype(float), pd.Series(b).dropna().astype(float)
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return np.nan
    sp = np.sqrt(((n1 - 1) * a.std(ddof=1) ** 2 + (n2 - 1) * b.std(ddof=1) ** 2) / (n1 + n2 - 2))
    if sp <= 0:
        return np.nan
    d = (b.mean() - a.mean()) / sp
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    return d * J


def kw_pairwise(df, feature):
    vals = {g: df.loc[df["group"] == g, feature].dropna().values for g in GROUP_ORDER}
    vals = {g: v for g, v in vals.items() if len(v) >= 3}
    if len(vals) < 2:
        return None, []
    H, p_kw = stats.kruskal(*vals.values())
    n_total = sum(len(v) for v in vals.values())
    k = len(vals)
    eta2 = max((H - k + 1) / (n_total - k), 0)
    omnibus = {"feature": feature, "n": n_total, "H": H, "p": p_kw, "eta2": eta2}
    pairs = []
    raw_ps = []
    for g1, g2 in combinations(GROUP_ORDER, 2):
        if g1 not in vals or g2 not in vals:
            continue
        u, p = stats.mannwhitneyu(vals[g1], vals[g2], alternative="two-sided")
        g = hedges_g(vals[g1], vals[g2])
        pairs.append(
            {
                "feature": feature,
                "g1": g1,
                "g2": g2,
                "n1": len(vals[g1]),
                "n2": len(vals[g2]),
                "U": u,
                "p": p,
                "hedges_g": g,
            }
        )
        raw_ps.append(p)
    p_fdr = fdr(raw_ps)
    for i in range(len(pairs)):
        pairs[i]["p_fdr"] = p_fdr[i]
    return omnibus, pairs


def run_ancova(df, feature):
    dat = df.dropna(subset=[feature, "group", "BDI_T1", "BAI_T1", "AGE", "SEX"]).copy()
    if len(dat) < 40 or dat["group"].nunique() < 2:
        return []
    formula = f"Q('{feature}') ~ C(group, Treatment('HC')) + BDI_T1 + BAI_T1 + AGE + C(SEX)"
    try:
        m = smf.ols(formula, data=dat).fit()
        aov = anova_lm(m, typ=2)
        sse = aov.loc["Residual", "sum_sq"]
        p_eta2 = aov.drop("Residual")["sum_sq"] / (aov.drop("Residual")["sum_sq"] + sse)
        rows = []
        for term in aov.index:
            if term == "Residual":
                continue
            rows.append(
                {
                    "feature": feature,
                    "term": term,
                    "n": len(dat),
                    "F": aov.loc[term, "F"],
                    "p": aov.loc[term, "PR(>F)"],
                    "partial_eta2": p_eta2[term],
                }
            )
        for label, key in [("MEW", "MEW"), ("MD", "MD"), ("SD", "SD")]:
            t = f"C(group, Treatment('HC'))[T.{key}]"
            if t in m.params.index:
                rows.append(
                    {
                        "feature": feature,
                        "term": f"{label}_vs_HC_beta",
                        "n": len(dat),
                        "beta": m.params[t],
                        "p": m.pvalues[t],
                    }
                )
        return rows
    except Exception as e:
        return [{"feature": feature, "term": "ERROR", "error": str(e)}]


def main():
    psd_files = sorted(PSD_DIR.glob("psd_*.mat"))
    psd_files = [p for p in psd_files if p.stem.replace("psd_", "") not in EXCLUDED]
    print(f"Found {len(psd_files)} PSD .mat files")

    long_rows = []
    failed = []
    for i, p in enumerate(psd_files, 1):
        sid = p.stem.replace("psd_", "")
        try:
            freqs, psd, ch_names, n_ep = load_psd_mat(p)
            rows = fit_subject(sid, freqs, psd, ch_names, n_ep)
            long_rows.extend(rows)
            if i % 25 == 0 or i == len(psd_files):
                print(f"  [{i:3d}/{len(psd_files)}] {sid} ok (n_ep={n_ep})")
        except Exception as e:
            failed.append((sid, str(e)))
            print(f"  [{i:3d}/{len(psd_files)}] {sid} FAILED: {e}")

    df_long = pd.DataFrame(long_rows)
    df_long.to_csv(OUT_DIR / "01_specparam_features.csv", index=False)
    print(f"\nLong-format: {df_long.shape}")

    fq = (
        df_long.groupby("subject")
        .agg(
            r2_mean=("r2", "mean"),
            r2_min=("r2", "min"),
            fit_error_mean=("fit_error", "mean"),
            n_epochs=("n_epochs", "first"),
        )
        .reset_index()
    )
    fq.to_csv(OUT_DIR / "02_fit_quality.csv", index=False)
    print(f"Fit quality: mean R² = {fq['r2_mean'].mean():.4f} | min R² = {fq['r2_min'].min():.4f}")
    print(
        f"Epochs used: mean = {fq['n_epochs'].mean():.1f}, min = {fq['n_epochs'].min()}, max = {fq['n_epochs'].max()}"
    )

    df_subj = aggregate_features(df_long)
    df_subj.to_csv(OUT_DIR / "03_subject_brain_avg.csv", index=False)
    print(f"Subject-level features: {df_subj.shape}")

    groups = load_groups()
    df_subj["NEW_NUMBER"] = df_subj["subject"]
    df = df_subj.merge(groups, on="NEW_NUMBER", how="inner")
    df["group"] = pd.Categorical(df["group"], categories=GROUP_ORDER, ordered=True)
    print(f"\nMerged with groups: n = {len(df)}")
    print(df["group"].value_counts())

    KEY_FEATURES = [
        "exponent_BRAIN",
        "offset_BRAIN",
        "exponent_FRONTAL",
        "exponent_CENTRAL",
        "exponent_PARIETAL",
        "exponent_OCCIPITAL",
        "exponent_TEMPORAL",
        "beta_pw_BRAIN",
        "beta_cf_BRAIN",
        "beta_bw_BRAIN",
        "beta_pw_FRONTAL",
        "beta_pw_CENTRAL",
        "alpha_pw_BRAIN",
        "alpha_cf_BRAIN",
        "alpha_bw_BRAIN",
        "theta_pw_BRAIN",
    ]

    omnibus_rows, pairwise_rows, ancova_rows = [], [], []
    for feat in KEY_FEATURES:
        if feat not in df.columns:
            continue
        o, ps = kw_pairwise(df, feat)
        if o is not None:
            omnibus_rows.append(o)
            pairwise_rows.extend(ps)
        ancova_rows.extend(run_ancova(df, feat))

    odf = pd.DataFrame(omnibus_rows)
    if not odf.empty:
        odf["p_fdr"] = fdr(odf["p"].values)
        odf = odf.sort_values("eta2", ascending=False)
    odf.to_csv(OUT_DIR / "04_group_kw.csv", index=False)
    pd.DataFrame(pairwise_rows).to_csv(OUT_DIR / "05_pairwise_fdr.csv", index=False)
    pd.DataFrame(ancova_rows).to_csv(OUT_DIR / "06_ancova_sex_age.csv", index=False)

    lines = []
    P = lines.append
    P("=" * 80)
    P(f"Task B Phase 2 — FOOOF (Full sample, n_merged={len(df)})")
    P("=" * 80)
    P(f"PSD source: {PSD_DIR}")
    P(
        f"FOOOF settings: freq={FREQ_RANGE} | aperiodic={APERIODIC_MODE} | "
        f"peaks=[{PEAK_WIDTH_LIMITS}, n≤{MAX_N_PEAKS}, height≥{MIN_PEAK_HEIGHT}, thresh={PEAK_THRESHOLD}σ]"
    )
    P("")
    P("Group composition:")
    P(df["group"].value_counts().to_string())
    P("")
    P(f"Fit quality: mean R² = {fq['r2_mean'].mean():.4f}, min R² = {fq['r2_min'].min():.4f}")
    P("")
    P("─" * 80)
    P("Omnibus KW (sorted by η²):")
    P("─" * 80)
    P(odf.to_string(index=False) if not odf.empty else "(empty)")
    P("")
    pdf = pd.DataFrame(pairwise_rows)
    P("─" * 80)
    P("SD vs HC pairwise (核心):")
    P("─" * 80)
    if not pdf.empty:
        sd_hc = pdf[(pdf["g1"] == "HC") & (pdf["g2"] == "SD")][
            ["feature", "n1", "n2", "hedges_g", "p", "p_fdr"]
        ]
        P(sd_hc.to_string(index=False))
    P("")
    P("─" * 80)
    P("MEW vs HC pairwise (解離型態核心):")
    P("─" * 80)
    if not pdf.empty:
        mew_hc = pdf[(pdf["g1"] == "HC") & (pdf["g2"] == "MEW")][
            ["feature", "n1", "n2", "hedges_g", "p", "p_fdr"]
        ]
        P(mew_hc.to_string(index=False))
    P("")
    P("─" * 80)
    P("ANCOVA group partial η² × p（控制 BDI/BAI/AGE/SEX）:")
    P("─" * 80)
    adf = pd.DataFrame(ancova_rows)
    if not adf.empty:
        gp = adf[adf["term"].astype(str).str.contains("group", na=False)][
            ["feature", "F", "p", "partial_eta2"]
        ]
        P(gp.to_string(index=False))
        P("")
        P("SD vs HC β（ANCOVA-adjusted）:")
        sd_b = adf[adf["term"] == "SD_vs_HC_beta"][["feature", "beta", "p"]]
        P(sd_b.to_string(index=False))
        P("")
        P("MEW vs HC β:")
        mew_b = adf[adf["term"] == "MEW_vs_HC_beta"][["feature", "beta", "p"]]
        P(mew_b.to_string(index=False))

    if failed:
        P("")
        P(f"Failed: {len(failed)}")
        for sid, e in failed[:10]:
            P(f"  {sid}: {e}")

    txt = "\n".join(lines)
    (OUT_DIR / "07_summary.txt").write_text(txt, encoding="utf-8")
    print("\n" + txt)


if __name__ == "__main__":
    main()
