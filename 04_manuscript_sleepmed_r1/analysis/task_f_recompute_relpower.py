"""從同一批 PSD 重算相對功率，修正舊分母重複加總子頻段的問題。

分母只加 Delta／Theta／Alpha／Beta／Gamma 五段（1–45 Hz），
Beta1 等子頻段只當分子。保留舊值作對照，輸出到 S17_relative_power。"""

from __future__ import annotations

import sqlite3
import warnings
from itertools import combinations
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import loadmat
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

BASE_DIR = Path(__file__).parent
PSD_DIR = BASE_DIR / "outputs" / "補充資料" / "analysis" / "S10_fooof_psd_input"
OUT_DIR = BASE_DIR / "outputs" / "補充資料" / "analysis" / "S17_relative_power"
OUT_DIR.mkdir(parents=True, exist_ok=True)

from _dynamic_n import DB_PATH as _DB, EXCLUDED

DB_PATH = Path(_DB)
SEED = 42
ISI_ITEMS = [f"ISI_T1_Q{i}" for i in range(1, 8)]

GROUP_ORDER = ["HC", "MEW", "MD", "SD"]
GROUP_ORDER_CN = ["非失眠對照", "睡眠維持/早醒型", "輕-日常型", "嚴重-高困擾型"]
CODE_MAP = dict(zip(GROUP_ORDER_CN, GROUP_ORDER))

# 文獻標準頻段（使用者拍板）
PRIMARY_BANDS = {  # 進分母（連續非重疊）
    "Delta": (1.0, 4.0),
    "Theta": (4.0, 8.0),
    "Alpha": (8.0, 13.0),
    "Beta": (13.0, 30.0),
    "Gamma": (30.0, 45.0),
}
DERIVED_BANDS = {  # 不進分母，僅算 子/total
    "Alpha1": (8.0, 10.0),
    "Alpha2": (10.0, 13.0),
    "Beta1": (13.0, 18.0),
    "Beta2": (18.0, 24.0),
    "Beta3": (24.0, 30.0),
    "HighBeta": (30.0, 40.0),
}
ALL_BANDS = {**PRIMARY_BANDS, **DERIVED_BANDS}

CH_RENAME = {"T7": "T3", "T8": "T4", "P7": "T5", "P8": "T6"}
BRAIN_REGIONS = {
    "FRONTAL": ["FP1", "FP2", "F3", "F4", "F7", "F8", "FZ"],
    "CENTRAL": ["C3", "C4", "CZ"],
    "PARIETAL": ["P3", "P4", "PZ"],
    "TEMPORAL": ["T3", "T4", "T5", "T6"],
    "OCCIPITAL": ["O1", "O2"],
}
ALL_19 = sorted({c for chs in BRAIN_REGIONS.values() for c in chs})


def band_abs_power(freqs, psd_1d, lo, hi):
    """梯形積分 PSD 於 [lo, hi]。"""
    m = (freqs >= lo) & (freqs <= hi)
    if m.sum() < 2:
        return np.nan
    return np.trapz(psd_1d[m], freqs[m])


def load_psd(path):
    m = loadmat(str(path), squeeze_me=True)
    freqs = np.asarray(m["freqs"], dtype=float).ravel()
    psd = np.asarray(m["psd"], dtype=float)  # [n_ch × n_freq]
    ch = m["ch_names"]
    ch = [str(c).strip() for c in ch] if isinstance(ch, np.ndarray) else [str(ch)]
    ch = [CH_RENAME.get(c.upper(), c.upper()) for c in ch]
    return freqs, psd, ch


def compute_subject(sid, freqs, psd, ch_names):
    """回傳 subject-level dict：各 band 之 BRAIN + region 相對功率 + 比值。"""
    ch_abs = {}  # ch -> {band: abs}
    for ci, ch in enumerate(ch_names):
        p = psd[ci]
        bp = {b: band_abs_power(freqs, p, lo, hi) for b, (lo, hi) in ALL_BANDS.items()}
        # 分母固定用五個主要頻帶的總和（1–45 Hz）；子頻段不重複計入。
        total = np.nansum([bp[b] for b in PRIMARY_BANDS])
        ch_abs[ch] = {"_total": total, **bp}

    def region_rel(chs, band):
        """region 平均之相對功率 = mean_ch(band_abs)/mean_ch(total)*100（pooled）。"""
        absvals, totvals = [], []
        for ch in chs:
            if ch in ch_abs and not np.isnan(ch_abs[ch][band]) and ch_abs[ch]["_total"] > 0:
                absvals.append(ch_abs[ch][band])
                totvals.append(ch_abs[ch]["_total"])
        if not absvals:
            return np.nan
        return float(np.mean(absvals) / np.mean(totvals) * 100)

    def region_abs(chs, band):
        vals = [ch_abs[ch][band] for ch in chs if ch in ch_abs and not np.isnan(ch_abs[ch][band])]
        return float(np.mean(vals)) if vals else np.nan

    row = {"subject": sid}
    for band in ALL_BANDS:
        row[f"{band}_rel_BRAIN"] = region_rel(ALL_19, band)
        row[f"{band}_abs_BRAIN"] = region_abs(ALL_19, band)
    for region, chs in BRAIN_REGIONS.items():
        row[f"Beta_rel_{region}"] = region_rel(chs, band="Beta")
    # 比值（用 abs，total 抵消）
    a = {b: region_abs(ALL_19, b) for b in ALL_BANDS}

    def safe(n, d):
        return float(n / d) if (d and not np.isnan(n) and not np.isnan(d) and d != 0) else np.nan

    row["TBR_BRAIN"] = safe(a["Theta"], a["Beta"])
    row["ABR_BRAIN"] = safe(a["Alpha"], a["Beta"])
    row["TAR_BRAIN"] = safe(a["Theta"], a["Alpha"])
    row["BAR_BRAIN"] = safe(a["Beta"], a["Alpha"])
    row["BDR_BRAIN"] = safe(a["Beta"], a["Delta"])
    row["DBR_BRAIN"] = safe(a["Delta"], a["Beta"])
    row["DTABR_BRAIN"] = safe(a["Delta"] + a["Theta"], a["Alpha"] + a["Beta"])
    return row


def load_groups(valid_sids=None):
    quoted = ", ".join(
        f'"{c}"' for c in ["NEW_NUMBER", "AGE", "SEX", "ISI_T1", "BDI_T1", "BAI_T1"] + ISI_ITEMS
    )
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql(
        f"SELECT {quoted} FROM isi_raw_data WHERE ISI_T1 IS NOT NULL AND ISI_T1_Q1 IS NOT NULL",
        conn,
    )
    # OLD 相對 Beta（DB 重疊分母版）供對照
    old = pd.read_sql(
        "SELECT NEW_NUMBER, EEG_PWR_REL_BETA_BRAIN_AVG AS Beta_rel_OLD "
        "FROM isi_raw_data_recalc_3s60 WHERE EEG_PWR_REL_BETA_BRAIN_AVG IS NOT NULL",
        conn,
    )
    conn.close()
    for c in df.columns:
        if c != "NEW_NUMBER":
            df[c] = pd.to_numeric(df[c], errors="coerce")
    old["Beta_rel_OLD"] = pd.to_numeric(old["Beta_rel_OLD"], errors="coerce")

    # 分群限於 EEG 與 ISI 都完整的樣本，不能直接用整份資料庫。
    if valid_sids is not None:
        before = len(df)
        df = df[df["NEW_NUMBER"].isin(valid_sids)].copy()
        print(f"  [v2.2 new methodology] K-means/HC 限制到 EEG∩DB: {before} → {len(df)}")
    df_ins = df[df["ISI_T1"] >= 9].dropna(subset=ISI_ITEMS).copy()
    df_ctrl = df[df["ISI_T1"] < 9].dropna(subset=["ISI_T1"]).copy()
    imp, scl = SimpleImputer(strategy="median"), StandardScaler()
    X = scl.fit_transform(imp.fit_transform(df_ins[ISI_ITEMS]))
    km = KMeans(n_clusters=3, random_state=SEED, n_init=100)
    df_ins["c"] = km.fit_predict(X)
    order = df_ins.groupby("c")["ISI_T1"].mean().sort_values().index.tolist()
    q23 = {
        ci: df_ins.loc[df_ins["c"] == ci, ["ISI_T1_Q2", "ISI_T1_Q3"]].mean().sum() for ci in order
    }
    maint = max(q23, key=q23.get)
    sev = order[2]
    if maint == sev:
        q57 = {
            ci: df_ins.loc[df_ins["c"] == ci, ["ISI_T1_Q5", "ISI_T1_Q7"]].mean().sum()
            for ci in order
        }
        sev = max(q57, key=q57.get)
    light = [ci for ci in order if ci not in (maint, sev)][0]
    cn = {maint: "睡眠維持/早醒型", light: "輕-日常型", sev: "嚴重-高困擾型"}
    df_ins["group_cn"] = df_ins["c"].map(cn)
    df_ctrl["group_cn"] = "非失眠對照"
    out = pd.concat([df_ctrl, df_ins], ignore_index=True)
    out["group"] = out["group_cn"].map(CODE_MAP)
    out = out.merge(old, on="NEW_NUMBER", how="left")
    return out[["NEW_NUMBER", "group", "AGE", "SEX", "BDI_T1", "BAI_T1", "ISI_T1", "Beta_rel_OLD"]]


def hedges_g(a, b):
    """Hedges' g (weighted pooled SD + J correction). Sign: positive means b > a.

    Note (v2.2): 與 make_v22_stats_package.py 同公式 — paper / 補充資料正式
    效果量來源是 S30 master，此處輸出僅作為診斷用途，請勿直接引用此處數字。
    """
    a, b = pd.Series(a).dropna(), pd.Series(b).dropna()
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return np.nan
    sp = np.sqrt(((n1 - 1) * a.std(ddof=1) ** 2 + (n2 - 1) * b.std(ddof=1) ** 2) / (n1 + n2 - 2))
    if sp <= 0:
        return np.nan
    d = (b.mean() - a.mean()) / sp
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    return d * J


def fdr(p):
    p = np.asarray(p, float)
    m = ~np.isnan(p)
    o = np.full_like(p, np.nan)
    if m.sum():
        o[m] = multipletests(p[m], method="fdr_bh")[1]
    return o


def main():
    psd_files = sorted(PSD_DIR.glob("psd_*.mat"))
    psd_files = [p for p in psd_files if p.stem.replace("psd_", "") not in EXCLUDED]
    print(f"PSD .mat: {len(psd_files)}")

    rows = []
    for i, p in enumerate(psd_files, 1):
        sid = p.stem.replace("psd_", "")
        try:
            f, psd, ch = load_psd(p)
            rows.append(compute_subject(sid, f, psd, ch))
        except Exception as e:
            print(f"  {sid} FAIL: {e}")
        if i % 50 == 0:
            print(f"  [{i}/{len(psd_files)}]")
    sub = pd.DataFrame(rows)
    sub.to_csv(OUT_DIR / "01_relpower_subject.csv", index=False)
    print(f"Subject-level: {sub.shape}")

    # 分群只用 EEG 有效樣本，和主分析一致。
    valid_sids = set(sub["subject"].tolist())
    g = load_groups(valid_sids=valid_sids)
    sub["NEW_NUMBER"] = sub["subject"]
    df = sub.merge(g, on="NEW_NUMBER", how="inner")
    df["group"] = pd.Categorical(df["group"], categories=GROUP_ORDER, ordered=True)
    print(f"Merged: n={len(df)}\n{df['group'].value_counts().to_string()}")

    # 4 群 KW + pairwise
    KEY = [
        f"{b}_rel_BRAIN"
        for b in [
            "Delta",
            "Theta",
            "Alpha",
            "Beta",
            "Gamma",
            "Beta1",
            "Beta2",
            "Beta3",
            "HighBeta",
            "Alpha1",
            "Alpha2",
        ]
    ]
    KEY += [
        "TBR_BRAIN",
        "ABR_BRAIN",
        "BDR_BRAIN",
        "Beta_rel_FRONTAL",
        "Beta_rel_CENTRAL",
        "Beta_rel_PARIETAL",
        "Beta_rel_TEMPORAL",
        "Beta_rel_OCCIPITAL",
    ]
    omni, pair = [], []
    for feat in KEY:
        if feat not in df.columns:
            continue
        vals = {gg: df.loc[df["group"] == gg, feat].dropna().values for gg in GROUP_ORDER}
        vals = {gg: v for gg, v in vals.items() if len(v) >= 3}
        if len(vals) < 2:
            continue
        H, pv = stats.kruskal(*vals.values())
        N = sum(len(v) for v in vals.values())
        k = len(vals)
        eta2 = max((H - k + 1) / (N - k), 0)
        omni.append({"feature": feat, "n": N, "H": H, "p": pv, "eta2": eta2})
        rp = []
        for g1, g2 in combinations(GROUP_ORDER, 2):
            if g1 in vals and g2 in vals:
                u, pp = stats.mannwhitneyu(vals[g1], vals[g2], alternative="two-sided")
                g = hedges_g(vals[g1], vals[g2])
                pair.append(
                    {
                        "feature": feat,
                        "g1": g1,
                        "g2": g2,
                        "n1": len(vals[g1]),
                        "n2": len(vals[g2]),
                        "hedges_g": g,
                        "p": pp,
                    }
                )
                rp.append(pp)
        pf = fdr(rp)
        for j, idx in enumerate([x for x in range(len(pair)) if pair[x]["feature"] == feat]):
            pair[idx]["p_fdr"] = pf[j]
    odf = pd.DataFrame(omni)
    if not odf.empty:
        odf["p_fdr"] = fdr(odf["p"].values)
        odf = odf.sort_values("eta2", ascending=False)
    odf.to_csv(OUT_DIR / "02_group_kw.csv", index=False)
    pd.DataFrame(pair).to_csv(OUT_DIR / "03_pairwise_fdr.csv", index=False)

    # OLD vs NEW Beta 對照
    cmp_rows = []
    for gg in GROUP_ORDER:
        sg = df[df["group"] == gg]
        cmp_rows.append(
            {
                "group": gg,
                "n": len(sg),
                "Beta_rel_OLD_mean": sg["Beta_rel_OLD"].mean(),
                "Beta_rel_OLD_sd": sg["Beta_rel_OLD"].std(ddof=1),
                "Beta_rel_NEW_mean": sg["Beta_rel_BRAIN"].mean(),
                "Beta_rel_NEW_sd": sg["Beta_rel_BRAIN"].std(ddof=1),
            }
        )
    cmpdf = pd.DataFrame(cmp_rows)
    # d vs HC for OLD & NEW
    hc = df[df["group"] == "HC"]
    for gg in ["MEW", "MD", "SD"]:
        sg = df[df["group"] == gg]
        cmpdf.loc[cmpdf["group"] == gg, "g_vs_HC_OLD"] = hedges_g(
            hc["Beta_rel_OLD"], sg["Beta_rel_OLD"]
        )
        cmpdf.loc[cmpdf["group"] == gg, "g_vs_HC_NEW"] = hedges_g(
            hc["Beta_rel_BRAIN"], sg["Beta_rel_BRAIN"]
        )
    cmpdf.to_csv(OUT_DIR / "04_old_vs_new.csv", index=False)

    # Summary
    L = []
    P = L.append
    P("=" * 80)
    P("Task F — 相對功率重算（文獻標準 Beta 13–30, total = 1–45 連續積分）")
    P("=" * 80)
    P(f"來源 PSD: {len(psd_files)} 個 .mat | merged n={len(df)}")
    P(f"主頻段(進分母): {', '.join(f'{b}{v}' for b, v in PRIMARY_BANDS.items())}")
    P("")
    P("─" * 80)
    P("OLD（DB 重疊分母）vs NEW（1–45 連續積分）Beta 相對功率對照")
    P("─" * 80)
    P(cmpdf.to_string(index=False))
    P("")
    P("─" * 80)
    P("NEW 四群 KW（依 η² 排序）")
    P("─" * 80)
    P(odf.to_string(index=False))
    P("")
    P("─" * 80)
    P("NEW Beta_rel_BRAIN pairwise（核心）")
    P("─" * 80)
    pdf = pd.DataFrame(pair)
    if not pdf.empty:
        bp = pdf[pdf["feature"] == "Beta_rel_BRAIN"][
            ["g1", "g2", "n1", "n2", "hedges_g", "p", "p_fdr"]
        ]
        P(bp.to_string(index=False))
    txt = "\n".join(L)
    (OUT_DIR / "05_summary.txt").write_text(txt, encoding="utf-8")
    print("\n" + txt)


if __name__ == "__main__":
    main()
