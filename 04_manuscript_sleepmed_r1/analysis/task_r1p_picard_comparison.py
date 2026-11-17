"""比較 PICARD＋CleanLine 與先前 Infomax＋notch 的結果。

逐人報相關、ICC(A,1) 與一致界限，再用同一批人重算主要組間比較。
分組沿用 ISI，不因前處理版本或 EEG 結果改動。"""

from __future__ import annotations
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from statsmodels.stats.weightstats import ttost_ind

from _r1_common import AN, COV, OUTCOME, fmt_p, hedges_g, load_all, icc_a1
from task_b_fooof_full import aggregate_features, fit_subject, load_psd_mat
from task_f_recompute_relpower import compute_subject

warnings.filterwarnings("ignore")
PIC = Path(
    os.environ.get("PSYML_ALT_PSD", str(AN / "_v22_infomax_N271_20261116" / "S10_fooof_psd_input"))
)  # 先前執行（v2.2）的頻譜檔
ORIG_PSD = AN / "S10_fooof_psd_input"
OUT = AN / "S_R1p_picard"
OUT.mkdir(parents=True, exist_ok=True)
MEAS = [
    (OUTCOME, "Periodic beta peak power, whole-scalp"),
    ("Beta_rel_BRAIN", "Relative beta power (13–30 Hz)"),
    ("exponent_BRAIN", "Aperiodic exponent"),
    ("Alpha_rel_BRAIN", "Relative alpha power (8–13 Hz)"),
]
SESOI = 0.37


def results(d, col, label):
    hc, lo, hi = (d.loc[d.k2 == g, col].dropna() for g in ("HC", "Lower", "Higher"))
    allins = d.loc[d.k2 != "HC", col].dropna()
    r = dict(
        data=label, n=int(d[col].notna().sum()), n_hc=len(hc), n_lower=len(lo), n_higher=len(hi)
    )
    for nm, a, b in [
        ("higher", hi, hc),
        ("lower", lo, hc),
        ("hi_vs_lo", hi, lo),
        ("pooled", allins, hc),
    ]:
        g, l, h, _, _ = hedges_g(a, b)
        r.update({f"g_{nm}": g, f"g_{nm}_lo": l, f"g_{nm}_hi": h})
    r["p_higher"] = float(mannwhitneyu(hi, hc).pvalue)
    r["p_hi_vs_lo"] = float(mannwhitneyu(hi, lo).pvalue)
    for tag, extra in [("M1", ""), ("M2", " + BDI_T1 + BAI_T1")]:
        m = smf.ols(f"{col} ~ C(k2, Treatment('HC')) + {COV}{extra}", d).fit()
        k = [c for c in m.params.index if "Higher" in c][0]
        r.update({f"{tag}_b": m.params[k], f"{tag}_p": m.pvalues[k]})
    dd = d.assign(thr=(d.ISI_T1 >= 18).astype(int))
    m = smf.ols(f"{col} ~ thr + {COV}", dd).fit()
    r.update(thr_b=m.params["thr"], thr_p=m.pvalues["thr"])
    sp = np.sqrt(
        ((len(lo) - 1) * lo.var(ddof=1) + (len(hc) - 1) * hc.var(ddof=1)) / (len(lo) + len(hc) - 2)
    )
    r["tost_p"] = float(ttost_ind(lo, hc, -SESOI * sp, SESOI * sp, usevar="pooled")[0])
    return r


def main():
    df = load_all()
    pic_ids = {p.stem.replace("psd_", "") for p in PIC.glob("psd_*.mat")}
    orig_ids = {p.stem.replace("psd_", "") for p in ORIG_PSD.glob("psd_*.mat")}
    samp = pd.DataFrame(
        [
            dict(
                recordings_passing_analysed=len(orig_ids),
                recordings_passing_earlier=len(pic_ids),
                in_both=len(orig_ids & pic_ids),
                only_analysed=len(orig_ids - pic_ids),
                only_earlier=len(pic_ids - orig_ids),
                analytic_sample=len(df),
                analytic_sample_with_earlier=int(df.subject.isin(pic_ids).sum()),
            )
        ]
    )
    samp.to_csv(OUT / "01_sample.csv", index=False)
    long, rel = [], []
    for sid in df.subject:
        if sid not in pic_ids:
            continue
        freqs, psd, ch, n_ep = load_psd_mat(PIC / f"psd_{sid}.mat")
        long += fit_subject(sid, freqs, psd, ch, n_ep)
        rel.append(compute_subject(sid, freqs, psd, ch))
    pic = aggregate_features(pd.DataFrame(long)).merge(pd.DataFrame(rel), on="subject")
    cols = [c for c, _ in MEAS]
    d = df.merge(
        pic[["subject"] + cols].rename(columns={c: c + "_earlier" for c in cols}),
        on="subject",
        how="left",
    )
    d[["subject"] + cols + [c + "_earlier" for c in cols]].to_csv(
        OUT / "_subject_level_internal.csv", index=False
    )
    ag = []
    for c, lab in MEAS:
        s = d[[c, c + "_earlier"]].dropna()
        x, y = s[c].to_numpy(), s[c + "_earlier"].to_numpy()
        diff = y - x
        ag.append(
            dict(
                measure=lab,
                n=len(s),
                pearson_r=pearsonr(x, y)[0],
                spearman_rho=spearmanr(x, y)[0],
                icc_a1=icc_a1(x, y),
                mean_analysed=x.mean(),
                mean_earlier=y.mean(),
                mean_diff=diff.mean(),
                loa_low=diff.mean() - 1.96 * diff.std(ddof=1),
                loa_high=diff.mean() + 1.96 * diff.std(ddof=1),
            )
        )
    ag = pd.DataFrame(ag)
    ag.to_csv(OUT / "02_agreement.csv", index=False)
    both = d.dropna(subset=[OUTCOME, OUTCOME + "_earlier"])
    res = pd.DataFrame(
        [
            results(both, OUTCOME, "PICARD + CleanLine (analysed)"),
            results(both, OUTCOME + "_earlier", "Extended Infomax + notch (earlier run)"),
        ]
    )
    res.to_csv(OUT / "03_main_results.csv", index=False)

    s0 = samp.iloc[0]
    L = [
        "# S_R1p — 兩次前處理的對照（PICARD＋CleanLine 對 extended Infomax＋notch）",
        "",
        f"- 通過 60 個 epoch 標準的記錄：分析資料 {s0.recordings_passing_analysed} 筆、先前的執行 {s0.recordings_passing_earlier} 筆；兩者皆通過 {s0.in_both} 筆（只在分析資料 {s0.only_analysed}、只在先前的執行 {s0.only_earlier}）。",
        f"- 分析樣本 {s0.analytic_sample} 人中，先前的執行有資料者 {s0.analytic_sample_with_earlier} 人。",
        "",
        "## 逐人一致性",
        "",
        "| 指標 | n | Pearson r | Spearman rho | ICC(A,1) | 分析資料平均 | 先前執行平均 | 平均差（先前 − 分析） | 95% 一致界限 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in ag.itertuples():
        L.append(
            f"| {r.measure} | {r.n} | {r.pearson_r:.3f} | {r.spearman_rho:.3f} | {r.icc_a1:.3f} | {r.mean_analysed:.3f} | {r.mean_earlier:.3f} | {r.mean_diff:+.3f} | [{r.loa_low:+.3f}, {r.loa_high:+.3f}] |"
        )
    L += [
        "",
        "## 主要結果（週期性 beta；同一批人、同樣的分組）",
        "",
        "| 項目 | " + " | ".join(res.data) + " |",
        "|---|---|---|",
    ]

    def row(name, f):
        L.append(f"| {name} | " + " | ".join(f(r) for r in res.itertuples()) + " |")

    row("n（健康／低／高）", lambda r: f"{r.n}（{r.n_hc}／{r.n_lower}／{r.n_higher}）")
    row(
        "高嚴重度對健康組 g [95% CI]",
        lambda r: f"{r.g_higher:+.2f} [{r.g_higher_lo:+.2f}, {r.g_higher_hi:+.2f}]",
    )
    row(
        "低嚴重度對健康組 g [95% CI]",
        lambda r: f"{r.g_lower:+.2f} [{r.g_lower_lo:+.2f}, {r.g_lower_hi:+.2f}]",
    )
    row("高對低 g（p）", lambda r: f"{r.g_hi_vs_lo:+.2f}（{fmt_p(r.p_hi_vs_lo)}）")
    row("合併對健康組 g", lambda r: f"{r.g_pooled:+.2f}")
    row("M1 高嚴重度 b（p）", lambda r: f"{r.M1_b:.3f}（{fmt_p(r.M1_p)}）")
    row("M2 高嚴重度 b（p）", lambda r: f"{r.M2_b:.3f}（{fmt_p(r.M2_p)}）")
    row("ISI ≥ 18 係數（p）", lambda r: f"{r.thr_b:.3f}（{fmt_p(r.thr_p)}）")
    row("等效檢定 p", lambda r: fmt_p(r.tost_p))
    (OUT / "04_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
