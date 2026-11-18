"""用 surface Laplacian 頻譜檢查 beta 差異是否仍在。

和主分析使用同一批 60 個 epoch、同一套 SpecParam 設定。
19 電極的外圈缺少外側鄰居，估計較不準，解讀以中央與內圈為主。"""

from __future__ import annotations
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import pearsonr, spearmanr

from _r1_common import AN, COV, COV_COLS, OUTCOME, fmt_p, hedges_g, load_all, mwu_p
from task_b_fooof_full import fit_subject, load_psd_mat

warnings.filterwarnings("ignore")

from _dynamic_n import local_path

PIPE = Path(
    local_path("eeg_pipeline_root", "EEG_PIPELINE_ROOT")
)
LAP_DIR = PIPE / "outputs" / "psd_mat_laplacian"
REF_DIR = AN / "S10_fooof_psd_input"
OUT = AN / "S_R1t_surface_laplacian"
OUT.mkdir(parents=True, exist_ok=True)
THR = 18

SETS = {
    "central 3 (C3, Cz, C4)": ["C3", "CZ", "C4"],
    "inner 9 (F3, F4, Fz, C3, C4, Cz, P3, P4, Pz)": [
        "F3",
        "F4",
        "FZ",
        "C3",
        "C4",
        "CZ",
        "P3",
        "P4",
        "PZ",
    ],
    "outer 10 (Fp1, Fp2, F7, F8, T3, T4, T5, T6, O1, O2)": [
        "FP1",
        "FP2",
        "F7",
        "F8",
        "T3",
        "T4",
        "T5",
        "T6",
        "O1",
        "O2",
    ],
    "all 19": None,
}
SHORT = {
    "central 3 (C3, Cz, C4)": "c3",
    "inner 9 (F3, F4, Fz, C3, C4, Cz, P3, P4, Pz)": "in9",
    "outer 10 (Fp1, Fp2, F7, F8, T3, T4, T5, T6, O1, O2)": "out10",
    "all 19": "all",
}
MEAS = [
    ("peak", "Periodic beta peak power"),
    ("adj_13_30", "Aperiodic-adjusted power, 13–30 Hz"),
    ("adj_13_20", "Aperiodic-adjusted power, 13–20 Hz"),
    ("adj_20_30", "Aperiodic-adjusted power, 20–30 Hz"),
    ("adj_30_45", "Aperiodic-adjusted power, 30–45 Hz"),
    ("raw_30_45", "Raw log power, 30–45 Hz"),
]
VERS = [("ref", "Average reference (main analysis)"), ("lap", "Surface Laplacian")]


def subject_features(sid, src):
    """一位受試者、一種頻譜 → 各電極組的各項指標。"""
    freqs, psd, ch, n_ep = load_psd_mat(src / f"psd_{sid}.mat")
    fit = pd.DataFrame(fit_subject(sid, freqs, psd, ch, n_ep)).set_index("channel")
    sel = (freqs >= 1.0) & (freqs <= 45.0)
    f = freqs[sel]
    lp = np.log10(psd[:, sel])
    ap = (
        fit.loc[ch, "offset"].to_numpy()[:, None]
        - fit.loc[ch, "exponent"].to_numpy()[:, None] * np.log10(f)[None, :]
    )
    adj = lp - ap
    band = lambda spec, lo, hi: spec[:, (f >= lo) & ((f < hi) if hi < 45 else (f <= hi))].mean(1)
    per_ch = pd.DataFrame(
        dict(
            peak=fit.loc[ch, "beta_pw"].to_numpy(),
            adj_13_30=band(adj, 13, 30),
            adj_13_20=band(adj, 13, 20),
            adj_20_30=band(adj, 20, 30),
            adj_30_45=band(adj, 30, 45),
            raw_30_45=band(lp, 30, 45),
            r2=fit.loc[ch, "r2"].to_numpy(),
        ),
        index=ch,
    )
    row = {"subject": sid, "r2": per_ch["r2"].mean()}
    for lab, chs in SETS.items():
        g = per_ch if chs is None else per_ch.loc[[c for c in chs if c in per_ch.index]]
        for m, _ in MEAS:
            row[f"{m}_{SHORT[lab]}"] = g[m].dropna().mean()
    return row


def main():
    df = load_all().dropna(subset=COV_COLS).reset_index(drop=True)
    have = [s for s in df["subject"] if (LAP_DIR / f"psd_{s}.mat").exists()]
    if len(have) < len(df):
        print(f"注意：{len(df) - len(have)} 人沒有 step6 的輸出，未納入")
    df = df[df["subject"].isin(have)].reset_index(drop=True)
    feats = {}
    for v, src in [("ref", REF_DIR), ("lap", LAP_DIR)]:
        t = pd.DataFrame([subject_features(s, src) for s in df["subject"]])
        feats[v] = t.rename(columns={c: f"{v}_{c}" for c in t.columns if c != "subject"})
    d = df.merge(feats["ref"], on="subject").merge(feats["lap"], on="subject")
    d.to_csv(OUT / "_subject_level_internal.csv", index=False, encoding="utf-8-sig")

    # 核對：用同一條程式路徑重算的平均參考版，應與主分析的週期性 beta 相同
    chk = d[["ref_peak_all", OUTCOME]].dropna()
    r_chk = pearsonr(chk["ref_peak_all"], chk[OUTCOME])[0]
    d_chk = float((chk["ref_peak_all"] - chk[OUTCOME]).abs().max())

    hc = d[d["k2"] == "HC"]
    d["step"] = (d["ISI_T1"] >= THR).astype(int)
    rows = []
    for v, vlab in VERS:
        for lab in SETS:
            for m, mlab in MEAS:
                col = f"{v}_{m}_{SHORT[lab]}"
                x = d.dropna(subset=[col])
                sd = hc[col].std(ddof=1)
                m1 = smf.ols(f"{col} ~ C(k2, Treatment('HC')) + {COV}", x).fit()
                m2 = smf.ols(f"{col} ~ C(k2, Treatment('HC')) + {COV} + BDI_T1 + BAI_T1", x).fit()
                ms = smf.ols(f"{col} ~ step + {COV}", x[x["k2"] != "HC"]).fit()
                for g in ["Higher", "Lower"]:
                    a = x.loc[x["k2"] == g, col]
                    b = hc[col].dropna()
                    gg, lo, hi, n1, n2 = hedges_g(a, b)
                    k = [c for c in m1.params.index if g in c][0]
                    rows.append(
                        dict(
                            version=vlab,
                            electrodes=lab,
                            measure=mlab,
                            contrast=f"{g} vs HC",
                            n=n1,
                            n_hc=n2,
                            g=gg,
                            ci_low=lo,
                            ci_high=hi,
                            p_mwu=mwu_p(a, b),
                            m1_diff_in_HC_SD=m1.params[k] / sd,
                            m1_p=m1.pvalues[k],
                            m2_diff_in_HC_SD=m2.params[k] / sd,
                            m2_p=m2.pvalues[k],
                            step18_in_HC_SD=ms.params["step"] / sd if g == "Higher" else np.nan,
                            step18_p=ms.pvalues["step"] if g == "Higher" else np.nan,
                        )
                    )
    C = pd.DataFrame(rows)
    C.round(4).to_csv(OUT / "01_contrasts.csv", index=False, encoding="utf-8-sig")

    ag = []
    for lab in SETS:
        for m, mlab in MEAS[:2]:
            s = d[[f"ref_{m}_{SHORT[lab]}", f"lap_{m}_{SHORT[lab]}"]].dropna()
            ag.append(
                dict(
                    electrodes=lab,
                    measure=mlab,
                    n=len(s),
                    pearson_r=pearsonr(s.iloc[:, 0], s.iloc[:, 1])[0],
                    spearman_rho=spearmanr(s.iloc[:, 0], s.iloc[:, 1])[0],
                )
            )
    A = pd.DataFrame(ag)
    A.round(4).to_csv(OUT / "02_agreement.csv", index=False, encoding="utf-8-sig")

    sg = lambda x: f"{x:+.2f}"

    def line(v, lab, mlab):
        h = C[
            (C.version == v)
            & (C.electrodes == lab)
            & (C.measure == mlab)
            & (C.contrast == "Higher vs HC")
        ].iloc[0]
        l = C[
            (C.version == v)
            & (C.electrodes == lab)
            & (C.measure == mlab)
            & (C.contrast == "Lower vs HC")
        ].iloc[0]
        return (
            f"| {mlab} | {v} | {int(h.n)}／{int(h.n_hc)} | {sg(h.g)} [{sg(h.ci_low)}, {sg(h.ci_high)}] | {fmt_p(h.p_mwu)} | {sg(h.m1_diff_in_HC_SD)}（{fmt_p(h.m1_p)}） | "
            f"{sg(h.m2_diff_in_HC_SD)}（{fmt_p(h.m2_p)}） | {sg(h.step18_in_HC_SD)}（{fmt_p(h.step18_p)}） | {sg(l.g)} [{sg(l.ci_low)}, {sg(l.ci_high)}] |"
        )

    L = [
        "# S_R1t — 表面拉普拉斯轉換後的 beta",
        "",
        f"- 分析樣本 {len(d)} 人（健康組 {int((d.k2 == 'HC').sum())}、低嚴重度 {int((d.k2 == 'Lower').sum())}、高嚴重度 {int((d.k2 == 'Higher').sum())}）；每人 60 個 epoch，與主分析相同。",
        f"- 核對：用同一條程式路徑重算的平均參考版，全頭皮週期性 beta 與主分析的相關 r = {r_chk:.4f}，最大差 {d_chk:.4f}。",
        f"- 擬合品質（R² 平均）：平均參考 {d['ref_r2'].mean():.3f}，表面拉普拉斯 {d['lap_r2'].mean():.3f}。",
        "- 19 個電極的表面拉普拉斯在外圈電極較不準，主要看中央 3 個與內圈 9 個電極。",
        "",
    ]
    for lab in SETS:
        L += [
            f"## {lab}",
            "",
            "| 指標 | 資料 | n（高／健康） | 高嚴重度 g [95% CI] | p | M1 差距（p） | M2 差距（p） | ISI ≥ 18（p） | 低嚴重度 g [95% CI] |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for m, mlab in MEAS:
            for v, vlab in VERS:
                L.append(line(vlab, lab, mlab))
        L.append("")
    L += [
        "M1、M2 與 ISI ≥ 18 的係數以健康組標準差為單位；M1 校正年齡、性別、教育，M2 再加 BDI-II、BAI；ISI ≥ 18 的落差在失眠症狀者中估計。",
        "",
        "## 兩種資料的逐人相關",
        "",
        "| 電極 | 指標 | n | Pearson r | Spearman rho |",
        "|---|---|---|---|---|",
    ]
    for r in A.itertuples():
        L.append(
            f"| {r.electrodes} | {r.measure} | {r.n} | {r.pearson_r:.2f} | {r.spearman_rho:.2f} |"
        )
    (OUT / "03_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
