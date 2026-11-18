"""比較低、高 beta 子頻段，並直接檢定高段減低段的組間差。

不能用「一段顯著、另一段不顯著」來推論兩段不同。
同時計算去背景與原始 log 功率，另列峰頻率與峰值功率。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from _r1_common import AN, COV, COV_COLS, fmt_p, hedges_g, load_all, mwu_p
from _r1_spectra import band_mean, load_spectra

OUT = AN / "S_R1n_beta_subbands"
OUT.mkdir(parents=True, exist_ok=True)
BANDS = [
    ("13–20 Hz", 13, 20),
    ("20–30 Hz", 20, 30),
    ("13–18 Hz", 13, 18),
    ("18–24 Hz", 18, 24),
    ("24–30 Hz", 24, 30),
    ("13–30 Hz", 13, 30),
]
THR = 18


def main():
    df = load_all().dropna(subset=COV_COLS).reset_index(drop=True)
    freq, chn, adj, raw = load_spectra(df["subject"])
    for kind, spec in [("adj", adj), ("raw", raw)]:
        for lab, lo, hi in BANDS:
            df[f"{kind}_{lo}_{hi}"] = band_mean(spec, freq, lo, hi)
    df["adj_hi_minus_lo"] = df["adj_20_30"] - df["adj_13_20"]
    df["raw_hi_minus_lo"] = df["raw_20_30"] - df["raw_13_20"]
    df["thr"] = (df["ISI_T1"] >= THR).astype(int)
    hc = df[df["k2"] == "HC"]

    def contrast(col):
        rows = []
        m = smf.ols(f"{col} ~ C(k2, Treatment('HC')) + {COV}", df).fit()
        sd = hc[col].std(ddof=1)
        m2 = smf.ols(
            f"{col} ~ C(k2, Treatment('HC')) + {COV} + BDI_T1 + BAI_T1", df
        ).fit()  # M2：再控制憂鬱與焦慮
        for g in ["Higher", "Lower"]:
            a = df.loc[df["k2"] == g, col]
            gg, lo, hi, n1, n2 = hedges_g(a, hc[col])
            k = [c for c in m.params.index if g in c][0]
            rows.append(
                dict(
                    contrast=f"{g} vs HC",
                    n=n1,
                    n_hc=n2,
                    g=gg,
                    ci_low=lo,
                    ci_high=hi,
                    p_mwu=mwu_p(a, hc[col]),
                    adj_diff_in_HC_SD=m.params[k] / sd,
                    p_adjusted=m.pvalues[k],
                    m2_diff_in_HC_SD=m2.params[k] / sd,
                    p_m2=m2.pvalues[k],
                    retained_pct=100 * m2.params[k] / m.params[k],
                )
            )
        return rows

    # A
    A = []
    for kind, klab in [("adj", "aperiodic-adjusted"), ("raw", "raw log power")]:
        for lab, lo, hi in BANDS:
            for r in contrast(f"{kind}_{lo}_{hi}"):
                A.append(dict(band=lab, measure=klab, **r))
    A = pd.DataFrame(A)
    A.to_csv(OUT / "01_subband_contrasts.csv", index=False)
    # B
    Bt = []
    for col, klab in [
        ("adj_hi_minus_lo", "aperiodic-adjusted"),
        ("raw_hi_minus_lo", "raw log power"),
    ]:
        for r in contrast(col):
            Bt.append(dict(measure=f"20–30 minus 13–20 Hz, {klab}", **r))
    gg, lo, hi, n1, n2 = hedges_g(
        df.loc[df["k2"] == "Higher", "adj_hi_minus_lo"],
        df.loc[df["k2"] == "Lower", "adj_hi_minus_lo"],
    )
    Bt.append(
        dict(
            measure="20–30 minus 13–20 Hz, aperiodic-adjusted",
            contrast="Higher vs Lower",
            n=n1,
            n_hc=n2,
            g=gg,
            ci_low=lo,
            ci_high=hi,
            p_mwu=mwu_p(
                df.loc[df["k2"] == "Higher", "adj_hi_minus_lo"],
                df.loc[df["k2"] == "Lower", "adj_hi_minus_lo"],
            ),
            adj_diff_in_HC_SD=np.nan,
            p_adjusted=np.nan,
            m2_diff_in_HC_SD=np.nan,
            p_m2=np.nan,
            retained_pct=np.nan,
        )
    )
    Bt = pd.DataFrame(Bt)
    Bt.to_csv(OUT / "02_high_minus_low.csv", index=False)
    # C
    Ct = []
    for lab, lo, hi in BANDS:
        col = f"adj_{lo}_{hi}"
        m = smf.ols(f"{col} ~ thr + {COV}", df).fit()
        ci = m.conf_int().loc["thr"]
        sd = hc[col].std(ddof=1)
        Ct.append(
            dict(
                band=lab,
                step_in_HC_SD=m.params["thr"] / sd,
                ci_low=ci[0] / sd,
                ci_high=ci[1] / sd,
                p=m.pvalues["thr"],
                n_above=int(df["thr"].sum()),
            )
        )
    col = "adj_hi_minus_lo"
    m = smf.ols(f"{col} ~ thr + {COV}", df).fit()
    ci = m.conf_int().loc["thr"]
    sd = hc[col].std(ddof=1)
    Ct.append(
        dict(
            band="20–30 minus 13–20 Hz",
            step_in_HC_SD=m.params["thr"] / sd,
            ci_low=ci[0] / sd,
            ci_high=ci[1] / sd,
            p=m.pvalues["thr"],
            n_above=int(df["thr"].sum()),
        )
    )
    Ct = pd.DataFrame(Ct)
    Ct.to_csv(OUT / "03_threshold.csv", index=False)
    # D
    ch = pd.read_csv(AN / "S10_fooof_full" / "01_specparam_features.csv")
    ch = ch[ch["subject"].isin(df["subject"])].merge(df[["subject", "k2"]], on="subject")
    ch["loc"] = np.where(
        ch["beta_cf"].isna(), "no peak", np.where(ch["beta_cf"] < 20, "13–20 Hz", "20–30 Hz")
    )
    per = (
        ch.groupby(["subject", "k2"])
        .apply(
            lambda g: pd.Series(
                {
                    "pct_low": 100 * (g["loc"] == "13–20 Hz").mean(),
                    "pct_high": 100 * (g["loc"] == "20–30 Hz").mean(),
                    "pw_low": g.loc[g["loc"] == "13–20 Hz", "beta_pw"].mean(),
                    "pw_high": g.loc[g["loc"] == "20–30 Hz", "beta_pw"].mean(),
                    "cf": g["beta_cf"].mean(),
                }
            )
        )
        .reset_index()
    )
    Dt = []
    ph = per[per["k2"] == "HC"]
    for v, lab in [
        ("pct_low", "Channels with dominant beta peak at 13–20 Hz, %"),
        ("pct_high", "Channels with dominant beta peak at 20–30 Hz, %"),
        ("pw_low", "Peak power when peak at 13–20 Hz"),
        ("pw_high", "Peak power when peak at 20–30 Hz"),
        ("cf", "Beta peak frequency, Hz"),
    ]:
        row = dict(measure=lab, HC=ph[v].mean())
        for g in ["Lower", "Higher"]:
            a = per.loc[per["k2"] == g, v]
            gg, lo, hi, n1, n2 = hedges_g(a.dropna(), ph[v].dropna())
            row.update({g: a.mean(), f"g_{g}": gg, f"p_{g}": mwu_p(a, ph[v]), f"n_{g}": n1})
        Dt.append(row)
    Dt = pd.DataFrame(Dt)
    Dt.to_csv(OUT / "04_peak_location.csv", index=False)

    L = [
        "# S_R1n — beta 子頻段",
        "",
        f"N = {len(df)}（HC {len(hc)}、低嚴重度 {(df.k2 == 'Lower').sum()}、高嚴重度 {(df.k2 == 'Higher').sum()}）。",
        "",
        "## A 各子頻段對健康組的效果量（去背景功率）",
        "",
        "M1 校正年齡、性別、教育；M2 再加 BDI-II、BAI；保留＝M2 係數／M1 係數。",
        "",
        "| 頻段 | 高嚴重度 g [95% CI] | M1 p | M2 p | 保留 | 低嚴重度 g [95% CI] | M1 p |",
        "|---|---|---|---|---|---|---|",
    ]
    for lab, lo, hi in BANDS:
        x = A[(A.band == lab) & (A.measure == "aperiodic-adjusted")].set_index("contrast")
        h, l = x.loc["Higher vs HC"], x.loc["Lower vs HC"]
        L.append(
            f"| {lab} | {h.g:+.2f} [{h.ci_low:+.2f}, {h.ci_high:+.2f}] | {fmt_p(h.p_adjusted)} | {fmt_p(h.p_m2)} | {h.retained_pct:.0f}% | {l.g:+.2f} [{l.ci_low:+.2f}, {l.ci_high:+.2f}] | {fmt_p(l.p_adjusted)} |"
        )
    L += [
        "",
        "原始 log 功率",
        "",
        "| 頻段 | 高嚴重度 g [95% CI] | M1 p | M2 p | 保留 | 低嚴重度 g | M1 p |",
        "|---|---|---|---|---|---|---|",
    ]
    for lab, lo, hi in BANDS:
        x = A[(A.band == lab) & (A.measure == "raw log power")].set_index("contrast")
        h, l = x.loc["Higher vs HC"], x.loc["Lower vs HC"]
        L.append(
            f"| {lab} | {h.g:+.2f} [{h.ci_low:+.2f}, {h.ci_high:+.2f}] | {fmt_p(h.p_adjusted)} | {fmt_p(h.p_m2)} | {h.retained_pct:.0f}% | {l.g:+.2f} | {fmt_p(l.p_adjusted)} |"
        )
    L += [
        "",
        "## B 高段減低段（20–30 減 13–20 Hz）",
        "",
        "| 指標 | 比較 | g [95% CI] | p（Mann–Whitney） | p（校正後） |",
        "|---|---|---|---|---|",
    ]
    for r in Bt.itertuples():
        L.append(
            f"| {r.measure} | {r.contrast} | {r.g:+.2f} [{r.ci_low:+.2f}, {r.ci_high:+.2f}] | {fmt_p(r.p_mwu)} | {fmt_p(r.p_adjusted) if not np.isnan(r.p_adjusted) else '–'} |"
        )
    L += [
        "",
        f"## C 門檻（ISI ≥ {THR}）在各子頻段的係數（去背景功率，健康組標準差單位）",
        "",
        "| 頻段 | 係數 [95% CI] | p |",
        "|---|---|---|",
    ]
    for r in Ct.itertuples():
        L.append(
            f"| {r.band} | {r.step_in_HC_SD:+.2f} [{r.ci_low:+.2f}, {r.ci_high:+.2f}] | {fmt_p(r.p)} |"
        )
    L += [
        "",
        "## D 偵測到的 beta 峰落在哪一段",
        "",
        "| 指標 | 健康組 | 低嚴重度（g, p） | 高嚴重度（g, p） |",
        "|---|---|---|---|",
    ]
    for r in Dt.itertuples():
        L.append(
            f"| {r.measure} | {r.HC:.2f} | {r.Lower:.2f}（{r.g_Lower:+.2f}, {fmt_p(r.p_Lower)}） | {r.Higher:.2f}（{r.g_Higher:+.2f}, {fmt_p(r.p_Higher)}） |"
        )
    (OUT / "05_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
