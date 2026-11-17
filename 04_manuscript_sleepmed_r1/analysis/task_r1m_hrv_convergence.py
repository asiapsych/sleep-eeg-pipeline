"""探索 periodic beta 與 HRV 的關聯及調整後的穩健性。

HRV 讀 task_r1q 通過品管的重算值；資料庫的軟體數值只作對照。
標準化係數同報一般與 HC3 標準誤，避免殘差變異不均時過度樂觀。"""

from __future__ import annotations
import sqlite3

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

from _dynamic_n import DB_PATH
from _r1_common import AN, COV, COV_COLS, OUTCOME, fmt_p, hedges_g, load_all
from _r1_spectra import band_mean, load_spectra

OUT = AN / "S_R1m_hrv"
OUT.mkdir(parents=True, exist_ok=True)
FIG = AN / "S_R1_figures"
FIG.mkdir(parents=True, exist_ok=True)
B = 5000
rng = np.random.default_rng(20261116)
HRV = [
    ("lnRMSSD", "ln RMSSD"),
    ("lnHF", "ln HF power"),
    ("lnSDNN", "ln SDNN"),
    ("HR", "Heart rate"),
    ("lnLFHF", "ln LF/HF"),
]
EEG = [
    (OUTCOME, "Periodic beta peak power"),
    ("adj_13_30", "Aperiodic-adjusted beta power, 13–30 Hz"),
    ("adj_13_20", "Aperiodic-adjusted power, 13–20 Hz"),
    ("adj_20_30", "Aperiodic-adjusted power, 20–30 Hz"),
    ("Beta_rel_BRAIN", "Relative beta power"),
    ("logBeta_abs", "Absolute beta power (log)"),
    ("exponent_BRAIN", "Aperiodic exponent"),
    ("offset_BRAIN", "Aperiodic offset"),
    ("alpha_pw_BRAIN", "Periodic alpha peak power"),
    ("Alpha_rel_BRAIN", "Relative alpha power"),
    ("Theta_rel_BRAIN", "Relative theta power"),
    ("adj_30_45", "Aperiodic-adjusted power, 30–45 Hz"),
    ("raw_30_45", "Raw log power, 30–45 Hz"),
]
COL = {"HC": "#898781", "Lower": "#2a78d6", "Higher": "#eb6834"}


def z(s):
    return (s - s.mean()) / s.std(ddof=1)


def std_beta(d, y, x, extra="", boot=0):
    """標準化係數。殘差變異不均（Breusch–Pagan 顯著）時一般標準誤會偏樂觀，所以另報 HC3 穩健標準誤的 p；
    boot > 0 時再給 bootstrap 百分位信賴區間（對受試者重抽）。"""
    dd = d.dropna(subset=[y, x]).copy()
    dd["_y"], dd["_x"] = z(dd[y]), z(dd[x])
    f = f"_y ~ _x + {COV}{extra}"
    m = smf.ols(f, dd).fit()
    ci = m.conf_int().loc["_x"]
    r = m.get_robustcov_results("HC3")
    i = list(m.params.index).index("_x")
    out = dict(
        n=int(m.nobs),
        beta=m.params["_x"],
        ci_low=ci[0],
        ci_high=ci[1],
        p=m.pvalues["_x"],
        p_hc3=float(r.pvalues[i]),
    )
    if boot:
        g = np.random.default_rng(20261117)
        bs = []
        for _ in range(boot):
            try:
                bs.append(smf.ols(f, dd.iloc[g.integers(0, len(dd), len(dd))]).fit().params["_x"])
            except Exception:
                pass
        out.update(boot_low=float(np.percentile(bs, 2.5)), boot_high=float(np.percentile(bs, 97.5)))
    return out


def partial_r(d, y, x, covs):
    dd = d.dropna(subset=[y, x] + covs)
    C = np.column_stack([np.ones(len(dd))] + [dd[c].to_numpy(float) for c in covs])
    ry = dd[y].to_numpy(float) - C @ np.linalg.lstsq(C, dd[y].to_numpy(float), rcond=None)[0]
    rx = dd[x].to_numpy(float) - C @ np.linalg.lstsq(C, dd[x].to_numpy(float), rcond=None)[0]
    return np.corrcoef(rx, ry)[0, 1], rx, ry, dd


def main():
    df = load_all()
    with sqlite3.connect(DB_PATH) as con:
        h = pd.read_sql(
            "SELECT NEW_NUMBER AS subject, HRV_RMSSD_MS, HRV_HF, HRV_SDNN_MS, HRV_EKG_HR, HRV_LF_HF, HRV_RESP_RATE, PSQI_T1_C6 FROM isi_raw_data WHERE NEW_NUMBER IS NOT NULL",
            con,
        )
    for c in h.columns[1:]:
        h[c] = pd.to_numeric(h[c], errors="coerce")
    pos = lambda s: np.log(s.where(s > 0))
    # HRV 的主要來源：task_r1q_hrv_recompute.py 由原始心電訊號計算、且通過訊號品管的數值。
    # 擷取軟體算出、登錄在總表的數值（HRV_*）只當對照（sw_lnRMSSD）。
    q = pd.read_csv(AN / "S_R1q_hrv_recompute" / "_subject_level_internal.csv")
    q = q[q["ok"].astype(bool)]
    h = h.merge(q[["subject", "RMSSD", "HF", "SDNN", "LFHF", "HR"]], on="subject", how="left")
    h["lnRMSSD"], h["lnHF"], h["lnSDNN"], h["lnLFHF"] = (
        pos(h.RMSSD),
        pos(h.HF),
        pos(h.SDNN),
        pos(h.LFHF),
    )
    h["sw_lnRMSSD"] = pos(h.HRV_RMSSD_MS)
    h["RESP"] = h.HRV_RESP_RATE.where(
        h.HRV_RESP_RATE.between(4, 40)
    )  # 資料庫有 2 筆呼吸速率明顯是壞值（342、1.6e41），設為缺值
    h["sleepmed"] = (h.PSQI_T1_C6 > 0).astype(float).where(h.PSQI_T1_C6.notna())
    d = df.merge(h, on="subject", how="left").dropna(subset=COV_COLS).reset_index(drop=True)
    freq, chn, adj, raw = load_spectra(d["subject"])
    d["adj_13_30"], d["adj_13_20"], d["adj_20_30"] = (
        band_mean(adj, freq, 13, 30),
        band_mean(adj, freq, 13, 20),
        band_mean(adj, freq, 20, 30),
    )
    d["adj_30_45"], d["raw_30_45"] = band_mean(adj, freq, 30, 45), band_mean(raw, freq, 30, 45)
    d["logBeta_abs"] = np.log10(d["Beta_abs_BRAIN"])
    d["thr"] = (d["ISI_T1"] >= 18).astype(int)
    dh = d.dropna(subset=["lnRMSSD"]).reset_index(drop=True)
    hc = dh[dh.k2 == "HC"]
    n_by = dh.k2.value_counts().to_dict()

    # A
    A = []
    for v, lab in HRV + [("RESP", "Respiration rate")]:
        m = smf.ols(f"{v} ~ C(k2, Treatment('HC')) + {COV}", dh).fit()
        row = dict(measure=lab, HC_mean=hc[v].mean(), HC_sd=hc[v].std(ddof=1))
        for g in ["Lower", "Higher"]:
            a = dh.loc[dh.k2 == g, v]
            gg, lo, hi, n1, _ = hedges_g(a.dropna(), hc[v].dropna())
            k = [c for c in m.params.index if g in c][0]
            row.update(
                {
                    f"{g}_mean": a.mean(),
                    f"{g}_g": gg,
                    f"{g}_ci_low": lo,
                    f"{g}_ci_high": hi,
                    f"{g}_p_adjusted": m.pvalues[k],
                }
            )
        ms = smf.ols(f"{v} ~ thr + {COV}", dh).fit()
        row.update(step18_p=ms.pvalues["thr"])
        A.append(row)
    A = pd.DataFrame(A)
    A.to_csv(OUT / "01_group.csv", index=False)

    # B
    Bt = []
    for hv, hlab in HRV:
        for ev, elab in EEG:
            r = std_beta(dh, ev, hv)
            s = dh[[ev, hv]].dropna()
            rho, prho = stats.spearmanr(s[ev], s[hv])
            Bt.append(dict(hrv=hlab, eeg=elab, **r, spearman_rho=rho, spearman_p=prho))
    Bt = pd.DataFrame(Bt)
    Bt["q_within_hrv"] = np.nan
    for hlab in Bt.hrv.unique():
        i = Bt.hrv == hlab
        Bt.loc[i, "q_within_hrv"] = multipletests(Bt.loc[i, "p"], method="fdr_bh")[1]
    Bt.to_csv(OUT / "02_associations.csv", index=False)

    # C
    y, x = OUTCOME, "lnRMSSD"
    Ct = []

    def add(label, dd, extra="", ylab=y, xlab=x):
        Ct.append(dict(model=label, **std_beta(dd, ylab, xlab, extra, boot=2000)))

    add("Age, sex, education", dh)
    add("+ ISI total", dh, " + ISI_T1")
    add("+ ISI total, BDI-II, BAI", dh, " + ISI_T1 + BDI_T1 + BAI_T1")
    add("+ severity cluster (HC / lower / higher)", dh, " + C(k2)")
    add("+ sleep-medication use (PSQI component 6 > 0)", dh, " + sleepmed")
    add("+ respiration rate", dh, " + RESP")
    add("+ heart rate", dh, " + HR")
    add("+ 30–45 Hz power (raw log)", dh, " + raw_30_45")
    add("+ all of the above", dh, " + ISI_T1 + BDI_T1 + BAI_T1 + sleepmed + RESP + raw_30_45")
    rk = dh.dropna(subset=[y, x]).copy()
    rk["ry"], rk["rx"] = rk[y].rank(), rk[x].rank()
    add("Rank-transformed", rk, ylab="ry", xlab="rx")
    zz = dh.dropna(subset=[y, x])
    keep = (np.abs(z(zz[y])) <= 3) & (np.abs(z(zz[x])) <= 3)
    add(f"Excluding |z| > 3 (n removed = {int((~keep).sum())})", zz[keep])
    add("No sleep-medication use only", dh[dh.sleepmed == 0])
    for g in ["HC", "Lower", "Higher"]:
        add(f"Within {g}", dh[dh.k2 == g])
    add("Acquisition-software values, recordings passing signal QC", dh, xlab="sw_lnRMSSD")
    add("Acquisition-software values, all with a recorded value", d, xlab="sw_lnRMSSD")
    Ct = pd.DataFrame(Ct)
    r0, rx, ry, dd = partial_r(dh, y, x, COV_COLS)
    bs = []
    for _ in range(B):
        i = rng.integers(0, len(dd), len(dd))
        bs.append(partial_r(dd.iloc[i], y, x, COV_COLS)[0])
    pr = dict(
        partial_r=r0,
        ci_low=float(np.percentile(bs, 2.5)),
        ci_high=float(np.percentile(bs, 97.5)),
        n=len(dd),
    )
    Ct.to_csv(OUT / "03_robustness.csv", index=False)

    # D
    Dt = []
    for tag, extra in [("without HRV", ""), ("with ln RMSSD", " + lnRMSSD")]:
        m = smf.ols(f"{y} ~ C(k2, Treatment('HC')) + {COV}{extra}", dh.dropna(subset=[y])).fit()
        sd = hc[y].std(ddof=1)
        for g in ["Higher", "Lower"]:
            k = [c for c in m.params.index if g in c][0]
            Dt.append(
                dict(
                    model=tag,
                    group=g,
                    diff_in_HC_SD=m.params[k] / sd,
                    p=m.pvalues[k],
                    n=int(m.nobs),
                )
            )
        if extra:
            Dt.append(
                dict(
                    model=tag,
                    group="ln RMSSD (per SD)",
                    diff_in_HC_SD=m.params["lnRMSSD"] * dh.lnRMSSD.std(ddof=1) / sd,
                    p=m.pvalues["lnRMSSD"],
                    n=int(m.nobs),
                )
            )
    Dt = pd.DataFrame(Dt)
    Dt.to_csv(OUT / "04_group_with_hrv.csv", index=False)

    # 圖
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw=dict(width_ratios=[1, 1.15]))
    for g in ["HC", "Lower", "Higher"]:
        i = (dd.k2 == g).to_numpy()
        ax[0].scatter(
            rx[i],
            ry[i],
            s=16,
            color=COL[g],
            alpha=0.75,
            label={"HC": "Controls", "Lower": "Lower severity", "Higher": "Higher severity"}[g],
            linewidths=0,
        )
    b1, b0 = np.polyfit(rx, ry, 1)
    xs = np.linspace(rx.min(), rx.max(), 50)
    ax[0].plot(xs, b0 + b1 * xs, color="#222", lw=1.3)
    ax[0].set_xlabel("ln RMSSD (adjusted for age, sex, education)")
    ax[0].set_ylabel("Periodic beta peak power (adjusted)")
    ax[0].set_title(
        f"A  Partial r = {r0:.2f} [{pr['ci_low']:.2f}, {pr['ci_high']:.2f}], n = {pr['n']}",
        loc="left",
        fontsize=10,
    )
    ax[0].legend(frameon=False, fontsize=8)
    sub = Bt[Bt.hrv == "ln RMSSD"].iloc[::-1].reset_index(drop=True)
    for i, r in sub.iterrows():
        c = "#eb6834" if r.q_within_hrv < 0.05 else "#898781"
        ax[1].plot([r.ci_low, r.ci_high], [i, i], color=c, lw=1.6)
        ax[1].plot(r.beta, i, "o", color=c, ms=5)
    ax[1].axvline(0, color="#222", lw=0.8)
    ax[1].set_yticks(range(len(sub)))
    ax[1].set_yticklabels(sub.eeg, fontsize=8)
    ax[1].set_xlabel("Standardised coefficient for ln RMSSD [95% CI]")
    ax[1].set_title("B  Association of ln RMSSD with each EEG measure", loc="left", fontsize=10)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig9_hrv.png", dpi=200)
    plt.close(fig)

    L = [
        "# S_R1m — 週期性 beta 與心率變異",
        "",
        f"HRV 由原始心電訊號計算並通過訊號品管的人：{len(dh)}（HC {n_by.get('HC', 0)}、低嚴重度 {n_by.get('Lower', 0)}、高嚴重度 {n_by.get('Higher', 0)}）；全樣本 {len(df)}。",
        "",
        "## A HRV 的組間差",
        "",
        "| 指標 | 健康組 M (SD) | 低嚴重度 g [95% CI]（p） | 高嚴重度 g [95% CI]（p） | ISI ≥ 18 階梯 p |",
        "|---|---|---|---|---|",
    ]
    for r in A.itertuples():
        L.append(
            f"| {r.measure} | {r.HC_mean:.2f} ({r.HC_sd:.2f}) | {r.Lower_g:+.2f} [{r.Lower_ci_low:+.2f}, {r.Lower_ci_high:+.2f}]（{fmt_p(r.Lower_p_adjusted)}） | {r.Higher_g:+.2f} [{r.Higher_ci_low:+.2f}, {r.Higher_ci_high:+.2f}]（{fmt_p(r.Higher_p_adjusted)}） | {fmt_p(r.step18_p)} |"
        )
    L += ["", "## B 各 HRV 指標對各腦波指標（標準化係數；校正年齡、性別、教育）", ""]
    for hlab in Bt.hrv.unique():
        L += [
            f"**{hlab}**",
            "",
            "| 腦波指標 | n | 係數 [95% CI] | p | q | Spearman rho |",
            "|---|---|---|---|---|---|",
        ]
        for r in Bt[Bt.hrv == hlab].itertuples():
            L.append(
                f"| {r.eeg} | {r.n} | {r.beta:+.2f} [{r.ci_low:+.2f}, {r.ci_high:+.2f}] | {fmt_p(r.p)} | {fmt_p(r.q_within_hrv)} | {r.spearman_rho:+.2f} |"
            )
        L.append("")
    L += [
        f"## C 週期性 beta × ln RMSSD 的穩健性（偏相關 r = {r0:.2f}，bootstrap 95% CI [{pr['ci_low']:.2f}, {pr['ci_high']:.2f}]，n = {pr['n']}）",
        "",
        "| 模型 | n | ln RMSSD 的標準化係數 [bootstrap 95% CI] | p（一般標準誤） | p（HC3 穩健標準誤） |",
        "|---|---|---|---|---|",
    ]
    for r in Ct.itertuples():
        L.append(
            f"| {r.model} | {r.n} | {r.beta:+.2f} [{r.boot_low:+.2f}, {r.boot_high:+.2f}] | {fmt_p(r.p)} | {fmt_p(r.p_hc3)} |"
        )
    L += ["", "殘差變異不均，一般標準誤的 p 偏樂觀；解讀以 bootstrap 信賴區間與 HC3 的 p 為準。"]
    L += [
        "",
        "## D 加入 HRV 前後的組間差（週期性 beta，健康組標準差單位）",
        "",
        "| 模型 | 項目 | 係數 | p |",
        "|---|---|---|---|",
    ]
    for r in Dt.itertuples():
        L.append(f"| {r.model} | {r.group} | {r.diff_in_HC_SD:+.2f} | {fmt_p(r.p)} |")
    (OUT / "05_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
