"""檢查 beta 升高的頻率、嚴重度與頭皮分布。

去背景頻譜用既有的 offset／exponent，不重新擬合。
H1–H3 先調整年齡、性別、教育；H4 的逐電極效果量用未調整值。
逐頻率比較以 cluster permutation 控制多重比較。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, Patch
from scipy import stats
from scipy.io import loadmat
from statsmodels.stats.multitest import multipletests

from _r1_common import AN, COV_COLS, hedges_g, load_all, mwu_p

THR = 18
B = 5000
OUT = AN / "S_R1h_beta_characterisation"
OUT.mkdir(parents=True, exist_ok=True)
FIG = AN / "S_R1_figures"
FIG.mkdir(parents=True, exist_ok=True)
PSD_DIR = AN / "S10_fooof_psd_input"
CH_RENAME = {"T7": "T3", "T8": "T4", "P7": "T5", "P8": "T6"}
CORE = ["F3", "F4", "FZ", "C3", "C4", "CZ", "P3", "P4", "PZ"]
BANDS = {
    "theta 4–8 Hz": (4, 8),
    "alpha 8–13 Hz": (8, 13),
    "low beta 13–20 Hz": (13, 20),
    "high beta 20–30 Hz": (20, 30),
    "beta 13–30 Hz": (13, 30),
    "30–45 Hz": (30, 45),
}
MD = []
rng = np.random.default_rng(20261116)


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    MD.append(s)


df = load_all().dropna(subset=COV_COLS).reset_index(drop=True)
chp = pd.read_csv(AN / "S10_fooof_full" / "01_specparam_features.csv").set_index(
    ["subject", "channel"]
)

# 去背景頻譜
flat_all, flat_core, raw_all, flat_ch, CHN, FREQ = [], [], [], [], None, None
for sid in df["subject"]:
    m = loadmat(str(PSD_DIR / f"psd_{sid}.mat"), squeeze_me=True)
    f = np.asarray(m["freqs"], float)
    sel = (f >= 1.0) & (f <= 45.0)
    names = [CH_RENAME.get(str(c).upper(), str(c).upper()) for c in np.atleast_1d(m["ch_names"])]
    lp = np.log10(np.asarray(m["psd"], float)[:, sel])
    ff = f[sel]
    ap = np.array([[chp.loc[(sid, c), "offset"], chp.loc[(sid, c), "exponent"]] for c in names])
    fl = lp - (ap[:, [0]] - ap[:, [1]] * np.log10(ff)[None, :])
    order = np.argsort(names)
    if CHN is None:
        CHN, FREQ = [names[i] for i in order], ff
    fl = fl[order]
    lp = lp[order]
    flat_all.append(fl.mean(0))
    raw_all.append(lp.mean(0))
    flat_ch.append(fl)
    flat_core.append(fl[[CHN.index(c) for c in CORE]].mean(0))
flat_all, flat_core, raw_all, flat_ch = map(np.asarray, (flat_all, flat_core, raw_all, flat_ch))

C = np.column_stack(
    [
        np.ones(len(df)),
        df["AGE"].to_numpy(float),
        df["SEX"].to_numpy(float),
        df["EDU"].to_numpy(float),
    ]
)


def resid(M):
    M2 = M.reshape(len(M), -1)
    out = M2 - C @ np.linalg.lstsq(C, M2, rcond=None)[0] + M2.mean(0)
    return out.reshape(M.shape)


R_all, R_core, R_raw, R_ch = resid(flat_all), resid(flat_core), resid(raw_all), resid(flat_ch)
isHC = (df["group"] == "HC").to_numpy()
isIns = ~isHC
MASK = {
    "Higher vs HC": (df["k2"] == "Higher").to_numpy(),
    "Lower vs HC": (df["k2"] == "Lower").to_numpy(),
    f"ISI >={THR} vs HC": isIns & (df["ISI_T1"] >= THR).to_numpy(),
    f"ISI 9-{THR - 1} vs HC": isIns & (df["ISI_T1"] < THR).to_numpy(),
}
P("# S_R1h — beta 升高的樣貌\n")
P(
    f"N = {len(df)}（HC {isHC.sum()}、低嚴重度 {MASK['Lower vs HC'].sum()}、高嚴重度 {MASK['Higher vs HC'].sum()}）；頻率 1–45 Hz，解析度 {FREQ[1] - FREQ[0]:.2f} Hz，共 {len(FREQ)} 點\n"
)


def tstat(A, Bm):
    return (A.mean(0) - Bm.mean(0)) / np.sqrt(
        A.var(0, ddof=1) / len(A) + Bm.var(0, ddof=1) / len(Bm)
    )


def gvec(A, Bm):
    n1, n2 = len(A), len(Bm)
    sp = np.sqrt(((n1 - 1) * A.var(0, ddof=1) + (n2 - 1) * Bm.var(0, ddof=1)) / (n1 + n2 - 2))
    return (A.mean(0) - Bm.mean(0)) / sp * (1 - 3 / (4 * (n1 + n2) - 9))


def clusters(t, thr):
    out, start = [], None
    sgn = np.where(t > thr, 1, np.where(t < -thr, -1, 0))
    for i in range(len(t) + 1):
        s = sgn[i] if i < len(t) else 0
        if start is not None and s != sgn[start]:
            out.append((start, i - 1, float(np.abs(t[start:i]).sum()), int(sgn[start])))
            start = None
        if start is None and s != 0:
            start = i
    return out


def cluster_test(M, mask, label, measure):
    X = np.vstack([M[mask], M[isHC]])
    n1 = int(mask.sum())
    thr = stats.t.ppf(0.975, len(X) - 2)
    t = tstat(X[:n1], X[n1:])
    obs = clusters(t, thr)
    null = np.empty(B)
    for b in range(B):
        p = rng.permutation(len(X))
        cl = clusters(tstat(X[p[:n1]], X[p[n1:]]), thr)
        null[b] = max([c[2] for c in cl], default=0.0)
    g = gvec(X[:n1], X[n1:])
    rows = []
    for s, e, mass, sg in obs:
        j = s + int(np.argmax(np.abs(g[s : e + 1])))
        rows.append(
            dict(
                measure=measure,
                contrast=label,
                direction="higher than HC" if sg > 0 else "lower than HC",
                from_Hz=round(FREQ[s], 2),
                to_Hz=round(FREQ[e], 2),
                n_bins=e - s + 1,
                cluster_mass=round(mass, 1),
                p_cluster=round((np.sum(null >= mass) + 1) / (B + 1), 4),
                peak_g=round(float(g[j]), 3),
                peak_Hz=round(FREQ[j], 2),
            )
        )
    return g, rows


# H1
G, rows = {}, []
for meas, M in [
    ("aperiodic-adjusted, whole-scalp", R_all),
    ("aperiodic-adjusted, core 9 electrodes", R_core),
    ("raw log power, whole-scalp", R_raw),
]:
    for lab, mk in MASK.items():
        g, r = cluster_test(M, mk, lab, meas)
        G[(meas, lab)] = g
        rows += r
cl = pd.DataFrame(rows).sort_values(["measure", "contrast", "from_Hz"])
cl.to_csv(OUT / "01_frequency_clusters.csv", index=False, encoding="utf-8-sig")
pd.DataFrame({"freq_Hz": FREQ, **{f"g | {m} | {l}": v for (m, l), v in G.items()}}).round(4).to_csv(
    OUT / "02_g_by_frequency.csv", index=False, encoding="utf-8-sig"
)
P("## H1. 逐頻率比較（cluster-based permutation，B = %d；只列通過逐點門檻的連續區段）\n" % B)
P(cl[cl["n_bins"] >= 3].to_markdown(index=False))
P("")

# H2
rows = []
for bn, (lo, hi) in BANDS.items():
    sel = (FREQ >= lo) & (FREQ <= hi)
    for meas, M in [("aperiodic-adjusted", R_all), ("raw log power", R_raw)]:
        v = M[:, sel].mean(1)
        for lab in ["Higher vs HC", "Lower vs HC"]:
            g, l, h, n1, n2 = hedges_g(v[MASK[lab]], v[isHC])
            rows.append(
                dict(
                    band=bn,
                    measure=meas,
                    contrast=lab,
                    n_a=n1,
                    n_HC=n2,
                    g=round(g, 3),
                    ci_low=round(l, 3),
                    ci_high=round(h, 3),
                    p=round(mwu_p(v[MASK[lab]], v[isHC]), 4),
                )
            )
bd = pd.DataFrame(rows)
bd.to_csv(OUT / "03_band_averages.csv", index=False, encoding="utf-8-sig")
P("## H2. 頻帶平均功率（全頭皮；aperiodic-adjusted＝扣掉背景後，不需偵測峰值）\n")
P(bd.to_markdown(index=False))
P("")

# H3
BINS = [(9, 10), (11, 12), (13, 14), (15, 16), (17, 18), (19, 28)]
hz = np.arange(1, 45)
M1 = np.column_stack([R_all[:, (FREQ >= a) & (FREQ < a + 1)].mean(1) for a in hz])
isi = df["ISI_T1"].to_numpy()
HM = np.vstack([gvec(M1[isIns & (isi >= lo) & (isi <= hi)], M1[isHC]) for lo, hi in BINS])
pd.DataFrame(
    HM,
    index=[f"{lo}-{hi}" if hi < 28 else f">={lo}" for lo, hi in BINS],
    columns=[f"{a}-{a + 1} Hz" for a in hz],
).round(3).to_csv(OUT / "04_frequency_by_severity.csv", encoding="utf-8-sig")

# H4
selb = (FREQ >= 13) & (FREQ <= 30)
chb = R_ch[:, :, selb].mean(2)  # 人 × 電極：13–30 Hz 去背景功率
pkw = chp["beta_pw"].unstack("channel").reindex(df["subject"])[CHN].to_numpy()
rows = []
for lab in ["Higher vs HC", "Lower vs HC"]:
    tmp = []
    for j, c in enumerate(CHN):
        g1, l1, h1, n1, _ = hedges_g(chb[MASK[lab], j], chb[isHC, j])
        g2, l2, h2, n2, _ = hedges_g(pkw[MASK[lab], j], pkw[isHC, j])
        tmp.append(
            dict(
                contrast=lab,
                channel=c,
                g_band=round(g1, 3),
                p_band=mwu_p(chb[MASK[lab], j], chb[isHC, j]),
                g_peak=round(g2, 3),
                n_with_peak=n2,
                p_peak=mwu_p(pkw[MASK[lab], j], pkw[isHC, j]),
            )
        )
    t = pd.DataFrame(tmp)
    t["q_band"] = multipletests(t["p_band"], method="fdr_bh")[1]
    t["q_peak"] = multipletests(t["p_peak"], method="fdr_bh")[1]
    rows.append(t)
tp = pd.concat(rows).round(4)
tp.to_csv(OUT / "05_topography.csv", index=False, encoding="utf-8-sig")
P(
    "## H4. 各電極的效果量（g_band＝13–30 Hz 去背景功率；g_peak＝偵測到的 beta 峰值功率；q 為 19 個電極內的 FDR）\n"
)
P(tp.to_markdown(index=False))
P("")
for lab in ["Higher vs HC", "Lower vs HC"]:
    t = tp[tp["contrast"] == lab]
    P(
        f"- {lab}：g_band 範圍 {t['g_band'].min():.2f} 到 {t['g_band'].max():.2f}，{int((t['q_band'] < 0.05).sum())}/19 個電極 q < .05；"
        f"最大在 {t.loc[t['g_band'].idxmax(), 'channel']}，最小在 {t.loc[t['g_band'].idxmin(), 'channel']}"
    )
P("\n數字全部由本腳本即時計算；尚未寫回 SSOT。")
(OUT / "06_summary.md").write_text("\n".join(MD), encoding="utf-8")

# 圖
SURFACE, INK, INK2, MUTED, GRID, AXIS = (
    "#fcfcfb",
    "#0b0b0b",
    "#52514e",
    "#898781",
    "#e1e0d9",
    "#c3c2b7",
)
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#898781"
plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "Helvetica",
        "font.size": 10,
        "text.color": INK,
        "axes.labelcolor": INK2,
        "axes.edgecolor": AXIS,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "axes.axisbelow": True,
        "legend.frameon": False,
        "legend.fontsize": 9,
    }
)
DIV = LinearSegmentedColormap.from_list(
    "div", ["#1c5cab", "#86b6ef", "#f0efec", "#f3a583", "#c4481a"]
)


def plabel(ax, s):
    ax.text(
        -0.02,
        1.07,
        s,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="right",
        va="bottom",
    )


fig = plt.figure(figsize=(13.6, 8.4))
gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.34)
a = fig.add_subplot(gs[0, 0])
for lab, mk, col in [
    ("Healthy controls", isHC, GRAY),
    ("Lower-severity cluster", MASK["Lower vs HC"], BLUE),
    ("Higher-severity cluster", MASK["Higher vs HC"], ORANGE),
]:
    mu, se = R_all[mk].mean(0), R_all[mk].std(0, ddof=1) / np.sqrt(mk.sum())
    a.fill_between(FREQ, mu - se, mu + se, color=col, alpha=0.18, lw=0)
    a.plot(FREQ, mu, color=col, lw=2, label=f"{lab} (n = {mk.sum()})")
a.axvspan(13, 30, color=GRID, alpha=0.45, lw=0)
a.text(
    21.5,
    a.get_ylim()[0] + 0.04,
    "beta 13–30 Hz",
    ha="center",
    va="bottom",
    fontsize=8.5,
    color=INK2,
)
a.set(
    xlabel="Frequency (Hz)",
    ylabel="Aperiodic-adjusted power (log10)",
    xlim=(1, 45),
    title="Spectra after removing the aperiodic component",
)
a.legend(loc="upper right")
plabel(a, "A")

b = fig.add_subplot(gs[0, 1])
meas = "aperiodic-adjusted, whole-scalp"
for lab, col in [("Higher vs HC", ORANGE), ("Lower vs HC", BLUE)]:
    b.plot(
        FREQ,
        G[(meas, lab)],
        color=col,
        lw=2,
        label=lab.replace(" vs HC", "-severity cluster vs HC"),
    )
    for r in cl[
        (cl["measure"] == meas) & (cl["contrast"] == lab) & (cl["p_cluster"] < 0.05)
    ].itertuples():
        b.axvspan(r.from_Hz, r.to_Hz, color=col, alpha=0.14, lw=0)
        b.text(
            (r.from_Hz + r.to_Hz) / 2,
            1.02 if col == ORANGE else -0.62,
            f"{r.from_Hz:.0f}–{r.to_Hz:.0f} Hz\n"
            + ("p < .001" if r.p_cluster < 0.001 else f"p = {r.p_cluster:.3f}".replace("0.", ".")),
            ha="center",
            va="top",
            fontsize=8.5,
            color=INK,
        )
b.axhline(0, color=INK2, lw=1)
b.set(
    xlabel="Frequency (Hz)",
    ylabel="Hedges' g vs healthy controls",
    xlim=(1, 45),
    ylim=(-0.75, 1.1),
    title="Group difference at each frequency (shaded: significant clusters)",
)
b.legend(loc="lower right")
plabel(b, "B")

c = fig.add_subplot(gs[1, 0])
im = c.imshow(
    HM,
    aspect="auto",
    cmap=DIV,
    vmin=-1,
    vmax=1,
    origin="lower",
    extent=[1, 45, -0.5, len(BINS) - 0.5],
)
c.set_yticks(range(len(BINS)))
c.set_yticklabels([f"{lo}–{hi}" if hi < 28 else f"≥ {lo}" for lo, hi in BINS])
for xx in (13, 30):
    c.axvline(xx, color=INK2, lw=0.8, ls=(0, (3, 3)))
c.set(
    xlabel="Frequency (Hz)",
    ylabel="ISI total score band",
    title="Difference from healthy controls by frequency and severity",
)
c.grid(False)
cb = fig.colorbar(im, ax=c, fraction=0.035, pad=0.02)
cb.outline.set_visible(False)
cb.set_label("Hedges' g", fontsize=9)
cb.ax.tick_params(length=0, labelsize=8)
plabel(c, "C")

d = fig.add_subplot(gs[1, 1])
mr = "raw log power, whole-scalp"
for m_, ls, nm in [(meas, "-", "aperiodic-adjusted"), (mr, (0, (4, 2)), "raw power")]:
    d.plot(
        FREQ,
        G[(m_, "Higher vs HC")],
        color=ORANGE,
        lw=2,
        ls=ls,
        label=f"Higher-severity cluster, {nm}",
    )
d.axhline(0, color=INK2, lw=1)
d.axvspan(13, 30, color=GRID, alpha=0.45, lw=0)
d.set(
    xlabel="Frequency (Hz)",
    ylabel="Hedges' g vs healthy controls",
    xlim=(1, 45),
    ylim=(-0.75, 1.1),
    title="With and without removing the aperiodic component",
)
d.legend(loc="lower right")
plabel(d, "D")
fig.savefig(FIG / "fig7_spectral_specificity.png", dpi=200, bbox_inches="tight")
plt.close(fig)

POS = {
    "FP1": (-0.25, 0.80),
    "FP2": (0.25, 0.80),
    "F7": (-0.70, 0.52),
    "F3": (-0.36, 0.45),
    "FZ": (0, 0.42),
    "F4": (0.36, 0.45),
    "F8": (0.70, 0.52),
    "T3": (-0.87, 0),
    "C3": (-0.43, 0),
    "CZ": (0, 0),
    "C4": (0.43, 0),
    "T4": (0.87, 0),
    "T5": (-0.70, -0.52),
    "P3": (-0.36, -0.45),
    "PZ": (0, -0.42),
    "P4": (0.36, -0.45),
    "T6": (0.70, -0.52),
    "O1": (-0.25, -0.80),
    "O2": (0.25, -0.80),
}
fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.4))
for ax, lab, ttl in zip(
    axes,
    ["Higher vs HC", "Lower vs HC"],
    ["Higher-severity cluster vs healthy controls", "Lower-severity cluster vs healthy controls"],
):
    t = tp[tp["contrast"] == lab].set_index("channel")
    ax.add_patch(Circle((0, 0), 1.0, fill=False, ec=AXIS, lw=1.5))
    ax.plot([-0.08, 0, 0.08], [0.995, 1.09, 0.995], color=AXIS, lw=1.5)
    for cname, (x, y) in POS.items():
        sig = t.loc[cname, "q_band"] < 0.05
        ax.scatter(
            x,
            y,
            s=900,
            c=[t.loc[cname, "g_band"]],
            cmap=DIV,
            vmin=-1,
            vmax=1,
            edgecolor=INK if sig else AXIS,
            lw=1.8 if sig else 0.8,
            zorder=3,
        )
        ax.text(
            x,
            y + 0.02,
            cname.capitalize().replace("Fp", "Fp"),
            ha="center",
            va="bottom",
            fontsize=7.5,
            color=INK2,
            zorder=4,
        )
        ax.text(
            x,
            y - 0.015,
            f"{t.loc[cname, 'g_band']:+.2f}",
            ha="center",
            va="top",
            fontsize=8.5,
            fontweight="bold" if sig else "normal",
            zorder=4,
        )
    ax.set(xlim=(-1.15, 1.15), ylim=(-1.1, 1.18))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(ttl, loc="center")
sm = plt.cm.ScalarMappable(cmap=DIV, norm=plt.Normalize(-1, 1))
cb = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02)
cb.outline.set_visible(False)
cb.set_label("Hedges' g, aperiodic-adjusted 13–30 Hz power", fontsize=9)
cb.ax.tick_params(length=0, labelsize=8)
fig.legend(
    handles=[
        Patch(fc=SURFACE, ec=INK, lw=1.8, label="FDR q < .05 across 19 electrodes"),
        Patch(fc=SURFACE, ec=AXIS, lw=0.8, label="not significant"),
    ],
    loc="lower center",
    ncol=2,
    bbox_to_anchor=(0.45, -0.02),
)
fig.savefig(FIG / "fig8_topography.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\n輸出 → {OUT}\n圖 → {FIG}")
