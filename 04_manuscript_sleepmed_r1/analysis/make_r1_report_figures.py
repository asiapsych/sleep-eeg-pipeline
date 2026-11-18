"""用 R1 分析結果畫報告用圖，存到 S_R1_figures。

送期刊的版本另外由 sub_fig*.py 產生；兩套圖讀同一批分析輸出。"""

from __future__ import annotations
import re
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from _r1_common import AN, ISI_ITEMS, OUTCOME, load_all

OUT = AN / "S_R1_figures"
OUT.mkdir(parents=True, exist_ok=True)

SURFACE, INK, INK2, MUTED, GRID, AXIS = (
    "#fcfcfb",
    "#0b0b0b",
    "#52514e",
    "#898781",
    "#e1e0d9",
    "#c3c2b7",
)
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
C_K3 = {"MDI": BLUE, "MEW": AQUA, "SCA": ORANGE}
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


def save(fig, name):
    fig.savefig(OUT / name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved", OUT / name)


def panel_label(ax, s):
    ax.text(
        -0.02,
        1.08,
        s,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="right",
        va="bottom",
    )


df = load_all()
hc = df[df["group"] == "HC"]
ins = df[df["group"] != "HC"].copy()
Y = OUTCOME

# Figure 1
bins = pd.read_csv(AN / "S_R1a_severity_form" / "03_isi_bins.csv")
b2 = bins[bins["scheme"] == "2-point bins"].reset_index(drop=True)
pooled = float(bins.loc[bins["scheme"].str.startswith("all insomnia"), "g_vs_HC"].iloc[0])
thr = pd.read_csv(AN / "S_R1a_severity_form" / "05_threshold.csv")
C_HAT = int(
    float(thr.loc[thr["quantity"].str.startswith("data-selected threshold"), "value"].iloc[0])
)

fig, (a, b) = plt.subplots(
    1, 2, figsize=(12.5, 4.6), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.24)
)
rng = np.random.default_rng(1)
jit = lambda n: rng.uniform(-0.28, 0.28, n)
d = hc.dropna(subset=[Y])
a.scatter(
    d["ISI_T1"] + jit(len(d)),
    d[Y],
    s=14,
    color=GRAY,
    alpha=0.45,
    lw=0,
    label=f"Healthy controls (n = {len(d)})",
)
for k, col, lab in [
    ("Lower", BLUE, "Lower-severity group"),
    ("Higher", ORANGE, "Higher-severity group"),
]:
    d = ins[ins["k2"] == k].dropna(subset=[Y])
    a.scatter(
        d["ISI_T1"] + jit(len(d)),
        d[Y],
        s=16,
        color=col,
        alpha=0.6,
        lw=0,
        label=f"{lab} (n = {len(d)})",
    )
# 連續的曲線：與正文的非線性檢定同一個模型、同一批人（task_r1u_shape_robust.py）：
#   具失眠症狀者（ISI ≥ 9）；自然三次樣條 cr(ISI, df=3, constraints='center')（相對於直線多 2 個非線性自由度）；
#   校正年齡、性別、教育年數；在共變項平均值上預測，95% CI 用 HC3 標準誤。
#   注意 patsy 的 cr(x, df=3) 不加 constraints='center' 時，含截距的模型只多 1 個非線性自由度，與檢定不是同一個模型。
# 健康組只畫平均值當參照（虛線）。不畫分段平均與 ISI 18 的垂直線：本資料分不出階梯與快速上升。
import statsmodels.formula.api as _smf

_d = ins.dropna(subset=[Y, "AGE", "SEX", "EDU"]).copy()
_d["SEXN"] = _d["SEX"].astype(float)
_m = _smf.ols(f"{Y} ~ AGE + SEXN + EDU + cr(ISI_T1, df=3, constraints='center')", _d).fit(
    cov_type="HC3"
)
assert int(_m.df_model) == 6, _m.df_model  # 3 個共變項 + 3 個樣條基底（= 直線 + 2 個非線性自由度）
_g = pd.DataFrame(
    dict(
        ISI_T1=np.linspace(_d["ISI_T1"].min(), _d["ISI_T1"].max(), 200),
        AGE=_d["AGE"].mean(),
        SEXN=_d["SEXN"].mean(),
        EDU=_d["EDU"].mean(),
    )
)
_p = _m.get_prediction(_g).summary_frame(alpha=0.05)
a.axhline(hc[Y].mean(), color=MUTED, lw=1.2, ls=(0, (4, 3)), label="Mean of healthy controls")
a.fill_between(_g["ISI_T1"], _p["mean_ci_lower"], _p["mean_ci_upper"], color=INK, alpha=0.13, lw=0)
a.plot(
    _g["ISI_T1"],
    _p["mean"],
    color=INK,
    lw=2.2,
    label=f"Fitted curve, insomnia symptoms (n = {len(_d)}; spline, 95% CI)",
)
# 23 分以上只有少數人：在圖上標出來，讓讀者知道曲線末端是靠幾個人估的
_n23 = int((_d["ISI_T1"] >= 23).sum())
a.axvspan(22.5, 26.6, color=MUTED, alpha=0.10, lw=0)
a.text(24.55, 0.105, f"ISI ≥ 23\nn = {_n23}", fontsize=8.5, color=INK2, ha="center", va="bottom")
a.axvline(8.5, color=AXIS, lw=1, ls=(0, (3, 3)))
a.text(8.7, 0.98, "ISI ≥ 9", fontsize=8.5, color=INK2, va="top")
a.set(
    xlabel="ISI total score",
    ylabel="Periodic beta peak power (whole-scalp)",
    xlim=(-1, 27),
    ylim=(0.08, 1.0),
    title="Periodic beta across insomnia severity",
)
a.legend(loc="upper left", bbox_to_anchor=(0, -0.16), ncol=2, handletextpad=0.4, columnspacing=1.2)
panel_label(a, "A")

# 19 分以上再細分（task_r1u 的 05_upper_range.csv），讓 B 與 A 的曲線講同一件事
_up = pd.read_csv(AN / "S_R1u_shape_robust" / "05_upper_range.csv", encoding="utf-8-sig")
_up = _up[_up["ISI_range"] != "17-18"][["ISI_range", "n", "g_vs_HC", "ci_low", "ci_high"]]
b2 = pd.concat(
    [
        b2[~b2["ISI_range"].str.startswith(">=")][
            ["ISI_range", "n", "g_vs_HC", "ci_low", "ci_high"]
        ],
        _up,
    ],
    ignore_index=True,
)
x = np.arange(len(b2))
cols = [INK2] * len(b2)  # 各分段同一個顏色，不預設哪裡是分界
for i, r in b2.iterrows():
    b.plot([i, i], [r["ci_low"], r["ci_high"]], color=cols[i], lw=2, solid_capstyle="round")
    b.scatter(i, r["g_vs_HC"], s=70, color=cols[i], zorder=3, edgecolor=SURFACE, lw=1.5)
    b.text(
        i, r["ci_high"] + 0.06, f"{r['g_vs_HC']:+.2f}", ha="center", fontsize=9, fontweight="bold"
    )
    b.text(i, -0.86, f"n = {int(r['n'])}", ha="center", fontsize=8.5, color=INK2)
b.axhline(0, color=INK2, lw=1)
b.axhline(pooled, color=MUTED, lw=1.2, ls=(0, (4, 3)))
b.text(
    2.0,
    pooled + 0.05,
    f"all insomnia symptoms pooled: {pooled:+.2f}",
    ha="center",
    fontsize=8.5,
    color=INK2,
)
b.set_xticks(x)
b.set_xticklabels([r.replace(">=", "≥ ").replace("-", "–") for r in b2["ISI_range"]])
b.set(
    xlabel="ISI total score band",
    ylabel="Hedges' g vs healthy controls (95% CI)",
    ylim=(-0.95, 2.5),
    title="Effect size by severity band",
)
b.grid(axis="x", visible=False)
panel_label(b, "B")
save(fig, "fig1_severity_shape.png")

# Figure 2
bc = pd.read_csv(AN / "S_R1b_k234" / "03_beta_by_cluster.csv")
rows = []
for sol in ["K=2", "K=3", "K=4"]:
    t = bc[bc["solution"] == sol].sort_values("ISI_mean")
    for _, r in t.iterrows():
        rows.append((sol, r["cluster"], r))
fig, ax = plt.subplots(figsize=(8.6, 5.2))
ypos, labels, y = [], [], 0
NAME = {
    "Lower": "Lower-severity",
    "Higher": "Higher-severity",
    "K4_C1": "Cluster 1",
    "K4_C2": "Cluster 2",
    "K4_C3": "Cluster 3",
    "K4_C4": "Cluster 4",
}
prev = None
for sol, cl, r in rows:
    if prev is not None and sol != prev:
        y -= 0.7
    if sol != prev:
        ax.text(
            -0.575,
            y + 0.62,
            sol.replace("=", " = ") + (" (selected by internal indices)" if sol == "K=2" else ""),
            fontsize=10,
            fontweight="bold",
            va="center",
            ha="right",
            clip_on=False,
        )
    prev = sol
    col = ORANGE if r["ISI_mean"] >= 17 else BLUE
    ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=col, lw=2, solid_capstyle="round")
    ax.scatter(r["g_vs_HC"], y, s=70, color=col, zorder=3, edgecolor=SURFACE, lw=1.5)
    ax.text(
        1.42,
        y,
        f"{r['g_vs_HC']:+.2f} [{r['ci_low']:+.2f}, {r['ci_high']:+.2f}]",
        va="center",
        fontsize=9,
    )
    ypos.append(y)
    labels.append(f"{NAME.get(cl, cl)}  (n = {int(r['n_beta'])}, ISI {r['ISI_mean']:.1f})")
    y -= 1
ax.axvline(0, color=INK2, lw=1)
ax.axvline(pooled, color=MUTED, lw=1.2, ls=(0, (4, 3)))
ax.text(pooled + 0.02, y + 0.35, f"pooled g = {pooled:+.2f}", fontsize=8.5, color=INK2)
ax.set_yticks(ypos)
ax.set_yticklabels(labels)
ax.set_ylim(y + 0.2, 1.1)
ax.set_xlim(-0.55, 1.4)
ax.set(
    xlabel="Hedges' g vs healthy controls (95% CI)",
    title="Periodic beta by cluster under K = 2, 3 and 4",
)
ax.scatter([], [], s=60, color=ORANGE, label="Cluster mean ISI ≥ 17")
ax.scatter([], [], s=60, color=BLUE, label="Cluster mean ISI < 17")
ax.legend(loc="upper left", bbox_to_anchor=(0, -0.1), ncol=2)
ax.grid(axis="y", visible=False)
ax.tick_params(axis="y", length=0)
save(fig, "fig2_k234_forest.png")

# Figure 3
note = thr.loc[thr["quantity"].str.startswith("bootstrap distribution"), "note"].iloc[0]
dist = {int(k): float(v) for k, v in re.findall(r"(\d+): ([\d.]+)%", note)}
mc = pd.read_csv(AN / "S_R1a_severity_form" / "04_model_comparison.csv")
mi = mc[mc["scope"] == "insomnia only"].copy()
lin = float(mi.loc[mi["model"] == "linear ISI total", "AIC"].iloc[0])
keep = [
    "linear ISI total",
    "quadratic ISI total",
    "spline df=3",
    "spline df=4",
    "spline df=5",
    "single threshold ISI>=15 (conventional)",
    f"single threshold ISI>={C_HAT} (data-selected)",
    "K=2 higher-severity cluster",
    "K=3 profile (3 levels)",
    "SCA indicator",
    "PC1",
]
mi = mi[mi["model"].isin(keep)].set_index("model").loc[keep].reset_index()
mi["dAIC"] = mi["AIC"] - lin
LAB = {
    "linear ISI total": "Linear ISI total (reference)",
    "quadratic ISI total": "Quadratic",
    "spline df=3": "Spline, 3 df",
    "spline df=4": "Spline, 4 df",
    "spline df=5": "Spline, 5 df",
    "single threshold ISI>=15 (conventional)": "Threshold at ISI ≥ 15 (conventional)",
    f"single threshold ISI>={C_HAT} (data-selected)": f"Threshold at ISI ≥ {C_HAT} (data-selected)",
    "K=2 higher-severity cluster": "K = 2 cluster membership",
    "K=3 profile (3 levels)": "K = 3 profile membership",
    "SCA indicator": "Highest-severity profile vs rest",
    "PC1": "First principal component",
}
fig, (a, b) = plt.subplots(
    1, 2, figsize=(12.5, 4.4), gridspec_kw=dict(width_ratios=[1, 1.15], wspace=0.75)
)
ks = sorted(dist)
a.bar(ks, [dist[k] for k in ks], width=0.72, color=[ORANGE if 17 <= k <= 20 else BLUE for k in ks])
for k in ks:
    if dist[k] >= 1:
        a.text(
            k,
            dist[k] + 1.2,
            f"{dist[k]:.0f}%",
            ha="center",
            fontsize=9,
            fontweight="bold" if k == C_HAT else "normal",
        )
a.set_xticks(ks)
a.set(
    xlabel="Selected ISI threshold",
    ylabel="Bootstrap resamples (%)",
    ylim=(0, 70),
    title="Where the step falls (5000 bootstrap resamples)",
)
from matplotlib.patches import Patch

a.legend(
    handles=[
        Patch(color=ORANGE, label="ISI 17–20: 94% of resamples"),
        Patch(color=BLUE, label="Other values"),
    ],
    loc="upper left",
)
a.grid(axis="x", visible=False)
panel_label(a, "A")
yy = np.arange(len(mi))[::-1]
b.barh(yy, mi["dAIC"], height=0.62, color=BLUE)
for yi, v in zip(yy, mi["dAIC"]):
    b.text(
        v - 0.25 if v < 0 else v + 0.25,
        yi,
        f"{v:+.1f}",
        va="center",
        ha="right" if v < 0 else "left",
        fontsize=9,
    )
b.axvline(0, color=INK2, lw=1)
b.set_yticks(yy)
b.set_yticklabels([LAB[m] for m in mi["model"]])
b.set_xlim(-10.5, 3.5)
b.set(
    xlabel="AIC difference from the linear model (lower is better)",
    title="How severity relates to periodic beta",
)
b.grid(axis="y", visible=False)
b.tick_params(axis="y", length=0)
panel_label(b, "B")
save(fig, "fig3_threshold_models.png")

# Figure 4
mt = pd.read_csv(AN / "S_R1e_clustering_robustness" / "01_methods_k2_k3.csv")
meth = [m for m in mt["method"].unique() if m != "K-means"]
fig, (a, b) = plt.subplots(1, 2, figsize=(12.8, 4.6), gridspec_kw=dict(wspace=0.62))
yy = np.arange(len(meth))[::-1]
for k, col, off in [(2, ORANGE, 0.16), (3, BLUE, -0.16)]:
    t = mt[(mt["K"] == k) & (mt["method"] != "K-means")].set_index("method").loc[meth]
    a.hlines(
        yy + off, t["ARI_ci_low"].astype(float), t["ARI_ci_high"].astype(float), color=col, lw=2
    )
    a.scatter(
        t["ARI_vs_Kmeans"],
        yy + off,
        s=60,
        color=col,
        zorder=3,
        edgecolor=SURFACE,
        lw=1.5,
        label=f"K = {k}",
    )
    for yi, v in zip(yy + off, t["ARI_vs_Kmeans"]):
        a.text(1.04, yi, f"{v:.2f}", va="center", fontsize=8.5, color=INK)
for v, s in [(0.6, "moderate"), (0.8, "high")]:
    a.axvline(v, color=AXIS, lw=1, ls=(0, (3, 3)))
    a.text(v + 0.01, len(meth) - 0.45, s, fontsize=8, color=INK2)
a.set_yticks(yy)
a.set_yticklabels(meth)
a.set_xlim(0, 1.02)
a.set_ylim(-0.6, len(meth) - 0.25)
a.set(
    xlabel="Adjusted Rand index vs K-means (95% bootstrap CI)",
    title="Agreement between clustering methods",
)
a.legend(loc="upper left", bbox_to_anchor=(0, -0.14), ncol=2)
a.grid(axis="y", visible=False)
a.tick_params(axis="y", length=0)
panel_label(a, "A")

t = mt[mt["K"] == 2].reset_index(drop=True)
yy = np.arange(len(t))[::-1]
ci = lambda s: [float(v) for v in re.findall(r"-?\d+\.\d+", s)]
for i, r in t.iterrows():
    for g, c, col, off in [
        (r["largest_g"], ci(r["largest_ci"]), ORANGE, 0.14),
        (r["smallest_g"], ci(r["smallest_ci"]), BLUE, -0.14),
    ]:
        b.plot(c, [yy[i] + off] * 2, color=col, lw=2, solid_capstyle="round")
        b.scatter(g, yy[i] + off, s=60, color=col, zorder=3, edgecolor=SURFACE, lw=1.5)
b.scatter([], [], s=60, color=ORANGE, label="Higher-severity cluster")
b.scatter([], [], s=60, color=BLUE, label="Lower-severity cluster")
b.axvline(0, color=INK2, lw=1)
b.set_yticks(yy)
b.set_yticklabels([f"{m}  ({s})" for m, s in zip(t["method"], t["sizes"])])
b.set_xlim(-0.3, 1.4)
b.set(
    xlabel="Hedges' g vs healthy controls (95% CI)", title="Two-cluster solution under nine methods"
)
b.legend(loc="upper left", bbox_to_anchor=(0, -0.14), ncol=2)
b.grid(axis="y", visible=False)
b.tick_params(axis="y", length=0)
panel_label(b, "B")
save(fig, "fig4_method_agreement.png")

# Figure 5
Z = StandardScaler().fit_transform(ins[ISI_ITEMS])
pca = PCA().fit(Z)
S = pca.transform(Z)
comp = pca.components_.copy()
if np.corrcoef(S[:, 0], ins["ISI_T1"])[0, 1] < 0:
    S[:, 0] *= -1
if comp[1][1] + comp[1][2] < 0:
    S[:, 1] *= -1
pb = pd.read_csv(AN / "S_R1f_pc_biomarkers" / "01_pc_by_biomarker_insomnia.csv")
BL = {
    "beta_pw_BRAIN": "Periodic beta, whole-brain",
    "beta_pw_CENTRAL": "Periodic beta, central",
    "beta_pw_MIDLINE": "Periodic beta, midline",
    "Beta_rel_BRAIN": "Relative beta, whole-brain",
    "Beta_rel_CENTRAL": "Relative beta, central",
    "exponent_BRAIN": "Aperiodic exponent, whole-brain",
    "exponent_CENTRAL": "Aperiodic exponent, central",
    "offset_BRAIN": "Aperiodic offset, whole-brain",
    "offset_CENTRAL": "Aperiodic offset, central",
    "alpha_pw_BRAIN": "Periodic alpha",
    "alpha_cf_BRAIN": "Alpha peak frequency",
    "theta_pw_BRAIN": "Periodic theta",
    "Theta_rel_BRAIN": "Relative theta",
    "Alpha_rel_BRAIN": "Relative alpha",
    "Delta_rel_BRAIN": "Relative delta",
    "Gamma_rel_BRAIN": "Relative 30–45 Hz",
}
fig, (a, b) = plt.subplots(
    1, 2, figsize=(12.8, 5.6), gridspec_kw=dict(width_ratios=[1, 0.62], wspace=0.95)
)
for g in ["MDI", "MEW", "SCA"]:
    m = (ins["group"] == g).to_numpy()
    a.scatter(
        S[m, 0],
        S[m, 1],
        s=26,
        color=C_K3[g],
        alpha=0.75,
        edgecolor=SURFACE,
        lw=0.6,
        label=f"{g} (n = {m.sum()})",
    )
    a.text(
        S[m, 0].mean(),
        S[m, 1].mean(),
        g,
        fontsize=11,
        fontweight="bold",
        ha="center",
        va="center",
        bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec=C_K3[g], lw=1.5),
    )
ev = pca.explained_variance_ratio_
a.set(
    xlabel=f"PC1 ({ev[0]:.0%} of variance; r = .97 with ISI total)",
    ylabel=f"PC2 ({ev[1]:.0%}; maintenance and early waking vs daytime items)",
    title="Symptom space of the K = 3 profiles",
)
a.legend(loc="upper left", bbox_to_anchor=(0, -0.12), ncol=3)
panel_label(a, "A")
M = pb[["PC1_beta", "PC2_beta", "PC3_beta"]].to_numpy()
Q = pb[["PC1_q", "PC2_q", "PC3_q"]].to_numpy()
cmap = LinearSegmentedColormap.from_list(
    "div", ["#1c5cab", "#86b6ef", "#f0efec", "#f3a583", "#c4481a"]
)
im = b.imshow(M, cmap=cmap, vmin=-0.3, vmax=0.3, aspect="auto")
for i in range(M.shape[0]):
    for j in range(3):
        b.text(
            j,
            i,
            f"{M[i, j]:+.2f}" + ("*" if Q[i, j] < 0.05 else ""),
            ha="center",
            va="center",
            fontsize=8.5,
            fontweight="bold" if Q[i, j] < 0.05 else "normal",
            color=INK,
        )
b.set_xticks(range(3))
b.set_xticklabels(["PC1", "PC2", "PC3"])
b.set_yticks(range(len(pb)))
b.set_yticklabels([BL[f] for f in pb["biomarker"]])
b.set_title("Standardised coefficients (* FDR q < .05)")
b.grid(False)
b.tick_params(length=0)
for s in b.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=b, fraction=0.05, pad=0.04)
cb.outline.set_visible(False)
cb.ax.tick_params(length=0, labelsize=8)
panel_label(b, "B")
save(fig, "fig5_pc_biomarkers.png")

# Figure 6
e = pd.read_csv(AN / "S_R1c_sensitivity" / "03_emg_covariate.csv").rename(
    columns={"added_covariate": "setting"}
)
v = pd.read_csv(AN / "S_R1c_sensitivity" / "05b_vigilance_covariate.csv").rename(
    columns={"added_covariate": "setting"}
)
q = pd.read_csv(AN / "S_R1d_preproc_qc" / "02_qc_covariate.csv")
CON = "K2 Higher vs HC"
SETS = [
    (e, "none (age, sex, education only)", "Age, sex, education only"),
    (e, "+ 30–45 Hz relative power", "+ 30–45 Hz relative power"),
    (e, "+ log10 30–45 Hz absolute power", "+ 30–45 Hz absolute power (log)"),
    (e, "+ aperiodic exponent (whole-brain)", "+ aperiodic exponent"),
    (v, "+ theta/alpha ratio", "+ theta/alpha ratio"),
    (v, "+ theta relative power", "+ relative theta"),
    (v, "+ alpha peak frequency", "+ alpha peak frequency"),
    (q, "+ n_ic_removed", "+ components removed"),
    (q, "+ n_ret_muscle_dom", "+ retained muscle components"),
    (q, "+ max_ret_muscle_p", "+ max retained muscle probability"),
    (q, "+ reject_rate", "+ epoch rejection rate"),
    (q, "+ all of the above QC metrics", "+ all preprocessing metrics"),
]
rs = pd.read_csv(AN / "S_R1d_preproc_qc" / "03_restricted_samples.csv")
SUB = [
    ("full sample", "Full sample"),
    ("no retained muscle-dominant component", "No retained muscle component"),
    ("max retained muscle probability < 0.50", "Max muscle probability < .50"),
    ("max retained muscle probability < 0.20", "Max muscle probability < .20"),
    ("no interpolated channel", "No interpolated channel"),
    ("epoch rejection rate <= 10%", "Epoch rejection ≤ 10%"),
]
fig, (a, b) = plt.subplots(1, 2, figsize=(12.8, 4.8), gridspec_kw=dict(wspace=0.72))
yy = np.arange(len(SETS))[::-1]
GC = [GRAY] + [ORANGE] * 3 + [AQUA] * 3 + [BLUE] * 5
for yi, (t, s, lab), col in zip(yy, SETS, GC):
    r = t[(t["setting"] == s) & (t["contrast"] == CON)].iloc[0]
    sd = r["b"] / r["b_in_HC_SD"]
    a.plot([r["ci_low"] / sd, r["ci_high"] / sd], [yi, yi], color=col, lw=2, solid_capstyle="round")
    a.scatter(r["b_in_HC_SD"], yi, s=60, color=col, zorder=3, edgecolor=SURFACE, lw=1.5)
    a.text(1.36, yi, f"{r['b_in_HC_SD']:.2f}", va="center", fontsize=9)
for col, lab in [
    (GRAY, "Base model"),
    (ORANGE, "Muscle proxies"),
    (AQUA, "Vigilance markers"),
    (BLUE, "Preprocessing metrics"),
]:
    a.scatter([], [], s=60, color=col, label=lab)
a.axvline(0, color=INK2, lw=1)
a.set_yticks(yy)
a.set_yticklabels([s[2] for s in SETS])
a.set_xlim(-0.05, 1.34)
a.set(
    xlabel="Adjusted difference, higher-severity cluster − HC (HC SD units, 95% CI)",
    title="Adding covariates",
)
a.legend(
    loc="upper left", bbox_to_anchor=(-0.55, -0.14), ncol=4, columnspacing=1.0, handletextpad=0.3
)
a.grid(axis="y", visible=False)
a.tick_params(axis="y", length=0)
panel_label(a, "A")
yy = np.arange(len(SUB))[::-1]
for yi, (s, lab) in zip(yy, SUB):
    for con, col, off in [("K2 Higher vs HC", ORANGE, 0.14), ("K2 Lower vs HC", BLUE, -0.14)]:
        r = rs[(rs["subsample"] == s) & (rs["contrast"] == con)].iloc[0]
        b.plot([r["ci_low"], r["ci_high"]], [yi + off] * 2, color=col, lw=2, solid_capstyle="round")
        b.scatter(r["g"], yi + off, s=60, color=col, zorder=3, edgecolor=SURFACE, lw=1.5)
        b.text(1.22, yi + off, f"{r['g']:+.2f}  (n = {int(r['n_a'])})", va="center", fontsize=8.5)
b.scatter([], [], s=60, color=ORANGE, label="Higher-severity cluster")
b.scatter([], [], s=60, color=BLUE, label="Lower-severity cluster")
b.axvline(0, color=INK2, lw=1)
b.set_yticks(yy)
b.set_yticklabels([s[1] for s in SUB])
b.set_xlim(-0.35, 1.2)
b.set(xlabel="Hedges' g vs healthy controls (95% CI)", title="Restricting to cleaner recordings")
b.legend(loc="upper left", bbox_to_anchor=(0, -0.14), ncol=2)
b.grid(axis="y", visible=False)
b.tick_params(axis="y", length=0)
panel_label(b, "B")
save(fig, "fig6_robustness.png")

# Figure 10：樣本的嚴重度組成與合併比較（主文圖 6）
cost = pd.read_csv(AN / "S_R1u_shape_robust" / "03_pooling_cost.csv", encoding="utf-8-sig")
cost["c"] = cost["cutoff"].str.extract(r"(\d+)").astype(int)
BOUND = 0.37  # 統合分析的清醒 beta 效果量（Zhao et al., 2021），也是等效檢定的界值
fig, a = plt.subplots(figsize=(8.4, 4.8))
a.fill_between(cost["c"], cost["ci_low"], cost["ci_high"], color=ORANGE, alpha=0.16, lw=0)
a.plot(
    cost["c"],
    cost["g"],
    color=ORANGE,
    lw=2.2,
    marker="o",
    ms=5.5,
    label="Included: ISI at or above the minimum",
)
m = cost["g_excluded"].notna()
a.fill_between(
    cost.loc[m, "c"],
    cost.loc[m, "ci_low_excluded"],
    cost.loc[m, "ci_high_excluded"],
    color=BLUE,
    alpha=0.13,
    lw=0,
)
a.plot(
    cost.loc[m, "c"],
    cost.loc[m, "g_excluded"],
    color=BLUE,
    lw=2.2,
    marker="o",
    ms=5.5,
    label="Not included: ISI from 9 to below the minimum",
)
a.axhline(0, color=INK2, lw=1)
a.axhline(BOUND, color=MUTED, lw=1.1, ls=(0, (4, 3)))
a.text(
    20.65,
    BOUND,
    "Meta-analytic\nestimate\n(g = 0.37)",
    fontsize=8.5,
    color=INK2,
    va="center",
    ha="left",
)
for c_, lab in [
    (9, "All participants with\ninsomnia symptoms"),
    (15, "Moderate clinical\ninsomnia (ISI ≥ 15)"),
]:
    r = cost[cost["c"] == c_].iloc[0]
    a.annotate(
        f"{lab}\ng = {r['g']:+.2f}",
        (c_, r["g"]),
        xytext=(c_ + (0.1 if c_ == 9 else -0.2), r["ci_high"] + 0.33),
        fontsize=8.5,
        color=INK2,
        ha="left" if c_ == 9 else "center",
        va="bottom",
        arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.9),
    )
a.set_xticks(cost["c"])
a.set_xlim(8.5, 20.5)
a.set_ylim(-0.35, 1.75)
a.set(
    xlabel="Minimum ISI total required for inclusion in the insomnia-symptom group",
    ylabel="Hedges' g vs healthy controls (95% CI)",
    title="Pooled difference in periodic beta by severity composition",
)
a.legend(loc="upper left", bbox_to_anchor=(0, -0.16), ncol=2, columnspacing=1.2, handletextpad=0.4)
save(fig, "fig10_pooling_cost.png")

# Figure 11：症狀的兩個向度（主文圖 1；讀 task_r1w 的輸出）
# 版面：左邊一格放負荷量；右邊上下兩格，上面是第二個向度的外部關聯（只畫校正後顯著的），下面是兩個向度對週期性 beta。
# 圖上的字盡量少，說明放圖說。
W_ = AN / "S_R1w_pc_interpretation"
ld = pd.read_csv(W_ / "01_loadings.csv", index_col=0, encoding="utf-8-sig")
ex = pd.read_csv(W_ / "03_external_correlates.csv", encoding="utf-8-sig")
pvb = pd.read_csv(
    W_ / "07_two_component_model.csv", encoding="utf-8-sig"
)  # 只放平行分析保留的兩個成分
items = [i for i in ld.index if i[0].isdigit()]
SHORT = {
    "Sex (0 = female, 1 = male)": "Male sex",
    "DBAS (dysfunctional beliefs about sleep)": "Dysfunctional beliefs (DBAS)",
    "PSQI sleep duration (higher = shorter)": "Shorter sleep duration (PSQI)",
}
fig = plt.figure(figsize=(11.6, 5.0))
gs = fig.add_gridspec(
    2, 2, width_ratios=[1.0, 0.95], height_ratios=[3, 2.2], wspace=0.62, hspace=0.75
)
a = fig.add_subplot(gs[:, 0])
b = fig.add_subplot(gs[0, 1])
c = fig.add_subplot(gs[1, 1])
yy = np.arange(len(items))[::-1]
a.barh(yy + 0.19, ld.loc[items, "PC1"], height=0.36, color=GRAY, label="Dimension 1 (severity)")
a.barh(
    yy - 0.19,
    ld.loc[items, "PC2"],
    height=0.36,
    color=ORANGE,
    label="Dimension 2 (night-time vs daytime)",
)
a.axvline(0, color=INK2, lw=1)
a.set_yticks(yy)
a.set_yticklabels([i[2:] for i in items])
a.set_xlim(-0.55, 0.85)
a.set(xlabel="Component loading", title="ISI item loadings")
a.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.13), ncol=1)
a.grid(axis="y", visible=False)
a.tick_params(axis="y", length=0)
panel_label(a, "A")
e2 = (
    ex[ex["q_PC2"] < 0.05].sort_values("rho_PC2").reset_index(drop=True)
)  # 主文只畫校正後顯著的；13 個變項的完整結果在補充表 S1
z = np.arctanh(e2["rho_PC2"])
se = 1.06 / np.sqrt(e2["n"] - 3)  # Spearman 的 Fisher z 近似標準誤
lo_, hi_ = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
for i, r in e2.iterrows():
    b.plot([lo_[i], hi_[i]], [i, i], color=ORANGE, lw=2.2, solid_capstyle="round")
    b.scatter(r["rho_PC2"], i, s=64, color=ORANGE, zorder=3, edgecolor=SURFACE, lw=1.5)
    if r["rho_PC2"] >= 0:
        b.text(hi_[i] + 0.03, i, f"{r['rho_PC2']:+.2f}", va="center", fontsize=9.5)
    else:
        b.text(
            lo_[i] - 0.03,
            i,
            f"{r['rho_PC2']:+.2f}".replace("-", "−"),
            va="center",
            ha="right",
            fontsize=9.5,
        )
b.axvline(0, color=INK2, lw=1)
b.set_yticks(range(len(e2)))
b.set_yticklabels([SHORT.get(v, v) for v in e2["variable"]])
b.set_ylim(-0.6, len(e2) - 0.4)
b.set_xlim(-0.55, 0.7)
b.set(
    xlabel="Spearman ρ with dimension 2 (95% CI)", title="Dimension 2 and variables outside the ISI"
)
b.grid(axis="y", visible=False)
b.tick_params(axis="y", length=0)
panel_label(b, "B")
pp = pvb[pvb["component"].isin(["PC1", "PC2"])].reset_index(drop=True)
for i, (r, col) in enumerate(zip(pp.itertuples(), [GRAY, ORANGE])):
    c.plot([r.ci_low, r.ci_high], [1 - i, 1 - i], color=col, lw=2.4, solid_capstyle="round")
    c.scatter(r.std_beta, 1 - i, s=80, color=col, zorder=3, edgecolor=SURFACE, lw=1.5)
    c.text(
        r.ci_high + 0.03, 1 - i, f"{r.std_beta:+.2f}".replace("-", "−"), va="center", fontsize=9.5
    )
c.axvline(0, color=INK2, lw=1)
c.set_yticks([1, 0])
c.set_yticklabels(["Dimension 1", "Dimension 2"])
c.set_ylim(-0.6, 1.6)
c.set_xlim(-0.3, 0.62)
c.set(xlabel="Standardised coefficient (95% CI)", title="Association with periodic beta")
c.grid(axis="y", visible=False)
c.tick_params(axis="y", length=0)
panel_label(c, "C")
save(fig, "fig11_symptom_dimensions.png")

# Figure 12（補充圖 S2）：主要與次要 beta 指標隨 ISI 的變化
# 三個面板用同一批人（失眠症狀者）、同一組共變項、同一套樣條（與 task_r1u 的檢定相同），各自保留縱軸；
# 都標出 ISI ≥ 23（n = 8）。不把三條線疊在同一個座標上：三個指標的量尺與意義不同。
from _r1_spectra import band_mean as _bm, load_spectra as _ls
import statsmodels.formula.api as _smf2

_fq, _cn, _ad, _rw = _ls(df["subject"])
_dd = df.copy()
_dd["adj_13_30"] = _bm(_ad, _fq, 13, 30)
_dd["SEXN"] = _dd["SEX"].astype(float)
_MEAS = [
    (Y, "Periodic beta peak power", "Primary outcome"),
    (
        "adj_13_30",
        "Aperiodic-adjusted beta power, 13–30 Hz (log10)",
        "Does not require peak detection",
    ),
    ("Beta_rel_BRAIN", "Relative beta power (%)", "Share of total 1–45 Hz power"),
]
fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.7), gridspec_kw=dict(wspace=0.30))
_rng = np.random.default_rng(1)
for ax, (col, ylab, sub), pl in zip(axes, _MEAS, "ABC"):
    _h = _dd[_dd["ISI_T1"] < 9].dropna(subset=[col])
    _i = _dd[_dd["ISI_T1"] >= 9].dropna(subset=[col, "AGE", "SEXN", "EDU"])
    ax.scatter(
        _h["ISI_T1"] + _rng.uniform(-0.28, 0.28, len(_h)),
        _h[col],
        s=12,
        color=GRAY,
        alpha=0.4,
        lw=0,
        label="Healthy controls",
    )
    ax.scatter(
        _i["ISI_T1"] + _rng.uniform(-0.28, 0.28, len(_i)),
        _i[col],
        s=14,
        color=BLUE,
        alpha=0.5,
        lw=0,
        label="Insomnia symptoms",
    )
    _i = _i.assign(YV=_i[col].astype(float))  # 不用 patsy 的 Q()：本腳本前面有同名的變數
    _mm = _smf2.ols("YV ~ AGE + SEXN + EDU + cr(ISI_T1, df=3, constraints='center')", _i).fit(
        cov_type="HC3"
    )
    _gg = pd.DataFrame(
        dict(
            ISI_T1=np.linspace(9, _i["ISI_T1"].max(), 200),
            AGE=_i["AGE"].mean(),
            SEXN=_i["SEXN"].mean(),
            EDU=_i["EDU"].mean(),
        )
    )
    _pp = _mm.get_prediction(_gg).summary_frame(alpha=0.05)
    ax.axhline(
        _h[col].mean(), color=MUTED, lw=1.2, ls=(0, (4, 3)), label="Mean of healthy controls"
    )
    ax.fill_between(
        _gg["ISI_T1"], _pp["mean_ci_lower"], _pp["mean_ci_upper"], color=INK, alpha=0.13, lw=0
    )
    ax.plot(_gg["ISI_T1"], _pp["mean"], color=INK, lw=2.2, label="Fitted curve (spline, 95% CI)")
    ax.axvspan(22.5, 26.6, color=MUTED, alpha=0.10, lw=0)
    lo_, hi_ = np.nanpercentile(_dd[col], [0.5, 99.5])
    pad = 0.08 * (hi_ - lo_)
    ax.set_ylim(lo_ - 2.2 * pad, hi_ + pad)
    ax.text(
        24.55,
        lo_ - 2.0 * pad,
        "ISI ≥ 23\nn = 8",
        fontsize=8.5,
        color=INK2,
        ha="center",
        va="bottom",
    )
    ax.set(xlabel="ISI total score", ylabel=ylab, xlim=(-1, 27), title=sub)
    panel_label(ax, pl)
axes[0].legend(
    loc="upper left", bbox_to_anchor=(0, -0.16), ncol=4, handletextpad=0.4, columnspacing=1.2
)
save(fig, "fig12_severity_by_measure.png")
