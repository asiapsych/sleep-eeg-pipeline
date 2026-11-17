"""Figure 2：beta 隨 ISI 的曲線，以及各分數帶對健康組的效果量。

樣條、共變項與 HC3 區間和 task_r1u 一致。
ISI ≥ 23 的人很少，圖上用灰底標出，避免把曲線末端看得太確定。"""

import statsmodels.formula.api as smf

from _r1_common import OUTCOME as Y, load_all
from _r1_submission import (
    AN,
    AXIS,
    BLUE,
    GRAY,
    INK,
    INK2,
    MUTED,
    ORANGE,
    SURFACE,
    np,
    panel_label,
    pd,
    plt,
    save_fig,
    sg,
)

df = load_all()
hc = df[df["ISI_T1"] < 9]
ins = df[df["ISI_T1"] >= 9].copy()
bins = pd.read_csv(AN / "S_R1a_severity_form" / "03_isi_bins.csv", encoding="utf-8-sig")
b2 = bins[bins["scheme"] == "2-point bins"]
pooled = float(bins.loc[bins["scheme"].str.startswith("all insomnia"), "g_vs_HC"].iloc[0])
up = pd.read_csv(AN / "S_R1u_shape_robust" / "05_upper_range.csv", encoding="utf-8-sig")
cols = ["ISI_range", "n", "g_vs_HC", "ci_low", "ci_high"]
b2 = pd.concat(
    [b2[~b2["ISI_range"].str.startswith(">=")][cols], up[up["ISI_range"] != "17-18"][cols]],
    ignore_index=True,
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
m_ = ins.dropna(subset=[Y, "AGE", "SEX", "EDU"]).copy()
m_["SEXN"] = m_["SEX"].astype(float)
m_["YV"] = m_[Y].astype(float)
fit = smf.ols("YV ~ AGE + SEXN + EDU + cr(ISI_T1, df=3, constraints='center')", m_).fit(
    cov_type="HC3"
)
assert (
    int(fit.df_model) == 6
)  # 3 個共變項 + 3 個樣條基底（直線 + 2 個非線性自由度），與 task_r1u 的檢定相同
g_ = pd.DataFrame(
    dict(
        ISI_T1=np.linspace(m_["ISI_T1"].min(), m_["ISI_T1"].max(), 200),
        AGE=m_["AGE"].mean(),
        SEXN=m_["SEXN"].mean(),
        EDU=m_["EDU"].mean(),
    )
)
pr = fit.get_prediction(g_).summary_frame(alpha=0.05)
a.axhline(hc[Y].mean(), color=MUTED, lw=1.2, ls=(0, (4, 3)), label="Mean of healthy controls")
a.fill_between(g_["ISI_T1"], pr["mean_ci_lower"], pr["mean_ci_upper"], color=INK, alpha=0.13, lw=0)
a.plot(
    g_["ISI_T1"],
    pr["mean"],
    color=INK,
    lw=2.2,
    label=f"Fitted curve, insomnia symptoms (n = {len(m_)}; spline, 95% CI)",
)
n23 = int((m_["ISI_T1"] >= 23).sum())
a.axvspan(22.5, 26.6, color=MUTED, alpha=0.10, lw=0)
a.text(24.55, 0.105, f"ISI ≥ 23\nn = {n23}", fontsize=8.5, color=INK2, ha="center", va="bottom")
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

for i, r in b2.iterrows():
    b.plot([i, i], [r["ci_low"], r["ci_high"]], color=INK2, lw=2, solid_capstyle="round")
    b.scatter(i, r["g_vs_HC"], s=70, color=INK2, zorder=3, edgecolor=SURFACE, lw=1.5)
    b.text(i, r["ci_high"] + 0.06, sg(r["g_vs_HC"]), ha="center", fontsize=9, fontweight="bold")
    b.text(i, -0.86, f"n = {int(r['n'])}", ha="center", fontsize=8.5, color=INK2)
b.axhline(0, color=INK2, lw=1)
b.axhline(
    pooled,
    color=MUTED,
    lw=1.2,
    ls=(0, (4, 3)),
    label=f"All insomnia symptoms pooled ({sg(pooled)})",
)
b.legend(loc="upper left", fontsize=8.5, handlelength=2.6)  # 寫在線旁會壓到 13–14 分的數值
b.set_xticks(range(len(b2)))
b.set_xticklabels([r.replace(">=", "≥ ").replace("-", "–") for r in b2["ISI_range"]])
b.set(
    xlabel="ISI total score band",
    ylabel="Hedges' g vs healthy controls (95% CI)",
    ylim=(-0.95, 2.85),
    title="Effect size by severity band",
)
b.grid(axis="x", visible=False)
panel_label(b, "B")
save_fig(fig, "Figure2")
