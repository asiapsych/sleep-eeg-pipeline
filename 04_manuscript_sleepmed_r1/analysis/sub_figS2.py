"""Figure S2：峰值、去背景與相對 beta 功率的 ISI 曲線。

三個面板用同一批人、同一套樣條，各自保留縱軸；灰底是 ISI ≥ 23 的 8 人。"""

import statsmodels.formula.api as smf

from _r1_common import OUTCOME as Y, load_all
from _r1_spectra import band_mean, load_spectra
from _r1_submission import (
    BLUE,
    GRAY,
    INK,
    INK2,
    MUTED,
    np,
    panel_label,
    pd,
    plt,
    save_fig,
)

df = load_all()
freq, chn, adj, raw = load_spectra(df["subject"])
df["adj_13_30"] = band_mean(adj, freq, 13, 30)
df["SEXN"] = df["SEX"].astype(float)
MEAS = [
    (Y, "Periodic beta peak power", "Primary outcome"),
    (
        "adj_13_30",
        "Aperiodic-adjusted beta power, 13–30 Hz (log10)",
        "Does not require peak detection",
    ),
    ("Beta_rel_BRAIN", "Relative beta power (%)", "Share of total 1–45 Hz power"),
]
fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.7), gridspec_kw=dict(wspace=0.30))
rng = np.random.default_rng(1)
for ax, (col, ylab, sub), pl in zip(axes, MEAS, "ABC"):
    h = df[df["ISI_T1"] < 9].dropna(subset=[col])
    i = df[df["ISI_T1"] >= 9].dropna(subset=[col, "AGE", "SEXN", "EDU"])
    i = i.assign(YV=i[col].astype(float))
    ax.scatter(
        h["ISI_T1"] + rng.uniform(-0.28, 0.28, len(h)),
        h[col],
        s=12,
        color=GRAY,
        alpha=0.4,
        lw=0,
        label="Healthy controls",
    )
    ax.scatter(
        i["ISI_T1"] + rng.uniform(-0.28, 0.28, len(i)),
        i[col],
        s=14,
        color=BLUE,
        alpha=0.5,
        lw=0,
        label="Insomnia symptoms",
    )
    fit = smf.ols("YV ~ AGE + SEXN + EDU + cr(ISI_T1, df=3, constraints='center')", i).fit(
        cov_type="HC3"
    )
    g_ = pd.DataFrame(
        dict(
            ISI_T1=np.linspace(9, i["ISI_T1"].max(), 200),
            AGE=i["AGE"].mean(),
            SEXN=i["SEXN"].mean(),
            EDU=i["EDU"].mean(),
        )
    )
    pr = fit.get_prediction(g_).summary_frame(alpha=0.05)
    ax.axhline(h[col].mean(), color=MUTED, lw=1.2, ls=(0, (4, 3)), label="Mean of healthy controls")
    ax.fill_between(
        g_["ISI_T1"], pr["mean_ci_lower"], pr["mean_ci_upper"], color=INK, alpha=0.13, lw=0
    )
    ax.plot(g_["ISI_T1"], pr["mean"], color=INK, lw=2.2, label="Fitted curve (spline, 95% CI)")
    ax.axvspan(22.5, 26.6, color=MUTED, alpha=0.10, lw=0)
    lo_, hi_ = np.nanpercentile(df[col], [0.5, 99.5])
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
save_fig(fig, "FigureS2")
