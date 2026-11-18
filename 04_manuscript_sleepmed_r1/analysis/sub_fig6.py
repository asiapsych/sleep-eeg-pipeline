"""Figure 6：提高 ISI 納入下限時，合併效果量怎麼變。

橘線是納入者，藍線是未納入的失眠症狀者；虛線為文獻的 g = 0.37。
數值讀 task_r1u 的 03_pooling_cost.csv。"""

from _r1_submission import (
    AN,
    AXIS,
    BLUE,
    INK2,
    MUTED,
    ORANGE,
    pd,
    plt,
    save_fig,
    sg,
)

cost = pd.read_csv(AN / "S_R1u_shape_robust" / "03_pooling_cost.csv", encoding="utf-8-sig")
cost["c"] = cost["cutoff"].str.extract(r"(\d+)").astype(int)
BOUND = 0.37
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
    (15, "At least moderate\nsymptom severity (ISI ≥ 15)"),
]:
    r = cost[cost["c"] == c_].iloc[0]
    a.annotate(
        f"{lab}\ng = {sg(r['g'])}",
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
)
a.legend(loc="upper left", bbox_to_anchor=(0, -0.16), ncol=2, columnspacing=1.2, handletextpad=0.4)
save_fig(fig, "Figure6")
