"""Figure 1：ISI 負荷量、第二向度的外部關聯、兩個向度對 beta 的係數。

讀 task_r1w 的結果；外部關聯面板只畫 BH 校正後顯著的項目。"""

from _r1_submission import (
    AN,
    GRAY,
    INK2,
    ORANGE,
    SURFACE,
    np,
    panel_label,
    pd,
    plt,
    save_fig,
    sg,
)

W = AN / "S_R1w_pc_interpretation"
ld = pd.read_csv(W / "01_loadings.csv", index_col=0, encoding="utf-8-sig")
ex = pd.read_csv(W / "03_external_correlates.csv", encoding="utf-8-sig")
pvb = pd.read_csv(W / "07_two_component_model.csv", encoding="utf-8-sig")
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

e2 = ex[ex["q_PC2"] < 0.05].sort_values("rho_PC2").reset_index(drop=True)
z = np.arctanh(e2["rho_PC2"])
se = 1.06 / np.sqrt(e2["n"] - 3)  # Spearman 的 Fisher z 近似標準誤
lo_, hi_ = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
for i, r in e2.iterrows():
    b.plot([lo_[i], hi_[i]], [i, i], color=ORANGE, lw=2.2, solid_capstyle="round")
    b.scatter(r["rho_PC2"], i, s=64, color=ORANGE, zorder=3, edgecolor=SURFACE, lw=1.5)
    if r["rho_PC2"] >= 0:
        b.text(hi_[i] + 0.03, i, sg(r["rho_PC2"]), va="center", fontsize=9.5)
    else:
        b.text(lo_[i] - 0.03, i, sg(r["rho_PC2"]), va="center", ha="right", fontsize=9.5)
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
    c.text(r.ci_high + 0.03, 1 - i, sg(r.std_beta), va="center", fontsize=9.5)
c.axvline(0, color=INK2, lw=1)
c.set_yticks([1, 0])
c.set_yticklabels(["Dimension 1", "Dimension 2"])
c.set_ylim(-0.6, 1.6)
c.set_xlim(-0.3, 0.62)
c.set(xlabel="Standardised coefficient (95% CI)", title="Association with periodic beta")
c.grid(axis="y", visible=False)
c.tick_params(axis="y", length=0)
panel_label(c, "C")
save_fig(fig, "Figure1")
