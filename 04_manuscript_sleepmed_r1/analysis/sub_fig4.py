"""Figure 4：頻率 × ISI 分數帶，以及 19 個電極的 beta 效果量。

頭皮圖的粗框表示電極內 FDR 校正後顯著；數值來自 task_r1h。"""

from matplotlib.patches import Circle, Patch

from _r1_submission import (
    AN,
    AXIS,
    DIV,
    INK,
    INK2,
    SURFACE,
    panel_label,
    pd,
    plt,
    save_fig,
    sg,
)

hm = pd.read_csv(
    AN / "S_R1h_beta_characterisation" / "04_frequency_by_severity.csv",
    index_col=0,
    encoding="utf-8-sig",
)
tp = pd.read_csv(AN / "S_R1h_beta_characterisation" / "05_topography.csv", encoding="utf-8-sig")
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

fig = plt.figure(figsize=(16.0, 4.9))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1], wspace=0.12)
a = fig.add_subplot(gs[0])
im = a.imshow(
    hm.to_numpy(float),
    aspect="auto",
    cmap=DIV,
    vmin=-1,
    vmax=1,
    origin="lower",
    extent=[1, 45, -0.5, len(hm) - 0.5],
)
a.set_yticks(range(len(hm)))
a.set_yticklabels([str(i).replace(">=", "≥ ").replace("-", "–") for i in hm.index])
for xx in (13, 30):
    a.axvline(xx, color=INK2, lw=0.8, ls=(0, (3, 3)))
a.set(
    xlabel="Frequency (Hz)",
    ylabel="ISI total score band",
    title="Difference from healthy controls by frequency and severity",
)
a.grid(False)
cb = fig.colorbar(im, ax=a, fraction=0.035, pad=0.02)
cb.outline.set_visible(False)
cb.set_label("Hedges' g", fontsize=9)
cb.ax.tick_params(length=0, labelsize=8)
panel_label(a, "A")
axes = [fig.add_subplot(gs[1]), fig.add_subplot(gs[2])]
for ax, lab, ttl in zip(
    axes,
    ["Higher vs HC", "Lower vs HC"],
    ["Higher-severity group vs HC", "Lower-severity group vs HC"],
):
    t = tp[tp["contrast"] == lab].set_index("channel")
    ax.add_patch(Circle((0, 0), 1.0, fill=False, ec=AXIS, lw=1.5))
    ax.plot([-0.08, 0, 0.08], [0.995, 1.09, 0.995], color=AXIS, lw=1.5)
    for cname, (x, y) in POS.items():
        sig = t.loc[cname, "q_band"] < 0.05
        ax.scatter(
            x,
            y,
            s=760,
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
            cname.capitalize(),
            ha="center",
            va="bottom",
            fontsize=7,
            color=INK2,
            zorder=4,
        )
        ax.text(
            x,
            y - 0.015,
            sg(t.loc[cname, "g_band"]),
            ha="center",
            va="top",
            fontsize=8,
            fontweight="bold" if sig else "normal",
            zorder=4,
        )
    ax.set(xlim=(-1.15, 1.15), ylim=(-1.22, 1.18))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(ttl, loc="center")
panel_label(axes[0], "B", x=0.02, y=1.0)
fig.legend(
    handles=[
        Patch(fc=SURFACE, ec=INK, lw=1.8, label="FDR q < .05 across 19 electrodes"),
        Patch(fc=SURFACE, ec=AXIS, lw=0.8, label="not significant"),
    ],
    loc="lower center",
    ncol=2,
    bbox_to_anchor=(0.70, -0.02),
)
save_fig(fig, "Figure4")
