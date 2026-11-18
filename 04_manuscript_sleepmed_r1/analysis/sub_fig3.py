"""Figure 3：去背景頻譜、逐頻率效果量與顯著區段。

區段直接讀 task_r1h 的 cluster permutation 結果，不在畫圖時重做檢定。"""

from _r1_common import COV_COLS, load_all
from _r1_spectra import load_spectra
from _r1_submission import (
    AN,
    BLUE,
    GRAY,
    GRID,
    INK,
    INK2,
    ORANGE,
    nl,
    np,
    panel_label,
    pd,
    plt,
    save_fig,
)

df = load_all().dropna(subset=COV_COLS).reset_index(drop=True)
freq, chn, adj, raw = load_spectra(df["subject"])
flat = adj.mean(1)  # 人 × 頻率（全頭皮平均）
C = np.column_stack(
    [
        np.ones(len(df)),
        df["AGE"].to_numpy(float),
        df["SEX"].to_numpy(float),
        df["EDU"].to_numpy(float),
    ]
)
R = flat - C @ np.linalg.lstsq(C, flat, rcond=None)[0] + flat.mean(0)
H = pd.read_csv(AN / "S_R1h_beta_characterisation" / "02_g_by_frequency.csv", encoding="utf-8-sig")
cl = pd.read_csv(
    AN / "S_R1h_beta_characterisation" / "01_frequency_clusters.csv", encoding="utf-8-sig"
)
MEAS, RAWM = "aperiodic-adjusted, whole-scalp", "raw log power, whole-scalp"
G = lambda m, c: H[f"g | {m} | {c}"].to_numpy()
F = H["freq_Hz"].to_numpy()

fig, (a, b, c) = plt.subplots(1, 3, figsize=(16.5, 4.5), gridspec_kw=dict(wspace=0.27))
for lab, mk, col in [
    ("Healthy controls", (df["k2"] == "HC").to_numpy(), GRAY),
    ("Lower-severity group", (df["k2"] == "Lower").to_numpy(), BLUE),
    ("Higher-severity group", (df["k2"] == "Higher").to_numpy(), ORANGE),
]:
    mu, se = R[mk].mean(0), R[mk].std(0, ddof=1) / np.sqrt(mk.sum())
    a.fill_between(freq, mu - se, mu + se, color=col, alpha=0.18, lw=0)
    a.plot(freq, mu, color=col, lw=2, label=f"{lab} (n = {mk.sum()})")
a.axvspan(13, 30, color=GRID, alpha=0.45, lw=0)
a.set(
    xlabel="Frequency (Hz)",
    ylabel="Aperiodic-adjusted power (log10)",
    xlim=(1, 45),
    title="Spectra after removing the aperiodic component",
)
a.text(
    21.5,
    a.get_ylim()[0] + 0.01,
    "beta 13–30 Hz",
    ha="center",
    va="bottom",
    fontsize=8.5,
    color=INK2,
)
a.legend(loc="upper right")
panel_label(a, "A")

for lab, col, nm in [
    ("Higher vs HC", ORANGE, "Higher-severity group vs HC"),
    ("Lower vs HC", BLUE, "Lower-severity group vs HC"),
]:
    b.plot(F, G(MEAS, lab), color=col, lw=2, label=nm)
    for r in cl[
        (cl["measure"] == MEAS) & (cl["contrast"] == lab) & (cl["p_cluster"] < 0.05)
    ].itertuples():
        b.axvspan(r.from_Hz, r.to_Hz, color=col, alpha=0.14, lw=0)
        b.text(
            (r.from_Hz + r.to_Hz) / 2,
            1.02,
            f"{r.from_Hz:.0f}–{r.to_Hz:.0f} Hz\n"
            + ("p < .001" if r.p_cluster < 0.001 else "p = " + nl(r.p_cluster, 3)),
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
    title="Group difference at each frequency",
)
b.legend(loc="lower right")
panel_label(b, "B")

for m_, ls, nm in [(MEAS, "-", "Aperiodic component removed"), (RAWM, (0, (4, 2)), "Not removed")]:
    c.plot(F, G(m_, "Higher vs HC"), color=ORANGE, lw=2, ls=ls, label=nm)
c.axhline(0, color=INK2, lw=1)
c.axvspan(13, 30, color=GRID, alpha=0.45, lw=0)
c.set(
    xlabel="Frequency (Hz)",
    ylabel="Hedges' g, higher-severity group vs HC",
    xlim=(1, 45),
    ylim=(-0.75, 1.1),
    title="With and without removing the aperiodic component",
)
c.legend(loc="lower right")
panel_label(c, "C")
save_fig(fig, "Figure3")
