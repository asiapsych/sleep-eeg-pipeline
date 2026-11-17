"""Figure 5：肌電、清醒度、前處理與 Laplacian 的敏感度分析。

A 畫調整後組間差，B 畫各子樣本對健康組的 Hedges' g，兩者量尺不同。"""

from _r1_submission import (
    AN,
    AQUA,
    BLUE,
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

e = pd.read_csv(AN / "S_R1c_sensitivity" / "03_emg_covariate.csv", encoding="utf-8-sig").rename(
    columns={"added_covariate": "setting"}
)
v = pd.read_csv(
    AN / "S_R1c_sensitivity" / "05b_vigilance_covariate.csv", encoding="utf-8-sig"
).rename(columns={"added_covariate": "setting"})
q = pd.read_csv(AN / "S_R1d_preproc_qc" / "02_qc_covariate.csv", encoding="utf-8-sig")
rs = pd.read_csv(AN / "S_R1d_preproc_qc" / "03_restricted_samples.csv", encoding="utf-8-sig")
lap = pd.read_csv(AN / "S_R1t_surface_laplacian" / "01_contrasts.csv", encoding="utf-8-sig")
ep = pd.read_csv(AN / "S_R1s_vigilance_epochs" / "02_sensitivity.csv").set_index("analysis")
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


# B 面板的列：(標籤, 高嚴重度的 g/lo/hi/n, 低嚴重度的 g/lo/hi/n)
def from_rs(s):
    out = []
    for con in ("K2 Higher vs HC", "K2 Lower vs HC"):
        r = rs[(rs["subsample"] == s) & (rs["contrast"] == con)].iloc[0]
        out.append((r["g"], r["ci_low"], r["ci_high"], int(r["n_a"])))
    return out


def from_lap():
    out = []
    for con in ("Higher vs HC", "Lower vs HC"):
        r = lap[
            (lap["version"] == "Surface Laplacian")
            & lap["electrodes"].str.startswith("central 3")
            & (lap["measure"] == "Periodic beta peak power")
            & (lap["contrast"] == con)
        ].iloc[0]
        out.append((r["g"], r["ci_low"], r["ci_high"], int(r["n"])))
    return out


def from_ep(name):
    r = ep.loc[name]
    return [
        (r["g_higher"], r["g_higher_lo"], r["g_higher_hi"], int(r["n_higher"])),
        (r["g_lower"], r["g_lower_lo"], r["g_lower_hi"], int(r["n_lower"])),
    ]


ROWS = [
    ("Full sample", from_rs("full sample")),
    ("No retained muscle component", from_rs("no retained muscle-dominant component")),
    ("Max muscle probability < .50", from_rs("max retained muscle probability < 0.50")),
    ("Max muscle probability < .20", from_rs("max retained muscle probability < 0.20")),
    ("No interpolated channel", from_rs("no interpolated channel")),
    ("Epoch rejection ≤ 10%", from_rs("epoch rejection rate <= 10%")),
    ("Surface Laplacian, central electrodes", from_lap()),
    ("Epochs with theta ≥ alpha removed", from_ep("Excluding epochs with theta ≥ alpha")),
    (
        "Participants with no such epoch",
        from_ep("Participants with no flagged epoch (criterion A)"),
    ),
]

fig, (a, b) = plt.subplots(1, 2, figsize=(13.4, 5.6), gridspec_kw=dict(wspace=0.82))
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
    xlabel="Adjusted difference, higher-severity group − HC (HC SD units, 95% CI)",
    title="Adding covariates",
)
a.legend(
    loc="upper left", bbox_to_anchor=(-0.55, -0.12), ncol=4, columnspacing=1.0, handletextpad=0.3
)
a.grid(axis="y", visible=False)
a.tick_params(axis="y", length=0)
panel_label(a, "A")
yy = np.arange(len(ROWS))[::-1]
for yi, (lab, vals) in zip(yy, ROWS):
    for (g, lo, hi, n), col, off in zip(vals, [ORANGE, BLUE], [0.16, -0.16]):
        b.plot([lo, hi], [yi + off] * 2, color=col, lw=2, solid_capstyle="round")
        b.scatter(g, yi + off, s=55, color=col, zorder=3, edgecolor=SURFACE, lw=1.5)
        b.text(1.72, yi + off, f"{sg(g)}  (n = {n})", va="center", fontsize=8.5)
b.scatter([], [], s=60, color=ORANGE, label="Higher-severity group")
b.scatter([], [], s=60, color=BLUE, label="Lower-severity group")
b.axvline(0, color=INK2, lw=1)
b.set_yticks(yy)
b.set_yticklabels([r[0] for r in ROWS])
b.set_xlim(-0.6, 1.7)
b.set(
    xlabel="Hedges' g vs healthy controls (95% CI)",
    title="Restricted samples, surface Laplacian and vigilance",
)
b.legend(loc="upper left", bbox_to_anchor=(0, -0.12), ncol=2)
b.grid(axis="y", visible=False)
b.tick_params(axis="y", length=0)
panel_label(b, "B")
save_fig(fig, "Figure5")
