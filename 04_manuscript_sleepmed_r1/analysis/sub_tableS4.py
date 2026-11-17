"""Table S4：前處理、肌電、清醒度與頻譜擬合的檢查。

彙整 R9–R11 及 task_r1g／k／p／s／t 的結果，只重新排表。"""

from _r1_submission import (
    AN,
    ci,
    ci_txt,
    csv_part,
    csv_sections,
    fp,
    gci,
    minus,
    nl,
    pd,
    sg,
    uz,
    write_table,
)

rd = lambda p: pd.read_csv(p, encoding="utf-8-sig")
T = AN / "S_R1_tables"

# (a)、(b)、(f)
part_a = csv_part(T / "R10.csv", "(a) Preprocessing summary by group, M (SD)")
h11, sec = csv_sections(T / "R11.csv")
part_b = ("(b) High-frequency power", h11, sec["High-frequency power (muscle-activity proxies)"])
part_f = ("(f) Vigilance markers, one value per participant", h11, sec["Vigilance markers"])

# (c)、(d)：表 R9 的兩半
r9 = pd.read_csv(T / "R9.csv", dtype=str, keep_default_na=False, encoding="utf-8-sig")
v9 = r9.values.tolist()
iA = next(i for i, r in enumerate(v9) if r[0].startswith("Panel A"))
iB = next(i for i, r in enumerate(v9) if r[0].startswith("Panel B"))
part_c = (
    "(c) Difference in periodic beta peak power from healthy controls after adding one covariate at a time, in control SD units [95% CI]",
    ["Covariate added", "Higher severity vs HC", "p", "Lower severity vs HC", "p"],
    [[minus(c) for c in r] for r in v9[iA + 1 : iB]],
)
part_d = (
    "(d) Restricted samples, Hedges' g relative to healthy controls [95% CI]",
    v9[iB + 1],
    v9[iB + 2 :],
)

# (e)
lp = rd(AN / "S_R1t_surface_laplacian" / "01_contrasts.csv")
EL = {"all 19": "All 19", "inner 9": "Inner 9", "central 3": "Central 3", "outer 10": "Outer 10"}
MS = ["Periodic beta peak power", "Aperiodic-adjusted power, 13–30 Hz", "Raw log power, 30–45 Hz"]
e = []
for ver in ("Average reference (main analysis)", "Surface Laplacian"):
    for k, klab in EL.items():
        for m in MS:
            s = lp[
                (lp["version"] == ver) & lp["electrodes"].str.startswith(k) & (lp["measure"] == m)
            ].set_index("contrast")
            hi, lo = s.loc["Higher vs HC"], s.loc["Lower vs HC"]
            e.append(
                [
                    ver.replace(" (main analysis)", ""),
                    klab,
                    m.replace("Raw log", "Log"),
                    f"{int(hi['n'])} / {int(lo['n'])} / {int(hi['n_hc'])}",
                    gci(hi["g"], hi["ci_low"], hi["ci_high"]),
                    fp(hi["p_mwu"]),
                    sg(hi["m2_diff_in_HC_SD"]),
                    fp(hi["m2_p"]),
                    gci(lo["g"], lo["ci_low"], lo["ci_high"]),
                    fp(lo["p_mwu"]),
                ]
            )
part_e = (
    "(e) Average reference and surface Laplacian",
    [
        "Reference",
        "Electrodes",
        "Measure",
        "n higher / lower / HC",
        "Higher severity, g [95% CI]",
        "p",
        "Higher severity, Model 2 difference",
        "p",
        "Lower severity, g [95% CI]",
        "p",
    ],
    e,
)

# (g)、(h)
S = AN / "S_R1s_vigilance_epochs"
fl = rd(S / "01_flagged_by_group.csv")
se = rd(S / "02_sensitivity.csv")
mm = lambda r, g: f"{uz(r[g + '_mean'])} ({uz(r[g + '_median'])})"
cov = lambda r: r["analysis"].startswith("+ covariate")  # 加共變項不改變未調整的 g，該格不重複列
part_g = (
    "(g) Vigilance markers computed epoch by epoch",
    [
        "Measure",
        "HC, M (median)",
        "Lower severity, M (median)",
        "g",
        "p",
        "Higher severity, M (median)",
        "g",
        "p",
        "Kruskal–Wallis p",
        "ρ with periodic beta",
        "p",
    ],
    [
        [
            r["measure"],
            mm(r, "HC"),
            mm(r, "Lower"),
            sg(r["Lower_g"]),
            fp(r["Lower_p"]),
            mm(r, "Higher"),
            sg(r["Higher_g"]),
            fp(r["Higher_p"]),
            fp(r["kw_p"]),
            sg(r["rho_with_beta"]),
            fp(r["p_with_beta"]),
        ]
        for _, r in fl.iterrows()
    ],
)
part_h = (
    "(h) Periodic beta peak power after controlling for or excluding flagged epochs",
    [
        "Analysis",
        "n higher / lower / HC",
        "Higher severity, g [95% CI]",
        "Lower severity, g [95% CI]",
        "Higher severity, Model 1 difference",
        "p",
    ],
    [
        [
            r["analysis"],
            f"{int(r['n_higher'])} / {int(r['n_lower'])} / {int(r['n_hc'])}",
            "–" if cov(r) else gci(r["g_higher"], r["g_higher_lo"], r["g_higher_hi"]),
            "–" if cov(r) else gci(r["g_lower"], r["g_lower_lo"], r["g_lower_hi"]),
            sg(r["adj_higher_in_HC_SD"]),
            fp(r["adj_higher_p"]),
        ]
        for _, r in se.iterrows()
    ],
)

# (i)、(j)
G = AN / "S_R1g_followup"
FE = {
    "prop_beta_peak": "Proportion of electrodes with a beta peak",
    "beta_cf_BRAIN": "Beta peak frequency, Hz",
    "beta_bw_BRAIN": "Beta peak bandwidth, Hz",
    "beta_pw_zero_imputed": "Peak power with zero entered where no peak was detected",
    "beta_pw_BRAIN": "Peak power (primary outcome)",
    "beta_pw_CENTRAL": "Peak power, central",
    "exponent_BRAIN": "Aperiodic exponent, whole scalp",
    "exponent_CENTRAL": "Aperiodic exponent, central",
    "offset_BRAIN": "Aperiodic offset, whole scalp",
    "offset_CENTRAL": "Aperiodic offset, central",
}
pk = rd(G / "02_beta_peak_detection.csv")
rf = rd(G / "04c_refit_frequency_range.csv")
gc = lambda g, c: f"{sg(g)} {ci_txt(c)}"
part_i = (
    "(i) Detection and parameters of the beta peak",
    [
        "Measure",
        "HC, M",
        "Lower severity, M",
        "Higher severity, M",
        "Kruskal–Wallis H",
        "p",
        "Lower severity, g [95% CI]",
        "p",
        "Higher severity, g [95% CI]",
        "p",
    ],
    [
        [
            FE[r["feature"]],
            f"{r['HC_mean']:.3f}",
            f"{r['Lower_mean']:.3f}",
            f"{r['Higher_mean']:.3f}",
            f"{r['KW_H']:.2f}",
            fp(r["KW_p"]),
            gc(r["g_Lower"], r["ci_Lower"]),
            fp(r["p_Lower"]),
            gc(r["g_Higher"], r["ci_Higher"]),
            fp(r["p_Higher"]),
        ]
        for _, r in pk.iterrows()
    ],
)
part_j = (
    "(j) Spectral parameterisation refitted with lower upper bounds",
    [
        "Fitted range",
        "Measure",
        "Mean fit R²",
        "r with the 1–45 Hz estimate",
        "Higher severity, g [95% CI]",
        "p",
        "Lower severity, g [95% CI]",
        "p",
    ],
    [
        [
            r["fit_range"],
            FE[r["feature"]]
            .replace("Peak power (primary outcome)", "Periodic beta peak power, whole scalp")
            .replace("Peak power, central", "Periodic beta peak power, central"),
            nl(r["mean_R2"], 3),
            "–" if r["fit_range"].startswith("1–45") else nl(r["r_with_original"], 2),
            gc(r["g_Higher"], r["ci_Higher"]),
            fp(r["p_Higher"]),
            gc(r["g_Lower"], r["ci_Lower"]),
            fp(r["p_Lower"]),
        ]
        for _, r in rf.iterrows()
    ],
)

# (k)、(l)、(m)
P = AN / "S_R1p_picard"
ag = rd(P / "02_agreement.csv")
mr = rd(P / "03_main_results.csv")
sa = rd(AN / "S_R1k_recruitment_flags" / "01_sensitivity.csv")
part_k = (
    "(k) Agreement between the two preprocessing runs",
    [
        "Measure",
        "n",
        "Pearson r",
        "Spearman ρ",
        "ICC(A,1)",
        "Mean difference",
        "95% limits of agreement",
    ],
    [
        [
            r["measure"].replace("whole-scalp", "whole scalp"),
            str(int(r["n"])),
            nl(r["pearson_r"], 3),
            nl(r["spearman_rho"], 3),
            nl(r["icc_a1"], 3),
            sg(r["mean_diff"], 4),
            ci(r["loa_low"], r["loa_high"], 3),
        ]
        for _, r in ag.iterrows()
    ],
)
part_l = (
    "(l) Main results under the two preprocessing runs",
    [
        "Preprocessing",
        "n (HC / lower / higher)",
        "Higher severity, g [95% CI]",
        "Lower severity, g [95% CI]",
        "Higher vs lower, g [95% CI]",
        "All with symptoms, g [95% CI]",
        "Higher severity, Model 2 b",
        "p",
    ],
    [
        [
            r["data"],
            f"{int(r['n'])} ({int(r['n_hc'])} / {int(r['n_lower'])} / {int(r['n_higher'])})",
            gci(r["g_higher"], r["g_higher_lo"], r["g_higher_hi"]),
            gci(r["g_lower"], r["g_lower_lo"], r["g_lower_hi"]),
            gci(r["g_hi_vs_lo"], r["g_hi_vs_lo_lo"], r["g_hi_vs_lo_hi"]),
            gci(r["g_pooled"], r["g_pooled_lo"], r["g_pooled_hi"]),
            uz(r["M2_b"], 3),
            fp(r["M2_p"]),
        ]
        for _, r in mr.iterrows()
    ],
)
part_m = (
    "(m) Main result before and after removal of four participants who met the exclusion criteria",
    [
        "Sample",
        "n (HC / lower / higher)",
        "Higher severity, g [95% CI]",
        "Lower severity, g [95% CI]",
        "Higher severity, Model 2 b",
        "p",
    ],
    [
        [
            r["sample"]
            .replace("Submitted sample", "Before removal")
            .replace("Revised sample (analytic sample)", "Analytic sample"),
            f"{int(r['N'])} ({int(r['n_hc'])} / {int(r['n_lower'])} / {int(r['n_higher'])})",
            gci(r["Higher_g"], r["Higher_lo"], r["Higher_hi"]),
            gci(r["Lower_g"], r["Lower_lo"], r["Lower_hi"]),
            uz(r["M2_higher_b"], 3),
            fp(r["M2_higher_p"]),
        ]
        for _, r in sa.iterrows()
    ],
)

write_table(
    "TableS4",
    "Supplementary Table S4",
    "Muscle activity, preprocessing, vigilance and spectral parameterisation: checks on the difference in periodic beta peak power.",
    [
        part_a,
        part_b,
        part_c,
        part_d,
        part_e,
        part_f,
        part_g,
        part_h,
        part_i,
        part_j,
        part_k,
        part_l,
        part_m,
    ],
    "g = Hedges' g relative to healthy controls (HC); p for g from two-sided Mann–Whitney U tests. Model 1 adjusts for age, sex and education; Model 2 additionally for BDI-II and BAI; adjusted differences are in control SD units unless b is given. "
    "In (a), H is the Kruskal–Wallis statistic across the three groups. In (c), the first row is Model 1 and each following row adds the covariate named. "
    "In (e), inner 9 = F3, Fz, F4, C3, Cz, C4, P3, Pz, P4; central 3 = C3, Cz, C4; outer 10 = the remaining electrodes, where Laplacian estimates from a 19-channel montage are least reliable; n is given for periodic beta peak power, which requires a detected peak. "
    "In (g), power was taken at posterior electrodes (O1, O2, P3, P4, Pz); an epoch was flagged when theta power was at least equal to alpha power, or when alpha power fell below half of the participant's own median; ρ = Spearman correlation across all participants. "
    "In (h), rows beginning with + add the covariate to Model 1 and leave the unadjusted g unchanged. "
    "In (k) and (l), the analysed data used PICARD and CleanLine and the earlier run used extended Infomax and a notch filter; comparisons include participants with an estimate in both runs. "
    "These analyses limit, but cannot exclude, a contribution of muscle activity or reduced vigilance.",
)
