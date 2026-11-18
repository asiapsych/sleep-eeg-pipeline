"""Table S2：beta 與嚴重度的關係型態及切點敏感度。

分數帶、函數形式、切點掃描與 HC3 對照來自 task_r1a／g／u。
樣條自由度按實際項數標示，和正文一致。"""

import re
from _r1_submission import (
    AN,
    ci,
    csv_part,
    fp,
    gci,
    minus,
    nl,
    pd,
    sg,
    uz,
    write_table,
)

U, G, A = AN / "S_R1u_shape_robust", AN / "S_R1g_followup", AN / "S_R1a_severity_form"
rd = lambda p: pd.read_csv(p, encoding="utf-8-sig")
SC = {"insomnia": "Insomnia symptoms", "full": "Full sample"}
NAMES = {
    "beta_pw_BRAIN": "Periodic beta peak power, whole scalp",
    "beta_pw_CENTRAL": "Periodic beta peak power, central",
    "Beta_rel_BRAIN": "Relative beta power, whole scalp",
    "Beta_rel_CENTRAL": "Relative beta power, central",
    "logBeta_abs_BRAIN": "Absolute beta power (log), whole scalp",
    "exponent_BRAIN": "Aperiodic exponent, whole scalp",
    "exponent_CENTRAL": "Aperiodic exponent, central",
    "offset_BRAIN": "Aperiodic offset, whole scalp",
    "offset_CENTRAL": "Aperiodic offset, central",
}

# (a)
sub, h5, r5 = csv_part(
    AN / "S_R1_tables" / "R5.csv",
    "(a) Periodic beta peak power relative to healthy controls by ISI band",
)
part_a = (sub, h5[:-1], [r[:-1] for r in r5])  # 最後一欄是初投三群的人數，不放

# (b)
up = rd(U / "07_upper_range_other_measures.csv")
part_b = (
    "(b) Upper range of severity: primary and secondary beta measures, Hedges' g relative to healthy controls [95% CI]",
    [
        "ISI band",
        "n",
        "Periodic beta peak power",
        "Aperiodic-adjusted beta power, 13–30 Hz",
        "Relative beta power",
    ],
    [
        [r["ISI_range"].replace(">=", "≥ ").replace("-", "–"), str(int(r["n"]))]
        + [
            gci(r[f"g_{k}"], r[f"ci_low_{k}"], r[f"ci_high_{k}"])
            for k in ("peak_power", "aperiodic_adjusted_13_30", "relative_beta")
        ]
        for _, r in up.iterrows()
    ],
)

# (c)
mc = rd(A / "04_model_comparison.csv")
th = rd(A / "05_threshold.csv")
c_hat = int(
    float(th.loc[th["quantity"].str.startswith("data-selected threshold"), "value"].iloc[0])
)
mi = mc[mc["scope"] == "insomnia only"].set_index("model")
lin = mi.loc["linear ISI total", "AIC"]
ML = [
    ("covariates only", "Covariates only"),
    ("linear ISI total", "Linear ISI total"),
    ("quadratic ISI total", "Quadratic ISI total"),
    ("spline df=3", "Natural cubic spline, 2 df"),
    ("spline df=4", "Natural cubic spline, 3 df (curve in Figure 2)"),
    ("spline df=5", "Natural cubic spline, 4 df"),
    ("single threshold ISI>=15 (conventional)", "Step at ISI ≥ 15 (conventional category)"),
    (
        f"single threshold ISI>={c_hat} (data-selected)",
        f"Step at ISI ≥ {c_hat} (location selected from the data)ᵃ",
    ),
    ("K=2 higher-severity cluster", "Two severity groups (K-means)"),
    ("PC1", "First symptom dimension (component 1)"),
]  # 三群的兩列在 sub_response_tables.py
assert int(mi["n"].iloc[0]) == 162
part_c = (
    "(c) Functional forms relating ISI total to periodic beta peak power (participants with insomnia symptoms, n = 162)",
    ["Severity term", "k", "R²", "Adjusted R²", "AIC", "BIC", "ΔAIC"],
    [
        [
            lab,
            str(int(mi.loc[m, "k_severity"])),
            nl(mi.loc[m, "R2"], 3),
            nl(mi.loc[m, "adj_R2"], 3),
            uz(mi.loc[m, "AIC"], 1),
            uz(mi.loc[m, "BIC"], 1),
            "Reference" if m == "linear ISI total" else sg(mi.loc[m, "AIC"] - lin, 1),
        ]
        for m, lab in ML
    ],
)

# (d)、(e)
st = rd(U / "02_shape_tests.csv")
bm = rd(U / "08_shape_by_measure.csv")
d = []
for _, r in st.iterrows():
    spl = r["test"].startswith("spline")
    cand = re.search(r"candidates (\d+)-?–?(\d+)", r["test"].replace("–", "-")) if not spl else None
    d.append(
        [
            SC[r["scope"]],
            "Model 2" if r["adjustment"] == "M2" else "Model 1",
            str(int(r["n"])),
            "Spline (3 df) vs linear"
            if spl
            else f"Largest step vs linear, candidates ISI {cand.group(1)}–{cand.group(2)}",
            "–" if spl else str(int(r["selected_cutoff"])),
            fp(r["p_hc3"]) if spl else fp(r["p_selection_corrected"]),
        ]
    )
part_d = (
    "(d) Tests against a linear association (periodic beta peak power)",
    ["Sample", "Adjustment", "n", "Test", "Location of the largest step", "p"],
    d,
)
part_e = (
    "(e) Spline versus linear association by beta measure (participants with insomnia symptoms)",
    ["Measure", "n", "p", "p excluding ISI ≥ 23"],
    [
        [
            r["measure"][0].upper() + r["measure"][1:].replace("13-30", "13–30"),
            str(int(r["n"])),
            fp(r["spline_vs_linear_p_hc3"]),
            fp(r["excluding_isi23plus_p_hc3"]),
        ]
        for _, r in bm.iterrows()
    ],
)

# (f)
tv = lambda key: str(th.loc[th["quantity"].str.startswith(key)].iloc[0]["value"])
tn = lambda key: str(th.loc[th["quantity"].str.startswith(key)].iloc[0]["note"])
n_ab, n_be = re.search(r"n above = (\d+), n below = (\d+)", tn("data-selected threshold")).groups()
b_, se_, _p = [x.strip() for x in tv("step coefficient").split(",")]
perm = float(re.search(r"= ([\d.]+)$", tn("max-scanned")).group(1))
share = {int(k): float(v) for k, v in re.findall(r"(\d+): ([\d.]+)%", tn("bootstrap distribution"))}
med, lo_, hi_ = re.search(
    r"median (\d+); 2\.5–97\.5 pct (\d+)–(\d+)", tv("bootstrap distribution")
).groups()
part_f = (
    "(f) Scan for the largest step (candidates ISI 11–22, at least 10 participants on each side)",
    ["Quantity", "Value", "Note"],
    [
        [
            "Location of the largest step",
            f"ISI ≥ {tv('data-selected threshold')}",
            f"{n_ab} participants at or above, {n_be} below",
        ],
        [
            "Step coefficient in a step-only model, b (SE)",
            f"{b_} ({se_})",
            "p not corrected for the selection of the location",
        ],
        [
            "Largest step F across candidates",
            tv("max-scanned"),
            f"Permutation p {'< .001' if perm < 0.001 else '= ' + nl(perm, 3)} against no association (5,000 permutations, scan repeated in each)",
        ],
        [
            "Step added to the linear model",
            minus(tn("step added")),
            f"p = {tv('step added')}, not corrected for selection",
        ],
        [
            "Linear term added to the step model",
            minus(tn("linear added")),
            f"p = {tv('linear added')}",
        ],
        [
            "Location selected in 5,000 bootstrap resamples",
            f"Median {med}; 95% interval {lo_}–{hi_}",
            f"ISI {c_hat} selected in {share[c_hat]:.1f}% of resamples; 17 to 20 in {sum(share[k] for k in (17, 18, 19, 20)):.1f}%",
        ],
        [
            "Full sample, candidates ISI 5–22",
            f"ISI ≥ {tv('full-sample data-selected')}",
            tn("full-sample data-selected"),
        ],
    ],
)

# (g)
SEG = {
    "all insomnia": "All participants with insomnia symptoms",
    "below threshold": "Below the largest step",
    "at/above threshold": "At or above the largest step",
    "K=2 lower-severity cluster": "Lower-severity group",
    "K=2 higher-severity cluster": "Higher-severity group",
}


def seg(s):
    for k, v in SEG.items():
        if s.startswith(k):
            rng = re.search(r"\((.*)\)", s)
            return v + (f" ({rng.group(1).replace('>=', '≥ ').replace('-', '–')})" if rng else "")
    raise ValueError(s)


sl = rd(G / "01a_segment_slopes.csv")
_lin = rs_ = (
    pd.read_csv(U / "01_robust_se.csv", encoding="utf-8-sig")
    .set_index("result")
    .loc["ISI total (per point), M1"]
)
assert abs(_lin["b"] - sl.loc[0, "slope_per_ISI_point"]) < 1e-4 and sl.loc[0, "segment"].startswith(
    "all insomnia"
)  # 同一個模型
sl.loc[0, "p"] = _lin["p_ols"]
part_g = (
    "(g) Slope of periodic beta peak power on ISI total within segments",
    [
        "Segment",
        "n",
        "Slope per ISI point [95% CI]",
        "p",
        "Slope per 5 points, in control SD [95% CI]",
    ],
    [
        [
            seg(r["segment"]),
            str(int(r["n"])),
            f"{uz(r['slope_per_ISI_point'], 4)} [{uz(r['ci_low'], 4)}, {uz(r['ci_high'], 4)}]",
            fp(r["p"]),
            f"{sg(r['slope_in_HC_SD_per_5_points'])} {ci(r['ci_low_SD5'], r['ci_high_SD5'])}",
        ]
        for _, r in sl.iterrows()
    ],
)

# (h)
rs = rd(U / "01_robust_se.csv")
lab_h = lambda s: (
    s.replace(">=", "≥")
    .replace(", M1", ", Model 1")
    .replace(", M2", ", Model 2")
    .replace("Higher vs HC", "Higher severity vs HC")
    .replace("Lower vs HC", "Lower severity vs HC")
)
part_h = (
    "(h) Main regression results with ordinary and heteroscedasticity-consistent (HC3) standard errors",
    ["Result", "Sample", "n", "b", "p, ordinary", "95% CI, HC3", "p, HC3"],
    [
        [
            lab_h(r["result"]),
            SC[r["scope"]],
            str(int(r["n"])),
            uz(r["b"], 4),
            fp(r["p_ols"]),
            f"[{uz(r['ci_low_hc3'], 4)}, {uz(r['ci_high_hc3'], 4)}]",
            fp(r["p_hc3"]),
        ]
        for _, r in rs.iterrows()
    ],
)

# (i)
rc = rd(G / "05b_recluster_by_cutoff.csv")
part_i = (
    "(i) Severity groups re-derived after raising the lower ISI bound for inclusion",
    ["Participants included", "Group", "ISI M", "n", "g vs HC [95% CI]", "p"],
    [
        [
            r["insomnia_definition"].replace(">=", "≥"),
            r["cluster"] + " severity",
            f"{r['ISI_mean']:.1f}",
            str(int(r["n_a"])),
            gci(r["g"], r["ci_low"], r["ci_high"]),
            fp(r["p"]),
        ]
        for _, r in rc.iterrows()
    ],
)

# (j)
sm = rd(G / "06c_shape_by_measure.csv")
part_j = (
    "(j) Location of the largest step by measure (full sample, candidates ISI 5–22)",
    ["Measure", "n", "Location of the largest step", "Direction at or above", "Permutation p"],
    [
        [
            NAMES[r["feature"]],
            str(int(r["n"])),
            f"ISI ≥ {int(r['selected_threshold'])}",
            "Higher" if r["direction"] == "+" else "Lower",
            fp(r["perm_p"]),
        ]
        for _, r in sm.iterrows()
    ],
)

write_table(
    "TableS2",
    "Supplementary Table S2",
    "Periodic beta across the range of insomnia severity: band-wise effect sizes, tests of the form of the association and sensitivity analyses.",
    [part_a, part_b, part_c, part_d, part_e, part_f, part_g, part_h, part_i, part_j],
    "g = Hedges' g relative to healthy controls (HC); p in (a) from two-sided Mann–Whitney U tests. Only eight participants had ISI totals of 23 or above. "
    "In (c), models are ordinary least squares with age, sex and education; k = number of severity parameters; ΔAIC is relative to the linear model. "
    "ᵃOne additional parameter is counted for the selected location. "
    "In (d) and (e), the spline is a natural cubic spline with 3 df and its test is the nested comparison with the linear model (2 additional df; HC3 standard errors). "
    "The largest-step test repeats the scan of candidate locations in each of 4,999 wild-bootstrap resamples generated under the linear model, so that its p value accounts for the selection. "
    "The permutation tests in (f) and (j) use the null hypothesis of no association and are therefore not tests against linearity. "
    "Model 1 adjusts for age, sex and education; Model 2 additionally for BDI-II and BAI.",
)
