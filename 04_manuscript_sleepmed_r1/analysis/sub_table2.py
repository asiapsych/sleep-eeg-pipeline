"""Table 2：症狀向度、嚴重度組及其他 beta／aperiodic 指標。

主要結果取自 R4；其他指標的未調整 g 在這裡計算，調整後結果讀 task_r1i。"""

from _r1_common import hedges_g, load_all, mwu_p
from _r1_submission import (
    AN,
    ci_txt,
    csv_part,
    fp,
    gci,
    np,
    pd,
    sg,
    write_table,
)

df = load_all()
df["logBeta_abs_BRAIN"] = np.log10(df["Beta_abs_BRAIN"])
hc = df[df["k2"] == "HC"]
mood = pd.read_csv(AN / "S_R1i_mood_adjustment" / "01_mood_adjustment.csv", encoding="utf-8-sig")
VARS = [
    ("Beta_rel_BRAIN", "Relative beta power (13–30 Hz), whole scalp"),
    ("logBeta_abs_BRAIN", "Absolute beta power (log), whole scalp"),
    ("exponent_BRAIN", "Aperiodic exponent, whole scalp"),
    ("exponent_CENTRAL", "Aperiodic exponent, central"),
    ("offset_BRAIN", "Aperiodic offset, whole scalp"),
    ("offset_CENTRAL", "Aperiodic offset, central"),
]
rows = []
for var, lab in VARS:
    for grp, k in [("Higher severity vs HC", "Higher"), ("Lower severity vs HC", "Lower")]:
        a = df.loc[df["k2"] == k, var]
        g, lo_, hi_, _, _ = hedges_g(a, hc[var])
        m = mood[(mood["variable"] == var) & (mood["group"] == grp)].iloc[0]
        rows.append(
            [
                lab if k == "Higher" else "",
                grp,
                gci(g, lo_, hi_),
                fp(mwu_p(a, hc[var])),
                f"{sg(m['M1_diff_in_HC_SD'])} {ci_txt(m['M1_ci'])}",
                fp(m["M1_p"]),
                f"{sg(m['M2_diff_in_HC_SD'])} {ci_txt(m['M2_ci'])}",
                fp(m["M2_p"]),
            ]
        )
import re

_sub, _h, _r = csv_part(AN / "S_R1_tables" / "R4.csv", "(a) Whole-scalp periodic beta peak power")
part_a = (
    _sub,
    _h,
    [[re.sub(r"(?<=\d)a$", "ᵃ", c) for c in r] for r in _r],
)  # 報告表的註記 a 改成上標
write_table(
    "Table2",
    "Table 2",
    "Periodic beta peak power in relation to the two symptom dimensions and the two severity groups, with conventional beta power and aperiodic parameters.",
    [
        part_a,
        (
            "(b) Conventional beta power and aperiodic parameters",
            [
                "Measure",
                "Contrast",
                "Unadjusted g [95% CI]",
                "p",
                "Model 1 difference [95% CI]",
                "p",
                "Model 2 difference [95% CI]",
                "p",
            ],
            rows,
        ),
    ],
    "Symptom dimensions are the two retained principal components of the seven ISI items, entered together and adjusted for age, sex and education (standardised coefficients). "
    "Unadjusted contrasts are Hedges' g with p values from two-sided Mann–Whitney U tests. In (a), adjusted estimates are unstandardised coefficients (b) for group membership; "
    "in (b), adjusted differences are in units of the control SD. Model 1 adjusts for age, sex and education; Model 2 additionally for BDI-II and BAI. "
    "The equivalence test is a two one-sided-tests procedure with bounds of ±0.37 pooled SD, the meta-analytic effect size for waking beta power in insomnia; p < .05 indicates that the difference lies within the bounds. "
    "Bayes factors are Jeffreys–Zellner–Siow Bayes factors for the two-sample comparison; BF₀₁ quantifies evidence for no difference and BF₁₀ evidence for a difference. "
    "ᵃRange of BF₀₁ for prior widths of 0.5 to 1.0. HC = healthy controls.",
)
