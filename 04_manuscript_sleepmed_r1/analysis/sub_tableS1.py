"""Table S1：ISI 的向度、外部關聯、EEG 關聯與分群依據。

兩個保留成分的 EEG 模型讀 task_r1w/09，不能誤用同放三個成分的 task_r1f。
三群與四群的表留在 sub_response_tables。"""

from _r1_submission import (
    AN,
    ci,
    ci_txt,
    csv_part,
    csv_sections,
    fp,
    nl,
    pd,
    sg,
    write_table,
)

W = AN / "S_R1w_pc_interpretation"
ld = pd.read_csv(W / "01_loadings.csv", index_col=0, encoding="utf-8-sig")
rb = pd.read_csv(W / "06_robustness.csv", encoding="utf-8-sig").set_index("check")
ex = pd.read_csv(W / "03_external_correlates.csv", encoding="utf-8-sig")
adj = pd.read_csv(W / "08_external_adjusted_sex_age.csv", encoding="utf-8-sig").set_index(
    "variable"
)
PCS = ["PC1", "PC2", "PC3"]

a = [[i[2:]] + [sg(ld.loc[i, pc]) for pc in PCS] for i in ld.index if i[0].isdigit()]
a.append(["Eigenvalue"] + [f"{rb.loc[f'{pc} eigenvalue', 'value']:.2f}" for pc in PCS])
a.append(
    ["95th percentile of eigenvalues in permuted data"]
    + [rb.loc[f"{pc} eigenvalue", "reference"].split("= ")[1] for pc in PCS]
)
a.append(["Proportion of variance"] + [nl(ld.loc["explained variance", pc], 3) for pc in PCS])
a.append(["Correlation with ISI total"] + [nl(ld.loc["r with ISI total", pc], 2) for pc in PCS])
a.append(
    ["Bootstrap congruence, median (5th percentile)"]
    + [
        f"{nl(rb.loc[f'{pc} bootstrap congruence (median)', 'value'], 2)} ({rb.loc[f'{pc} bootstrap congruence (median)', 'reference'].split('= ')[1].replace('0.', '.')})"
        for pc in PCS
    ]
)

b = []
for _, r in ex.iterrows():
    v = r["variable"]
    extra = (
        [
            f"{sg(adj.loc[v, 'beta_adjusted_sex_age'])} {ci(adj.loc[v, 'ci_low'], adj.loc[v, 'ci_high'])}",
            fp(adj.loc[v, "p_adjusted"]),
        ]
        if v in adj.index
        else ["–", "–"]
    )
    b.append(
        [v, str(int(r["n"])), sg(r["rho_PC1"]), fp(r["q_PC1"]), sg(r["rho_PC2"]), fp(r["q_PC2"])]
        + extra
    )

NAMES = {
    "beta_pw_BRAIN": "Periodic beta peak power, whole scalp",
    "beta_pw_CENTRAL": "Periodic beta peak power, central",
    "beta_pw_MIDLINE": "Periodic beta peak power, midline",
    "Beta_rel_BRAIN": "Relative beta power, whole scalp",
    "Beta_rel_CENTRAL": "Relative beta power, central",
    "exponent_BRAIN": "Aperiodic exponent, whole scalp",
    "exponent_CENTRAL": "Aperiodic exponent, central",
    "offset_BRAIN": "Aperiodic offset, whole scalp",
    "offset_CENTRAL": "Aperiodic offset, central",
    "alpha_pw_BRAIN": "Periodic alpha peak power",
    "alpha_cf_BRAIN": "Alpha peak frequency",
    "theta_pw_BRAIN": "Periodic theta peak power",
    "Theta_rel_BRAIN": "Relative theta power",
    "Alpha_rel_BRAIN": "Relative alpha power",
    "Delta_rel_BRAIN": "Relative delta power",
    "Gamma_rel_BRAIN": "Relative 30–45 Hz power",
}
t2c = pd.read_csv(W / "09_two_component_by_eeg.csv", encoding="utf-8-sig")
c = [
    [NAMES.get(r["biomarker"], r["biomarker"]), str(int(r["n"]))]
    + sum(
        [
            [f"{sg(r[f'{pc}_beta'])} {ci(r[f'{pc}_ci_low'], r[f'{pc}_ci_high'])}", fp(r[f"{pc}_q"])]
            for pc in ("PC1", "PC2")
        ],
        [],
    )
    for _, r in t2c.iterrows()
]

_ari = (
    pd.read_csv(AN / "S_R1e_clustering_robustness" / "03_ari_summary.csv", encoding="utf-8-sig")
    .set_index("set")
    .loc["all 36 method pairs"]
)
_h8, _s8 = csv_sections(
    AN / "S_R1_tables" / "R8.csv"
)  # 表 R8 的 K = 3 那一半只給審稿人看（sub_response_tables.py）

write_table(
    "TableS1",
    "Supplementary Table S1",
    "Structure of insomnia symptoms: principal components of the ISI items, their external correlates and EEG associations, and the basis for the severity groups.",
    [
        (
            "(a) Principal components of the seven ISI items (participants with insomnia symptoms, n = 162)",
            ["Item or quantity", "Component 1", "Component 2", "Component 3"],
            a,
        ),
        (
            "(b) Correlates of the two retained dimensions among variables outside the ISI",
            [
                "Variable",
                "n",
                "Dimension 1 ρ",
                "q",
                "Dimension 2 ρ",
                "q",
                "Dimension 2, adjusted for sex and age: β [95% CI]",
                "p",
            ],
            b,
        ),
        (
            "(c) The two retained dimensions and 16 EEG measures (standardised coefficients)",
            ["EEG measure", "n", "Dimension 1 β [95% CI]", "q", "Dimension 2 β [95% CI]", "q"],
            c,
        ),
        csv_part(
            AN / "S_R1_tables" / "R1.csv",
            "(d) Validity indices and bootstrap stability of K-means solutions with two to six clusters",
        ),
        (
            "(e) Agreement among nine clustering methods for two clusters",
            _h8[:3],
            [r[:3] for r in _s8["K = 2"]],
        ),
    ],  # 各方法下的效果量不列：不拿腦波結果當分群的依據
    "Component signs are arbitrary; component 1 is oriented towards a higher ISI total and component 2 towards higher scores on difficulty staying asleep and early-morning awakening. "
    "Parallel analysis retained components whose eigenvalue exceeded the 95th percentile of eigenvalues from 2,000 data sets with each item permuted independently; components 1 and 2 were retained. "
    "In (b), ρ = Spearman correlation; q = Benjamini–Hochberg-adjusted p across the 13 variables within each dimension; adjusted coefficients are standardised and were estimated for the variables listed. "
    "PSQI component scores are coded so that higher values indicate worse sleep. In (c), the two retained dimensions were entered together with age, sex and education, as in the main analysis (Table 2); q is adjusted across the 16 measures within each dimension. "
    f"In (e), the mean ARI over all 36 pairs of methods was {nl(_ari['mean_ARI_K2'], 2)} {ci_txt(_ari['K2_ci']).replace('+', '').replace('0.', '.')}. "
    "ARI = adjusted Rand index; DBAS = Dysfunctional Beliefs and Attitudes about Sleep scale.",
)
