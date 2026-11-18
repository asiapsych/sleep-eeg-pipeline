"""Table 1：健康組、全部失眠症狀者、低與高嚴重度組的特徵。

描述統計直接計算。三組用 Kruskal–Wallis（13 項內 BH 校正），
健康組對全部失眠症狀者用未校正的 Mann–Whitney U；性別用卡方檢定。"""

from statsmodels.stats.multitest import multipletests

from _r1_common import load_all, mwu_p
from _r1_submission import (
    fp,
    kw3,
    msd,
    nl,
    pd,
    stats,
    write_table,
)

df = load_all()
hc = df[df["k2"] == "HC"]
ins = df[df["k2"] != "HC"]
lo = df[df["k2"] == "Lower"]
hi = df[df["k2"] == "Higher"]
LAB = [
    ("AGE", "Age, years"),
    ("EDU", "Education, years"),
    ("ISI_T1", "ISI total"),
    ("ISI_T1_Q1", "  1. Difficulty falling asleep"),
    ("ISI_T1_Q2", "  2. Difficulty staying asleep"),
    ("ISI_T1_Q3", "  3. Early-morning awakening"),
    ("ISI_T1_Q4", "  4. Dissatisfaction with sleep"),
    ("ISI_T1_Q5", "  5. Interference with daytime functioning"),
    ("ISI_T1_Q6", "  6. Noticeability to others"),
    ("ISI_T1_Q7", "  7. Worry or distress"),
    ("PSQI_T1", "PSQI global score"),
    ("BDI_T1", "BDI-II"),
    ("BAI_T1", "BAI"),
]
kw = [kw3(hc[c], lo[c], hi[c]) for c, _ in LAB]
q = multipletests([k[1] for k in kw], method="fdr_bh")[1]
fem = lambda d: f"{int((d['SEX'] == 0).sum())} ({100 * (d['SEX'] == 0).mean():.1f})"
chi2, pchi, _, _ = stats.chi2_contingency(pd.crosstab(df["k2"], df["SEX"]))
pchi_all = stats.chi2_contingency(pd.crosstab(df["k2"] == "HC", df["SEX"]))[1]
rows = [
    [
        "Female, n (%)",
        fem(hc),
        fem(ins),
        fp(pchi_all),
        fem(lo),
        fem(hi),
        f"χ² = {chi2:.2f}",
        fp(pchi),
        "–",
    ]
]
for (c, lab), (H, p, e), qq in zip(LAB, kw, q):
    rows.append(
        [
            lab,
            msd(hc[c]),
            msd(ins[c]),
            fp(mwu_p(hc[c], ins[c])),
            msd(lo[c]),
            msd(hi[c]),
            f"{H:.2f}",
            fp(qq),
            nl(e, 2),
        ]
    )
headers = [
    "Variable",
    f"HC (n = {len(hc)})",
    f"Insomnia symptoms (n = {len(ins)})",
    "pᵃ",
    f"Lower severity (n = {len(lo)})",
    f"Higher severity (n = {len(hi)})",
    "H",
    "q",
    "η²",
]
n_miss = int(df["PSQI_T1"].isna().sum())
write_table(
    "Table1",
    "Table 1",
    "Demographic and clinical characteristics of healthy controls, all participants with insomnia symptoms and the two severity groups.",
    [(None, headers, rows)],
    "Values are M (SD) unless stated otherwise. ᵃHealthy controls versus all participants with insomnia symptoms (Mann–Whitney U test; χ² test for sex), unadjusted. "
    "Severity groups were obtained by K-means clustering of the seven ISI items among participants with ISI ≥ 9 and are used for presentation. "
    "H = Kruskal–Wallis statistic across healthy controls and the two severity groups; q = Benjamini–Hochberg-adjusted p across the 13 continuous variables; "
    "η² = (H − 2)/(N − 3). Because the ISI items define the groups, group differences on the items and the ISI total are descriptive. "
    f"The PSQI was missing for {n_miss} participant{'s' if n_miss != 1 else ''}. ISI = Insomnia Severity Index; PSQI = Pittsburgh Sleep Quality Index; "
    "BDI-II = Beck Depression Inventory-II; BAI = Beck Anxiety Inventory; HC = healthy controls.",
)
