"""Table S5：探索性的 EEG–HRV 關聯。

HRV 只用 task_r1q 通過訊號品管的重算值；task_r1m 的 BH 校正在各 HRV 指標內進行。"""

from _r1_submission import (
    AN,
    ci,
    fp,
    gci,
    nl,
    pd,
    sg,
    write_table,
)

rd = lambda p: pd.read_csv(p, encoding="utf-8-sig")
Q, M = AN / "S_R1q_hrv_recompute", AN / "S_R1m_hrv"

qc = rd(Q / "01_qc.csv")
GL = {
    "HC": "Healthy controls",
    "Lower": "Lower severity",
    "Higher": "Higher severity",
    "All": "All",
}
REASONS = [
    "no usable ECG channel",
    "inconsistent beat morphology",
    "corrected beats > 5%",
    "mean heart rate outside 40–120 bpm",
    "RR interval outside 300–2000 ms",
]
part_a = (
    "(a) Signal quality control of the electrocardiogram",
    ["Group", "n in the sample", "n with a raw recording", "n retained (%)"]
    + [r[0].upper() + r[1:] for r in REASONS],
    [
        [
            GL[r["group"]],
            str(int(r["n_sample"])),
            str(int(r["n_with_raw_file"])),
            f"{int(r['n_pass'])} ({r['pct_pass_of_raw']:.1f})",
        ]
        + [str(int(r[k])) for k in REASONS]
        for _, r in qc.iterrows()
    ],
)

asc = rd(M / "02_associations.csv")
HV = ["ln RMSSD", "ln HF power", "ln SDNN", "ln LF/HF", "Heart rate"]
eegs = list(dict.fromkeys(asc["eeg"]))
ix = asc.set_index(["eeg", "hrv"])
part_b = (
    "(b) EEG measures and heart rate variability, standardised β (q)",
    ["EEG measure", "n"] + HV,
    [
        [e_.replace("Raw log", "Log"), str(int(ix.loc[(e_, HV[0]), "n"]))]
        + [f"{sg(ix.loc[(e_, h), 'beta'])} ({fp(ix.loc[(e_, h), 'q_within_hrv'])})" for h in HV]
        for e_ in eegs
    ],
)

rb = rd(M / "03_robustness.csv")
rb = rb[
    ~rb["model"].str.contains("medication|all of the above", case=False)
]  # 含用藥的模型不放（task_r1o 備用）
part_c = (
    "(c) Periodic beta peak power and ln RMSSD under alternative models",
    ["Model", "n", "β [95% CI]", "p", "p, HC3", "Bootstrap 95% CI"],
    [
        [
            r["model"]
            .replace("Within HC", "Within healthy controls")
            .replace("Within Lower", "Within lower severity")
            .replace("Within Higher", "Within higher severity")
            .replace("severity cluster", "severity group"),
            str(int(r["n"])),
            gci(r["beta"], r["ci_low"], r["ci_high"]),
            fp(r["p"]),
            fp(r["p_hc3"]),
            ci(r["boot_low"], r["boot_high"]),
        ]
        for _, r in rb.iterrows()
    ],
)

ag = rd(Q / "02_agreement.csv")
part_d = (
    "(d) Agreement between recomputed values and those exported by the acquisition software",
    [
        "Measure",
        "n",
        "Spearman ρ",
        "Median, software",
        "Median, recomputed",
        "n within 10%",
        "n differing by more than 10%",
    ],
    [
        [
            r["measure"],
            str(int(r["n"])),
            nl(r["spearman_rho"], 2),
            f"{r['median_software']:.1f}",
            f"{r['median_computed']:.1f}",
            str(int(r["n_within_10pct"])),
            str(int(r["n_diff_gt_10pct"])),
        ]
        for _, r in ag.iterrows()
    ],
)

write_table(
    "TableS5",
    "Supplementary Table S5",
    "Exploratory analysis of resting heart rate variability and EEG measures.",
    [part_a, part_b, part_c, part_d],
    "The electrocardiogram was recorded concurrently with the EEG in a subset of participants; a recording could fail more than one quality criterion. "
    "In (b), each cell is the standardised coefficient of the EEG measure in a model of the heart rate variability index adjusted for age, sex and education, with q the Benjamini–Hochberg-adjusted p across the 13 EEG measures within each index. "
    "In (c), rows beginning with + add the named covariate to the first model; the bootstrap interval is from 2,000 resamples of participants. "
    "These analyses were not planned before data collection, the subset was determined by the availability and quality of recordings, and the results are reported as exploratory.",
)
