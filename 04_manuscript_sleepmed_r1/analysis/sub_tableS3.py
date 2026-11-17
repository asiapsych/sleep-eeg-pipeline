"""Table S3：頻率範圍、頭皮分布及 aperiodic 參數。

頻率區段列 p < .20 的結果；若沒有，列 p 最小的一段。
這是呈現規則，不代表把顯著門檻改成 .20。"""

from _r1_submission import (
    AN,
    ci_txt,
    csv_sections,
    fp,
    gci,
    pd,
    sg,
    write_table,
)

H, N = AN / "S_R1h_beta_characterisation", AN / "S_R1n_beta_subbands"
rd = lambda p: pd.read_csv(p, encoding="utf-8-sig")
GRP = {"Higher vs HC": "Higher severity vs HC", "Lower vs HC": "Lower severity vs HC"}
GRP3 = {**GRP, "Higher vs Lower": "Higher vs lower severity"}
cap = lambda s: s[0].upper() + s[1:]

# (a)
fc = rd(H / "01_frequency_clusters.csv")
MEAS = {
    "aperiodic-adjusted, whole-scalp": "Aperiodic-adjusted power, whole scalp",
    "aperiodic-adjusted, core 9 electrodes": "Aperiodic-adjusted power, nine fronto-centro-parietal electrodes",
    "raw log power, whole-scalp": "Log power, whole scalp",
}
a = []
for m, lab in MEAS.items():
    for c, clab in GRP.items():
        s = fc[(fc["measure"] == m) & (fc["contrast"] == c)]
        if s.empty:
            a.append(
                [lab, clab, "No bin exceeded the cluster-forming threshold", "–", "–", "–", "–"]
            )
            continue
        keep = s[s["p_cluster"] < 0.20]
        keep = keep if len(keep) else s.nsmallest(1, "p_cluster")
        for _, r in keep.sort_values("from_Hz").iterrows():
            rng = (
                f"{r['from_Hz']:.1f}–{r['to_Hz']:.1f}" if r["n_bins"] > 1 else f"{r['from_Hz']:.1f}"
            )
            a.append(
                [
                    lab,
                    clab,
                    rng,
                    cap(r["direction"]),
                    str(int(r["n_bins"])),
                    fp(r["p_cluster"]),
                    f"{sg(r['peak_g'])} at {r['peak_Hz']:.1f} Hz",
                ]
            )
part_a = (
    "(a) Frequency-wise comparison: clusters of adjacent frequency bins",
    ["Spectrum", "Contrast", "Range, Hz", "Direction", "Bins", "Cluster p", "Largest g"],
    a,
)


# (b)、(c)
def contrast_rows(t, first):
    return [
        [
            *first(r),
            GRP3[r["contrast"]],
            f"{int(r['n'])} / {int(r['n_hc'])}",
            gci(r["g"], r["ci_low"], r["ci_high"]),
            fp(r["p_mwu"]),
            sg(r["adj_diff_in_HC_SD"]),
            fp(r["p_adjusted"]),
            sg(r["m2_diff_in_HC_SD"]),
            fp(r["p_m2"]),
        ]
        for _, r in t.iterrows()
    ]


H9 = [
    "Contrast",
    "n / n comparison",
    "g [95% CI]",
    "p",
    "Model 1 difference",
    "p",
    "Model 2 difference",
    "p",
]
sb = rd(N / "01_subband_contrasts.csv")
hl = rd(N / "02_high_minus_low.csv")
ML = {"aperiodic-adjusted": "Aperiodic-adjusted", "raw log power": "Log power"}
part_b = (
    "(b) Power in beta sub-bands (whole scalp)",
    ["Band", "Power"] + H9,
    contrast_rows(
        sb.sort_values(["measure", "band", "contrast"]), lambda r: [r["band"], ML[r["measure"]]]
    ),
)
part_c = (
    "(c) Difference between 20–30 Hz and 13–20 Hz power",
    ["Measure"] + H9,
    contrast_rows(hl, lambda r: [cap(r["measure"]).replace("raw log power", "log power")]),
)

# (d)
tp = rd(H / "05_topography.csv")
hi, lo = (tp[tp["contrast"] == c].set_index("channel") for c in ("Higher vs HC", "Lower vs HC"))
ORDER = [
    "FP1",
    "FP2",
    "F7",
    "F3",
    "FZ",
    "F4",
    "F8",
    "T3",
    "C3",
    "CZ",
    "C4",
    "T4",
    "T5",
    "P3",
    "PZ",
    "P4",
    "T6",
    "O1",
    "O2",
]
assert set(ORDER) == set(hi.index) == set(lo.index)
ch = lambda c: c[0] + c[1:].lower() if c[1:] in ("P1", "P2", "Z") else c
part_d = (
    "(d) Effect sizes at each electrode",
    [
        "Electrode",
        "Higher severity: adjusted 13–30 Hz power, g",
        "q",
        "Higher severity: peak power, g",
        "q",
        "n with a peak",
        "Lower severity: adjusted 13–30 Hz power, g",
        "q",
        "Lower severity: peak power, g",
        "q",
    ],
    [
        [
            ch(c),
            sg(hi.loc[c, "g_band"]),
            fp(hi.loc[c, "q_band"]),
            sg(hi.loc[c, "g_peak"]),
            fp(hi.loc[c, "q_peak"]),
            str(int(hi.loc[c, "n_with_peak"])),
            sg(lo.loc[c, "g_band"]),
            fp(lo.loc[c, "q_band"]),
            sg(lo.loc[c, "g_peak"]),
            fp(lo.loc[c, "q_peak"]),
        ]
        for c in ORDER
    ],
)

# (e)、(f)
h11, sec = csv_sections(AN / "S_R1_tables" / "R11.csv")
fix = lambda rows: [[r[0].replace("whole-brain", "whole scalp")] + r[1:] for r in rows]
part_e = (
    "(e) Periodic beta peak power by region",
    ["Region"] + h11[1:],
    sec["Periodic beta peak power by region"],
)
part_f = ("(f) Aperiodic component, unadjusted", h11, fix(sec["Aperiodic component"]))

# (g)
mo = rd(AN / "S_R1i_mood_adjustment" / "01_mood_adjustment.csv")
ap = mo[mo["variable"].str.startswith(("exponent", "offset"))]
part_g = (
    "(g) Aperiodic component, adjusted difference from healthy controls in control SD units [95% CI]",
    ["Measure", "Contrast", "n", "Model 1", "p", "Model 2", "p"],
    [
        [
            r["measure"].replace("whole-scalp", "whole scalp"),
            r["group"],
            str(int(r["n"])),
            f"{sg(r['M1_diff_in_HC_SD'])} {ci_txt(r['M1_ci'])}",
            fp(r["M1_p"]),
            f"{sg(r['M2_diff_in_HC_SD'])} {ci_txt(r['M2_ci'])}",
            fp(r["M2_p"]),
        ]
        for _, r in ap.iterrows()
    ],
)

write_table(
    "TableS3",
    "Supplementary Table S3",
    "Frequency range and scalp distribution of the difference in beta activity, and the aperiodic component.",
    [part_a, part_b, part_c, part_d, part_e, part_f, part_g],
    "g = Hedges' g relative to healthy controls (HC; n = 109); higher severity n = 60, lower severity n = 102. "
    "In (a), spectra were residualised for age, sex and education and compared at each frequency bin (0.33-Hz resolution); clusters were formed at two-sided p < .05 and tested with 5,000 permutations of cluster mass. "
    "Clusters with p < .20 are listed; where none reached this level, the cluster with the smallest p is shown. The limits of a cluster depend on the cluster-forming threshold and are not estimates of where an effect begins or ends. "
    "In (b) and (c), p for g is from two-sided Mann–Whitney U tests; Model 1 adjusts for age, sex and education and Model 2 additionally for BDI-II and BAI, with differences in control SD units. "
    "In (d), peak power is periodic beta peak power among participants with a detected beta peak at that electrode; q is Benjamini–Hochberg-adjusted across the 19 electrodes within each contrast and measure.",
)
