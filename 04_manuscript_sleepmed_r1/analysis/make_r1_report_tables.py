"""把 R1 分析結果排成報告用的三線表。

Word、CSV、Markdown 放在 S_R1_tables；PNG 經 LibreOffice 與 pdftoppm 轉出。
描述統計及 R2、R10 的 Kruskal–Wallis 在這裡計算，其餘檢定讀既有 CSV。"""

from __future__ import annotations
import re
import subprocess

import numpy as np
import pandas as pd
from PIL import Image, ImageChops
from scipy import stats
from statsmodels.stats.multitest import multipletests

from apa_docx_helper import section as _section, subheader as _subheader, write_all


def section(label):
    return _section(str(label).replace("*", ""))


def subheader(cells):
    return _subheader([str(c).replace("*", "") for c in cells])


from _r1_common import AN, ISI_ITEMS, OUTCOME, hedges_g, load_all, mwu_p

OUT = AN / "S_R1_tables"
OUT.mkdir(parents=True, exist_ok=True)
SOFFICE = "/Applications/LibreOffice.app/Contents/MacOS/soffice"
Y = OUTCOME


# 格式
def nl(x, d=2):
    """去前導 0（相關、p、比例這類上限為 1 的數）。"""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    s = f"{x:.{d}f}"
    return s.replace("-0.", "−.").replace("0.", ".", 1) if abs(x) < 1 else s.replace("-", "−")


def fp(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "–"
    return "< .001" if p < 0.001 else nl(p, 3)


def sg(x, d=2):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    if round(x, d) == 0:
        return f"{0:.{d}f}"
    return f"{x:+.{d}f}".replace("-", "−")


def uz(x, d=2):
    """不帶正號；避免 −0.00。"""
    return f"{0:.{d}f}" if round(x, d) == 0 else f"{x:.{d}f}".replace("-", "−")


def ci(lo, hi, d=2):
    return f"[{sg(lo, d)}, {sg(hi, d)}]"


def msd(s, d=2):
    s = pd.Series(s).dropna()
    return f"{s.mean():.{d}f} ({s.std(ddof=1):.{d}f})"


def kw3(a, b, c):
    a, b, c = (pd.Series(v).dropna() for v in (a, b, c))
    H, p = stats.kruskal(a, b, c)
    n = len(a) + len(b) + len(c)
    return H, p, max((H - 2) / (n - 3), 0)


def to_png(num):
    docx = OUT / f"{num}.docx"
    subprocess.run(
        [SOFFICE, "--headless", "--convert-to", "pdf", "--outdir", str(OUT), str(docx)],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            "pdftoppm",
            "-png",
            "-r",
            "200",
            "-f",
            "1",
            "-l",
            "1",
            "-singlefile",
            str(OUT / f"{num}.pdf"),
            str(OUT / num),
        ],
        check=True,
        capture_output=True,
    )
    im = Image.open(OUT / f"{num}.png").convert("RGB")
    box = ImageChops.difference(im, Image.new("RGB", im.size, "white")).getbbox()
    pad = 40
    im.crop(
        (
            max(box[0] - pad, 0),
            max(box[1] - pad, 0),
            min(box[2] + pad, im.width),
            min(box[3] + pad, im.height),
        )
    ).save(OUT / f"{num}.png")
    (OUT / f"{num}.pdf").unlink()
    print("saved", OUT / f"{num}.png")


def emit(num, title, headers, rows, note, widths=None, landscape=False):
    num = num.replace("Table", "")  # helper 會自動加上 "Supplementary Table "
    write_all(OUT, num, title, headers, rows, note=note, col_widths=widths, landscape=landscape)
    to_png(num)


df = load_all()
hc = df[df["group"] == "HC"]
ins = df[df["group"] != "HC"]
lo = ins[ins["k2"] == "Lower"]
hi = ins[ins["k2"] == "Higher"]
hcv = hc[Y].dropna()
R1a, R1b, R1c, R1d, R1e, R1f = (
    AN / d
    for d in [
        "S_R1a_severity_form",
        "S_R1b_k234",
        "S_R1c_sensitivity",
        "S_R1d_preproc_qc",
        "S_R1e_clustering_robustness",
        "S_R1f_pc_biomarkers",
    ]
)

# Table R1：群數的效度指標
v = pd.read_csv(R1b / "01_validity_indices.csv")
gap = pd.read_csv(R1e / "04_gap_statistic.csv").set_index("K")
rows = [
    [
        str(int(r.K)),
        r.sizes.replace("/", " / "),
        nl(r.silhouette, 3),
        f"{r.calinski_harabasz:.1f}",
        f"{r.davies_bouldin:.3f}",
        nl(gap.loc[r.K, "gap"], 3),
        nl(r.bootstrap_ARI_median, 2),
        f"[{nl(r.bootstrap_ARI_2_5, 2)}, {nl(r.bootstrap_ARI_97_5, 2)}]",
    ]
    for r in v.itertuples()
]
emit(
    "TableR1",
    "Cluster-Validity Indices and Bootstrap Stability for K-Means Solutions With Two to Six Clusters",
    [
        "*K*",
        "Cluster sizes",
        "Silhouette",
        "Calinski–Harabasz",
        "Davies–Bouldin",
        "Gap",
        "Bootstrap ARI",
        "95% interval",
    ],
    rows,
    "K-means on the seven standardised ISI items in participants with clinically significant insomnia symptoms (*n* = "
    + str(len(ins))
    + "). "
    "Higher silhouette, Calinski–Harabasz and gap values and lower Davies–Bouldin values indicate better-separated clusters. "
    "Bootstrap ARI = median adjusted Rand index between the original solution and solutions refitted in 1,000 bootstrap resamples. "
    "All four internal indices and the gap statistic (Tibshirani criterion) selected *K* = 2.",
    widths=[0.4, 1.7, 0.8, 1.2, 1.1, 0.6, 0.9, 1.0],
    landscape=True,
)

# Table R2：兩群與 HC 的特徵
LAB = [
    ("AGE", "Age, years"),
    ("EDU", "Education, years"),
    ("ISI_T1", "ISI total"),
    ("ISI_T1_Q1", "  1. Sleep-onset difficulty"),
    ("ISI_T1_Q2", "  2. Sleep-maintenance difficulty"),
    ("ISI_T1_Q3", "  3. Early-morning awakening"),
    ("ISI_T1_Q4", "  4. Dissatisfaction with sleep"),
    ("ISI_T1_Q5", "  5. Interference with daytime functioning"),
    ("ISI_T1_Q6", "  6. Noticeability to others"),
    ("ISI_T1_Q7", "  7. Worry about sleep"),
    ("PSQI_T1", "PSQI global score"),
    ("BDI_T1", "BDI-II"),
    ("BAI_T1", "BAI"),
]
tmp = [(lab, msd(hc[c]), msd(lo[c]), msd(hi[c]), *kw3(hc[c], lo[c], hi[c])) for c, lab in LAB]
q = multipletests([t[5] for t in tmp], method="fdr_bh")[1]
allc = [
    (msd(ins[c]), fp(mwu_p(hc[c], ins[c]))) for c, lab in LAB
]  # 全部症狀者，與健康組的 Mann–Whitney p
pchi_all = stats.chi2_contingency(pd.crosstab(df["k2"] == "HC", df["SEX"]))[1]
fem = lambda d: f"{int((d['SEX'] == 0).sum())} ({100 * (d['SEX'] == 0).mean():.1f})"
chi2, pchi, _, _ = stats.chi2_contingency(pd.crosstab(df["k2"], df["SEX"]))
rows = [
    [
        "Female, *n* (%)",
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
rows += [
    [lab, a, al, pa, b, c, f"{H:.2f}", fp(qq), nl(e, 2)]
    for (lab, a, b, c, H, p, e), qq, (al, pa) in zip(tmp, q, allc)
]
emit(
    "TableR2",
    "Demographic and Clinical Characteristics of Healthy Controls, All Participants With Insomnia Symptoms and the Two Severity Groups",
    [
        "Variable",
        f"HC (*n* = {len(hc)})",
        f"Insomnia symptoms (*n* = {len(ins)})",
        "*p*^{a}",
        f"Lower severity (*n* = {len(lo)})",
        f"Higher severity (*n* = {len(hi)})",
        "*H*",
        "*q*",
        "η²",
    ],
    rows,
    "Values are *M* (*SD*) unless stated otherwise. ^{a}HC versus all participants with insomnia symptoms (Mann–Whitney *U* test; χ² test for sex), unadjusted. Severity groups are the K-means *K* = 2 solution on the seven ISI items among participants "
    "with ISI ≥ 9. *H* = Kruskal–Wallis statistic across the three groups; *q* = Benjamini–Hochberg false-discovery-rate-adjusted *p* "
    "across the 13 continuous variables; η² = (*H* − 2)/(*N* − 3). ISI = Insomnia Severity Index; PSQI = Pittsburgh Sleep Quality Index; "
    "BDI-II = Beck Depression Inventory-II; BAI = Beck Anxiety Inventory. The PSQI was missing for one participant in the lower-severity "
    "cluster. ISI items and the ISI total define the clusters, so their group differences are descriptive.",
    widths=[2.5, 1.1, 1.3, 0.6, 1.3, 1.3, 0.8, 0.6, 0.5],
    landscape=True,
)

# Table R3：K = 2／3／4 各群的 periodic beta
bc = pd.read_csv(R1b / "03_beta_by_cluster.csv")
NAME = {
    "Lower": "Lower severity",
    "Higher": "Higher severity",
    "K4_C1": "Cluster 1",
    "K4_C2": "Cluster 2",
    "K4_C3": "Cluster 3",
    "K4_C4": "Cluster 4",
    "MDI": "Cluster 1 (MDI)",
    "MEW": "Cluster 2 (MEW)",
    "SCA": "Cluster 3 (SCA)",
}
COL = {"K=2": "k2", "K=3": "group", "K=4": "k4"}
rows = [["Healthy controls", str(len(hcv)), "–", msd(hcv, 3), "–", "–", "–", "–"]]
for sol in ["K=2", "K=3", "K=4"]:
    t = bc[bc["solution"] == sol].sort_values("ISI_mean")
    rows.append(
        section(
            f"{sol.replace('=', ' = ')} (Kruskal–Wallis *H* = {t['KW_H'].iloc[0]:.2f}, *p* {'= ' + nl(t['KW_p'].iloc[0], 3) if t['KW_p'].iloc[0] >= 0.001 else '< .001'})"
        )
    )
    for r in t.itertuples():
        d = ins.loc[ins[COL[sol]] == r.cluster, Y]
        rows.append(
            [
                NAME[r.cluster],
                str(int(r.n_beta)),
                f"{r.ISI_mean:.1f}",
                msd(d, 3),
                sg(r.g_vs_HC),
                ci(r.ci_low, r.ci_high),
                fp(r.MWU_p),
                fp(r.FDR_q),
            ]
        )
g, l, h, n1, _ = hedges_g(ins[Y], hcv)
rows.append(section("All participants with insomnia symptoms pooled"))
rows.append(
    [
        "ISI ≥ 9",
        str(n1),
        f"{ins['ISI_T1'].mean():.1f}",
        msd(ins[Y], 3),
        sg(g),
        ci(l, h),
        fp(mwu_p(ins[Y], hcv)),
        "–",
    ]
)
emit(
    "TableR3",
    "Whole-Brain Periodic Beta Peak Power by Cluster Under Two-, Three- and Four-Cluster Solutions",
    ["Group", "*n*", "ISI *M*", "Periodic beta, *M* (*SD*)", "*g*", "95% CI", "*p*", "*q*"],
    rows,
    "*g* = Hedges' *g* relative to healthy controls (positive values indicate higher periodic beta than controls); *p* from two-sided "
    "Mann–Whitney *U* tests; *q* = Benjamini–Hochberg-adjusted *p* within each solution. Kruskal–Wallis tests compare the clusters of a "
    "solution and healthy controls. "
    + (
        ""
        if ins[Y].notna().all()
        else f"{int(ins[Y].isna().sum())} participant(s) with insomnia symptoms had no detectable in-band beta peak and are excluded. "
    )
    + "The three-cluster labels (MDI, MEW, SCA) are those used in the original submission.",
    widths=[2.0, 0.5, 0.7, 1.7, 0.7, 1.3, 0.7, 0.7],
    landscape=True,
)

# Table R4：K = 2 主分析
k2 = pd.read_csv(R1b / "04_k2_primary.csv")


def pick(an, con):
    return k2[(k2["analysis"] == an) & (k2["contrast"] == con)].iloc[0]


def nums(s):
    return [float(v) for v in re.findall(r"-?\d+\.\d+", str(s))]


def est(r):
    v = nums(r["estimate"])[0]
    return sg(v) if "g =" in r["estimate"] else uz(v, 3)


def cif(r):
    v = nums(r["ci"])
    if len(v) < 2:
        return "–"
    return ci(v[0], v[1]) if "g =" in r["estimate"] else f"[{uz(v[0], 3)}, {uz(v[1], 3)}]"


def pf(r):
    s = str(r["p"])
    return "< .001" if s == "<.001" else s


_dim = pd.read_csv(
    AN / "S_R1w_pc_interpretation" / "07_two_component_model.csv", encoding="utf-8-sig"
).set_index("component")
rows = [section("Symptom dimensions among participants with insomnia symptoms (standardised β)")]
for _c, _lab in [
    ("PC1", "Dimension 1 (severity)"),
    ("PC2", "Dimension 2 (sleep maintenance and early awakening vs daytime impact)"),
]:
    _r = _dim.loc[_c]
    rows.append(
        [_lab, str(int(_r["n"])), sg(_r["std_beta"]), ci(_r["ci_low"], _r["ci_high"]), fp(_r["p"])]
    )
rows.append(section("Unadjusted (Hedges' *g*)"))
for c, lab in [
    ("Higher vs HC", "Higher severity vs HC"),
    ("Lower vs HC", "Lower severity vs HC"),
    ("Higher vs Lower", "Higher vs lower severity"),
    ("All insomnia (pooled) vs HC", "All insomnia symptoms pooled vs HC"),
]:
    r = pick("unadjusted", c)
    rows.append([lab, f"{int(r['n_a'])} / {int(r['n_b'])}", est(r), cif(r), pf(r)])
for an, lab in [
    ("ANCOVA (age, sex, education)", "Adjusted for age, sex and education (*b*)"),
    ("ANCOVA + BDI + BAI", "Additionally adjusted for BDI-II and BAI (*b*)"),
]:
    rows.append(section(lab))
    for c, cl in [
        ("Higher vs HC", "Higher severity vs HC"),
        ("Lower vs HC", "Lower severity vs HC"),
        ("Higher vs Lower", "Higher vs lower severity"),
    ]:
        r = pick(an, c)
        rows.append([cl, str(int(r["n_a"])), est(r), cif(r), pf(r)])
r = k2[k2["analysis"].str.startswith("equivalence")].iloc[0]
rows.append(section("Equivalence test, bounds ±0.37 *SD* (Hedges' *g*, 90% CI)"))
rows.append(["Lower severity vs HC", f"{int(r['n_a'])} / {int(r['n_b'])}", est(r), cif(r), pf(r)])
bf = k2[k2["analysis"].str.startswith("Bayes factor")].set_index("contrast")
rows.append(section("Bayes factors (JZS prior, *r* = 0.707)"))
_lo, _hi = bf.loc["Lower vs HC"], bf.loc["Higher vs HC"]
_rng = nums(_lo["ci"])
rows.append(
    [
        "Lower severity vs HC, BF₀₁",
        f"{int(_lo['n_a'])} / {int(_lo['n_b'])}",
        f"{nums(_lo['estimate'])[0]:.2f}",
        f"{_rng[2]:.2f}–{_rng[3]:.2f}" + "^{a}",
        "–",
    ]
)
rows.append(
    [
        "Higher severity vs HC, BF₁₀",
        f"{int(_hi['n_a'])} / {int(_hi['n_b'])}",
        f"{int(re.findall(r'= (\d+)', _hi['estimate'])[0]):,}",
        "–",
        "–",
    ]
)
emit(
    "TableR4",
    "Periodic Beta Peak Power: Associations With the Two Symptom Dimensions and Contrasts Between the Two Severity Groups and Healthy Controls",
    ["Contrast", "*n*", "Estimate", "95% CI", "*p*"],
    rows,
    "Symptom dimensions are the two retained principal components of the seven ISI items, entered together and adjusted for age, sex and education. Unadjusted *p* values are from two-sided Mann–Whitney *U* tests. Adjusted estimates are unstandardised regression coefficients "
    f"(*b*) for group membership from ordinary least-squares models of whole-brain periodic beta peak power (*N* = {int(df[Y].notna().sum())}). "
    "The equivalence test is a two one-sided-tests procedure with bounds of ±0.37 pooled *SD*, the meta-analytic effect size for waking beta; "
    "*p* < .05 indicates that the difference lies within the bounds. Bayes factors are Jeffreys–Zellner–Siow Bayes factors for the two-sample "
    "comparison (Rouder et al., 2009); BF₀₁ quantifies evidence for no difference and BF₁₀ evidence for a difference. "
    "^{a}Range of BF₀₁ for prior widths of 0.5 to 1.0. BDI-II = Beck Depression Inventory-II; BAI = Beck Anxiety Inventory.",
    widths=[3.1, 0.9, 0.9, 1.5, 0.8],
)

# Table R5：嚴重度分段
b = pd.read_csv(R1a / "03_isi_bins.csv")


def band_rows(scheme, bounds):
    out = []
    for r, (a_, b_) in zip(b[b["scheme"] == scheme].itertuples(), bounds):
        d = ins[(ins["ISI_T1"] >= a_) & (ins["ISI_T1"] <= b_)][Y]
        out.append(
            [
                r.ISI_range.replace(">=", "≥ ").replace("-", "–"),
                str(int(r.n)),
                msd(d, 3),
                sg(r.g_vs_HC),
                ci(r.ci_low, r.ci_high),
                fp(r.MWU_p),
                f"{int(r.n_MEW)} / {int(r.n_MDI)} / {int(r.n_SCA)}",
            ]
        )
    return out


rows = [section("Two-point bands")] + band_rows(
    "2-point bins", [(9, 10), (11, 12), (13, 14), (15, 16), (17, 18), (19, 28)]
)
rows += [section("Conventional ISI categories")] + band_rows(
    "conventional ISI bands", [(9, 14), (15, 21), (22, 28)]
)
rows += [section("All participants with insomnia symptoms")] + band_rows(
    "all insomnia (pooled)", [(9, 28)]
)
emit(
    "TableR5",
    "Whole-Brain Periodic Beta Peak Power Relative to Healthy Controls by Insomnia Severity Index Band",
    ["ISI band", "*n*", "Periodic beta, *M* (*SD*)", "*g*", "95% CI", "*p*", "Original profiles"],
    rows,
    f"*g* = Hedges' *g* relative to healthy controls (*n* = {len(hcv)}; *M* = {hcv.mean():.3f}, *SD* = {hcv.std(ddof=1):.3f}); *p* from two-sided "
    "Mann–Whitney *U* tests, uncorrected. Original profiles = number of participants from the three-cluster solution of the original "
    "submission in each band (MEW / MDI / SCA).",
    widths=[1.1, 0.5, 1.8, 0.7, 1.4, 0.7, 1.5],
    landscape=True,
)

# Table R6：函數形式比較與門檻
mc = pd.read_csv(R1a / "04_model_comparison.csv")
thr = pd.read_csv(R1a / "05_threshold.csv")
C_HAT = int(
    float(thr.loc[thr["quantity"].str.startswith("data-selected threshold"), "value"].iloc[0])
)
ML = [
    ("covariates only", "Covariates only"),
    ("linear ISI total", "Linear ISI total"),
    ("quadratic ISI total", "Quadratic ISI total"),
    # patsy 的 cr(x, df=k) 不加 center 時含常數項，樣條項只有 k−1 個；這裡照樣條項的個數標示，與正文一致
    ("spline df=3", "Natural cubic spline, 2 *df*"),
    ("spline df=4", "Natural cubic spline, 3 *df*"),
    ("spline df=5", "Natural cubic spline, 4 *df*"),
    ("single threshold ISI>=15 (conventional)", "Single threshold, ISI ≥ 15 (conventional)"),
    (
        f"single threshold ISI>={C_HAT} (data-selected)",
        f"Single threshold, ISI ≥ {C_HAT} (data-selected)^{{a}}",
    ),
    ("K=2 higher-severity cluster", "Two-cluster membership"),
    ("K=3 profile (3 levels)", "Three-cluster membership"),
    ("SCA indicator", "Highest-severity cluster of three vs rest"),
    ("PC1", "First principal component"),
]
mi = mc[mc["scope"] == "insomnia only"].set_index("model")
lin = mi.loc["linear ISI total", "AIC"]
rows = [
    [
        lab,
        str(int(mi.loc[m, "k_severity"])),
        nl(mi.loc[m, "R2"], 3),
        nl(mi.loc[m, "adj_R2"], 3),
        f"{mi.loc[m, 'AIC']:.1f}".replace("-", "−"),
        f"{mi.loc[m, 'BIC']:.1f}".replace("-", "−"),
        "Reference" if m == "linear ISI total" else sg(mi.loc[m, "AIC"] - lin, 1),
    ]
    for m, lab in ML
]
tv = lambda key: str(thr.loc[thr["quantity"].str.startswith(key)].iloc[0]["value"])
tn = lambda key: str(thr.loc[thr["quantity"].str.startswith(key)].iloc[0]["note"])
perm = re.search(r"= ([\d.]+)$", tn("max-scanned")).group(1)
emit(
    "TableR6",
    "Comparison of Functional Forms Relating Insomnia Severity to Whole-Brain Periodic Beta Peak Power",
    ["Severity term", "*k*", "*R*²", "Adjusted *R*²", "AIC", "BIC", "ΔAIC"],
    rows,
    f"Ordinary least-squares models in participants with insomnia symptoms (*n* = 162), each including age, sex and education. *k* = number of "
    f"severity parameters. ΔAIC is relative to the linear model; lower AIC and BIC indicate better fit. ^{{a}}Threshold selected by scanning "
    f"ISI 11–22; one additional parameter is counted for the selected cut-point. The selected threshold was ISI ≥ {C_HAT} "
    f"(permutation *p* = {nl(float(perm), 4)} for the maximally selected step, 5,000 permutations; 94% of 5,000 bootstrap resamples selected "
    f"a threshold between 17 and 20). Added to the linear model, the step remained significant (*p* = {tv('step added').lstrip('0') if tv('step added')[0] != '<' else tv('step added')}); "
    f"added to the step model, the linear term did not (*p* = {tv('linear added')}).",
    widths=[3.6, 0.4, 0.6, 1.0, 0.8, 0.8, 0.7],
    landscape=True,
)

# Table R7：PCA 與 PC 迴歸
pl = pd.read_csv(R1a / "01_pca_loadings.csv", index_col=0)
pr = pd.read_csv(R1a / "02_pc_regression.csv")
ITL = dict(zip(ISI_ITEMS, [l.strip() for _, l in LAB[3:10]]))
rows = [section("Panel A. Component loadings")]
rows += [
    [ITL[i], nl(pl.loc[i, "PC1"]), nl(pl.loc[i, "PC2"]), nl(pl.loc[i, "PC3"]), ""]
    for i in ISI_ITEMS
]
rows += [
    [
        "Variance explained",
        *[f"{100 * pl.loc['explained_variance_ratio', c]:.1f}%" for c in ["PC1", "PC2", "PC3"]],
        "",
    ],
    [
        "Correlation with ISI total",
        *[nl(pl.loc["r_with_ISI_total", c]) for c in ["PC1", "PC2", "PC3"]],
        "",
    ],
]
rows += [
    [f"Mean score, {NAME[g_]}", *[sg(pl.loc[f"mean_{g_}", c]) for c in ["PC1", "PC2", "PC3"]], ""]
    for g_ in ["MDI", "MEW", "SCA"]
]
rows.append(section("Panel B. Regression of periodic beta peak power on the components"))
rows.append(subheader(["Term", "*b* (*SE*)", "β", "*t*", "*p*"]))
adj = pr[pr["model"] == "with covariates"]
rows += [
    [
        r.term,
        f"{float(r.b):.3f} ({float(r.se):.3f})".replace("-", "−"),
        sg(float(r.std_beta)),
        sg(float(r.t)),
        fp(float(r.p)),
    ]
    for r in adj.itertuples()
]
nf = pr[pr["model"].str.startswith("PC2+PC3")].iloc[0]
rows.append(
    [
        "PC2 and PC3 beyond PC1",
        f"Δ*R*² = {nl(float(nf['R2_model']), 3)}",
        "",
        nf["t"].replace("F(", "*F*(").replace("=", " = "),
        fp(float(nf["p"])),
    ]
)
emit(
    "TableR7",
    "Principal Components of the Seven Insomnia Severity Index Items and Their Association With Periodic Beta Peak Power",
    ["ISI item", "PC1", "PC2", "PC3", ""],
    rows,
    "Principal component analysis of the standardised ISI items in participants with insomnia symptoms (*n* = "
    + str(len(ins))
    + "). Panel B: ordinary "
    "least-squares model with the three components entered together, adjusted for age, sex and education (*n* = 162); β = standardised "
    "coefficient. Cluster labels are those of the three-cluster solution in the original submission.",
    widths=[2.7, 1.2, 0.7, 1.2, 0.6],
)

# Table R8：九種分群方法
mt = pd.read_csv(R1e / "01_methods_k2_k3.csv")
sm = pd.read_csv(R1e / "03_ari_summary.csv")


def arif(r):
    return (
        "Reference"
        if r.method == "K-means"
        else f"{nl(r.ARI_vs_Kmeans)} [{nl(float(r.ARI_ci_low))}, {nl(float(r.ARI_ci_high))}]"
    )


rows = []
for k in (2, 3):
    rows.append(section(f"*K* = {k}"))
    for r in mt[mt["K"] == k].itertuples():
        rows.append(
            [
                r.method,
                r.sizes.replace("/", " / "),
                arif(r),
                f"{sg(r.largest_g)} {r.largest_ci.replace('-', '−')}",
                f"{r.largest_ISI_mean:.1f}",
                f"{sg(r.smallest_g)} {r.smallest_ci.replace('-', '−')}",
                f"{r.smallest_ISI_mean:.1f}",
            ]
        )
a = sm.iloc[0]
c2 = lambda s: "[" + ", ".join(nl(float(v)) for v in re.findall(r"-?\d+\.\d+", s)) + "]"
emit(
    "TableR8",
    "Agreement Among Nine Clustering Methods and Periodic Beta Effect Sizes for Two- and Three-Cluster Solutions",
    [
        "Method",
        "Cluster sizes",
        "ARI vs K-means [95% CI]",
        "Largest *g* [95% CI]",
        "ISI *M*",
        "Smallest *g* [95% CI]",
        "ISI *M*",
    ],
    rows,
    "ARI = adjusted Rand index, with paired-bootstrap 95% confidence intervals (2,000 resamples). Largest and smallest *g* = Hedges' *g* "
    "relative to healthy controls for the cluster with the largest and the smallest effect; ISI *M* = mean ISI total of that cluster. "
    "In all 18 specifications the cluster with the largest effect was the cluster with the highest mean ISI. "
    f"Mean pairwise ARI across the 36 method pairs was {nl(a['mean_ARI_K2'])} {c2(a['K2_ci'])} for *K* = 2 and {nl(a['mean_ARI_K3'])} "
    f"{c2(a['K3_ci'])} for *K* = 3 (difference {nl(a['difference_K2_minus_K3'])} {c2(a['diff_ci'])}). "
    "GMM = Gaussian mixture model (scikit-learn; 30 initialisations); LCA = latent class analysis with a categorical measurement model "
    "(StepMix 3.0; 20 initialisations); hierarchical methods used Euclidean distance.",
    widths=[1.7, 1.1, 1.7, 1.6, 0.6, 1.6, 0.6],
    landscape=True,
)

# Table R9：敏感度分析
e = pd.read_csv(R1c / "03_emg_covariate.csv").rename(columns={"added_covariate": "setting"})
vg = pd.read_csv(R1c / "05b_vigilance_covariate.csv").rename(columns={"added_covariate": "setting"})
qc = pd.read_csv(R1d / "02_qc_covariate.csv")
rs = pd.read_csv(R1d / "03_restricted_samples.csv")
HC_SD = float(hc["beta_pw_BRAIN"].std(ddof=1))


def cov_row(t, s, lab):
    out = [lab]
    for con in ["K2 Higher vs HC", "K2 Lower vs HC"]:
        r = t[(t["setting"] == s) & (t["contrast"] == con)].iloc[0]
        sd = HC_SD  # 直接用健康組 SD；不從四捨五入後的係數倒推（2026-11-16 修）
        out += [
            f"{uz(r['b_in_HC_SD'])} [{uz(r['ci_low'] / sd)}, {uz(r['ci_high'] / sd)}]",
            fp(r["p"]),
        ]
    return out


rows = [
    section("Panel A. Adjusted difference from healthy controls, in control *SD* units [95% CI]")
]
rows.append(cov_row(e, "none (age, sex, education only)", "Age, sex and education only"))
rows += [
    cov_row(e, s, l)
    for s, l in [
        ("+ 30–45 Hz relative power", "+ 30–45 Hz relative power"),
        ("+ log10 30–45 Hz absolute power", "+ 30–45 Hz absolute power (log)"),
        ("+ aperiodic exponent (whole-brain)", "+ aperiodic exponent"),
    ]
]
rows += [
    cov_row(vg, s, l)
    for s, l in [
        ("+ theta/alpha ratio", "+ theta/alpha ratio"),
        ("+ theta relative power", "+ relative theta power"),
        ("+ alpha peak frequency", "+ alpha peak frequency"),
    ]
]
rows += [
    cov_row(qc, s, l)
    for s, l in [
        ("+ n_ic_removed", "+ number of components removed"),
        ("+ n_ret_muscle_dom", "+ number of retained muscle-dominant components"),
        ("+ max_ret_muscle_p", "+ maximum retained muscle probability"),
        ("+ reject_rate", "+ epoch rejection rate"),
        ("+ all of the above QC metrics", "+ all preprocessing metrics"),
    ]
]
rows.append(section("Panel B. Hedges' *g* vs healthy controls in restricted samples [95% CI]"))
rows.append(
    subheader(["Sample", "Higher severity", "*n* / *n* HC", "Lower severity", "*n* / *n* HC"])
)
for s, l in [
    ("full sample", "Full sample"),
    ("no retained muscle-dominant component", "No retained muscle-dominant component"),
    ("max retained muscle probability < 0.50", "Maximum retained muscle probability < .50"),
    ("max retained muscle probability < 0.20", "Maximum retained muscle probability < .20"),
    ("no interpolated channel", "No interpolated channel"),
    ("epoch rejection rate <= 10%", "Epoch rejection rate ≤ 10%"),
]:
    out = [l]
    for con in ["K2 Higher vs HC", "K2 Lower vs HC"]:
        r = rs[(rs["subsample"] == s) & (rs["contrast"] == con)].iloc[0]
        out += [
            f"{sg(r['g'])} {ci(r['ci_low'], r['ci_high'])}",
            f"{int(r['n_a'])} / {int(r['n_HC'])}",
        ]
    rows.append(out)
emit(
    "TableR9",
    "Sensitivity of the Periodic Beta Difference to Muscle-Activity Proxies, Vigilance Markers and Preprocessing Quality",
    ["Covariate added", "Higher severity vs HC", "*p*", "Lower severity vs HC", "*p*"],
    rows,
    f"Panel A: group coefficients from ordinary least-squares models of whole-brain periodic beta peak power (*N* = {int(df[Y].notna().sum())}), divided by the "
    "control-group *SD*; each row adds the listed covariate to age, sex and education. Panel B: unadjusted Hedges' *g* after restricting all "
    "groups to recordings meeting the criterion. Muscle-dominant component = independent component whose highest ICLabel class probability "
    "was muscle and that was retained because the probability was below the removal threshold of .80.",
    widths=[3.4, 1.9, 0.7, 1.9, 0.9],
    landscape=True,
)

# Table R10：前處理品管
sub = pd.read_csv(R1d / "_subject_level_internal.csv")
fq = pd.read_csv(AN / "S10_fooof_full" / "02_fit_quality.csv")[
    ["subject", "r2_mean", "fit_error_mean"]
]
sub = sub.merge(fq, on="subject", how="left")
sub["any_ret"] = (sub["n_ret_muscle_dom"] > 0).astype(float)
G3 = [sub[sub["k2"] == g_] for g_ in ["HC", "Lower", "Higher"]]
QL = [
    ("n_bad_ch", "Interpolated channels", 2),
    ("n_ic_removed", "Components removed", 2),
    ("n_removed_eye", "  Ocular", 2),
    ("n_removed_muscle", "  Muscle", 2),
    ("n_ret_muscle_dom", "Retained muscle-dominant components", 2),
    ("max_ret_muscle_p", "Maximum retained muscle probability", 2),
    ("reject_rate", "Epoch rejection rate", 3),
    ("r2_mean", "Spectral model fit, *R*²", 3),
    ("fit_error_mean", "Spectral model fit error", 3),
]
rows = []
for c, lab, d in QL:
    H, p, _ = kw3(*[g_[c] for g_ in G3])
    gh, l_, h_, _, _ = hedges_g(G3[2][c], G3[0][c])
    rows.append([lab, *[msd(g_[c], d) for g_ in G3], f"{H:.2f}", fp(p), f"{sg(gh)} {ci(l_, h_)}"])
pct = lambda g_: f"{int(g_['any_ret'].sum())} ({100 * g_['any_ret'].mean():.1f})"
chi2, pchi, _, _ = stats.chi2_contingency(pd.crosstab(sub["k2"], sub["any_ret"]))
rows.insert(
    5,
    [
        "  Participants with ≥ 1, *n* (%)",
        pct(G3[0]),
        pct(G3[1]),
        pct(G3[2]),
        f"χ² = {chi2:.2f}",
        fp(pchi),
        "–",
    ],
)
emit(
    "TableR10",
    "Preprocessing Quality Metrics in Healthy Controls and the Two Insomnia-Severity Clusters",
    [
        "Metric",
        f"HC (*n* = {len(G3[0])})",
        f"Lower severity (*n* = {len(G3[1])})",
        f"Higher severity (*n* = {len(G3[2])})",
        "*H*",
        "*p*",
        "*g*, higher vs HC [95% CI]",
    ],
    rows,
    "Values are *M* (*SD*) unless stated otherwise. Components were removed automatically when the ICLabel probability was ≥ .70 for ocular "
    "or ≥ .80 for muscle activity. Retained muscle-dominant components had muscle as their most probable class but fell below the removal "
    "threshold. Epoch rejection rate = proportion of 3-s epochs rejected before 60 artefact-free epochs were sampled per participant. "
    "*H* = Kruskal–Wallis statistic across the three groups, uncorrected.",
    widths=[2.9, 1.2, 1.5, 1.5, 0.8, 0.6, 1.9],
    landscape=True,
)

# Table R11：其他腦波指標
src = {
    "emg": pd.read_csv(R1c / "02_emg_indices.csv"),
    "vig": pd.read_csv(R1c / "05_vigilance.csv"),
    "ap": pd.read_csv(R1c / "06_aperiodic.csv"),
    "reg": pd.read_csv(R1c / "04_regional.csv"),
}
BLK = [
    (
        "Periodic beta peak power by region",
        "reg",
        [
            ("beta_pw_FRONTAL", "Frontal"),
            ("beta_pw_CENTRAL", "Central"),
            ("beta_pw_PARIETAL", "Parietal"),
            ("beta_pw_TEMPORAL", "Temporal"),
            ("beta_pw_OCCIPITAL", "Occipital"),
            ("beta_pw_MIDLINE", "Midline (Fz, Cz, Pz)"),
        ],
    ),
    (
        "Aperiodic component",
        "ap",
        [
            ("exponent_BRAIN", "Exponent, whole-brain"),
            ("offset_BRAIN", "Offset, whole-brain"),
            ("exponent_CENTRAL", "Exponent, central"),
            ("offset_CENTRAL", "Offset, central"),
        ],
    ),
    (
        "High-frequency power (muscle-activity proxies)",
        "emg",
        [
            ("Gamma_rel_BRAIN", "30–45 Hz, relative"),
            ("logGamma_abs", "30–45 Hz, absolute (log)"),
            ("HighBeta_rel_BRAIN", "30–40 Hz, relative"),
        ],
    ),
    (
        "Vigilance markers",
        "vig",
        [
            ("Theta_rel_BRAIN", "Relative theta"),
            ("Alpha_rel_BRAIN", "Relative alpha"),
            ("TAR_BRAIN", "Theta/alpha ratio"),
            ("alpha_cf_BRAIN", "Alpha peak frequency"),
            ("theta_pw_BRAIN", "Periodic theta peak power"),
        ],
    ),
]
rows = []
for title, key, feats in BLK:
    rows.append(section(title))
    t = src[key]
    for f, lab in feats:
        out = [lab]
        for con in ["K2 Higher vs HC", "K2 Lower vs HC"]:
            r = t[(t["feature"] == f) & (t["contrast"] == con)].iloc[0]
            out += [f"{sg(r['g'])} {ci(r['ci_low'], r['ci_high'])}", fp(r["MWU_p"])]
        rows.append(out)
emit(
    "TableR11",
    "Regional Periodic Beta, Aperiodic Parameters, High-Frequency Power and Vigilance Markers in the Two Insomnia-Severity Clusters",
    [
        "Measure",
        "Higher severity vs HC, *g* [95% CI]",
        "*p*",
        "Lower severity vs HC, *g* [95% CI]",
        "*p*",
    ],
    rows,
    "*g* = Hedges' *g* relative to healthy controls; *p* from two-sided Mann–Whitney *U* tests, uncorrected. Drowsiness is typically "
    "accompanied by higher theta, lower alpha, a higher theta/alpha ratio and a slower alpha peak frequency.",
    widths=[3.0, 2.3, 0.7, 2.3, 0.7],
    landscape=True,
)
