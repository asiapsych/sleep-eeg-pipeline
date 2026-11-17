"""並列 K = 2、3、4 的分群依據，重算兩群下的主要結果。

分群只用 ISI 七題；EEG 效果量不能當成選 K 的依據。
低嚴重度組對 HC 另報 TOST（±0.37 SD）與 JZS Bayes factor。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from scipy.integrate import quad
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ttost_ind

from _r1_common import AN, COV, COV_COLS, ISI_ITEMS, OUTCOME, SEED, fmt_p, hedges_g, load_all, mwu_p

Y = OUTCOME
SESOI = 0.37
OUT = AN / "S_R1b_k234"
OUT.mkdir(parents=True, exist_ok=True)
MD = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    MD.append(s)


df = load_all()
df["k3"] = df["group"]
hc = df[df["group"] == "HC"]
ins = df[df["group"] != "HC"]
hcv = hc[Y].dropna()
P("# S_R1b — K = 2／3／4 並列\n")
P(f"N = {len(df)}；HC = {len(hc)}；失眠 = {len(ins)}\n")

# 1. 內部效度指標
X = StandardScaler().fit_transform(SimpleImputer(strategy="median").fit_transform(ins[ISI_ITEMS]))
ari = pd.read_csv(AN / "S16_bootstrap_stability" / "05_selectK_bootstrap_ARI.csv").set_index("K")
rows = []
for k in range(2, 7):
    lab = KMeans(n_clusters=k, random_state=SEED, n_init=100).fit_predict(X)
    rows.append(
        dict(
            K=k,
            sizes="/".join(map(str, sorted(pd.Series(lab).value_counts().tolist(), reverse=True))),
            silhouette=round(silhouette_score(X, lab), 3),
            calinski_harabasz=round(calinski_harabasz_score(X, lab), 1),
            davies_bouldin=round(davies_bouldin_score(X, lab), 3),
            bootstrap_ARI_median=ari.loc[k, "ari_median"] if k in ari.index else np.nan,
            bootstrap_ARI_2_5=ari.loc[k, "ari_2.5"] if k in ari.index else np.nan,
            bootstrap_ARI_97_5=ari.loc[k, "ari_97.5"] if k in ari.index else np.nan,
        )
    )
val = pd.DataFrame(rows)
val.to_csv(OUT / "01_validity_indices.csv", index=False, encoding="utf-8-sig")
P("## 1. 內部效度指標與 bootstrap 穩定度（silhouette、CH 越大越好；DB 越小越好）\n")
P(val.to_markdown(index=False))
P("")

# 2–3. 各群特徵與 beta
SOL = {"K=2": "k2", "K=3": "k3", "K=4": "k4"}
sex_codes = sorted(df["SEX"].dropna().unique().tolist())
char, beta = [], []
for sol, col in SOL.items():
    names = ins.groupby(col)["ISI_T1"].mean().sort_values().index.tolist()
    ps, tmp = [], []
    for nm in names:
        d = ins[ins[col] == nm]
        r = dict(
            solution=sol,
            cluster=nm,
            n=len(d),
            ISI_mean=round(d["ISI_T1"].mean(), 2),
            ISI_sd=round(d["ISI_T1"].std(ddof=1), 2),
            ISI_min=int(d["ISI_T1"].min()),
            ISI_max=int(d["ISI_T1"].max()),
        )
        for q in ISI_ITEMS:
            r[q.replace("ISI_T1_", "")] = round(d[q].mean(), 2)
        r.update(
            BDI_mean=round(d["BDI_T1"].mean(), 2),
            BAI_mean=round(d["BAI_T1"].mean(), 2),
            AGE_mean=round(d["AGE"].mean(), 2),
        )
        for sc in sex_codes:
            r[f"SEX={sc:g}_pct"] = round(100 * (d["SEX"] == sc).mean(), 1)
        char.append(r)
        g, lo, hi, n1, n2 = hedges_g(d[Y], hcv)
        p = mwu_p(d[Y], hcv)
        ps.append(p)
        tmp.append(
            dict(
                solution=sol,
                cluster=nm,
                n_beta=n1,
                ISI_mean=round(d["ISI_T1"].mean(), 2),
                mean_beta=round(d[Y].mean(), 3),
                g_vs_HC=round(g, 3),
                ci_low=round(lo, 3),
                ci_high=round(hi, 3),
                MWU_p=round(p, 4),
            )
        )
    q = multipletests(ps, method="fdr_bh")[1]
    H, pk = stats.kruskal(hcv, *[ins.loc[ins[col] == nm, Y].dropna() for nm in names])
    for t, qq in zip(tmp, q):
        t.update(FDR_q=round(qq, 4), KW_H=round(H, 2), KW_p=round(pk, 5))
    beta += tmp
char = pd.DataFrame(char)
char.to_csv(OUT / "02_cluster_characteristics.csv", index=False, encoding="utf-8-sig")
beta = pd.DataFrame(beta)
beta.to_csv(OUT / "03_beta_by_cluster.csv", index=False, encoding="utf-8-sig")
P("## 2. 各解的群集特徵（依 ISI 總分由低到高）\n")
P(char.to_markdown(index=False))
P("")
gp, _, _, _, _ = hedges_g(ins[Y], hcv)
P(f"## 3. 各群對 HC 的 periodic beta（HC n = {len(hcv)}；全部失眠合併 g = {gp:.3f}）\n")
P(beta.to_markdown(index=False))
P("")

# 4. K = 2 當主分析
hi_ = ins[ins["k2"] == "Higher"]
lo_ = ins[ins["k2"] == "Lower"]
k2 = []
for label, a, b in [
    ("Higher vs HC", hi_[Y], hcv),
    ("Lower vs HC", lo_[Y], hcv),
    ("Higher vs Lower", hi_[Y], lo_[Y]),
    ("All insomnia (pooled) vs HC", ins[Y], hcv),
]:
    g, l, h, n1, n2 = hedges_g(a, b)
    k2.append(
        dict(
            analysis="unadjusted",
            contrast=label,
            n_a=n1,
            n_b=n2,
            estimate=f"g = {g:.3f}",
            ci=f"[{l:.3f}, {h:.3f}]",
            p=fmt_p(mwu_p(a, b)),
            note="Mann–Whitney U",
        )
    )
gh = hedges_g(hi_[Y], hcv)[0]
k2.append(
    dict(
        analysis="dilution",
        contrast="Higher g / pooled g",
        estimate=f"{gh / gp:.2f}",
        note="高分群效果量是合併效果量的幾倍",
    )
)

d = df.dropna(subset=[Y] + COV_COLS).copy()
d["k2"] = pd.Categorical(d["k2"], categories=["HC", "Lower", "Higher"])
for label, f, need in [
    ("ANCOVA (age, sex, education)", f"{Y} ~ {COV} + C(k2)", []),
    ("ANCOVA + BDI + BAI", f"{Y} ~ {COV} + BDI_T1 + BAI_T1 + C(k2)", ["BDI_T1", "BAI_T1"]),
]:
    dd = d.dropna(subset=need)
    m = smf.ols(f, dd).fit()
    for t in ["C(k2)[T.Lower]", "C(k2)[T.Higher]"]:
        ci = m.conf_int().loc[t]
        k2.append(
            dict(
                analysis=label,
                contrast=t.replace("C(k2)[T.", "").replace("]", "") + " vs HC",
                n_a=int(m.nobs),
                estimate=f"b = {m.params[t]:.3f}",
                ci=f"[{ci[0]:.3f}, {ci[1]:.3f}]",
                p=fmt_p(m.pvalues[t]),
                note=f"model R2 = {m.rsquared:.3f}",
            )
        )
    w = m.t_test("C(k2)[T.Higher] - C(k2)[T.Lower] = 0")
    k2.append(
        dict(
            analysis=label,
            contrast="Higher vs Lower",
            n_a=int(m.nobs),
            estimate=f"b = {float(w.effect):.3f}",
            p=fmt_p(float(w.pvalue)),
            note="",
        )
    )

a = lo_[Y].dropna().to_numpy()
b = hcv.to_numpy()
sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
p_tost, _, _ = ttost_ind(a, b, -SESOI * sp, SESOI * sp, usevar="pooled")
g, _, _, n1, n2 = hedges_g(a, b)
se = np.sqrt((n1 + n2) / (n1 * n2) + g**2 / (2 * (n1 + n2)))
k2.append(
    dict(
        analysis=f"equivalence (TOST, bound ±{SESOI} SD)",
        contrast="Lower vs HC",
        n_a=n1,
        n_b=n2,
        estimate=f"g = {g:.3f}",
        ci=f"90% CI [{g - 1.645 * se:.3f}, {g + 1.645 * se:.3f}]",
        p=fmt_p(p_tost),
        note="equivalent" if p_tost < 0.05 else "not established",
    )
)


def jzs_bf01(x, y, r=0.707):
    """兩獨立樣本 t 檢定的 JZS 貝氏因子 BF01（Rouder et al., 2009，式 5）；r 為 Cauchy 先驗的寬度。"""
    t = stats.ttest_ind(x, y).statistic
    nx, ny = len(x), len(y)
    n = nx * ny / (nx + ny)
    v = nx + ny - 2
    f = lambda gg: (
        (1 + n * gg * r**2) ** -0.5
        * (1 + t**2 / ((1 + n * gg * r**2) * v)) ** (-(v + 1) / 2)
        * (2 * np.pi) ** -0.5
        * gg**-1.5
        * np.exp(-1 / (2 * gg))
    )
    return (1 + t**2 / v) ** (-(v + 1) / 2) / quad(f, 0, np.inf)[0]


bf_lo = {r: jzs_bf01(a, b, r) for r in (0.5, 0.707, 1.0)}
k2.append(
    dict(
        analysis="Bayes factor (JZS, r = 0.707)",
        contrast="Lower vs HC",
        n_a=n1,
        n_b=n2,
        estimate=f"BF01 = {bf_lo[0.707]:.2f}",
        ci=f"r = 0.5 to 1.0: {bf_lo[0.5]:.2f} to {bf_lo[1.0]:.2f}",
        p="",
        note="BF01 > 3: moderate evidence for no difference",
    )
)
ah = hi_[Y].dropna().to_numpy()
k2.append(
    dict(
        analysis="Bayes factor (JZS, r = 0.707)",
        contrast="Higher vs HC",
        n_a=len(ah),
        n_b=len(b),
        estimate=f"BF10 = {1 / jzs_bf01(ah, b):.0f}",
        ci="",
        p="",
        note="",
    )
)
k2 = pd.DataFrame(k2)
k2.to_csv(OUT / "04_k2_primary.csv", index=False, encoding="utf-8-sig")
P("## 4. K = 2 當主分析時的結果\n")
P(k2.to_markdown(index=False))
P("")

# 5. 交叉表
ct23 = pd.crosstab(ins["k3"], ins["k2"])
ct24 = pd.crosstab(ins["k4"], ins["k2"])
ct = pd.concat([ct23.assign(table="K=3 profile × K=2"), ct24.assign(table="K=4 × K=2")])
ct.to_csv(OUT / "05_crosstabs.csv", encoding="utf-8-sig")
P("## 5. 交叉表\n")
P("K = 3 輪廓 × K = 2\n")
P(ct23.to_markdown())
P("")
P("K = 4 × K = 2\n")
P(ct24.to_markdown())
P("")

# 6. K = 2 下其他腦波指標
FEATS = [
    "beta_pw_BRAIN",
    "beta_pw_CENTRAL",
    "beta_pw_MIDLINE",
    "Beta_rel_BRAIN",
    "Beta_rel_CENTRAL",
    "exponent_BRAIN",
    "exponent_CENTRAL",
    "offset_BRAIN",
    "offset_CENTRAL",
    "alpha_pw_BRAIN",
    "theta_pw_BRAIN",
]
oe = []
for f in FEATS:
    for label, dd in [("Higher", hi_), ("Lower", lo_)]:
        g, l, h, n1, n2 = hedges_g(dd[f], hc[f])
        oe.append(
            dict(
                feature=f,
                contrast=f"{label} vs HC",
                n_a=n1,
                n_b=n2,
                g=round(g, 3),
                ci_low=round(l, 3),
                ci_high=round(h, 3),
                MWU_p=round(mwu_p(dd[f], hc[f]), 4),
            )
        )
oe = pd.DataFrame(oe)
oe["FDR_q"] = multipletests(oe["MWU_p"], method="fdr_bh")[1].round(4)
oe.to_csv(OUT / "06_k2_other_eeg.csv", index=False, encoding="utf-8-sig")
P("## 6. K = 2 下的其他腦波指標（FDR 為本表 22 個對比內的 BH 校正）\n")
P(oe.to_markdown(index=False))
P("")
P("數字全部由本腳本即時計算；CSV 為機讀來源。尚未寫回 SSOT。")

(OUT / "07_summary.md").write_text("\n".join(MD), encoding="utf-8")
print(f"\n輸出 → {OUT}")
