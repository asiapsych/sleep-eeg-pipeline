"""檢查 beta 的差異來自症狀嚴重度，還是嚴重度以外的症狀組態。

ISI 題項 PCA、分數帶、函數形式及同分數帶比較都從同一批樣本出發。
切點是掃描後挑的，置換檢定與 bootstrap 也要重做選點，不能當成預先指定。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from _r1_common import (
    AN,
    COV,
    COV_COLS,
    ISI_ITEMS,
    OUTCOME,
    fmt_p,
    g_only,
    hedges_g,
    load_all,
    mwu_p,
)

Y = OUTCOME
SEED = 20261112
B = 5000
OUT = AN / "S_R1a_severity_form"
OUT.mkdir(parents=True, exist_ok=True)
MD = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    MD.append(s)


df = load_all()
hc = df[df["group"] == "HC"].copy()
ins_all = df[df["group"] != "HC"].copy()  # n = 165（分群樣本）
ins = ins_all.dropna(subset=[Y] + COV_COLS).copy()  # n = 164（有 beta 峰）
hcv = hc[Y].dropna()
P(f"# S_R1a — 嚴重度形狀與症狀組態\n")
P(
    f"N = {len(df)}；HC = {len(hc)}（beta 有值 {len(hcv)}）；失眠 = {len(ins_all)}（beta 有值 {len(ins)}）\n"
)

# 1. PCA
Z = StandardScaler().fit_transform(ins_all[ISI_ITEMS])
pca = PCA().fit(Z)
S = pca.transform(Z)
comp = pca.components_.copy()
if np.corrcoef(S[:, 0], ins_all["ISI_T1"])[0, 1] < 0:  # PC1 方向：分數高＝嚴重
    S[:, 0] *= -1
    comp[0] *= -1
if comp[1][1] + comp[1][2] < 0:  # PC2 方向：正＝維持／早醒（Q2、Q3）
    S[:, 1] *= -1
    comp[1] *= -1
for i in range(3):
    ins_all[f"PC{i + 1}"] = S[:, i]
load = pd.DataFrame(comp[:3].T, index=ISI_ITEMS, columns=["PC1", "PC2", "PC3"]).round(3)
load.loc["explained_variance_ratio"] = pca.explained_variance_ratio_[:3].round(3)
load.loc["r_with_ISI_total"] = [
    round(float(np.corrcoef(ins_all[f"PC{i}"], ins_all["ISI_T1"])[0, 1]), 3) for i in (1, 2, 3)
]
for gname in ["MEW", "MDI", "SCA"]:
    load.loc[f"mean_{gname}"] = (
        ins_all.loc[ins_all["group"] == gname, ["PC1", "PC2", "PC3"]].mean().round(2).values
    )
load.to_csv(OUT / "01_pca_loadings.csv", encoding="utf-8-sig")
P("## 1. PCA（失眠組 n = %d，ISI 七題 z 分數）\n" % len(ins_all))
P(load.to_markdown())
P("")

ins = ins.merge(ins_all[["subject", "PC1", "PC2", "PC3"]], on="subject", how="left")

# 2. beta ~ PC
rows = []
for label, f in [
    ("with covariates", f"{Y} ~ {COV} + PC1 + PC2 + PC3"),
    ("unadjusted", f"{Y} ~ PC1 + PC2 + PC3"),
]:
    m = smf.ols(f, ins).fit()
    sdy = ins[Y].std(ddof=1)
    for pc in ["PC1", "PC2", "PC3"]:
        rows.append(
            dict(
                model=label,
                n=int(m.nobs),
                term=pc,
                b=round(m.params[pc], 4),
                se=round(m.bse[pc], 4),
                std_beta=round(m.params[pc] * ins[pc].std(ddof=1) / sdy, 3),
                t=round(m.tvalues[pc], 3),
                p=round(m.pvalues[pc], 4),
                R2_model=round(m.rsquared, 3),
            )
        )
m1 = smf.ols(f"{Y} ~ {COV} + PC1", ins).fit()
m3 = smf.ols(f"{Y} ~ {COV} + PC1 + PC2 + PC3", ins).fit()
a = anova_lm(m1, m3)
rows.append(
    dict(
        model="PC2+PC3 beyond PC1 (nested F)",
        n=int(m3.nobs),
        term="PC2+PC3",
        b="",
        se="",
        std_beta="",
        t=f"F({int(a['df_diff'].iloc[1])},{int(a['df_resid'].iloc[1])})={a['F'].iloc[1]:.2f}",
        p=round(float(a["Pr(>F)"].iloc[1]), 4),
        R2_model=round(m3.rsquared - m1.rsquared, 4),
    )
)
pcreg = pd.DataFrame(rows)
pcreg.to_csv(OUT / "02_pc_regression.csv", index=False, encoding="utf-8-sig")
P("## 2. beta 對主成分的迴歸（失眠組 n = %d）\n" % len(ins))
P(pcreg.to_markdown(index=False))
P("")


# 3. ISI 分段
def bin_table(bins, scheme):
    out = []
    for lo, hi in bins:
        s = ins[(ins["ISI_T1"] >= lo) & (ins["ISI_T1"] <= hi)]
        g, l, h, n1, _ = hedges_g(s[Y], hcv)
        comp_ = s["group"].value_counts().to_dict()
        out.append(
            dict(
                scheme=scheme,
                ISI_range=f"{lo}-{hi}" if hi < 28 else f">={lo}",
                n=n1,
                mean_beta=round(s[Y].mean(), 3),
                g_vs_HC=round(g, 3),
                ci_low=round(l, 3),
                ci_high=round(h, 3),
                MWU_p=round(mwu_p(s[Y], hcv), 4),
                n_MEW=comp_.get("MEW", 0),
                n_MDI=comp_.get("MDI", 0),
                n_SCA=comp_.get("SCA", 0),
            )
        )
    return out


bins = pd.DataFrame(
    bin_table([(9, 10), (11, 12), (13, 14), (15, 16), (17, 18), (19, 28)], "2-point bins")
    + bin_table([(9, 14), (15, 21), (22, 28)], "conventional ISI bands")
    + bin_table([(9, 28)], "all insomnia (pooled)")
)
bins.to_csv(OUT / "03_isi_bins.csv", index=False, encoding="utf-8-sig")
P(
    "## 3. 各 ISI 分段對 HC 的 Hedges' g（HC n = %d，HC 平均 beta = %.3f）\n"
    % (len(hcv), hcv.mean())
)
P(bins.to_markdown(index=False))
P("")


# 4. 函數形式比較
def scan_threshold(d, cands, min_n=10):
    best = None
    for c in cands:
        hi = d["ISI_T1"] >= c
        if hi.sum() < min_n or (~hi).sum() < min_n:
            continue
        m = smf.ols(f"{Y} ~ {COV} + thr", d.assign(thr=hi.astype(int))).fit()
        F = float(m.tvalues["thr"] ** 2)
        if best is None or F > best[1]:
            best = (c, F, m)
    return best


CANDS = list(range(11, 23))
c_hat, F_hat, m_thr = scan_threshold(ins, CANDS)
ins["thr_hat"] = (ins["ISI_T1"] >= c_hat).astype(int)
ins["thr15"] = (ins["ISI_T1"] >= 15).astype(int)
ins["k2_high"] = (ins["k2"] == "Higher").astype(int)
ins["sca"] = (ins["group"] == "SCA").astype(int)

forms = {
    "covariates only": COV,
    "linear ISI total": f"{COV} + ISI_T1",
    "quadratic ISI total": f"{COV} + ISI_T1 + I(ISI_T1**2)",
    "spline df=3": f"{COV} + cr(ISI_T1, df=3)",
    "spline df=4": f"{COV} + cr(ISI_T1, df=4)",
    "spline df=5": f"{COV} + cr(ISI_T1, df=5)",
    f"single threshold ISI>={c_hat} (data-selected)": f"{COV} + thr_hat",
    "single threshold ISI>=15 (conventional)": f"{COV} + thr15",
    "K=2 higher-severity cluster": f"{COV} + k2_high",
    "K=3 profile (3 levels)": f"{COV} + C(group)",
    "SCA indicator": f"{COV} + sca",
    "PC1": f"{COV} + PC1",
    "PC1 + PC2": f"{COV} + PC1 + PC2",
}
mc = []
for k, f in forms.items():
    m = smf.ols(f"{Y} ~ {f}", ins).fit()
    extra = 1 if "data-selected" in k else 0  # 掃描得到的切點算一個額外參數
    mc.append(
        dict(
            scope="insomnia only",
            model=k,
            n=int(m.nobs),
            k_severity=int(m.df_model) - 3 + extra,
            R2=round(m.rsquared, 4),
            adj_R2=round(m.rsquared_adj, 4),
            AIC=round(m.aic + 2 * extra, 2),
            BIC=round(m.bic + np.log(m.nobs) * extra, 2),
        )
    )

# 全樣本（含 HC，ISI 0–28）：同一組比較，回應「連續向度」的看法
full = df.dropna(subset=[Y] + COV_COLS).copy()
cF, FF, _ = scan_threshold(full, list(range(5, 23)))
full["thr_hat"] = (full["ISI_T1"] >= cF).astype(int)
full["ins"] = (full["ISI_T1"] >= 9).astype(int)
full["k2_high"] = (full["k2"] == "Higher").astype(int)
full["sca"] = (full["group"] == "SCA").astype(int)
forms_full = {
    "covariates only": COV,
    "linear ISI total": f"{COV} + ISI_T1",
    "spline df=3": f"{COV} + cr(ISI_T1, df=3)",
    "spline df=4": f"{COV} + cr(ISI_T1, df=4)",
    "spline df=5": f"{COV} + cr(ISI_T1, df=5)",
    "insomnia vs HC (ISI>=9; pooled contrast)": f"{COV} + ins",
    f"single threshold ISI>={cF} (data-selected)": f"{COV} + thr_hat",
    "K=2 higher-severity cluster vs rest": f"{COV} + k2_high",
    "SCA vs rest": f"{COV} + sca",
    "4 groups (HC/MEW/MDI/SCA)": f"{COV} + C(group)",
}
for k, f in forms_full.items():
    m = smf.ols(f"{Y} ~ {f}", full).fit()
    extra = 1 if "data-selected" in k else 0
    mc.append(
        dict(
            scope="full sample (HC included)",
            model=k,
            n=int(m.nobs),
            k_severity=int(m.df_model) - 3 + extra,
            R2=round(m.rsquared, 4),
            adj_R2=round(m.rsquared_adj, 4),
            AIC=round(m.aic + 2 * extra, 2),
            BIC=round(m.bic + np.log(m.nobs) * extra, 2),
        )
    )
mc = pd.DataFrame(mc)
mc.to_csv(OUT / "04_model_comparison.csv", index=False, encoding="utf-8-sig")
P("## 4. 嚴重度函數形式比較（AIC／BIC 越小越好；資料選出的切點已多算一個參數）\n")
P(mc.to_markdown(index=False))
P("")

# 門檻：互相加進去還剩多少、置換檢定、bootstrap 位置
thr_rows = []
m_lin = smf.ols(f"{Y} ~ {COV} + ISI_T1", ins).fit()
m_stp = smf.ols(f"{Y} ~ {COV} + thr_hat", ins).fit()
m_both = smf.ols(f"{Y} ~ {COV} + ISI_T1 + thr_hat", ins).fit()
thr_rows.append(
    dict(
        quantity=f"data-selected threshold (scan ISI {CANDS[0]}–{CANDS[-1]}, >=10 per side)",
        value=c_hat,
        note=f"n above = {int(ins['thr_hat'].sum())}, n below = {int((1 - ins['thr_hat']).sum())}",
    )
)
thr_rows.append(
    dict(
        quantity="step coefficient (b, SE, p) in step-only model",
        value=f"{m_stp.params['thr_hat']:.3f}, {m_stp.bse['thr_hat']:.3f}, {fmt_p(m_stp.pvalues['thr_hat'])}",
        note="p 未校正切點選擇，見下方置換檢定",
    )
)
thr_rows.append(
    dict(
        quantity="step added to linear model: p for step",
        value=fmt_p(m_both.pvalues["thr_hat"]),
        note=f"b = {m_both.params['thr_hat']:.3f}",
    )
)
thr_rows.append(
    dict(
        quantity="linear added to step model: p for linear ISI",
        value=fmt_p(m_both.pvalues["ISI_T1"]),
        note=f"b = {m_both.params['ISI_T1']:.4f}",
    )
)

rng = np.random.default_rng(SEED)
isi = ins["ISI_T1"].to_numpy()
null_F = np.empty(B)
for b in range(B):  # 置換 ISI 總分：虛無＝嚴重度與 beta 無關
    d = ins.assign(ISI_T1=rng.permutation(isi))
    r = scan_threshold(d, CANDS)
    null_F[b] = r[1] if r else np.nan
p_perm = (np.nansum(null_F >= F_hat) + 1) / (np.sum(~np.isnan(null_F)) + 1)
thr_rows.append(
    dict(
        quantity="max-scanned step F: observed",
        value=round(F_hat, 3),
        note=f"permutation p (B = {B}, corrects for cut-point selection) = {p_perm:.4f}",
    )
)

boot_c = np.empty(B)
n = len(ins)
for b in range(B):
    d = ins.iloc[rng.integers(0, n, n)]
    r = scan_threshold(d, CANDS)
    boot_c[b] = r[0] if r else np.nan
bc = pd.Series(boot_c).dropna()
dist = bc.value_counts(normalize=True).sort_index()
thr_rows.append(
    dict(
        quantity="bootstrap distribution of selected threshold",
        value=f"median {bc.median():.0f}; 2.5–97.5 pct {bc.quantile(0.025):.0f}–{bc.quantile(0.975):.0f}",
        note="; ".join(f"{int(k)}: {v:.1%}" for k, v in dist.items()),
    )
)
thr_rows.append(
    dict(
        quantity="full-sample data-selected threshold (scan ISI 5–22)",
        value=cF,
        note=f"F = {FF:.2f}",
    )
)
thr = pd.DataFrame(thr_rows)
thr.to_csv(OUT / "05_threshold.csv", index=False, encoding="utf-8-sig")
P("### 門檻位置與檢定\n")
P(thr.to_markdown(index=False))
P("")

# 5. 區段內相關
seg = []


def seg_row(label, d):
    if len(d) < 8:
        return
    r, p = stats.spearmanr(d["ISI_T1"], d[Y])
    seg.append(
        dict(
            segment=label,
            n=len(d),
            ISI_min=int(d["ISI_T1"].min()),
            ISI_max=int(d["ISI_T1"].max()),
            spearman_r=round(r, 3),
            p=round(p, 4),
        )
    )


seg_row("all insomnia", ins)
seg_row(f"below data-selected threshold (ISI < {c_hat})", ins[ins["thr_hat"] == 0])
seg_row(f"at/above data-selected threshold (ISI >= {c_hat})", ins[ins["thr_hat"] == 1])
seg_row("K=2 lower-severity cluster", ins[ins["k2"] == "Lower"])
seg_row("K=2 higher-severity cluster", ins[ins["k2"] == "Higher"])
seg_row("non-SCA (MEW + MDI)", ins[ins["group"] != "SCA"])
seg_row("SCA", ins[ins["group"] == "SCA"])
seg_row("full sample (HC + insomnia)", full)
seg_row("HC only (ISI < 9)", full[full["group"] == "HC"])
seg = pd.DataFrame(seg)
seg.to_csv(OUT / "06_within_segment.csv", index=False, encoding="utf-8-sig")
P("## 5. 各區段內 ISI 總分與 beta 的 Spearman 相關\n")
P(seg.to_markdown(index=False))
P("")

# 6. 同分數帶：SCA 對非 SCA
ov = []


def overlap(label, mask_a, name_a, name_b):
    A = ins[mask_a]
    Bd = ins[~mask_a]
    lo, hi = int(A["ISI_T1"].min()), int(Bd["ISI_T1"].max())
    band = ins[(ins["ISI_T1"] >= lo) & (ins["ISI_T1"] <= hi)].copy()
    a = band[mask_a.loc[band.index]]
    b = band[~mask_a.loc[band.index]]
    if len(a) < 5 or len(b) < 5:
        ov.append(dict(comparison=label, band=f"{lo}-{hi}", note="重疊範圍內人數不足"))
        return
    g, l, h, n1, n2 = hedges_g(a[Y], b[Y])
    band["ind"] = mask_a.loc[band.index].astype(int)
    m = smf.ols(f"{Y} ~ {COV} + ISI_T1 + ind", band).fit()
    ov.append(
        dict(
            comparison=label,
            band=f"{lo}-{hi}",
            n_a=n1,
            n_b=n2,
            ISI_mean_a=round(a["ISI_T1"].mean(), 2),
            ISI_mean_b=round(b["ISI_T1"].mean(), 2),
            g=round(g, 3),
            ci_low=round(l, 3),
            ci_high=round(h, 3),
            MWU_p=round(mwu_p(a[Y], b[Y]), 4),
            b_adj_for_ISI=round(m.params["ind"], 3),
            p_adj_for_ISI=round(m.pvalues["ind"], 4),
            note=f"{name_a} vs {name_b}; band = [{name_a} 最低分, {name_b} 最高分]",
        )
    )


overlap("SCA vs non-SCA at overlapping ISI totals", ins["group"] == "SCA", "SCA", "non-SCA")
overlap("K=2 higher vs lower at overlapping ISI totals", ins["k2"] == "Higher", "Higher", "Lower")
ov = pd.DataFrame(ov)
ov.to_csv(OUT / "07_overlap_band.csv", index=False, encoding="utf-8-sig")
P("## 6. 同分數帶比較（分數帶＝兩邊實際重疊的 ISI 範圍）\n")
P(ov.to_markdown(index=False))
P("")

# 7. spread 差值 bootstrap（S11）
MAIN_SIZES = [68, 50, 47]  # 與 task_za 相同：低／中／高
sev = ins_all.sort_values(["ISI_T1", "subject"], kind="mergesort").reset_index(drop=True)
bounds = np.cumsum([0] + MAIN_SIZES)
sev["tert"] = ""
for i, nm in enumerate(["low", "mid", "high"]):
    sev.loc[bounds[i] : bounds[i + 1] - 1, "tert"] = nm
ia = ins_all.merge(sev[["subject", "tert"]], on="subject")


def spreads(d_ins, d_hc):
    gp = [g_only(d_ins.loc[d_ins["group"] == k, Y], d_hc) for k in ["MEW", "MDI", "SCA"]]
    gs = [g_only(d_ins.loc[d_ins["tert"] == k, Y], d_hc) for k in ["low", "mid", "high"]]
    return (np.nanmax(gp) - np.nanmin(gp), np.nanmax(gs) - np.nanmin(gs), gp[2], gs[2])


obs = spreads(ia, hc[Y])
bs = np.empty((B, 4))
hy = hc[Y].to_numpy()
for b in range(B):
    bs[b] = spreads(
        ia.iloc[rng.integers(0, len(ia), len(ia))], hy[rng.integers(0, len(hy), len(hy))]
    )
d_spread = bs[:, 0] - bs[:, 1]
d_top = bs[:, 2] - bs[:, 3]
sp = pd.DataFrame(
    [
        dict(quantity="spread, K=3 profiles (max g − min g)", observed=round(obs[0], 3)),
        dict(quantity="spread, severity tertiles of equal sizes", observed=round(obs[1], 3)),
        dict(
            quantity="difference in spread (profiles − tertiles)",
            observed=round(obs[0] - obs[1], 3),
            boot_ci_low=round(np.percentile(d_spread, 2.5), 3),
            boot_ci_high=round(np.percentile(d_spread, 97.5), 3),
            prop_boot_gt0=round(float(np.mean(d_spread > 0)), 3),
        ),
        dict(
            quantity="g of SCA − g of highest-severity tertile",
            observed=round(obs[2] - obs[3], 3),
            boot_ci_low=round(np.percentile(d_top, 2.5), 3),
            boot_ci_high=round(np.percentile(d_top, 97.5), 3),
            prop_boot_gt0=round(float(np.mean(d_top > 0)), 3),
        ),
    ]
)
sp.to_csv(OUT / "08_spread_bootstrap.csv", index=False, encoding="utf-8-sig")
P(f"## 7. 兩種切法的 spread 差值（bootstrap B = {B}，受試者層級重抽，標籤固定）\n")
P(sp.to_markdown(index=False))
P("")
P(f"種子 {SEED}。數字全部由本腳本即時計算；CSV 為機讀來源。尚未寫回 SSOT。")

(OUT / "09_summary.md").write_text("\n".join(MD), encoding="utf-8")
print(f"\n輸出 → {OUT}")
