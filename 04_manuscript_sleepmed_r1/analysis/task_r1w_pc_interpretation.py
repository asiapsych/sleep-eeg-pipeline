"""用題項負荷量與外部變項解讀 ISI 的第二、第三主成分。

13 個外部變項全部列出，不按顯著性挑選；各成分內做 BH 校正。
PC1 朝總分增加、PC2 朝第 2／3 題增加、PC3 朝第 6 題增加。
成分正負號任意，比較係數大小時用絕對值，不能直接相減。"""

from __future__ import annotations
import sqlite3
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

from _dynamic_n import DB_PATH, ISI_ITEMS
from _r1_common import AN, COV, OUTCOME, load_all, fmt_p

OUT = AN / "S_R1w_pc_interpretation"
OUT.mkdir(parents=True, exist_ok=True)
ITEM = [
    "1 Difficulty falling asleep",
    "2 Difficulty staying asleep",
    "3 Early-morning awakening",
    "4 Dissatisfaction with sleep",
    "5 Interference with daytime functioning",
    "6 Noticeability to others",
    "7 Worry or distress",
]
PSQI = [f"PSQI_T1_C{i}" for i in range(1, 8)]
EXT = {
    "AGE": "Age",
    "SEX": "Sex (0 = female, 1 = male)",
    "EDU": "Education (years)",
    "BDI_T1": "BDI-II",
    "BAI_T1": "BAI",
    "DBAS": "DBAS (dysfunctional beliefs about sleep)",
    "PSQI_T1_C1": "PSQI subjective sleep quality",
    "PSQI_T1_C2": "PSQI sleep latency",
    "PSQI_T1_C3": "PSQI sleep duration (higher = shorter)",
    "PSQI_T1_C4": "PSQI habitual sleep efficiency",
    "PSQI_T1_C5": "PSQI sleep disturbances",
    "PSQI_T1_C6": "PSQI use of sleep medication",
    "PSQI_T1_C7": "PSQI daytime dysfunction",
}

df = load_all()
with sqlite3.connect(DB_PATH) as con:
    x = pd.read_sql(
        "SELECT NEW_NUMBER AS subject, DBAS, "
        + ", ".join(PSQI)
        + " FROM isi_raw_data WHERE NEW_NUMBER IS NOT NULL",
        con,
    )
df = df.merge(x, on="subject", how="left")
ins = df[df["ISI_T1"] >= 9].copy()
X = StandardScaler().fit_transform(SimpleImputer(strategy="median").fit_transform(ins[ISI_ITEMS]))
pca = PCA(3).fit(X)
S = pca.transform(X)
W = pca.components_.copy()
for j, ref in enumerate([ins["ISI_T1"].to_numpy(float), X[:, 1] + X[:, 2], X[:, 5]]):
    if np.corrcoef(S[:, j], ref)[0, 1] < 0:
        S[:, j] *= -1
        W[j] *= -1
PCS = ["PC1", "PC2", "PC3"]
ins[PCS] = S

t1 = pd.DataFrame(W.T, index=ITEM, columns=PCS).round(3)
t1.loc["explained variance"] = pca.explained_variance_ratio_.round(3)
t1.loc["r with ISI total"] = [
    round(float(np.corrcoef(S[:, j], ins["ISI_T1"])[0, 1]), 3) for j in range(3)
]
t1.to_csv(OUT / "01_loadings.csv", encoding="utf-8-sig")

ins["pc2_tertile"] = pd.qcut(
    ins["PC2"], 3, labels=["lowest third", "middle third", "highest third"]
)
t2 = ins.groupby("pc2_tertile", observed=True)[ISI_ITEMS + ["ISI_T1"]].mean().round(2)
t2.columns = ITEM + ["ISI total"]
t2.insert(0, "n", ins["pc2_tertile"].value_counts().reindex(t2.index).to_numpy())
t2.to_csv(OUT / "02_items_by_pc2_tertile.csv", encoding="utf-8-sig")

rows = []
for v, name in EXT.items():
    d = ins.dropna(subset=[v])
    r = dict(variable=name, n=len(d))
    for pc in PCS:
        rho, p = stats.spearmanr(d[pc], pd.to_numeric(d[v]))
        r[f"rho_{pc}"] = rho
        r[f"p_{pc}"] = p
    rows.append(r)
t3 = pd.DataFrame(rows)
for pc in PCS:
    t3[f"q_{pc}"] = multipletests(t3[f"p_{pc}"], method="fdr_bh")[1]
t3.to_csv(OUT / "03_external_correlates.csv", index=False, encoding="utf-8-sig")

z = ins.assign(**{c: (ins[c] - ins[c].mean()) / ins[c].std(ddof=1) for c in PCS})
m = smf.ols(f"{OUTCOME} ~ {COV} + PC1 + PC2 + PC3", z).fit()
sd = ins[OUTCOME].std(ddof=1)
t4 = pd.DataFrame(
    [
        dict(
            component=pc,
            n=int(m.nobs),
            std_beta=m.params[pc] / sd,
            ci_low=m.conf_int().loc[pc, 0] / sd,
            ci_high=m.conf_int().loc[pc, 1] / sd,
            p=m.pvalues[pc],
        )
        for pc in PCS
    ]
)
t4.to_csv(OUT / "04_pc_vs_beta.csv", index=False, encoding="utf-8-sig")

L = [
    "# ISI 題項的第二、第三主成分是什麼（task_r1w；描述性）",
    "",
    f"樣本：失眠症狀者 {len(ins)} 人。解釋變異：PC1 {pca.explained_variance_ratio_[0]:.1%}、PC2 {pca.explained_variance_ratio_[1]:.1%}、PC3 {pca.explained_variance_ratio_[2]:.1%}。",
    "",
    "## 負荷量",
    "",
    "| 題項 | PC1 | PC2 | PC3 |",
    "|---|---|---|---|",
]
for it in ITEM:
    L.append(
        f"| {it} | {t1.loc[it, 'PC1']:+.2f} | {t1.loc[it, 'PC2']:+.2f} | {t1.loc[it, 'PC3']:+.2f} |"
    )
L += [
    "",
    "## 與 ISI 以外變項的 Spearman 相關（q 為各成分內 13 個變項的 BH 校正）",
    "",
    "| 變項 | n | PC1 ρ | PC2 ρ | PC2 q | PC3 ρ | PC3 q |",
    "|---|---|---|---|---|---|---|",
]
for _, r in t3.iterrows():
    L.append(
        f"| {r['variable']} | {r['n']} | {r['rho_PC1']:+.2f} | {r['rho_PC2']:+.2f} | {fmt_p(r['q_PC2'])} | {r['rho_PC3']:+.2f} | {fmt_p(r['q_PC3'])} |"
    )
L += [
    "",
    "## 三個主成分對週期性 beta（同時放入；校正年齡、性別、教育年數）",
    "",
    "| 成分 | 標準化係數 [95% CI] | p |",
    "|---|---|---|",
]
for _, r in t4.iterrows():
    L.append(
        f"| {r['component']} | {r['std_beta']:+.2f} [{r['ci_low']:+.2f}, {r['ci_high']:+.2f}] | {fmt_p(r['p'])} |"
    )
# 穩健性
rng = np.random.default_rng(20261117)
n, k = X.shape
ev = PCA(k).fit(X).explained_variance_
perm = np.array(
    [
        PCA(k)
        .fit(np.column_stack([rng.permutation(X[:, j]) for j in range(k)]))
        .explained_variance_
        for _ in range(2000)
    ]
)
thr95 = np.percentile(perm, 95, axis=0)
cong = {j: [] for j in range(3)}
for _ in range(2000):
    i = rng.integers(0, n, n)
    Wb = PCA(3).fit(StandardScaler().fit_transform(X[i])).components_
    for j in range(3):
        cong[j].append(abs(W[j] @ Wb[j]) / np.sqrt((W[j] @ W[j]) * (Wb[j] @ Wb[j])))
dif = m.t_test("PC1 - PC2 = 0")
dif_r = m.get_robustcov_results(cov_type="HC3").t_test("PC1 - PC2 = 0")
ins["night_minus_day"] = ins[ISI_ITEMS[1:3]].mean(axis=1) - ins[ISI_ITEMS[4:7]].mean(axis=1)
m_nd = smf.ols(f"{OUTCOME} ~ {COV} + ISI_T1 + night_minus_day", ins).fit()
rob = [
    dict(
        check=f"PC{j + 1} eigenvalue",
        value=ev[j],
        reference=f"95th percentile of permuted data = {thr95[j]:.2f}",
    )
    for j in range(3)
]
rob += [
    dict(
        check=f"PC{j + 1} bootstrap congruence (median)",
        value=float(np.median(cong[j])),
        reference=f"5th percentile = {np.percentile(cong[j], 5):.2f}",
    )
    for j in range(3)
]
rob += [
    dict(
        check="PC1 minus PC2 (standardised coefficients on periodic beta)",
        value=float(dif.effect[0]) / sd,
        reference=f"p = {float(dif.pvalue):.4f}; HC3 p = {float(dif_r.pvalue):.4f}",
    ),
    dict(
        check="simple index (items 2-3 minus items 5-7): r with PC2",
        value=float(np.corrcoef(ins["night_minus_day"], ins["PC2"])[0, 1]),
        reference="",
    ),
    dict(
        check="simple index -> periodic beta, given ISI total (standardised)",
        value=m_nd.params["night_minus_day"] * ins["night_minus_day"].std(ddof=1) / sd,
        reference=f"p = {m_nd.pvalues['night_minus_day']:.3f}",
    ),
]
for v in ["SEX", "PSQI_T1_C3", "DBAS", "PSQI_T1_C7"]:
    d = ins.dropna(subset=[v])
    rho, pv = stats.spearmanr(d["night_minus_day"], pd.to_numeric(d[v]))
    rob.append(
        dict(check=f"simple index vs {EXT[v]} (Spearman)", value=rho, reference=f"p = {pv:.4f}")
    )
t6 = pd.DataFrame(rob)
t6.to_csv(OUT / "06_robustness.csv", index=False, encoding="utf-8-sig")
L += ["", "## 這兩個向度穩不穩", "", "| 檢查 | 數值 | 參照 |", "|---|---|---|"]
for _, r in t6.iterrows():
    L.append(f"| {r['check']} | {r['value']:+.2f} | {r['reference']} |")
(OUT / "05_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))


# 2026-11-17 補：獨立審視後加的三項（方向不變的比較、只放保留的兩個成分、校正性別與年齡）
# 主成分的正負號是任意的，「PC1 的係數減 PC2 的係數」會隨 PC2 的方向改變，不能當成兩者大小的檢定。
# 改報：(a) 只放平行分析保留的兩個成分的模型；(b) 兩個係數的絕對值相減，以 bootstrap 估計區間（不受方向影響）；
#       (c) 第二個成分與外部變項的關係，校正性別與年齡後是否還在（性別與第二個成分相關 .36）。
m2c = smf.ols(f"{OUTCOME} ~ {COV} + PC1 + PC2", z).fit()
r2c = m2c.get_robustcov_results(cov_type="HC3", use_t=True)
nm2 = list(m2c.params.index)
t8 = pd.DataFrame(
    [
        dict(
            component=pc,
            model="two retained components",
            n=int(m2c.nobs),
            std_beta=m2c.params[pc] / sd,
            ci_low=m2c.conf_int().loc[pc, 0] / sd,
            ci_high=m2c.conf_int().loc[pc, 1] / sd,
            p=m2c.pvalues[pc],
            p_hc3=r2c.pvalues[nm2.index(pc)],
        )
        for pc in ["PC1", "PC2"]
    ]
)
rng2 = np.random.default_rng(20261117)
ab = []
for _ in range(4000):
    bs = z.iloc[rng2.integers(0, len(z), len(z))]
    mb = smf.ols(f"{OUTCOME} ~ {COV} + PC1 + PC2", bs).fit()
    ab.append((abs(mb.params["PC1"]) - abs(mb.params["PC2"])) / sd)
ab = np.array(ab)
absdiff = dict(
    component="|PC1| minus |PC2|",
    model="two retained components; bootstrap (4,000), sign-invariant",
    n=int(m2c.nobs),
    std_beta=float(np.median(ab)),
    ci_low=float(np.percentile(ab, 2.5)),
    ci_high=float(np.percentile(ab, 97.5)),
    p=np.nan,
    p_hc3=np.nan,
)
t8 = pd.concat([t8, pd.DataFrame([absdiff])], ignore_index=True)
t8.to_csv(OUT / "07_two_component_model.csv", index=False, encoding="utf-8-sig")

adj = []
for v in ["PSQI_T1_C3", "DBAS", "PSQI_T1_C7", "BDI_T1", "PSQI_T1_C2"]:
    d = ins.dropna(subset=[v]).copy()
    d["v"] = (pd.to_numeric(d[v]) - pd.to_numeric(d[v]).mean()) / pd.to_numeric(d[v]).std(ddof=1)
    d["pc2z"] = (d["PC2"] - d["PC2"].mean()) / d["PC2"].std(ddof=1)
    m0 = smf.ols("pc2z ~ v", d).fit()
    m1 = smf.ols("pc2z ~ v + C(SEX) + AGE", d).fit()
    adj.append(
        dict(
            variable=EXT[v],
            n=len(d),
            beta_unadjusted=m0.params["v"],
            p_unadjusted=m0.pvalues["v"],
            beta_adjusted_sex_age=m1.params["v"],
            ci_low=m1.conf_int().loc["v", 0],
            ci_high=m1.conf_int().loc["v", 1],
            p_adjusted=m1.pvalues["v"],
        )
    )
t9 = pd.DataFrame(adj)
t9.to_csv(OUT / "08_external_adjusted_sex_age.csv", index=False, encoding="utf-8-sig")

L += [
    "",
    "## 只放平行分析保留的兩個成分（主文用這一組）",
    "",
    "| 項目 | 標準化係數 [95% CI] | p | HC3 p |",
    "|---|---|---|---|",
]
for _, r in t8.iterrows():
    pp_ = "" if pd.isna(r["p"]) else fmt_p(r["p"])
    ph_ = "" if pd.isna(r["p_hc3"]) else fmt_p(r["p_hc3"])
    L.append(
        f"| {r['component']} | {r['std_beta']:+.2f} [{r['ci_low']:+.2f}, {r['ci_high']:+.2f}] | {pp_} | {ph_} |"
    )
L += [
    "",
    "主成分的正負號是任意的，所以「PC1 減 PC2」的檢定會隨方向改變（上表「這兩個向度穩不穩」裡的那一列不要引用）。",
    f"不受方向影響的比較是兩個係數絕對值的差：bootstrap 中 {np.mean(ab > 0):.1%} 的重抽樣為正。",
    "",
    "## 第二個成分與外部變項：校正性別與年齡後",
    "",
    "| 變項 | n | 未校正 β | 校正後 β [95% CI] | 校正後 p |",
    "|---|---|---|---|---|",
]
for _, r in t9.iterrows():
    L.append(
        f"| {r['variable']} | {r['n']} | {r['beta_unadjusted']:+.2f} | {r['beta_adjusted_sex_age']:+.2f} [{r['ci_low']:+.2f}, {r['ci_high']:+.2f}] | {fmt_p(r['p_adjusted'])} |"
    )
(OUT / "05_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L[-14:]))


# 2026-11-17 補二：兩個保留的成分對 16 個腦波指標
# task_r1f 的同一張表同時放了三個成分；平行分析只保留兩個，主文的模型（Table 2、Figure 1C）也只放兩個。
# 這裡用同一組指標、同樣的標準化與共變項，只放 PC1 與 PC2，並在各成分內重新做 BH 校正，補充表 S1(c) 用這一份。
from _r1_common import COV_COLS

BIO = [
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
    "alpha_cf_BRAIN",
    "theta_pw_BRAIN",
    "Theta_rel_BRAIN",
    "Alpha_rel_BRAIN",
    "Delta_rel_BRAIN",
    "Gamma_rel_BRAIN",
]
rows = []
for f in BIO:
    dd = z.dropna(subset=[f] + COV_COLS).copy()
    dd["yz"] = (dd[f] - dd[f].mean()) / dd[f].std(ddof=1)
    mm = smf.ols(f"yz ~ {COV} + PC1 + PC2", dd).fit()
    ci_ = mm.conf_int()
    r = dict(biomarker=f, n=int(mm.nobs))
    for pc in ("PC1", "PC2"):
        r.update(
            {
                f"{pc}_beta": mm.params[pc],
                f"{pc}_ci_low": ci_.loc[pc, 0],
                f"{pc}_ci_high": ci_.loc[pc, 1],
                f"{pc}_p": mm.pvalues[pc],
            }
        )
    rows.append(r)
t10 = pd.DataFrame(rows)
for pc in ("PC1", "PC2"):
    t10[f"{pc}_q"] = multipletests(t10[f"{pc}_p"], method="fdr_bh")[1]
t10.to_csv(OUT / "09_two_component_by_eeg.csv", index=False, encoding="utf-8-sig")
sig1 = t10[t10["PC1_q"] < 0.05]
L2 = [
    "",
    "## 兩個保留的成分對 16 個腦波指標（只放 PC1 與 PC2；各成分內 BH 校正）",
    "",
    f"- PC1：q < .05 的指標 {len(sig1)} 個（{'、'.join(sig1['biomarker'])}；q 最大 {sig1['PC1_q'].max():.3f}）。",
    f"- PC2：q 最小 {t10['PC2_q'].min():.3f}，沒有通過校正的指標；未校正 p 最小 {t10['PC2_p'].min():.3f}（{t10.loc[t10['PC2_p'].idxmin(), 'biomarker']}）。",
]
with open(OUT / "05_summary.md", "a", encoding="utf-8") as fh:
    fh.write("\n".join(L2) + "\n")
print("\n".join(L2))
