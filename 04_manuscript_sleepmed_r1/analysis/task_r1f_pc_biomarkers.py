"""把 PC1–PC3 同時放進模型，對照 16 個 EEG 指標。

各 PC 內做 BH 校正，另列全樣本、MEW 對 MDI 及七題同時入模的結果。
修訂稿只放兩個保留成分的版本在 task_r1w/09。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

from _r1_common import AN, COV, COV_COLS, ISI_ITEMS, OUTCOME, hedges_g, load_all, mwu_p

OUT = AN / "S_R1f_pc_biomarkers"
OUT.mkdir(parents=True, exist_ok=True)
MD = []
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


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    MD.append(s)


def add_pcs(d):
    d = d.copy()
    Z = StandardScaler().fit_transform(d[ISI_ITEMS])
    pca = PCA().fit(Z)
    S, comp = pca.transform(Z), pca.components_.copy()
    if np.corrcoef(S[:, 0], d["ISI_T1"])[0, 1] < 0:
        S[:, 0] *= -1
        comp[0] *= -1
    if comp[1][1] + comp[1][2] < 0:
        S[:, 1] *= -1
        comp[1] *= -1
    for i in range(3):
        d[f"PC{i + 1}"] = S[:, i] / S[:, i].std(ddof=1)  # 標準化，係數可跨 PC 比較
    info = dict(
        evr=pca.explained_variance_ratio_[:3].round(3).tolist(),
        r_total=[round(float(np.corrcoef(d[f"PC{i}"], d["ISI_T1"])[0, 1]), 3) for i in (1, 2, 3)],
        loadings=pd.DataFrame(
            comp[:3].T,
            index=[q.replace("ISI_T1_", "") for q in ISI_ITEMS],
            columns=["PC1", "PC2", "PC3"],
        ).round(2),
    )
    return d, info


def pc_table(d):
    rows = []
    for f in BIO:
        dd = d.dropna(subset=[f] + COV_COLS).copy()
        dd["yz"] = (dd[f] - dd[f].mean()) / dd[f].std(ddof=1)
        m = smf.ols(f"yz ~ {COV} + PC1 + PC2 + PC3", dd).fit()
        r = dict(biomarker=f, n=int(m.nobs))
        for pc in ["PC1", "PC2", "PC3"]:
            r[f"{pc}_beta"] = round(m.params[pc], 3)
            r[f"{pc}_p"] = round(m.pvalues[pc], 4)
        rows.append(r)
    t = pd.DataFrame(rows)
    for pc in ["PC1", "PC2", "PC3"]:
        t[f"{pc}_q"] = multipletests(t[f"{pc}_p"], method="fdr_bh")[1].round(4)
    return t[
        [
            "biomarker",
            "n",
            "PC1_beta",
            "PC1_p",
            "PC1_q",
            "PC2_beta",
            "PC2_p",
            "PC2_q",
            "PC3_beta",
            "PC3_p",
            "PC3_q",
        ]
    ]


df = load_all()
ins, info_i = add_pcs(df[df["group"] != "HC"])
full, info_f = add_pcs(df)
P("# S_R1f — 症狀主成分 × 腦波指標\n")
for label, d, info, fn in [
    ("失眠組 PCA", ins, info_i, "01_pc_by_biomarker_insomnia.csv"),
    ("全樣本 PCA（含 HC）", full, info_f, "02_pc_by_biomarker_full.csv"),
]:
    t = pc_table(d)
    t.to_csv(OUT / fn, index=False, encoding="utf-8-sig")
    P(f"## {label}（n = {len(d)}）\n")
    P(f"解釋變異 PC1–PC3 = {info['evr']}；與 ISI 總分的相關 = {info['r_total']}\n")
    P("負荷量\n")
    P(info["loadings"].to_markdown())
    P("")
    P(
        "標準化係數（指標與 PC 皆為 z 分數；三個 PC 同時進模型；q 為每個 PC 在 16 個指標內的 BH-FDR）\n"
    )
    P(t.to_markdown(index=False))
    P("")

# MEW 對 MDI
mew, mdi = ins[ins["group"] == "MEW"], ins[ins["group"] == "MDI"]
rows = [
    dict(
        feature="ISI total",
        n_MEW=len(mew),
        n_MDI=len(mdi),
        mean_MEW=round(mew["ISI_T1"].mean(), 2),
        mean_MDI=round(mdi["ISI_T1"].mean(), 2),
        g=round(hedges_g(mew["ISI_T1"], mdi["ISI_T1"])[0], 3),
        MWU_p=round(mwu_p(mew["ISI_T1"], mdi["ISI_T1"]), 4),
    )
]
for f in ["PC1", "PC2", "PC3"] + BIO:
    g, lo, hi, n1, n2 = hedges_g(mew[f], mdi[f])
    rows.append(
        dict(
            feature=f,
            n_MEW=n1,
            n_MDI=n2,
            mean_MEW=round(mew[f].mean(), 3),
            mean_MDI=round(mdi[f].mean(), 3),
            g=round(g, 3),
            ci_low=round(lo, 3),
            ci_high=round(hi, 3),
            MWU_p=round(mwu_p(mew[f], mdi[f]), 4),
        )
    )
mm = pd.DataFrame(rows)
mm.to_csv(OUT / "03_mew_vs_mdi.csv", index=False, encoding="utf-8-sig")
P("## MEW 對 MDI 的直接比較（g 為 MEW − MDI）\n")
P(mm.to_markdown(index=False))
P("")

# 七題同時
dd = ins.dropna(subset=[OUTCOME] + COV_COLS)
m = smf.ols(f"{OUTCOME} ~ {COV} + " + " + ".join(ISI_ITEMS), dd).fit()
m0 = smf.ols(f"{OUTCOME} ~ {COV} + ISI_T1", dd).fit()
it = pd.DataFrame(
    [
        dict(
            item=q.replace("ISI_T1_", ""),
            b=round(m.params[q], 4),
            se=round(m.bse[q], 4),
            t=round(m.tvalues[q], 3),
            p=round(m.pvalues[q], 4),
        )
        for q in ISI_ITEMS
    ]
)
it.to_csv(OUT / "04_items_simultaneous.csv", index=False, encoding="utf-8-sig")
P(f"## 七個題項同時預測 periodic beta（失眠組 n = {int(m.nobs)}）\n")
P(
    f"七題模型 R² = {m.rsquared:.3f}、AIC = {m.aic:.1f}；只放總分 R² = {m0.rsquared:.3f}、AIC = {m0.aic:.1f}\n"
)
P(it.to_markdown(index=False))
P("")
P("數字全部由本腳本即時計算；尚未寫回 SSOT。")
(OUT / "05_summary.md").write_text("\n".join(MD), encoding="utf-8")
print(f"\n輸出 → {OUT}")
