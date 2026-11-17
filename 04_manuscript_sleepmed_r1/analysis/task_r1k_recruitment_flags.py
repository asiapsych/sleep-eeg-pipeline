"""比較初投與核對收案紀錄後的主要結果。

樣本排除另作幾種加嚴定義的對照；逐人名單只留內部，不進補充資料。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from _dynamic_n import ISI_ITEMS
from _r1_common import AN, COV, OUTCOME, SEED, hedges_g, load_all

OUT = AN / "S_R1k_recruitment_flags"


def one(d, label):
    d = d.copy()
    ins = d[d["ISI_T1"] >= 9].sort_values("subject")
    lab = KMeans(2, n_init=100, random_state=SEED).fit_predict(
        StandardScaler().fit_transform(ins[ISI_ITEMS])
    )
    hi = int(ins.groupby(lab)["ISI_T1"].mean().idxmax())
    d["k2"] = "HC"
    d.loc[ins.index, "k2"] = np.where(lab == hi, "Higher", "Lower")
    hc = d.loc[d["k2"] == "HC", OUTCOME].dropna()
    row = dict(
        sample=label,
        N=len(d),
        n_hc=int((d["k2"] == "HC").sum()),
        n_lower=int((d["k2"] == "Lower").sum()),
        n_higher=int((d["k2"] == "Higher").sum()),
    )
    for g in ["Higher", "Lower"]:
        gg, lo, hi_, _, _ = hedges_g(d.loc[d["k2"] == g, OUTCOME].dropna(), hc)
        row.update({f"{g}_g": gg, f"{g}_lo": lo, f"{g}_hi": hi_})
    d["thr18"] = (d["ISI_T1"] >= 18).astype(int)
    m = smf.ols(f"{OUTCOME} ~ thr18 + {COV}", d).fit()
    row.update(thr18_b=m.params["thr18"], thr18_p=m.pvalues["thr18"])
    m2 = smf.ols(f"{OUTCOME} ~ C(k2, Treatment('HC')) + {COV} + BDI_T1 + BAI_T1", d).fit()
    k = [c for c in m2.params.index if "Higher" in c][0]
    row.update(M2_higher_b=m2.params[k], M2_higher_p=m2.pvalues[k])
    return row


def main():
    rows = [
        one(load_all("submitted"), "Submitted sample"),
        one(load_all("audited"), "Revised sample (analytic sample)"),
    ]
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "01_sensitivity.csv", index=False)
    L = [
        "# 樣本定義改變前後的主結果（R1k）",
        "",
        "| 樣本 | N（健康／低／高） | 高嚴重度 g [95% CI] | 低嚴重度 g [95% CI] | ISI ≥ 18 係數（p） | M2 高嚴重度係數（p） |",
        "|---|---|---|---|---|---|",
    ]
    fp = lambda p: "< .001" if p < 0.001 else f"{p:.3f}".replace("0.", ".", 1)
    for r in t.itertuples():
        L.append(
            f"| {r.sample} | {r.N}（{r.n_hc}／{r.n_lower}／{r.n_higher}） | {r.Higher_g:+.2f} [{r.Higher_lo:+.2f}, {r.Higher_hi:+.2f}] | {r.Lower_g:+.2f} [{r.Lower_lo:+.2f}, {r.Lower_hi:+.2f}] "
            f"| {r.thr18_b:.3f}（{fp(r.thr18_p)}） | {r.M2_higher_b:.3f}（{fp(r.M2_higher_p)}） |"
        )
    (OUT / "02_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
