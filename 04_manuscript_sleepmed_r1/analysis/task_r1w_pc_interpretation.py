"""ISI 七題的主成分分析：負荷量、平行分析與重抽樣穩定度。

樣本為 ISI 總分 ≥ 9 者。成分正負號任意：
PC1 朝總分增加、PC2 朝第 2／3 題增加、PC3 朝第 6 題增加。"""

from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

from _dynamic_n import ISI_ITEMS
from _r1_common import AN, load_all

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
df = load_all()
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

L = [
    "# ISI 七題的主成分",
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
# 平行分析與重抽樣穩定度
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
t6 = pd.DataFrame(rob)
t6.to_csv(OUT / "06_robustness.csv", index=False, encoding="utf-8-sig")
L += ["", "## 平行分析與重抽樣穩定度", "", "| 檢查 | 數值 | 參照 |", "|---|---|---|"]
for _, r in t6.iterrows():
    L.append(f"| {r['check']} | {r['value']:+.2f} | {r['reference']} |")
(OUT / "05_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))
