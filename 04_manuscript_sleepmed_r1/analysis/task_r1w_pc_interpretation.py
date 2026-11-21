"""Principal component analysis of the seven ISI items: loadings, parallel analysis
and bootstrap stability.

Sample: ISI >= 9. Component signs are arbitrary; PC1 is oriented towards a higher ISI
total, PC2 towards items 2 and 3, and PC3 towards item 6."""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from _dynamic_n import ISI_ITEMS
from _r1_common import AN, load_all, prepare_isi_items

OUT = AN / "S_R1w_pc_interpretation"
OUT.mkdir(parents=True, exist_ok=True)
ITEM_NAMES = [
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
_, X = prepare_isi_items(ins[ISI_ITEMS])
pca = PCA(3).fit(X)
scores = pca.transform(X)
loadings = pca.components_.copy()
for j, ref in enumerate([ins["ISI_T1"].to_numpy(float), X[:, 1] + X[:, 2], X[:, 5]]):
    if np.corrcoef(scores[:, j], ref)[0, 1] < 0:
        scores[:, j] *= -1
        loadings[j] *= -1
PCS = ["PC1", "PC2", "PC3"]
ins[PCS] = scores

loading_table = pd.DataFrame(loadings.T, index=ITEM_NAMES, columns=PCS).round(3)
loading_table.loc["explained variance"] = pca.explained_variance_ratio_.round(3)
loading_table.loc["r with ISI total"] = [
    round(float(np.corrcoef(scores[:, j], ins["ISI_T1"])[0, 1]), 3) for j in range(3)
]
loading_table.to_csv(OUT / "01_loadings.csv", encoding="utf-8-sig")

ins["pc2_tertile"] = pd.qcut(
    ins["PC2"], 3, labels=["lowest third", "middle third", "highest third"]
)
tertile_table = ins.groupby("pc2_tertile", observed=True)[ISI_ITEMS + ["ISI_T1"]].mean().round(2)
tertile_table.columns = ITEM_NAMES + ["ISI total"]
tertile_table.insert(
    0, "n", ins["pc2_tertile"].value_counts().reindex(tertile_table.index).to_numpy()
)
tertile_table.to_csv(OUT / "02_items_by_pc2_tertile.csv", encoding="utf-8-sig")

summary_lines = [
    "# Principal components of the seven ISI items",
    "",
    f"n = {len(ins)} participants with insomnia symptoms. Explained variance: PC1 {pca.explained_variance_ratio_[0]:.1%}, PC2 {pca.explained_variance_ratio_[1]:.1%}, PC3 {pca.explained_variance_ratio_[2]:.1%}.",
    "",
    "## Loadings",
    "",
    "| Item | PC1 | PC2 | PC3 |",
    "|---|---|---|---|",
]
for it in ITEM_NAMES:
    summary_lines.append(
        f"| {it} | {loading_table.loc[it, 'PC1']:+.2f} | {loading_table.loc[it, 'PC2']:+.2f} | {loading_table.loc[it, 'PC3']:+.2f} |"
    )
# Parallel analysis: permuting each item separately keeps its distribution but
# removes the correlations between items.
rng = np.random.default_rng(20261117)
n, k = X.shape
eigenvalues = PCA(k).fit(X).explained_variance_
permuted_eigenvalues = np.array(
    [
        PCA(k)
        .fit(np.column_stack([rng.permutation(X[:, j]) for j in range(k)]))
        .explained_variance_
        for _ in range(2000)
    ]
)
null_95 = np.percentile(permuted_eigenvalues, 95, axis=0)
congruence = {j: [] for j in range(3)}
for _ in range(2000):
    i = rng.integers(0, n, n)
    resampled_loadings = PCA(3).fit(StandardScaler().fit_transform(X[i])).components_
    for j in range(3):
        congruence[j].append(
            abs(loadings[j] @ resampled_loadings[j])
            / np.sqrt((loadings[j] @ loadings[j]) * (resampled_loadings[j] @ resampled_loadings[j]))
        )
rob = [
    dict(
        check=f"PC{j + 1} eigenvalue",
        value=eigenvalues[j],
        reference=f"95th percentile of permuted data = {null_95[j]:.2f}",
    )
    for j in range(3)
]
rob += [
    dict(
        check=f"PC{j + 1} bootstrap congruence (median)",
        value=float(np.median(congruence[j])),
        reference=f"5th percentile = {np.percentile(congruence[j], 5):.2f}",
    )
    for j in range(3)
]
stability_table = pd.DataFrame(rob)
stability_table.to_csv(OUT / "03_robustness.csv", index=False, encoding="utf-8-sig")
summary_lines += ["", "## Parallel analysis and bootstrap stability", "", "| Check | Value | Reference |", "|---|---|---|"]
for _, r in stability_table.iterrows():
    summary_lines.append(f"| {r['check']} | {r['value']:+.2f} | {r['reference']} |")
(OUT / "04_summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
print("\n".join(summary_lines))
