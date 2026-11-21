"""Compare nine clustering methods and compute the gap statistic."""

from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture
from stepmix.stepmix import StepMix

from _r1_common import AN, ISI_ITEMS, SEED, load_all, prepare_isi_items

B = 2000
# Keep K = 2 for the manuscript. Set K = 3 to reproduce Table C in the response letter.
# Some output files are overwritten when K changes; keep the previous run first.
K = 2
OUT = AN / "S_R1e_clustering_robustness"
OUT.mkdir(parents=True, exist_ok=True)
summary_lines = []


def report(line):
    print(line)
    summary_lines.append(line)


df = load_all()
ins = df[df["group"] != "HC"].reset_index(drop=True)
filled_items, X = prepare_isi_items(ins[ISI_ITEMS])
# LCA treats the imputed 0-4 item scores as categories; the other methods use z scores.
categories = np.rint(filled_items).astype(int)

rng = np.random.default_rng(20261112)
n = len(ins)
boot_idx = rng.integers(0, n, (B, n))


def ari_ci(a, b):
    v = np.array([adjusted_rand_score(a[i], b[i]) for i in boot_idx])
    return np.percentile(v, 2.5), np.percentile(v, 97.5)


report("# Agreement across nine clustering methods\n")
report(f"participants with insomnia symptoms: n = {n}; paired bootstrap B = {B}\n")
labels = {"K-means": KMeans(n_clusters=K, random_state=SEED, n_init=100).fit_predict(X)}
for name, covariance in (
    ("spherical", "spherical"),
    ("diagonal", "diag"),
    ("tied", "tied"),
    ("full", "full"),
):
    labels[f"GMM {name}"] = GaussianMixture(
        n_components=K,
        covariance_type=covariance,
        random_state=SEED,
        n_init=30,
        max_iter=500,
        reg_covar=1e-4,
    ).fit_predict(X)
for name, linkage in (("Ward", "ward"), ("average", "average"), ("complete", "complete")):
    labels[f"Hierarchical {name}"] = AgglomerativeClustering(
        n_clusters=K, linkage=linkage
    ).fit_predict(X)
lca = StepMix(
    n_components=K,
    measurement="categorical",
    n_init=20,
    max_iter=300,
    random_state=SEED,
    progress_bar=0,
)
lca.fit(categories)
labels["LCA (categorical)"] = lca.predict(categories)

reference = labels["K-means"]
rows = []
for name, assigned in labels.items():
    sizes = pd.Series(assigned).value_counts()
    agreement = adjusted_rand_score(reference, assigned)
    low, high = (np.nan, np.nan) if name == "K-means" else ari_ci(reference, assigned)
    rows.append(
        dict(
            K=K,
            method=name,
            sizes="/".join(map(str, sorted(sizes.tolist(), reverse=True))),
            min_cluster_n=int(sizes.min()),
            ARI_vs_Kmeans=round(agreement, 3),
            ARI_ci_low="" if np.isnan(low) else round(low, 3),
            ARI_ci_high="" if np.isnan(high) else round(high, 3),
        )
    )

method_pairs = list(combinations(labels, 2))
pair_rows = [
    dict(K=K, method_1=a, method_2=b, ARI=round(adjusted_rand_score(labels[a], labels[b]), 3))
    for a, b in method_pairs
]
method_table = pd.DataFrame(rows)
method_table.to_csv(OUT / f"01_methods_k{K}.csv", index=False, encoding="utf-8-sig")
pair_table = pd.DataFrame(pair_rows)
pair_table.to_csv(OUT / "02_pairwise_ari.csv", index=False, encoding="utf-8-sig")
report(f"## K = {K}, nine methods\n")
report(method_table.drop(columns="K").to_markdown(index=False))
report("")

kmeans_pairs = [pair for pair in method_pairs if "K-means" in pair]


def mean_ari(indices, pairs):
    # Each pair uses the same resampled participants for both sets of labels.
    return np.mean([adjusted_rand_score(labels[a][indices], labels[b][indices]) for a, b in pairs])


full = np.arange(n)
summ = []
for label, pairs in [
    ("all 36 method pairs", method_pairs),
    ("8 pairs involving K-means", kmeans_pairs),
]:
    observed_ari = mean_ari(full, pairs)
    bootstrap_ari = np.array([mean_ari(indices, pairs) for indices in boot_idx[:1000]])
    summ.append(
        dict(
            set=label,
            **{
                f"mean_ARI_K{K}": round(observed_ari, 3),
                f"K{K}_ci": f"[{np.percentile(bootstrap_ari, 2.5):.3f}, {np.percentile(bootstrap_ari, 97.5):.3f}]",
            },
        )
    )
ar = pd.read_csv(AN / "S16_bootstrap_stability" / "05_selectK_bootstrap_ARI.csv").set_index("K")
summ.append(
    dict(
        set="K-means bootstrap re-clustering stability (B = 1000)",
        **{
            f"mean_ARI_K{K}": ar.loc[K, "ari_median"],
            f"K{K}_ci": f"[{ar.loc[K, 'ari_2.5']}, {ar.loc[K, 'ari_97.5']}]",
        },
    )
)
summ = pd.DataFrame(summ)
summ.to_csv(OUT / "03_ari_summary.csv", index=False, encoding="utf-8-sig")
report(f"## Agreement between methods (K = {K})\n")
report(summ.to_markdown(index=False))
report("")


def log_within_cluster_dispersion(items, k, n_init=20):
    if k == 1:
        return np.log(((items - items.mean(0)) ** 2).sum())
    return np.log(KMeans(n_clusters=k, random_state=SEED, n_init=n_init).fit(items).inertia_)


# Draw uniform reference samples over the observed range of each PC score.
pca = PCA().fit(X)
rotated = pca.transform(X)
lower, upper = rotated.min(0), rotated.max(0)
BG = 200
gap = []
for k in range(1, 7):
    ref = np.array(
        [
            log_within_cluster_dispersion(rng.uniform(lower, upper, rotated.shape), k, n_init=5)
            for _ in range(BG)
        ]
    )
    observed_logw = log_within_cluster_dispersion(X, k)
    gap.append(
        dict(
            K=k,
            logW=round(observed_logw, 4),
            gap=round(ref.mean() - observed_logw, 4),
            s_k=round(ref.std(ddof=0) * np.sqrt(1 + 1 / BG), 4),
        )
    )
gap = pd.DataFrame(gap)
gap["gap_K_ge_gap_K+1_minus_s"] = [
    bool(gap.loc[i, "gap"] >= gap.loc[i + 1, "gap"] - gap.loc[i + 1, "s_k"])
    if i + 1 < len(gap)
    else ""
    for i in range(len(gap))
]
sel = next(
    (
        int(gap.loc[i, "K"])
        for i in range(len(gap) - 1)
        if gap.loc[i, "gap_K_ge_gap_K+1_minus_s"] is True
    ),
    None,
)
gap.to_csv(OUT / "04_gap_statistic.csv", index=False, encoding="utf-8-sig")
report(f"## Gap statistic (Tibshirani rule: K = {sel})\n")
report(gap.to_markdown(index=False))
report("")
(OUT / "05_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")
print(f"\nresults in {OUT}")
