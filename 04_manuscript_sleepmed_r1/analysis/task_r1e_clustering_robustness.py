"""用固定的九種分群方法檢查兩群解的一致性，並計算 gap statistic。

分群只用 ISI 七題（ISI 總分 ≥ 9 者）。LCA 吃原始 0–4 分，其餘方法吃 z 分數。"""

from __future__ import annotations
import itertools
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from stepmix.stepmix import StepMix

from _r1_common import AN, ISI_ITEMS, SEED, load_all

B = 2000
OUT = AN / "S_R1e_clustering_robustness"
OUT.mkdir(parents=True, exist_ok=True)
MD = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    MD.append(s)


df = load_all()
ins = df[df["group"] != "HC"].reset_index(drop=True)
Ximp = SimpleImputer(strategy="median").fit_transform(ins[ISI_ITEMS])
X = StandardScaler().fit_transform(Ximp)
Xint = np.rint(Ximp).astype(int)


def gmm(k, cov):
    return GaussianMixture(
        n_components=k,
        covariance_type=cov,
        random_state=SEED,
        n_init=30,
        max_iter=500,
        reg_covar=1e-4,
    ).fit_predict(X)


def lca(k):
    m = StepMix(
        n_components=k,
        measurement="categorical",
        n_init=20,
        max_iter=300,
        random_state=SEED,
        progress_bar=0,
    )
    m.fit(Xint)
    return m.predict(Xint)


METHODS = {
    "K-means": lambda k: KMeans(n_clusters=k, random_state=SEED, n_init=100).fit_predict(X),
    "GMM spherical": lambda k: gmm(k, "spherical"),
    "GMM diagonal": lambda k: gmm(k, "diag"),
    "GMM tied": lambda k: gmm(k, "tied"),
    "GMM full": lambda k: gmm(k, "full"),
    "Hierarchical Ward": lambda k: AgglomerativeClustering(
        n_clusters=k, linkage="ward"
    ).fit_predict(X),
    "Hierarchical average": lambda k: AgglomerativeClustering(
        n_clusters=k, linkage="average"
    ).fit_predict(X),
    "Hierarchical complete": lambda k: AgglomerativeClustering(
        n_clusters=k, linkage="complete"
    ).fit_predict(X),
    "LCA (categorical)": lca,
}
rng = np.random.default_rng(20261112)
n = len(ins)
boot_idx = rng.integers(0, n, (B, n))


def ari_ci(a, b):
    v = np.array([adjusted_rand_score(a[i], b[i]) for i in boot_idx])
    return np.percentile(v, 2.5), np.percentile(v, 97.5)


P("# S_R1e — 九種分群方法的一致性\n")
P(f"失眠症狀者 n = {n}；配對 bootstrap B = {B}\n")
labels, rows, pair_rows = {}, [], []
for k in (2,):
    labels[k] = {name: np.asarray(f(k)) for name, f in METHODS.items()}
    ref = labels[k]["K-means"]
    for name, lab in labels[k].items():
        s = pd.Series(lab)
        sizes = s.value_counts()
        a = adjusted_rand_score(ref, lab)
        lo_, hi_ = (np.nan, np.nan) if name == "K-means" else ari_ci(ref, lab)
        rows.append(
            dict(
                K=k,
                method=name,
                sizes="/".join(map(str, sorted(sizes.tolist(), reverse=True))),
                min_cluster_n=int(sizes.min()),
                ARI_vs_Kmeans=round(a, 3),
                ARI_ci_low=round(lo_, 3) if lo_ == lo_ else "",
                ARI_ci_high=round(hi_, 3) if hi_ == hi_ else "",
            )
        )
    for m1, m2 in itertools.combinations(METHODS, 2):
        pair_rows.append(
            dict(
                K=k,
                method_1=m1,
                method_2=m2,
                ARI=round(adjusted_rand_score(labels[k][m1], labels[k][m2]), 3),
            )
        )
tab = pd.DataFrame(rows)
tab.to_csv(OUT / "01_methods_k2.csv", index=False, encoding="utf-8-sig")
pairs = pd.DataFrame(pair_rows)
pairs.to_csv(OUT / "02_pairwise_ari.csv", index=False, encoding="utf-8-sig")
for k in (2,):
    P(f"## K = {k}：九種方法\n")
    P(tab[tab["K"] == k].drop(columns="K").to_markdown(index=False))
    P("")

# 兩兩 ARI 平均（K = 2）
names = list(METHODS)
combos = list(itertools.combinations(names, 2))
km_combos = [(a, b) for a, b in combos if "K-means" in (a, b)]


def mean_ari(k, idx, cs):
    return np.mean([adjusted_rand_score(labels[k][a][idx], labels[k][b][idx]) for a, b in cs])


full = np.arange(n)
summ = []
for label, cs in [("all 36 method pairs", combos), ("8 pairs involving K-means", km_combos)]:
    o2 = mean_ari(2, full, cs)
    bs = np.array([mean_ari(2, i, cs) for i in boot_idx[:1000]])
    summ.append(
        dict(
            set=label,
            mean_ARI_K2=round(o2, 3),
            K2_ci=f"[{np.percentile(bs, 2.5):.3f}, {np.percentile(bs, 97.5):.3f}]",
        )
    )
ar = pd.read_csv(AN / "S16_bootstrap_stability" / "05_selectK_bootstrap_ARI.csv").set_index("K")
summ.append(
    dict(
        set="K-means bootstrap re-clustering stability (B = 1000)",
        mean_ARI_K2=ar.loc[2, "ari_median"],
        K2_ci=f"[{ar.loc[2, 'ari_2.5']}, {ar.loc[2, 'ari_97.5']}]",
    )
)
summ = pd.DataFrame(summ)
summ.to_csv(OUT / "03_ari_summary.csv", index=False, encoding="utf-8-sig")
P("## 跨方法一致性：K = 2（ARI 慣用讀法：> .80 高、.60–.80 中等、< .60 低）\n")
P(summ.to_markdown(index=False))
P("")


# Gap statistic
def logW(Xd, k, n_init=20):
    if k == 1:
        return np.log(((Xd - Xd.mean(0)) ** 2).sum())
    return np.log(KMeans(n_clusters=k, random_state=SEED, n_init=n_init).fit(Xd).inertia_)


pca = PCA().fit(X)
Xp = pca.transform(X)
mn, mx = Xp.min(0), Xp.max(0)
BG = 200
gap = []
for k in range(1, 7):
    ref = np.array([logW(rng.uniform(mn, mx, Xp.shape), k, n_init=5) for _ in range(BG)])
    gap.append(
        dict(
            K=k,
            logW=round(logW(X, k), 4),
            gap=round(ref.mean() - logW(X, k), 4),
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
P(
    f"## Gap statistic（Tibshirani 準則選出的 K = {sel}；選出 1 表示沒有超過單一連續分布的群集結構）\n"
)
P(gap.to_markdown(index=False))
P("")
P("數字全部由本腳本即時計算；尚未寫回 SSOT。")
(OUT / "05_summary.md").write_text("\n".join(MD), encoding="utf-8")
print(f"\n輸出 → {OUT}")
