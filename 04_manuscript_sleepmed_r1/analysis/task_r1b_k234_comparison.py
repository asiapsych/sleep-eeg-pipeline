"""K-means 二至六群的內部效度指標與重抽樣穩定度。

分群只用 ISI 七題（ISI 總分 ≥ 9 者）。重抽樣 ARI 讀 task_e2 的輸出。"""

from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

from _r1_common import AN, ISI_ITEMS, SEED, load_all

OUT = AN / "S_R1b_k234"
OUT.mkdir(parents=True, exist_ok=True)
MD = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    MD.append(s)


df = load_all()
ins = df[df["group"] != "HC"]
P("# S_R1b — K-means 二至六群\n")
P(f"失眠症狀者 n = {len(ins)}\n")

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
(OUT / "07_summary.md").write_text("\n".join(MD), encoding="utf-8")
