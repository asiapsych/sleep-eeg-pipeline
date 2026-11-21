"""Internal validity indices and bootstrap stability of K-means with 2 to 6 clusters.

Only the seven ISI items are clustered (ISI >= 9). The bootstrap ARIs come from task_e2."""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import (
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)

from _r1_common import AN, ISI_ITEMS, SEED, load_all, prepare_isi_items

OUT = AN / "S_R1b_k234"
OUT.mkdir(parents=True, exist_ok=True)
summary_lines = []


def report(line):
    print(line)
    summary_lines.append(line)


df = load_all()
ins = df[df["group"] != "HC"]
report("# K-means, 2 to 6 clusters\n")
report(f"participants with insomnia symptoms: n = {len(ins)}\n")

_, X = prepare_isi_items(ins[ISI_ITEMS])
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
report("## Internal validity and bootstrap stability\n")
report(val.to_markdown(index=False))
report("")
(OUT / "02_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")
