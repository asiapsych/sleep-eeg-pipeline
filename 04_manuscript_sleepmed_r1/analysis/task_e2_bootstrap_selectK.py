"""Predict-back bootstrap stability of K-means with 2 to 6 clusters.

Each resample is clustered, the fit is used to assign the original sample, and the
assignment is compared with the reference solution (ARI). The RNG is reset for each K."""

import sqlite3
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

from _dynamic_n import DB_PATH, ISI_ITEMS, S17_RELPOWER_CSV, SEED, excluded_ids
from _r1_common import prepare_isi_items

OUT_DIR = Path(__file__).parent / "outputs" / "analysis" / "S16_bootstrap_stability"
OUT_DIR.mkdir(parents=True, exist_ok=True)

K_RANGE = list(range(2, 7))
N_BOOT = 1000


def load_clustering_sample():
    """Clustering sample: usable EEG, not excluded, ISI >= 9 and all seven items answered."""
    quoted = ", ".join(f'"{c}"' for c in ["NEW_NUMBER", "ISI_T1"] + ISI_ITEMS)
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql(
        f"SELECT {quoted} FROM isi_raw_data WHERE ISI_T1 IS NOT NULL AND ISI_T1_Q1 IS NOT NULL",
        conn,
    )
    conn.close()
    for c in df.columns:
        if c != "NEW_NUMBER":
            df[c] = pd.to_numeric(df[c], errors="coerce")
    # Keep the database row order so the same seed draws the same participants.
    if S17_RELPOWER_CSV.exists():
        eeg_subjects = pd.read_csv(S17_RELPOWER_CSV)["subject"]
        df = df[df["NEW_NUMBER"].isin(eeg_subjects) & ~df["NEW_NUMBER"].isin(excluded_ids())].copy()
    df_ins = df[df["ISI_T1"] >= 9].dropna(subset=ISI_ITEMS).copy()
    return df_ins


def main():
    df_ins = load_clustering_sample()
    n = len(df_ins)
    print(f"clustering sample: n = {n}")

    _, X = prepare_isi_items(df_ins[ISI_ITEMS].values)
    n_idx = np.arange(n)

    rows = []
    for k in K_RANGE:
        km_ref = KMeans(n_clusters=k, random_state=SEED, n_init=100)
        lab_ref = km_ref.fit_predict(X)

        # same draws for every K
        rng = np.random.default_rng(SEED)
        aris = []
        for _ in range(N_BOOT):
            boot_idx = rng.choice(n_idx, size=n, replace=True)
            Xb = X[boot_idx]
            if len(np.unique(Xb, axis=0)) <= k:
                continue
            km_b = KMeans(n_clusters=k, random_state=SEED, n_init=10)
            km_b.fit(Xb)
            pred_full = km_b.predict(X)
            aris.append(adjusted_rand_score(lab_ref, pred_full))

        aris = np.asarray(aris)
        med = float(np.median(aris))
        lo, hi = np.percentile(aris, [2.5, 97.5])
        rows.append(
            {
                "K": k,
                "ari_median": round(med, 4),
                "ari_2.5": round(float(lo), 4),
                "ari_97.5": round(float(hi), 4),
                "n_boot": len(aris),
            }
        )
        print(f"  K={k}: ARI median={med:.4f}  95% [{lo:.4f}, {hi:.4f}]  (n={len(aris)})")

    out = pd.DataFrame(rows)
    output_csv = OUT_DIR / "05_selectK_bootstrap_ARI.csv"
    out.to_csv(output_csv, index=False)
    print(f"\nwrote {output_csv}")


if __name__ == "__main__":
    main()
