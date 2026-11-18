"""K-means 二至六群的 predict-back bootstrap 穩定度。

每次以重抽樣本擬合，再預測回原樣本，計算與參考解的 ARI；每個 K 重設同一個 RNG。"""

from __future__ import annotations

from pathlib import Path
import sqlite3
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler

from _dynamic_n import DB_PATH as _DB

DB_PATH = Path(_DB)
SEED = 42
ISI_ITEMS = [f"ISI_T1_Q{i}" for i in range(1, 8)]

OUT_DIR = Path(__file__).parent / "outputs" / "補充資料" / "analysis" / "S16_bootstrap_stability"
OUT_DIR.mkdir(parents=True, exist_ok=True)

K_RANGE = list(range(2, 7))
N_BOOT = 1000


def load_clustering_sample():
    """分群樣本：具可用 EEG、未排除、ISI ≥ 9 且七題完整。"""
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
    # bootstrap 只用 EEG 有效樣本，和主分析的納入條件一致。
    from pathlib import Path as _P

    relp_csv = (
        _P(__file__).parent
        / "outputs"
        / "補充資料"
        / "analysis"
        / "S17_relative_power"
        / "01_relpower_subject.csv"
    )
    if relp_csv.exists():
        valid_sids = set(pd.read_csv(relp_csv)["subject"].drop_duplicates().tolist())
        from _dynamic_n import excluded_ids

        df = df[df["NEW_NUMBER"].isin(valid_sids) & ~df["NEW_NUMBER"].isin(excluded_ids())].copy()
    df_ins = df[df["ISI_T1"] >= 9].dropna(subset=ISI_ITEMS).copy()
    return df_ins


def main():
    df_ins = load_clustering_sample()
    n = len(df_ins)
    print(f"分群樣本 n = {n}")

    imp = SimpleImputer(strategy="median")
    scl = StandardScaler()
    X = scl.fit_transform(imp.fit_transform(df_ins[ISI_ITEMS].values))
    n_idx = np.arange(n)

    rows = []
    for k in K_RANGE:
        # 參考解（n_init=100, seed=42）
        km_ref = KMeans(n_clusters=k, random_state=SEED, n_init=100)
        lab_ref = km_ref.fit_predict(X)

        # 每個 K 各自重設 RNG，每輪只抽一次樣本。
        rng = np.random.default_rng(SEED)
        aris = []
        for b in range(N_BOOT):
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
    out.to_csv(OUT_DIR / "05_selectK_bootstrap_ARI.csv", index=False)
    print(f"\n✓ 寫入 {OUT_DIR / '05_selectK_bootstrap_ARI.csv'}")


if __name__ == "__main__":
    main()
