"""把 task_e 的 predict-back bootstrap 延伸到 K = 2–6。

每個 K 重設同一個 RNG；K = 3 要能逐次重現 task_e。
不要混用舊的 in-bag ARI，兩種重抽程序的數字不能直接並列。"""

from __future__ import annotations

from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler

from task_e_bootstrap_stability import load_clustering_sample, ISI_ITEMS, SEED

OUT_DIR = Path(__file__).parent / "outputs" / "補充資料" / "analysis" / "S16_bootstrap_stability"
OUT_DIR.mkdir(parents=True, exist_ok=True)

K_RANGE = list(range(2, 7))
N_BOOT = 1000  # 與 task_e 相同


def main():
    df_ins = load_clustering_sample()
    n = len(df_ins)
    print(f"分群樣本 n = {n}（與 task_e 同一 loader）")

    imp = SimpleImputer(strategy="median")
    scl = StandardScaler()
    X = scl.fit_transform(imp.fit_transform(df_ins[ISI_ITEMS].values))
    n_idx = np.arange(n)

    rows = []
    for k in K_RANGE:
        # 參考解：與 task_e 相同設定（n_init=100, seed=42）
        km_ref = KMeans(n_clusters=k, random_state=SEED, n_init=100)
        lab_ref = km_ref.fit_predict(X)

        # 每個 K 各自重設 RNG，每輪只抽一次樣本，才可和 task_e 對照。
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

    # 自我驗證：K=3 必須重現 task_e 的已發表數字
    k3 = out.loc[out["K"] == 3].iloc[0]
    ok = abs(k3["ari_median"] - 0.6192) < 5e-5
    print(
        f"K=3 重現 task_e 0.6192: {'✓ 一致' if ok else '✗ 不一致——請檢查'} "
        f"(本次 {k3['ari_median']}, [{k3['ari_2.5']}, {k3['ari_97.5']}])"
    )


if __name__ == "__main__":
    main()
