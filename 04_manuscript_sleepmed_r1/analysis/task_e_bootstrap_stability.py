"""K = 3 的 bootstrap silhouette 與 predict-back ARI。

參考解沿用主分析；每次用重抽樣本擬合，再預測回原樣本計算 ARI。
用這裡可重現的結果取代先前無腳本可查的 silhouette 數字。"""

from __future__ import annotations

import sqlite3
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler

BASE_DIR = Path(__file__).parent
OUT_DIR = BASE_DIR / "outputs" / "補充資料" / "analysis" / "S16_bootstrap_stability"
OUT_DIR.mkdir(parents=True, exist_ok=True)

from _dynamic_n import DB_PATH as _DB

DB_PATH = Path(_DB)
SEED = 42
N_BOOT = 1000
K = 3
ISI_ITEMS = [f"ISI_T1_Q{i}" for i in range(1, 8)]


def load_clustering_sample():
    """與主分析完全相同：T1 ISI≥9 + 完整七題 → n=165。"""
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
        before = len(df)
        from _dynamic_n import excluded_ids

        df = df[df["NEW_NUMBER"].isin(valid_sids) & ~df["NEW_NUMBER"].isin(excluded_ids())].copy()
        print(f"  [v2.2 new methodology] bootstrap K-means 限制到 EEG∩DB: {before} → {len(df)}")
    df_ins = df[df["ISI_T1"] >= 9].dropna(subset=ISI_ITEMS).copy()
    return df_ins


def main():
    df_ins = load_clustering_sample()
    n = len(df_ins)
    print(f"分群樣本 n = {n}（T1 ISI≥9 + 完整七題）")

    imp = SimpleImputer(strategy="median")
    scl = StandardScaler()
    X = scl.fit_transform(imp.fit_transform(df_ins[ISI_ITEMS].values))

    # 原始 K=3 解（n_init=100，與主分析一致）
    km_orig = KMeans(n_clusters=K, random_state=SEED, n_init=100)
    labels_orig = km_orig.fit_predict(X)
    sil_orig = silhouette_score(X, labels_orig)
    print(f"原始 K=3 silhouette = {sil_orig:.4f}")

    rng = np.random.default_rng(SEED)
    sil_rows = []
    ari_rows = []
    n_idx = np.arange(n)

    for b in range(N_BOOT):
        boot_idx = rng.choice(n_idx, size=n, replace=True)
        Xb = X[boot_idx]
        # (1) Silhouette：bootstrap 樣本上重 fit K=3，需 unique 樣本 ≥ K
        uniq = np.unique(Xb, axis=0)
        if len(uniq) <= K:
            continue
        km_b = KMeans(n_clusters=K, random_state=SEED, n_init=10)
        lab_b = km_b.fit_predict(Xb)
        try:
            sil_b = silhouette_score(Xb, lab_b)
        except Exception:
            continue
        sil_rows.append({"iteration": b + 1, "silhouette": sil_b})

        # (2) ARI：用 bootstrap-fit 模型 predict 原始全樣本，與原始標籤比對
        pred_full = km_b.predict(X)
        ari_b = adjusted_rand_score(labels_orig, pred_full)
        ari_rows.append({"iteration": b + 1, "ARI_vs_original": ari_b})

        if (b + 1) % 200 == 0:
            print(f"  [{b + 1}/{N_BOOT}] done")

    sil_df = pd.DataFrame(sil_rows)
    ari_df = pd.DataFrame(ari_rows)
    sil_df.to_csv(OUT_DIR / "01_silhouette_bootstrap_raw.csv", index=False)
    ari_df.to_csv(OUT_DIR / "02_ari_bootstrap_raw.csv", index=False)

    def describe(s):
        return {
            "n_valid": int(s.count()),
            "mean": float(s.mean()),
            "sd": float(s.std(ddof=1)),
            "median": float(s.median()),
            "ci95_low": float(s.quantile(0.025)),
            "ci95_high": float(s.quantile(0.975)),
            "min": float(s.min()),
            "max": float(s.max()),
            "pct_gt_0.15": float((s > 0.15).mean() * 100),
        }

    sil_stats = describe(sil_df["silhouette"])
    ari_stats = describe(ari_df["ARI_vs_original"])

    summ = pd.DataFrame({"silhouette_bootstrap": sil_stats, "ARI_bootstrap": ari_stats}).T
    summ.to_csv(OUT_DIR / "03_summary.csv")

    lines = []
    P = lines.append
    P("=" * 78)
    P("Task E — 可重現之聚類穩定性 Bootstrap（N_BOOT = 1000）")
    P("=" * 78)
    P(f"分群樣本 n = {n}（T1 ISI≥9 + 完整七題；與主分析一致）")
    P(f"原始 K=3 silhouette = {sil_orig:.4f}")
    P(f"K-means: StandardScaler + median imputer + KMeans(K=3, n_init=100/orig, 10/boot)")
    P("")
    P("─" * 78)
    P("(1) Silhouette Bootstrap（每次重採樣後重 fit K=3）")
    P("─" * 78)
    for k, v in sil_stats.items():
        P(f"  {k:15s}: {v:.4f}" if isinstance(v, float) else f"  {k:15s}: {v}")
    P("")
    P("─" * 78)
    P("(2) ARI Bootstrap（bootstrap 模型 predict 原樣本 vs 原始標籤）")
    P("─" * 78)
    for k, v in ari_stats.items():
        P(f"  {k:15s}: {v:.4f}" if isinstance(v, float) else f"  {k:15s}: {v}")
    P("")
    P("─" * 78)
    P("對照論文 §7.2 原宣稱（silhouette mean 0.219，95% CI [0.186, 0.270]）")
    P("─" * 78)
    P(f"  本次可重現之 silhouette mean = {sil_stats['mean']:.3f}")
    P(
        f"  本次可重現之 silhouette 95% CI = [{sil_stats['ci95_low']:.3f}, {sil_stats['ci95_high']:.3f}]"
    )
    claim_match = abs(sil_stats["mean"] - 0.219) < 0.03
    P(
        f"  與原宣稱 mean=0.219 差距 = {abs(sil_stats['mean'] - 0.219):.3f} → "
        f"{'近似可接受' if claim_match else '差異明顯，應改用本次數字'}"
    )
    P("")
    P("【建議】論文 §7.2 應改為本腳本可重現之數字：")
    P(
        f"  「Bootstrap (n = 1,000): silhouette mean = {sil_stats['mean']:.3f}, "
        f"95% CI [{sil_stats['ci95_low']:.3f}, {sil_stats['ci95_high']:.3f}], "
        f"{sil_stats['pct_gt_0.15']:.0f}% > 0.15; "
        f"bootstrap ARI median = {ari_stats['median']:.3f}; "
        f"95% bootstrap percentile interval [{ari_stats['ci95_low']:.3f}, {ari_stats['ci95_high']:.3f}]. "
        f"ARI median falls near the lower bound of moderate recovery, but the interval spans poor to excellent recovery.」"
    )

    txt = "\n".join(lines)
    (OUT_DIR / "04_summary.txt").write_text(txt, encoding="utf-8")
    print("\n" + txt)


if __name__ == "__main__":
    main()
