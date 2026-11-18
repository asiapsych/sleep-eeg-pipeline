"""核對 58 人的人工審查與最終自動版：ICA、epoch 及 EEG 指標。

Welch 設定要和 MATLAB 相同；組間差只用兩版都有的同一批人。
樣本較小，這裡主要用來看一致性與方向。"""

from __future__ import annotations
import os
import warnings
from pathlib import Path

import mne
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.signal import welch
from scipy.signal.windows import hamming
from scipy.stats import pearsonr, spearmanr

from _r1_common import AN, hedges_g, load_all, icc_a1
from task_b_fooof_full import CH_RENAME, aggregate_features, fit_subject
from task_f_recompute_relpower import compute_subject

warnings.filterwarnings("ignore")
mne.set_log_level("ERROR")

from _dynamic_n import local_path

V1_DIR = Path(local_path("v1_labdata", "PSYML_V1_LABDATA"))
OUT = AN / "S_R1j_manual_vs_auto"
OUT.mkdir(parents=True, exist_ok=True)
SFREQ, NWIN = 300, 900  # 與 v2.2 相同：Hamming 900 點、nfft 900、不重疊
MEAS = [
    ("beta_pw_BRAIN", "Periodic beta peak power, whole-scalp"),
    ("Beta_rel_BRAIN", "Relative beta power (13–30 Hz), whole-scalp"),
    ("exponent_BRAIN", "Aperiodic exponent, whole-scalp"),
    ("Alpha_rel_BRAIN", "Relative alpha power (8–13 Hz), whole-scalp"),
]


def psd_from_set(path):
    ep = mne.read_epochs_eeglab(str(path))
    x = ep.get_data() * 1e6  # µV；[片段 × 電極 × 時間]
    f, p = welch(
        x,
        fs=SFREQ,
        window=hamming(NWIN, sym=True),
        nperseg=NWIN,
        noverlap=0,
        nfft=NWIN,
        detrend=False,
        axis=-1,
    )  # 與 MATLAB pwelch(hamming(900)) 一致
    keep = (f >= 0.5) & (f <= 70)
    ch = [CH_RENAME.get(c.upper(), c.upper()) for c in ep.ch_names]
    return f[keep], p.mean(0)[:, keep], ch, x.shape[0]


def main():
    # A 成分
    rows = []
    for f in sorted(V1_DIR.glob("*_ica_decisions.mat")):
        m = sio.loadmat(str(f), squeeze_me=True)
        a = set(np.atleast_1d(m["auto_remove"]).astype(int))
        fr = set(np.atleast_1d(m["final_remove"]).astype(int))
        rows.append(
            dict(
                subject=f.name[:7],
                n_auto_flagged=len(a),
                n_removed_after_review=len(fr),
                n_overridden_keep=len(a - fr),
                n_added_by_reviewer=len(fr - a),
            )
        )
    ic = pd.DataFrame(rows)
    ic_sum = pd.DataFrame(
        [
            dict(
                n_subjects=len(ic),
                n_auto_flagged=ic.n_auto_flagged.sum(),
                n_confirmed=ic.n_removed_after_review.sum() - ic.n_added_by_reviewer.sum(),
                n_overridden_keep=ic.n_overridden_keep.sum(),
                n_added_by_reviewer=ic.n_added_by_reviewer.sum(),
                n_subjects_with_override=int(
                    ((ic.n_overridden_keep + ic.n_added_by_reviewer) > 0).sum()
                ),
                agreement=1
                - (ic.n_overridden_keep.sum() + ic.n_added_by_reviewer.sum())
                / ic.n_auto_flagged.sum(),
            )
        ]
    )
    ic_sum.to_csv(OUT / "01_ic_review.csv", index=False)

    # B 片段
    er = pd.read_csv(V1_DIR / "epoch_review_log.csv")
    ep_sum = pd.DataFrame(
        [
            dict(
                n_subjects=len(er),
                mean_before=er.n_before_review.mean(),
                mean_removed=er.n_manual_removed.mean(),
                median_removed=er.n_manual_removed.median(),
                min_removed=er.n_manual_removed.min(),
                max_removed=er.n_manual_removed.max(),
                pct_removed=100 * er.n_manual_removed.sum() / er.n_before_review.sum(),
            )
        ]
    )
    ep_sum.to_csv(OUT / "02_epoch_review.csv", index=False)

    # C 指標
    auto = load_all()
    long, rel = [], []
    for sid in ic.subject:
        p = V1_DIR / f"{sid}_step3.set"
        if not p.exists():
            continue
        freqs, psd, ch, n_ep = psd_from_set(p)
        if n_ep < 60:
            continue  # 與最終規則一致：不足 60 個片段不納入
        long += fit_subject(sid, freqs, psd, ch, n_ep)
        rel.append(compute_subject(sid, freqs, psd, ch))
    man = aggregate_features(pd.DataFrame(long)).merge(pd.DataFrame(rel), on="subject")
    cols = [c for c, _ in MEAS]
    d = man[["subject"] + cols].merge(
        auto[["subject", "k2"] + cols], on="subject", suffixes=("_manual", "_auto")
    )
    d.to_csv(OUT / "_subject_level_internal.csv", index=False)

    out = []
    for c, lab in MEAS:
        s = d[[f"{c}_manual", f"{c}_auto"]].dropna()
        x, y = s.iloc[:, 0].values, s.iloc[:, 1].values
        diff = y - x
        out.append(
            dict(
                measure=lab,
                n=len(s),
                pearson_r=pearsonr(x, y)[0],
                spearman_rho=spearmanr(x, y)[0],
                icc_a1=icc_a1(x, y),
                mean_manual=x.mean(),
                mean_auto=y.mean(),
                mean_diff_auto_minus_manual=diff.mean(),
                loa_low=diff.mean() - 1.96 * diff.std(ddof=1),
                loa_high=diff.mean() + 1.96 * diff.std(ddof=1),
            )
        )
    ag = pd.DataFrame(out)
    ag.to_csv(OUT / "03_agreement.csv", index=False)

    # 同一批人的組間差（人數少，僅供對照方向）
    grp = []
    hc = d[d.k2 == "HC"]
    for g in ["Lower", "Higher"]:
        a = d[d.k2 == g]
        for v in ["manual", "auto"]:
            col = f"beta_pw_BRAIN_{v}"
            gg, lo, hi, n1, n2 = hedges_g(a[col].dropna(), hc[col].dropna())
            grp.append(dict(group=g, version=v, g=gg, lo=lo, hi=hi, n_group=n1, n_hc=n2))
    grp = pd.DataFrame(grp)

    L = [
        "# 人工審查版與全自動版的一致性（R1j）",
        "",
        f"- 目視審查的受試者：{len(ic)} 人；其中 {len(d)} 人在最終分析樣本內（其餘因不足 60 個片段或沒有完整的 ISI 未納入）。",
        "",
        "## A 成分",
        "",
        f"ICLabel 標記 {int(ic_sum.n_auto_flagged[0])} 個成分；研究者確認移除 {int(ic_sum.n_confirmed[0])} 個，改為保留 {int(ic_sum.n_overridden_keep[0])} 個（{int(ic_sum.n_subjects_with_override[0])} 人），"
        f"另外加移除 {int(ic_sum.n_added_by_reviewer[0])} 個。一致率 {100 * ic_sum.agreement[0]:.1f}%。",
        "",
        "## B 片段",
        "",
        f"三層自動剔除後每人平均剩 {ep_sum.mean_before[0]:.1f} 個片段；研究者目視再剔除中位數 {ep_sum.median_removed[0]:.0f} 個"
        f"（平均 {ep_sum.mean_removed[0]:.1f}，範圍 {int(ep_sum.min_removed[0])}–{int(ep_sum.max_removed[0])}），佔 {ep_sum.pct_removed[0]:.1f}%。",
        "",
        "## C 指標（目視審查後的資料 對 最終流程）",
        "",
        "| 指標 | n | Pearson r | Spearman rho | ICC(A,1) | 人工版平均 | 自動版平均 | 平均差 | 95% 一致界限 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in ag.itertuples():
        L.append(
            f"| {r.measure} | {r.n} | {r.pearson_r:.2f} | {r.spearman_rho:.2f} | {r.icc_a1:.2f} | {r.mean_manual:.3f} | {r.mean_auto:.3f} | {r.mean_diff_auto_minus_manual:+.3f} | [{r.loa_low:+.3f}, {r.loa_high:+.3f}] |"
        )
    L += [
        "",
        "## 同一批人的週期性 beta 組間差（人數少，只看方向）",
        "",
        "| 組 | 版本 | g | 95% CI | n（組／健康組） |",
        "|---|---|---|---|---|",
    ]
    for r in grp.itertuples():
        L.append(
            f"| {r.group} | {r.version} | {r.g:+.2f} | [{r.lo:+.2f}, {r.hi:+.2f}] | {r.n_group}／{r.n_hc} |"
        )
    (OUT / "04_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
