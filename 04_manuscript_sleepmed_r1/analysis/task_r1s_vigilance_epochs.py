"""以後腦區 theta／alpha 標記疑似嗜睡 epoch，再重算敏感度結果。

兩種標記：theta ≥ alpha，或 alpha 低於此人中位數的一半。
排除後至少留 30 個 epoch。沒有眼動電極，這些指標不能當成睡眠分期，
3 秒的 epoch 也無法排除更短的微睡眠。"""

from __future__ import annotations
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.io import loadmat
from scipy.stats import kruskal, mannwhitneyu, spearmanr

from _r1_common import AN, COV, COV_COLS, OUTCOME, fmt_p, hedges_g, load_all
from task_b_fooof_full import CH_RENAME, aggregate_features, fit_subject

warnings.filterwarnings("ignore")
from _dynamic_n import local_path

PIPE = Path(
    local_path("eeg_pipeline_root", "EEG_PIPELINE_ROOT")
)
EPD = PIPE / "outputs" / "psd_epochs"
OUT = AN / "S_R1s_vigilance_epochs"
OUT.mkdir(parents=True, exist_ok=True)
POST = ["O1", "O2", "P3", "P4", "PZ"]
MIN_KEEP = 30
THR = 18


def band(psd, f, lo, hi):
    sel = (f >= lo) & (f < hi)
    return np.trapz(psd[..., sel], f[sel], axis=-1)


def beta_from(sid, f, psd_mean, ch, n_ep):
    return aggregate_features(pd.DataFrame(fit_subject(sid, f, psd_mean, ch, n_ep)))[OUTCOME].iloc[
        0
    ]


def main():
    df = load_all()
    rows = []
    for sid in df["subject"]:
        p = EPD / f"psdep_{sid}.mat"
        if not p.exists():
            continue
        m = loadmat(str(p), squeeze_me=True)
        f = np.asarray(m["freqs"], float)
        P = np.asarray(m["psd_ep"], float)  # [epoch × 電極 × 頻率]
        ch = [
            CH_RENAME.get(str(c).strip().upper(), str(c).strip().upper())
            for c in np.atleast_1d(m["ch_names"])
        ]
        pi = [ch.index(c) for c in POST]
        th = band(P[:, pi, :], f, 4, 8).mean(1)
        al = band(P[:, pi, :], f, 8, 13).mean(1)  # 每個 epoch 一個值
        flagA = th >= al
        flagB = al < 0.5 * np.median(al)
        n = len(th)
        r = dict(
            subject=sid,
            ep_total=n,
            pctA=100 * flagA.mean(),
            pctB=100 * flagB.mean(),
            pctAB=100 * (flagA | flagB).mean(),
            ratio_median=float(np.median(th / al)),
            ratio_slope=float(np.polyfit(np.arange(n), np.log(th / al), 1)[0]) * n,
        )
        for tag, keep in [
            ("all", np.ones(n, bool)),
            ("exclA", ~flagA),
            ("exclB", ~flagB),
            ("exclAB", ~(flagA | flagB)),
            ("first", np.arange(n) < n // 2),
            ("second", np.arange(n) >= n // 2),
        ]:
            r[f"nkeep_{tag}"] = int(keep.sum())
            r[f"beta_{tag}"] = (
                beta_from(sid, f, P[keep].mean(0), ch, int(keep.sum()))
                if keep.sum() >= MIN_KEEP
                else np.nan
            )
        rows.append(r)
    ep = pd.DataFrame(rows)
    d = df.merge(ep, on="subject", how="left")
    d[
        ["subject", "k2", "pctA", "pctB", "pctAB", "ratio_median", "ratio_slope"]
        + [c for c in ep.columns if c.startswith(("nkeep_", "beta_"))]
    ].to_csv(OUT / "_subject_level_internal.csv", index=False)
    have = d.dropna(subset=["pctA"])
    hc = have[have.k2 == "HC"]
    chk = have[[OUTCOME, "beta_all"]].dropna()
    r_chk = np.corrcoef(chk[OUTCOME], chk["beta_all"])[0, 1]
    max_dev = float((chk[OUTCOME] - chk["beta_all"]).abs().max())

    # 1 各組比例
    A = []
    for v, lab in [
        ("pctA", "Epochs with theta ≥ alpha, %"),
        ("pctB", "Epochs with alpha < 50% of own median, %"),
        ("pctAB", "Epochs meeting either criterion, %"),
        ("ratio_median", "Median posterior theta/alpha ratio"),
        ("ratio_slope", "Change in log theta/alpha across the recording"),
    ]:
        row = dict(
            measure=lab,
            kw_p=float(
                kruskal(*[have.loc[have.k2 == g, v] for g in ("HC", "Lower", "Higher")]).pvalue
            ),
        )
        for g in ("HC", "Lower", "Higher"):
            x = have.loc[have.k2 == g, v]
            row.update({f"{g}_mean": x.mean(), f"{g}_median": x.median(), f"{g}_n": len(x)})
            if g != "HC":
                gg, lo, hi, _, _ = hedges_g(x, hc[v])
                row.update({f"{g}_g": gg, f"{g}_p": float(mannwhitneyu(x, hc[v]).pvalue)})
        rho, prho = spearmanr(have[v], have[OUTCOME], nan_policy="omit")
        row.update(rho_with_beta=rho, p_with_beta=prho)
        A.append(row)
    A = pd.DataFrame(A)
    A.to_csv(OUT / "01_flagged_by_group.csv", index=False)
    anyA = {g: 100 * (have.loc[have.k2 == g, "pctA"] > 0).mean() for g in ("HC", "Lower", "Higher")}
    allA = {g: int((have.loc[have.k2 == g, "pctA"] > 50).sum()) for g in ("HC", "Lower", "Higher")}

    # 2 敏感度
    S = []

    def res(col, label, dd=None, extra=""):
        dd = (have if dd is None else dd).dropna(subset=[col] + COV_COLS)
        h0 = dd.loc[dd.k2 == "HC", col]
        sd = h0.std(ddof=1)
        gh = hedges_g(dd.loc[dd.k2 == "Higher", col], h0)
        gl = hedges_g(dd.loc[dd.k2 == "Lower", col], h0)
        m = smf.ols(f"{col} ~ C(k2, Treatment('HC')) + {COV}{extra}", dd).fit()
        k = [c for c in m.params.index if "Higher" in c][0]
        mt = smf.ols(
            f"{col} ~ thr + {COV}{extra}", dd.assign(thr=(dd.ISI_T1 >= THR).astype(int))
        ).fit()
        S.append(
            dict(
                analysis=label,
                n=len(dd),
                n_hc=len(h0),
                n_lower=gl[3],
                n_higher=gh[3],
                g_higher=gh[0],
                g_higher_lo=gh[1],
                g_higher_hi=gh[2],
                g_lower=gl[0],
                g_lower_lo=gl[1],
                g_lower_hi=gl[2],
                adj_higher_in_HC_SD=m.params[k] / sd,
                adj_higher_p=m.pvalues[k],
                thr_in_HC_SD=mt.params["thr"] / sd,
                thr_p=mt.pvalues["thr"],
            )
        )

    res(OUTCOME, "All 60 epochs (primary analysis)")
    res(OUTCOME, "+ covariate: % epochs with theta ≥ alpha", extra=" + pctA")
    res(OUTCOME, "+ covariate: % epochs with alpha attenuation", extra=" + pctB")
    res("beta_exclA", "Excluding epochs with theta ≥ alpha")
    res("beta_exclB", "Excluding epochs with alpha attenuation")
    res("beta_exclAB", "Excluding epochs meeting either criterion")
    res("beta_first", "First half of the recording (epochs 1–30)")
    res("beta_second", "Second half of the recording (epochs 31–60)")
    res(OUTCOME, "Participants with no flagged epoch (criterion A)", dd=have[have.pctA == 0])
    S = pd.DataFrame(S)
    S.to_csv(OUT / "02_sensitivity.csv", index=False)

    L = [
        "# S_R1s — 逐 epoch 的清醒度指標",
        "",
        f"- 分析樣本 {len(df)} 人；有逐 epoch 頻譜 {len(have)} 人（每人 {int(have.ep_total.median())} 個 epoch）。",
        f"- 核對：60 個 epoch 取平均後重算的週期性 beta，與主分析的值相關 r = {r_chk:.4f}，最大差 {max_dev:.4f}。",
        "",
        "## 1 各組的疑似嗜睡 epoch",
        "",
        "| 指標 | 健康組 平均（中位數） | 低嚴重度 平均（中位數）；g，p | 高嚴重度 平均（中位數）；g，p | 三組 KW p | 與週期性 beta 的 rho（p） |",
        "|---|---|---|---|---|---|",
    ]
    for r in A.itertuples():
        L.append(
            f"| {r.measure} | {r.HC_mean:.2f}（{r.HC_median:.2f}） | {r.Lower_mean:.2f}（{r.Lower_median:.2f}）；{r.Lower_g:+.2f}，{fmt_p(r.Lower_p)} | {r.Higher_mean:.2f}（{r.Higher_median:.2f}）；{r.Higher_g:+.2f}，{fmt_p(r.Higher_p)} | {fmt_p(r.kw_p)} | {r.rho_with_beta:+.2f}（{fmt_p(r.p_with_beta)}） |"
        )
    L += [
        "",
        "至少有 1 個 epoch 符合定義 A 的人："
        + "、".join(f"{g} {anyA[g]:.1f}%" for g in ("HC", "Lower", "Higher"))
        + "；超過一半 epoch 符合定義 A 的人數："
        + "、".join(f"{g} {allA[g]}" for g in ("HC", "Lower", "Higher"))
        + "。",
        "",
        "## 2 敏感度分析（週期性 beta）",
        "",
        "| 分析 | n（健康／低／高） | 高嚴重度 g [95% CI] | 低嚴重度 g [95% CI] | 校正後高嚴重度（p） | ISI ≥ 18（p） |",
        "|---|---|---|---|---|---|",
    ]
    for r in S.itertuples():
        L.append(
            f"| {r.analysis} | {r.n}（{r.n_hc}／{r.n_lower}／{r.n_higher}） | {r.g_higher:+.2f} [{r.g_higher_lo:+.2f}, {r.g_higher_hi:+.2f}] | {r.g_lower:+.2f} [{r.g_lower_lo:+.2f}, {r.g_lower_hi:+.2f}] | {r.adj_higher_in_HC_SD:+.2f}（{fmt_p(r.adj_higher_p)}） | {r.thr_in_HC_SD:+.2f}（{fmt_p(r.thr_p)}） |"
        )
    L += [
        "",
        "校正後與門檻的係數以健康組標準差為單位；模型校正年齡、性別、教育。重算時至少要剩 30 個 epoch，不足者不納入該列。",
    ]
    (OUT / "03_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
