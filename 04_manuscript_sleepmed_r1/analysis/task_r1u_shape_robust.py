"""補做 HC3、非線性檢定與合併效果量的人數換算。

樣條與 wild bootstrap 都以「關係是直線」為虛無假設。
task_r1a 的置換檢定問的是「有沒有關係」，兩者不能互相替代。
掃描切點限 11–22 分、兩側至少各 10 人；bootstrap 每次重新選點。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import dmatrix
from statsmodels.stats.power import TTestIndPower

from _r1_common import AN, COV, OUTCOME, load_all, fmt_p, hedges_g

OUT = AN / "S_R1u_shape_robust"
OUT.mkdir(parents=True, exist_ok=True)
Y = OUTCOME
SEED = 20261117
B = 4999
CAND = range(11, 23)
MIN_SIDE = 10
THR = 18  # task_r1a 選出的切點

df = load_all()
df["k2"] = pd.Categorical(df["k2"], ["HC", "Lower", "Higher"])
ins_all = df[df["ISI_T1"] >= 9].copy()
MOOD = "BDI_T1 + BAI_T1"


# 1. 穩健標準誤
def both_se(formula, data, term, label, scope):
    m = smf.ols(formula, data).fit()
    r = m.get_robustcov_results(cov_type="HC3", use_t=True)
    i = list(m.params.index).index(term)
    lo, hi = r.conf_int()[i]
    return dict(
        result=label,
        scope=scope,
        n=int(m.nobs),
        b=m.params[term],
        se_ols=m.bse[term],
        p_ols=m.pvalues[term],
        se_hc3=r.bse[i],
        ci_low_hc3=lo,
        ci_high_hc3=hi,
        p_hc3=r.pvalues[i],
    )


rows = []
full_m2 = df.dropna(subset=["BDI_T1", "BAI_T1"])
ins_m2 = ins_all.dropna(subset=["BDI_T1", "BAI_T1"])
for g in ["Higher", "Lower"]:
    rows.append(both_se(f"{Y} ~ {COV} + C(k2)", df, f"C(k2)[T.{g}]", f"{g} vs HC, M1", "full"))
    rows.append(
        both_se(f"{Y} ~ {COV} + {MOOD} + C(k2)", full_m2, f"C(k2)[T.{g}]", f"{g} vs HC, M2", "full")
    )
ins_all["thr"] = (ins_all["ISI_T1"] >= THR).astype(int)
ins_m2 = ins_m2.assign(thr=(ins_m2["ISI_T1"] >= THR).astype(int))
rows.append(
    both_se(f"{Y} ~ {COV} + ISI_T1", ins_all, "ISI_T1", "ISI total (per point), M1", "insomnia")
)
rows.append(both_se(f"{Y} ~ {COV} + thr", ins_all, "thr", f"Step at ISI >= {THR}, M1", "insomnia"))
rows.append(
    both_se(f"{Y} ~ {COV} + {MOOD} + thr", ins_m2, "thr", f"Step at ISI >= {THR}, M2", "insomnia")
)
rows.append(
    both_se(
        f"{Y} ~ {COV} + ISI_T1 + thr",
        ins_all,
        "thr",
        f"Step at ISI >= {THR} added to linear, M1 (uncorrected for selection)",
        "insomnia",
    )
)
rows.append(
    both_se(
        f"{Y} ~ {COV} + ISI_T1 + thr",
        ins_all,
        "ISI_T1",
        f"ISI total added to step at {THR}, M1",
        "insomnia",
    )
)
rows.append(
    both_se(
        f"{Y} ~ {COV} + {MOOD} + ISI_T1", ins_m2, "ISI_T1", "ISI total (per point), M2", "insomnia"
    )
)
# 慣用臨床切點（ISI >= 15）的分組，只作敏感度分析
df["g15"] = pd.Categorical(
    np.where(df["ISI_T1"] < 9, "HC", np.where(df["ISI_T1"] >= 15, "Higher", "Lower")),
    ["HC", "Lower", "Higher"],
)
full_m2 = df.dropna(subset=["BDI_T1", "BAI_T1"])
for g in ["Higher", "Lower"]:
    rows.append(
        both_se(
            f"{Y} ~ {COV} + C(g15)",
            df,
            f"C(g15)[T.{g}]",
            f"ISI >= 15 grouping: {g} vs HC, M1",
            "full",
        )
    )
    rows.append(
        both_se(
            f"{Y} ~ {COV} + {MOOD} + C(g15)",
            full_m2,
            f"C(g15)[T.{g}]",
            f"ISI >= 15 grouping: {g} vs HC, M2",
            "full",
        )
    )
t1 = pd.DataFrame(rows)
t1.to_csv(OUT / "01_robust_se.csv", index=False, encoding="utf-8-sig")


# 2. 型態檢定
def design(data, mood):
    f = f"{COV}" + (f" + {MOOD}" if mood else "")
    C = np.asarray(dmatrix(f, data, return_type="dataframe"))
    return C


def perp(Z, X):
    """把 Z 對 X 正交化，並去掉已經落在 X 裡的方向（樣條基底包含直線）。"""
    Q = np.linalg.qr(X, mode="reduced")[0]
    R = Z - Q @ (Q.T @ Z)
    u, s, _ = np.linalg.svd(R, full_matrices=False)
    keep = s > 1e-8 * s.max()
    return u[:, keep]


def spline_vs_linear(data, mood):
    import statsmodels.api as sm

    y = data[Y].to_numpy()
    isi = data["ISI_T1"].to_numpy(float)
    X = np.column_stack([design(data, mood), isi])
    S = np.asarray(dmatrix("cr(x, df=3, constraints='center') - 1", {"x": isi}))
    Z = perp(S, X)
    m0 = sm.OLS(y, X).fit()
    m1 = sm.OLS(y, np.column_stack([X, Z])).fit()
    k = Z.shape[1]
    R = np.zeros((k, X.shape[1] + k))
    R[:, -k:] = np.eye(k)
    w = m1.get_robustcov_results(cov_type="HC3", use_t=True).f_test(R)
    return dict(
        test="spline (3 df) vs linear",
        added_df=k,
        p_hc3=float(w.pvalue),
        p_ols=float(m1.compare_f_test(m0)[1]),
        delta_r2=m1.rsquared - m0.rsquared,
        aic_linear=m0.aic,
        aic_alternative=m1.aic,
    )


def step_vs_linear(data, mood, cand=CAND):
    """挑出最大的階梯，並在「直線」的虛無假設下用 wild bootstrap 校正挑選。"""
    import statsmodels.api as sm

    y = data[Y].to_numpy()
    isi = data["ISI_T1"].to_numpy(float)
    X = np.column_stack([design(data, mood), isi])
    n, p = X.shape
    Q = np.linalg.qr(X, mode="reduced")[0]
    cs, Z = [], []
    for c in cand:
        z = (isi >= c).astype(float)
        if min(z.sum(), n - z.sum()) < MIN_SIDE:
            continue
        z = z - Q @ (Q.T @ z)
        Z.append(z / np.linalg.norm(z))
        cs.append(c)
    Z = np.column_stack(Z)

    def fmax(Ymat):
        E = Ymat - Q @ (Q.T @ Ymat)
        ss = (Z.T @ E) ** 2
        mx = ss.max(axis=0)
        return mx / ((np.sum(E * E, axis=0) - mx) / (n - p - 1)), ss.argmax(axis=0)

    obs, arg = fmax(y[:, None])
    m0 = sm.OLS(y, X).fit()
    h = np.sum(Q * Q, axis=1)
    e = m0.resid / np.sqrt(1 - h)
    rng = np.random.default_rng(SEED)
    sims = []
    for start in range(0, B, 500):
        k = min(500, B - start)
        Yb = m0.fittedvalues[:, None] + e[:, None] * rng.choice([-1.0, 1.0], size=(n, k))
        sims.append(fmax(Yb)[0])
    sims = np.concatenate(sims)
    return dict(
        test=f"selected step vs linear (candidates {min(cs)}-{max(cs)}; wild bootstrap, B = {B})",
        selected_cutoff=int(cs[int(arg[0])]),
        f_max=float(obs[0]),
        p_selection_corrected=float((1 + np.sum(sims >= obs[0])) / (B + 1)),
    )


rows = []
for scope, d0 in [("insomnia", ins_all), ("full", df)]:
    for mood in (False, True):
        d = d0.dropna(subset=["BDI_T1", "BAI_T1"]) if mood else d0
        tag = "M2" if mood else "M1"
        rows.append(dict(scope=scope, adjustment=tag, n=len(d), **spline_vs_linear(d, mood)))
        rows.append(dict(scope=scope, adjustment=tag, n=len(d), **step_vs_linear(d, mood)))
        if scope == "full":  # 全樣本另報較寬的候選範圍（5–22），供對照
            rows.append(
                dict(scope=scope, adjustment=tag, n=len(d), **step_vs_linear(d, mood, range(5, 23)))
            )
t2 = pd.DataFrame(rows)
t2.to_csv(OUT / "02_shape_tests.csv", index=False, encoding="utf-8-sig")


# 3. 合併的代價
curve = pd.read_csv(AN / "S_R1g_followup" / "05a_cutoff_dilution_curve.csv", encoding="utf-8-sig")
# 效果量在這裡直接重算：CSV 裡的數值已四捨五入到小數三位，再取兩位會多進一位（ISI >= 15 是 0.6446，不是 0.65）
hc_y = df.loc[df["ISI_T1"] < 9, Y]
cuts = curve["cutoff"].str.extract(r"(\d+)")[0].astype(int)
est = [hedges_g(df.loc[df["ISI_T1"] >= c, Y], hc_y) for c in cuts]
assert np.allclose([e[0] for e in est], curve["g"], atol=6e-4), "與 task_r1g 的效果量不一致"
curve["g"], curve["ci_low"], curve["ci_high"] = (
    [e[0] for e in est],
    [e[1] for e in est],
    [e[2] for e in est],
)
exc = [
    hedges_g(df.loc[(df["ISI_T1"] >= 9) & (df["ISI_T1"] < c), Y], hc_y)
    if c > 9
    else (np.nan, np.nan, np.nan, 0, 0)
    for c in cuts
]
curve["n_excluded"], curve["g_excluded"] = [e[3] for e in exc], [e[0] for e in exc]
curve["ci_low_excluded"], curve["ci_high_excluded"] = [e[1] for e in exc], [e[2] for e in exc]
pw = TTestIndPower()
curve["n_per_group_80pct_power"] = [
    int(np.ceil(pw.solve_power(effect_size=g, alpha=0.05, power=0.80))) for g in curve["g"]
]
curve["share_of_symptom_group"] = (curve["n_a"] / curve["n_a"].iloc[0]).round(3)
t3 = curve[
    [
        "cutoff",
        "n_a",
        "share_of_symptom_group",
        "g",
        "ci_low",
        "ci_high",
        "p",
        "n_excluded",
        "g_excluded",
        "ci_low_excluded",
        "ci_high_excluded",
        "n_per_group_80pct_power",
    ]
]
t3.to_csv(OUT / "03_pooling_cost.csv", index=False, encoding="utf-8-sig")


# 4. 摘要
def row(scope, adj, key):
    return t2[
        (t2["scope"] == scope) & (t2["adjustment"] == adj) & t2["test"].str.startswith(key)
    ].iloc[0]


L = [
    "# 型態檢定、穩健標準誤與合併的代價（task_r1u）",
    "",
    "## 1. 主要結果在穩健標準誤（HC3）下",
    "",
    "| 結果 | n | b | 一般 p | HC3 95% CI | HC3 p |",
    "|---|---|---|---|---|---|",
]
for _, r in t1.iterrows():
    L.append(
        f"| {r['result']} | {r['n']} | {r['b']:.4f} | {fmt_p(r['p_ols'])} | [{r['ci_low_hc3']:.4f}, {r['ci_high_hc3']:.4f}] | {fmt_p(r['p_hc3'])} |"
    )
L += [
    "",
    "## 2. 關係是不是直線",
    "",
    "| 範圍 | 模型 | n | 樣條對直線（HC3 p） | 挑出的切點 | 階梯對直線（校正挑切點的 p） |",
    "|---|---|---|---|---|---|",
]
for scope in ("insomnia", "full"):
    for adj in ("M1", "M2"):
        a = row(scope, adj, "spline")
        b = row(scope, adj, "selected step vs linear (candidates 11")
        L.append(
            f"| {scope} | {adj} | {a['n']} | {fmt_p(a['p_hc3'])} | {int(b['selected_cutoff'])} | {fmt_p(b['p_selection_corrected'])} |"
        )
L += [
    "",
    "全樣本另以 5–22 分為候選範圍的結果見 02_shape_tests.csv。",
    "",
    "## 3. 合併的代價",
    "",
    "| 失眠症狀的定義 | 人數 | 對健康組的 g [95% CI] | 每組所需人數（檢定力 .80） |",
    "|---|---|---|---|",
]
for _, r in t3.iterrows():
    L.append(
        f"| {r['cutoff']} | {r['n_a']} | {r['g']:+.2f} [{r['ci_low']:+.2f}, {r['ci_high']:+.2f}] | {r['n_per_group_80pct_power']} |"
    )
L += [
    "",
    "ISI ≥ 18 是由本樣本選出的切點，該列的效果量偏樂觀；ISI ≥ 15 是慣用的臨床切點，沒有這個問題。",
]
(OUT / "04_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))


# 2026-11-17 補：高分端的細分，以及拿掉最高分之後非線性是否還在
# 圖 2 的樣條曲線在 ISI 22 分以上下彎。這裡把 19 分以上再細分，並在排除 ISI ≥ 23 的參與者後重做樣條對直線的檢定，
# 用來確認「關係不是直線」不是只靠最高分那幾個人。
upper = []
for lab, mask in [
    ("17-18", (df["ISI_T1"] >= 17) & (df["ISI_T1"] <= 18)),
    ("19-20", (df["ISI_T1"] >= 19) & (df["ISI_T1"] <= 20)),
    ("21-22", (df["ISI_T1"] >= 21) & (df["ISI_T1"] <= 22)),
    (">=23", df["ISI_T1"] >= 23),
]:
    g_, lo_, hi_, n_, _ = hedges_g(df.loc[mask, Y], hc_y)
    upper.append(
        dict(
            ISI_range=lab,
            n=n_,
            mean_beta=df.loc[mask, Y].mean(),
            g_vs_HC=g_,
            ci_low=lo_,
            ci_high=hi_,
        )
    )
t5 = pd.DataFrame(upper)
t5.to_csv(OUT / "05_upper_range.csv", index=False, encoding="utf-8-sig")
sub22 = ins_all[ins_all["ISI_T1"] <= 22]
ex23 = dict(
    scope="insomnia, ISI <= 22", adjustment="M1", n=len(sub22), **spline_vs_linear(sub22, False)
)
pd.DataFrame([ex23]).to_csv(
    OUT / "06_spline_excluding_isi23plus.csv", index=False, encoding="utf-8-sig"
)
L += [
    "",
    "## 4. 高分端的細分",
    "",
    "| ISI | 人數 | 平均 beta | 對健康組的 g [95% CI] |",
    "|---|---|---|---|",
]
for _, r in t5.iterrows():
    L.append(
        f"| {r['ISI_range']} | {r['n']} | {r['mean_beta']:.3f} | {r['g_vs_HC']:+.2f} [{r['ci_low']:+.2f}, {r['ci_high']:+.2f}] |"
    )
L += [
    "",
    f"排除 ISI ≥ 23 的參與者後（{ex23['n']} 人），樣條對直線：HC3 p = {fmt_p(ex23['p_hc3'])}（一般 p = {fmt_p(ex23['p_ols'])}）。",
    "ISI 23 分以上只有 8 人，估計不精確；不能據此說最高分端沒有升高，也不能說有。",
]
(OUT / "04_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L[-10:]))


# 2026-11-17 補：最高分端在「不依賴峰值偵測」的指標上是否也下彎
# 圖 2 用的是峰值功率（有偵測到峰的電極取最大峰的高度）。這裡用同樣的分段，改看去背景 beta 功率（13–30 Hz 的頻帶平均）
# 與傳統相對 beta 功率，確認 ISI ≥ 23 的下彎是不是峰值功率這個指標特有的。
from _r1_spectra import band_mean, load_spectra

_freq, _chn, _adj, _raw = load_spectra(df["subject"])
df["adj_13_30"] = band_mean(_adj, _freq, 13, 30)
alt = []
for lab, mask in [
    ("9-16", (df["ISI_T1"] >= 9) & (df["ISI_T1"] <= 16)),
    ("17-18", (df["ISI_T1"] >= 17) & (df["ISI_T1"] <= 18)),
    ("19-22", (df["ISI_T1"] >= 19) & (df["ISI_T1"] <= 22)),
    (">=23", df["ISI_T1"] >= 23),
]:
    r = dict(ISI_range=lab, n=int(mask.sum()))
    for col, name in [
        (Y, "peak_power"),
        ("adj_13_30", "aperiodic_adjusted_13_30"),
        ("Beta_rel_BRAIN", "relative_beta"),
    ]:
        g_, lo_, hi_, _, _ = hedges_g(df.loc[mask, col], df.loc[df["ISI_T1"] < 9, col])
        r[f"g_{name}"], r[f"ci_low_{name}"], r[f"ci_high_{name}"] = g_, lo_, hi_
    alt.append(r)
t7 = pd.DataFrame(alt)
t7.to_csv(OUT / "07_upper_range_other_measures.csv", index=False, encoding="utf-8-sig")
L += [
    "",
    "## 5. 高分端在其他 beta 指標上的效果量（對健康組的 g [95% CI]）",
    "",
    "| ISI | 人數 | 峰值功率 | 去背景功率 13–30 Hz | 相對 beta |",
    "|---|---|---|---|---|",
]
for _, r in t7.iterrows():
    L.append(
        f"| {r['ISI_range']} | {r['n']} | "
        + " | ".join(
            f"{r['g_' + k]:+.2f} [{r['ci_low_' + k]:+.2f}, {r['ci_high_' + k]:+.2f}]"
            for k in ["peak_power", "aperiodic_adjusted_13_30", "relative_beta"]
        )
        + " |"
    )
L += [
    "",
    "ISI ≥ 23 的下彎主要出現在峰值功率；以不依賴峰值偵測的指標來看，這 8 人仍高於健康組，但區間很寬。",
]
(OUT / "04_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L[-9:]))


# 2026-11-17 補：主要與次要 beta 指標各自的型態檢定（同一套樣條、同一批人）
# Methods 2.6 的次要指標（去背景 beta 功率、傳統相對 beta 功率）也做同樣的樣條對直線檢定，並列出最高分端的峰寬，
# 供補充圖 S2（三個指標並排的嚴重度曲線）與補充表 S2 使用。
df["raw_13_30"] = band_mean(_raw, _freq, 13, 30)
ins_m = df[df["ISI_T1"] >= 9]


def _spl(d, col):
    d = d.dropna(subset=[col])
    return spline_vs_linear(d.rename(columns={col: "__y"}).assign(**{Y: d[col].to_numpy()}), False)


by_measure = []
for col, name in [
    (Y, "periodic beta peak power (primary)"),
    ("adj_13_30", "aperiodic-adjusted beta power, 13-30 Hz"),
    ("Beta_rel_BRAIN", "relative beta power"),
    ("raw_13_30", "absolute beta power (log), 13-30 Hz"),
]:
    a_ = _spl(ins_m, col)
    b_ = _spl(ins_m[ins_m["ISI_T1"] <= 22], col)
    by_measure.append(
        dict(
            measure=name,
            n=int(ins_m[col].notna().sum()),
            spline_vs_linear_p_hc3=a_["p_hc3"],
            excluding_isi23plus_p_hc3=b_["p_hc3"],
        )
    )
t8m = pd.DataFrame(by_measure)
t8m.to_csv(OUT / "08_shape_by_measure.csv", index=False, encoding="utf-8-sig")
bw = {
    lab: hedges_g(df.loc[m_, "beta_bw_BRAIN"], df.loc[df["ISI_T1"] < 9, "beta_bw_BRAIN"])[0]
    for lab, m_ in [
        ("19-22", (df["ISI_T1"] >= 19) & (df["ISI_T1"] <= 22)),
        (">=23", df["ISI_T1"] >= 23),
    ]
}
al = {
    lab: hedges_g(df.loc[m_, "Alpha_rel_BRAIN"], df.loc[df["ISI_T1"] < 9, "Alpha_rel_BRAIN"])[0]
    for lab, m_ in [
        ("19-22", (df["ISI_T1"] >= 19) & (df["ISI_T1"] <= 22)),
        (">=23", df["ISI_T1"] >= 23),
    ]
}
L += [
    "",
    "## 6. 各 beta 指標的樣條對直線（失眠症狀者；HC3 p）",
    "",
    "| 指標 | n | 全部 | 排除 ISI ≥ 23 |",
    "|---|---|---|---|",
]
for _, r in t8m.iterrows():
    L.append(
        f"| {r['measure']} | {r['n']} | {fmt_p(r['spline_vs_linear_p_hc3'])} | {fmt_p(r['excluding_isi23plus_p_hc3'])} |"
    )
L += [
    "",
    f"beta 峰寬對健康組的 g：19–22 分 {bw['19-22']:+.2f}、23 分以上 {bw['>=23']:+.2f}（最高分端的峰較寬、較低）。",
    f"相對 alpha 對健康組的 g：19–22 分 {al['19-22']:+.2f}、23 分以上 {al['>=23']:+.2f}（相對 beta 持平，部分來自 alpha 較低，不能當成 beta 沒有回落的證據）。",
]
(OUT / "04_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L[-9:]))
