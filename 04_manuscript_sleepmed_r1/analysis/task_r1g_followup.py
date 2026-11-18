"""補查斜率、峰偵測、電極子集、aperiodic、納入切點與傳統 beta。

分群只用 ISI。aperiodic 的重新擬合沿用 task_b 的設定，只改擬合頻率上限。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from fooof import FOOOFGroup
from scipy import stats
from scipy.io import loadmat
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

from _r1_common import AN, COV, COV_COLS, ISI_ITEMS, OUTCOME, SEED, hedges_g, load_all, mwu_p

Y = OUTCOME
THR = 18
B = 5000
OUT = AN / "S_R1g_followup"
OUT.mkdir(parents=True, exist_ok=True)
MD = []
rng = np.random.default_rng(20261116)


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    MD.append(s)


def gci(a, b):
    g, lo, hi, n1, n2 = hedges_g(a, b)
    return dict(
        n_a=n1,
        n_HC=n2,
        g=round(g, 3),
        ci_low=round(lo, 3),
        ci_high=round(hi, 3),
        p=round(mwu_p(a, b), 4),
    )


df = load_all()
df["logBeta_abs_BRAIN"] = np.log10(df["Beta_abs_BRAIN"])
ch = pd.read_csv(AN / "S10_fooof_full" / "01_specparam_features.csv")
qc = pd.read_csv(AN / "S_R1d_preproc_qc" / "_subject_level_internal.csv")[
    ["subject", "n_ret_muscle_dom", "max_ret_muscle_p"]
]
df = df.merge(qc, on="subject", how="left")
hc = df[df["group"] == "HC"]
ins = df[df["group"] != "HC"]
GR = {
    "Lower": ins[ins["k2"] == "Lower"],
    "Higher": ins[ins["k2"] == "Higher"],
    f"ISI 9-{THR - 1}": ins[ins["ISI_T1"] < THR],
    f"ISI >={THR}": ins[ins["ISI_T1"] >= THR],
}
d_all = df.dropna(subset=COV_COLS).copy()
d_all["k2c"] = pd.Categorical(d_all["k2"], categories=["HC", "Lower", "Higher"])
P("# S_R1g — 文獻追查後的補分析\n")
P(f"N = {len(df)}；HC = {len(hc)}；低嚴重度 = {len(GR['Lower'])}；高嚴重度 = {len(GR['Higher'])}\n")


# 快速門檻掃描（Frisch–Waugh：先把共變項的影響拿掉，再算相關）
def cmat(d):
    return np.column_stack(
        [
            np.ones(len(d)),
            d["AGE"].to_numpy(float),
            d["SEX"].to_numpy(float),
            d["EDU"].to_numpy(float),
        ]
    )


def resid(M, C):
    return M - C @ np.linalg.lstsq(C, M, rcond=None)[0]


def scan_fast(y, isi, C, cands, min_n=10):
    D = np.column_stack([(isi >= c).astype(float) for c in cands])
    ok = (D.sum(0) >= min_n) & ((1 - D).sum(0) >= min_n)
    ry = resid(y[:, None], C)[:, 0]
    rD = resid(D, C)
    den = np.sqrt((rD**2).sum(0) * (ry**2).sum())
    r = np.where(ok & (den > 0), (rD.T @ ry) / np.where(den > 0, den, 1), 0.0)
    F = r**2 * (len(y) - C.shape[1] - 1) / (1 - r**2)
    j = int(np.argmax(F))
    return cands[j], float(F[j]), float(np.sign(r[j]))


# G1. 斜率與折線模型
di = ins.dropna(subset=[Y] + COV_COLS).copy()
sd_hc = hc[Y].std(ddof=1)
rows = []
for label, d in [
    ("all insomnia (ISI 9-26)", di),
    (f"below threshold (ISI 9-{THR - 1})", di[di["ISI_T1"] < THR]),
    (f"at/above threshold (ISI >={THR})", di[di["ISI_T1"] >= THR]),
    ("K=2 lower-severity cluster", di[di["k2"] == "Lower"]),
    ("K=2 higher-severity cluster", di[di["k2"] == "Higher"]),
]:
    m = smf.ols(f"{Y} ~ {COV} + ISI_T1", d).fit()
    c = m.conf_int().loc["ISI_T1"]
    rows.append(
        dict(
            segment=label,
            n=int(m.nobs),
            slope_per_ISI_point=round(m.params["ISI_T1"], 4),
            ci_low=round(c[0], 4),
            ci_high=round(c[1], 4),
            p=round(m.pvalues["ISI_T1"], 4),
            slope_in_HC_SD_per_5_points=round(5 * m.params["ISI_T1"] / sd_hc, 3),
            ci_low_SD5=round(5 * c[0] / sd_hc, 3),
            ci_high_SD5=round(5 * c[1] / sd_hc, 3),
        )
    )
slopes = pd.DataFrame(rows)


def best(d, make, cands):
    out = None
    for c in cands:
        dd = d.assign(**make(d, c))
        if (dd["ISI_T1"] >= c).sum() < 10 or (dd["ISI_T1"] < c).sum() < 10:
            continue
        m = smf.ols(f"{Y} ~ {COV} + " + " + ".join(make(d, c).keys()), dd).fit()
        if out is None or m.aic < out[1].aic:
            out = (c, m)
    return out


def forms(d, cands, scope):
    r = []

    def add(name, m, extra=0, c=""):
        r.append(
            dict(
                scope=scope,
                model=name,
                selected_c=c,
                n=int(m.nobs),
                k=int(m.df_model) - 3 + extra,
                R2=round(m.rsquared, 4),
                AIC=round(m.aic + 2 * extra, 2),
                BIC=round(m.bic + np.log(m.nobs) * extra, 2),
            )
        )

    add("linear", smf.ols(f"{Y} ~ {COV} + ISI_T1", d).fit())
    add("spline df=4", smf.ols(f"{Y} ~ {COV} + cr(ISI_T1, df=4)", d).fit())
    c, m = best(d, lambda x, c: {"step": (x["ISI_T1"] >= c).astype(float)}, cands)
    add("step (flat, jump, flat)", m, 1, c)
    c, m = best(d, lambda x, c: {"hinge": np.clip(x["ISI_T1"] - c, 0, None)}, cands)
    add("hinge (flat, then linear rise)", m, 1, c)
    c, m = best(
        d, lambda x, c: {"lin": x["ISI_T1"], "hinge": np.clip(x["ISI_T1"] - c, 0, None)}, cands
    )
    add("broken stick (two slopes, continuous)", m, 1, c)
    c, m = best(
        d, lambda x, c: {"lin": x["ISI_T1"], "step": (x["ISI_T1"] >= c).astype(float)}, cands
    )
    add("linear + step", m, 1, c)
    c, m = best(
        d,
        lambda x, c: {
            "step": (x["ISI_T1"] >= c).astype(float),
            "lo": np.where(x["ISI_T1"] < c, x["ISI_T1"] - c, 0.0),
            "hi": np.where(x["ISI_T1"] >= c, x["ISI_T1"] - c, 0.0),
        },
        cands,
    )
    add("step + separate slopes", m, 1, c)
    return r


fm = pd.DataFrame(
    forms(di, range(11, 23), "insomnia only")
    + forms(df.dropna(subset=[Y] + COV_COLS), range(5, 23), "full sample (HC included)")
)
slopes.to_csv(OUT / "01a_segment_slopes.csv", index=False, encoding="utf-8-sig")
fm.to_csv(OUT / "01b_functional_forms.csv", index=False, encoding="utf-8-sig")
P("## G1. 門檻兩側的斜率，與折線模型的比較\n")
P("slope_in_HC_SD_per_5_points＝ISI 每增加 5 分，beta 變化幾個健康組標準差。\n")
P(slopes.to_markdown(index=False))
P("")
P("各種函數形式（資料選出的轉折點多算一個參數；AIC 越小越好）\n")
P(fm.to_markdown(index=False))
P("")

# G2. beta 峰的偵測
pk = (
    ch.groupby("subject")
    .agg(
        n_ch=("channel", "size"),
        n_beta_peak=("beta_pw", lambda s: s.notna().sum()),
        beta_pw_zero_imputed=("beta_pw", lambda s: s.fillna(0).mean()),
    )
    .reset_index()
)
pk["prop_beta_peak"] = pk["n_beta_peak"] / pk["n_ch"]
d2 = df.merge(pk, on="subject", how="left")
h2 = d2[d2["group"] == "HC"]
rows = []
for f in ["prop_beta_peak", "beta_cf_BRAIN", "beta_bw_BRAIN", "beta_pw_zero_imputed", Y]:
    H, pkw = stats.kruskal(
        *[d2.loc[d2["k2"] == g_, f].dropna() for g_ in ["HC", "Lower", "Higher"]]
    )
    r = dict(
        feature=f,
        HC_mean=round(h2[f].mean(), 3),
        Lower_mean=round(d2.loc[d2["k2"] == "Lower", f].mean(), 3),
        Higher_mean=round(d2.loc[d2["k2"] == "Higher", f].mean(), 3),
        KW_H=round(H, 2),
        KW_p=round(pkw, 4),
    )
    for g_ in ["Lower", "Higher"]:
        x = gci(d2.loc[d2["k2"] == g_, f], h2[f])
        r.update(
            {f"g_{g_}": x["g"], f"ci_{g_}": f"[{x['ci_low']}, {x['ci_high']}]", f"p_{g_}": x["p"]}
        )
    rows.append(r)
g2 = pd.DataFrame(rows)
dd = d2.dropna(subset=[Y] + COV_COLS).copy()
dd["k2c"] = pd.Categorical(dd["k2"], categories=["HC", "Lower", "Higher"])
m0 = smf.ols(f"{Y} ~ {COV} + C(k2c)", dd).fit()
m1 = smf.ols(f"{Y} ~ {COV} + prop_beta_peak + beta_cf_BRAIN + beta_bw_BRAIN + C(k2c)", dd).fit()
full19 = dd[dd["n_beta_peak"] == dd["n_ch"]]
g2.to_csv(OUT / "02_beta_peak_detection.csv", index=False, encoding="utf-8-sig")
P("## G2. beta 峰的偵測比例、峰頻率、峰寬\n")
P(g2.to_markdown(index=False))
P("")
P(
    f"- 高嚴重度群對 HC 的調整後差距：原本 b = {m0.params['C(k2c)[T.Higher]']:.4f}；加入偵測比例、峰頻率、峰寬後 b = {m1.params['C(k2c)[T.Higher]']:.4f}（p = {m1.pvalues['C(k2c)[T.Higher]']:.4f}）"
)
x = gci(full19.loc[full19["k2"] == "Higher", Y], full19.loc[full19["k2"] == "HC", Y])
y_ = gci(full19.loc[full19["k2"] == "Lower", Y], full19.loc[full19["k2"] == "HC", Y])
P(
    f"- 只留 19 個電極全部偵測到 beta 峰的人：高嚴重度 g = {x['g']} [{x['ci_low']}, {x['ci_high']}]（{x['n_a']} 對 {x['n_HC']}）；低嚴重度 g = {y_['g']} [{y_['ci_low']}, {y_['ci_high']}]（{y_['n_a']} 人）\n"
)

# G3. 電極子集
SETS = {
    "all 19 electrodes": None,
    "core 9 (F3 F4 Fz C3 C4 Cz P3 P4 Pz)": ["F3", "F4", "FZ", "C3", "C4", "CZ", "P3", "P4", "PZ"],
    "edge 10 (Fp1 Fp2 F7 F8 T3 T4 T5 T6 O1 O2)": [
        "FP1",
        "FP2",
        "F7",
        "F8",
        "T3",
        "T4",
        "T5",
        "T6",
        "O1",
        "O2",
    ],
    "midline 3 (Fz Cz Pz)": ["FZ", "CZ", "PZ"],
    "vertex (Cz)": ["CZ"],
}
rows = []
for name, chs in SETS.items():
    sub = ch if chs is None else ch[ch["channel"].isin(chs)]
    v = (
        sub.groupby("subject")["beta_pw"]
        .apply(lambda s: s.dropna().mean())
        .rename("v")
        .reset_index()
    )
    t = df[["subject", "k2", "ISI_T1", "group"]].merge(v, on="subject", how="left")
    h = t.loc[t["group"] == "HC", "v"]
    for lab, mask in [
        ("Higher vs HC", t["k2"] == "Higher"),
        ("Lower vs HC", t["k2"] == "Lower"),
        (f"ISI >={THR} vs HC", (t["group"] != "HC") & (t["ISI_T1"] >= THR)),
    ]:
        rows.append(dict(electrodes=name, contrast=lab, **gci(t.loc[mask, "v"], h)))
g3 = pd.DataFrame(rows)
g3.to_csv(OUT / "03_electrode_subsets.csv", index=False, encoding="utf-8-sig")
P("## G3. 電極子集的 periodic beta\n")
P(g3.to_markdown(index=False))
P("")

# G4. aperiodic
AP = ["exponent_BRAIN", "offset_BRAIN", "exponent_CENTRAL", "offset_CENTRAL"]
rows = []
for f in AP:
    sdh = hc[f].std(ddof=1)
    for name, extra, d in [
        ("age, sex, education", "", d_all),
        ("+ BDI-II + BAI", " + BDI_T1 + BAI_T1", d_all.dropna(subset=["BDI_T1", "BAI_T1"])),
        (
            "+ max retained muscle probability",
            " + max_ret_muscle_p",
            d_all.dropna(subset=["max_ret_muscle_p"]),
        ),
        (
            "+ BDI-II + BAI + muscle probability",
            " + BDI_T1 + BAI_T1 + max_ret_muscle_p",
            d_all.dropna(subset=["BDI_T1", "BAI_T1", "max_ret_muscle_p"]),
        ),
    ]:
        m = smf.ols(f"{f} ~ {COV}{extra} + C(k2c)", d).fit()
        r = dict(feature=f, covariates=name, n=int(m.nobs))
        for g_ in ["Higher", "Lower"]:
            t = f"C(k2c)[T.{g_}]"
            r[f"{g_}_b_in_HC_SD"] = round(m.params[t] / sdh, 3)
            r[f"{g_}_p"] = round(m.pvalues[t], 4)
        rows.append(r)
g4a = pd.DataFrame(rows)
clean = df[df["n_ret_muscle_dom"] == 0]
rows = []
for f in AP:
    for g_ in ["Higher", "Lower"]:
        rows.append(
            dict(
                feature=f,
                sample="no retained muscle-dominant component",
                contrast=f"{g_} vs HC",
                **gci(clean.loc[clean["k2"] == g_, f], clean.loc[clean["k2"] == "HC", f]),
            )
        )
g4b = pd.DataFrame(rows)

# 重新擬合：只改頻率範圍
PSD_DIR = AN / "S10_fooof_psd_input"
CH_RENAME = {"T7": "T3", "T8": "T4", "P7": "T5", "P8": "T6"}
CENTRAL = ["C3", "C4", "CZ"]


def refit(freq_range):
    out = []
    for sid in df["subject"]:
        m = loadmat(str(PSD_DIR / f"psd_{sid}.mat"), squeeze_me=True)
        freqs = np.asarray(m["freqs"], float).ravel()
        psd = np.asarray(m["psd"], float)
        names = [
            CH_RENAME.get(str(c).strip().upper(), str(c).strip().upper())
            for c in np.atleast_1d(m["ch_names"])
        ]
        fg = FOOOFGroup(
            peak_width_limits=(1.0, 12.0),
            max_n_peaks=8,
            min_peak_height=0.1,
            peak_threshold=2.0,
            aperiodic_mode="fixed",
            verbose=False,
        )
        fg.fit(freqs, psd, freq_range)
        rec = []
        for i, c in enumerate(names):
            fmod = fg.get_fooof(i)
            pks = (
                [(cf, pw) for cf, pw, bw in fmod.peak_params_ if 13.0 <= cf <= 30.0]
                if len(fmod.peak_params_)
                else []
            )
            rec.append(
                (
                    c,
                    fmod.aperiodic_params_[1],
                    fmod.aperiodic_params_[0],
                    max(pks, key=lambda t: t[1])[1] if pks else np.nan,
                    fmod.r_squared_,
                )
            )
        t = pd.DataFrame(rec, columns=["channel", "exponent", "offset", "beta_pw", "r2"])
        tc = t[t["channel"].isin(CENTRAL)]
        out.append(
            dict(
                subject=sid,
                exponent_BRAIN=t["exponent"].mean(),
                offset_BRAIN=t["offset"].mean(),
                beta_pw_BRAIN=t["beta_pw"].dropna().mean(),
                exponent_CENTRAL=tc["exponent"].mean(),
                offset_CENTRAL=tc["offset"].mean(),
                beta_pw_CENTRAL=tc["beta_pw"].dropna().mean(),
                r2=t["r2"].mean(),
            )
        )
    return pd.DataFrame(out)


rows = []
for rngHz in [(1.0, 45.0), (1.0, 35.0), (1.0, 30.0)]:
    rf = refit(rngHz).merge(df[["subject", "k2"]], on="subject")
    lab = f"{rngHz[0]:.0f}–{rngHz[1]:.0f} Hz"
    for f in ["beta_pw_BRAIN", "beta_pw_CENTRAL"] + AP:
        pair = pd.DataFrame(
            {"new": rf.set_index("subject")[f], "orig": df.set_index("subject")[f]}
        ).dropna()
        r = dict(
            fit_range=lab,
            feature=f,
            mean_R2=round(rf["r2"].mean(), 3),
            r_with_original=round(float(np.corrcoef(pair["new"], pair["orig"])[0, 1]), 3),
        )
        for g_ in ["Higher", "Lower"]:
            x = gci(rf.loc[rf["k2"] == g_, f], rf.loc[rf["k2"] == "HC", f])
            r.update(
                {
                    f"g_{g_}": x["g"],
                    f"ci_{g_}": f"[{x['ci_low']}, {x['ci_high']}]",
                    f"p_{g_}": x["p"],
                }
            )
        rows.append(r)
g4c = pd.DataFrame(rows)
g4a.to_csv(OUT / "04a_aperiodic_covariates.csv", index=False, encoding="utf-8-sig")
g4b.to_csv(OUT / "04b_aperiodic_no_muscle.csv", index=False, encoding="utf-8-sig")
g4c.to_csv(OUT / "04c_refit_frequency_range.csv", index=False, encoding="utf-8-sig")
P("## G4. aperiodic exponent 與 offset\n")
P("(a) 加入共變項後，兩群對 HC 的差距（以 HC 標準差為單位）\n")
P(g4a.to_markdown(index=False))
P("")
P("(b) 只留沒有殘留肌肉成分的人\n")
P(g4b.to_markdown(index=False))
P("")
P(
    "(c) 換擬合頻率範圍後重算（1–45 Hz 是原設定，用來確認重現；r_with_original＝與投稿版數值的相關）\n"
)
P(g4c.to_markdown(index=False))
P("")

# G5. 失眠切點
hcv = hc[Y].dropna()
rows = []
for c in range(9, 21):
    a = ins[ins["ISI_T1"] >= c]
    b_ = ins[ins["ISI_T1"] < c]
    r = dict(cutoff=f"ISI >= {c}", **{k: v for k, v in gci(a[Y], hcv).items() if k != "n_HC"})
    if len(b_) >= 10:
        x = gci(b_[Y], hcv)
        r.update(n_below=x["n_a"], g_below=x["g"], ci_below=f"[{x['ci_low']}, {x['ci_high']}]")
    rows.append(r)
g5a = pd.DataFrame(rows)
rows = []
for c in [9, 10, 11, 12]:
    s = ins[ins["ISI_T1"] >= c]
    X = StandardScaler().fit_transform(SimpleImputer(strategy="median").fit_transform(s[ISI_ITEMS]))
    lab = KMeans(n_clusters=2, random_state=SEED, n_init=100).fit_predict(X)
    hi_lab = s.groupby(lab)["ISI_T1"].mean().idxmax()
    for nm, mask in [("Higher", lab == hi_lab), ("Lower", lab != hi_lab)]:
        rows.append(
            dict(
                insomnia_definition=f"ISI >= {c}",
                cluster=nm,
                ISI_mean=round(s.loc[mask, "ISI_T1"].mean(), 2),
                **gci(s.loc[mask, Y], hcv),
            )
        )
g5b = pd.DataFrame(rows)
g5a.to_csv(OUT / "05a_cutoff_dilution_curve.csv", index=False, encoding="utf-8-sig")
g5b.to_csv(OUT / "05b_recluster_by_cutoff.csv", index=False, encoding="utf-8-sig")
P("## G5. 失眠切點的敏感度（健康組固定為 ISI < 9）\n")
P("(a) 以 ISI ≥ c 定義失眠症狀者時的合併效果量；右邊是 9 分到 c−1 分的人\n")
P(g5a.to_markdown(index=False))
P("")
P("(b) 換切點後重新分兩群\n")
P(g5b.to_markdown(index=False))
P("")

# G6. 傳統 beta 對 periodic beta；aperiodic 的形狀
FE = [
    "beta_pw_BRAIN",
    "Beta_rel_BRAIN",
    "logBeta_abs_BRAIN",
    "beta_pw_CENTRAL",
    "Beta_rel_CENTRAL",
    "exponent_BRAIN",
    "offset_BRAIN",
    "exponent_CENTRAL",
    "offset_CENTRAL",
]
rows = [
    dict(feature=f, contrast=f"{g_} vs HC", **gci(d[f], hc[f])) for f in FE for g_, d in GR.items()
]
g6a = pd.DataFrame(rows)


def boot_diff(trad, peri, grp):
    a = grp[[trad, peri]].dropna().to_numpy()
    h = hc[[trad, peri]].dropna().to_numpy()
    obs = hedges_g(a[:, 0], h[:, 0])[0] - hedges_g(a[:, 1], h[:, 1])[0]
    bs = np.empty(B)
    for i in range(B):
        x = a[rng.integers(0, len(a), len(a))]
        z = h[rng.integers(0, len(h), len(h))]
        bs[i] = hedges_g(x[:, 0], z[:, 0])[0] - hedges_g(x[:, 1], z[:, 1])[0]
    return dict(
        n=len(a),
        g_traditional=round(hedges_g(a[:, 0], h[:, 0])[0], 3),
        g_periodic=round(hedges_g(a[:, 1], h[:, 1])[0], 3),
        difference=round(obs, 3),
        ci_low=round(np.percentile(bs, 2.5), 3),
        ci_high=round(np.percentile(bs, 97.5), 3),
    )


rows = []
for trad, peri, reg in [
    ("Beta_rel_BRAIN", "beta_pw_BRAIN", "whole-scalp"),
    ("Beta_rel_CENTRAL", "beta_pw_CENTRAL", "central"),
    ("logBeta_abs_BRAIN", "beta_pw_BRAIN", "whole-scalp (absolute power)"),
]:
    for g_ in ["Lower", "Higher"]:
        rows.append(dict(region=reg, traditional=trad, group=g_, **boot_diff(trad, peri, GR[g_])))
g6b = pd.DataFrame(rows)

rows = []
for f in FE:
    d = df.dropna(subset=[f] + COV_COLS).copy()
    d["z"] = (d[f] - d[f].mean()) / d[f].std(ddof=1)
    d["s9"] = (d["ISI_T1"] >= 9).astype(float)
    d["s18"] = (d["ISI_T1"] >= THR).astype(float)
    aic = {
        k: smf.ols(f"z ~ {COV}{v}", d).fit().aic
        for k, v in [
            ("covariates", ""),
            ("linear", " + ISI_T1"),
            ("step at 9", " + s9"),
            ("step at 18", " + s18"),
            ("steps at 9 and 18", " + s9 + s18"),
        ]
    }
    mb = smf.ols(f"z ~ {COV} + s9 + s18", d).fit()
    y_, isi, Cm = d["z"].to_numpy(), d["ISI_T1"].to_numpy(float), cmat(d)
    c_hat, F_hat, sgn = scan_fast(y_, isi, Cm, list(range(5, 23)))
    nullF = np.array(
        [scan_fast(y_, rng.permutation(isi), Cm, list(range(5, 23)))[1] for _ in range(2000)]
    )
    rows.append(
        dict(
            feature=f,
            n=len(d),
            AIC_linear=round(aic["linear"] - aic["covariates"], 2),
            AIC_step9=round(aic["step at 9"] - aic["covariates"], 2),
            AIC_step18=round(aic["step at 18"] - aic["covariates"], 2),
            AIC_both=round(aic["steps at 9 and 18"] - aic["covariates"], 2),
            b_step9=round(mb.params["s9"], 3),
            p_step9=round(mb.pvalues["s9"], 4),
            b_step18=round(mb.params["s18"], 3),
            p_step18=round(mb.pvalues["s18"], 4),
            selected_threshold=c_hat,
            direction="+" if sgn > 0 else "−",
            perm_p=round((np.sum(nullF >= F_hat) + 1) / 2001, 4),
        )
    )
g6c = pd.DataFrame(rows)

rows = []
for trad, exps, reg in [
    ("Beta_rel_BRAIN", " + exponent_BRAIN + offset_BRAIN", "whole-scalp"),
    ("Beta_rel_CENTRAL", " + exponent_CENTRAL + offset_CENTRAL", "central"),
    ("beta_pw_BRAIN", " + exponent_BRAIN + offset_BRAIN", "whole-scalp"),
    ("beta_pw_CENTRAL", " + exponent_CENTRAL + offset_CENTRAL", "central"),
]:
    d = d_all.dropna(subset=[trad])
    sdh = hc[trad].std(ddof=1)
    for name, extra in [("age, sex, education", ""), ("+ aperiodic exponent and offset", exps)]:
        m = smf.ols(f"{trad} ~ {COV}{extra} + C(k2c)", d).fit()
        r = dict(outcome=trad, region=reg, covariates=name, n=int(m.nobs))
        for g_ in ["Lower", "Higher"]:
            t = f"C(k2c)[T.{g_}]"
            r[f"{g_}_b_in_HC_SD"] = round(m.params[t] / sdh, 3)
            r[f"{g_}_p"] = round(m.pvalues[t], 4)
        rows.append(r)
g6d = pd.DataFrame(rows)

BINS = [(9, 10), (11, 12), (13, 14), (15, 16), (17, 18), (19, 28)]
g6e = pd.DataFrame(
    [
        dict(
            feature=f,
            ISI_band=f"{lo}-{hi}" if hi < 28 else f">={lo}",
            **gci(ins.loc[(ins["ISI_T1"] >= lo) & (ins["ISI_T1"] <= hi), f], hc[f]),
        )
        for f in FE
        for lo, hi in BINS
    ]
)
for nm, t in [
    ("06a_groups.csv", g6a),
    ("06b_traditional_minus_periodic.csv", g6b),
    ("06c_shape_by_measure.csv", g6c),
    ("06d_adjust_for_aperiodic.csv", g6d),
    ("06e_bands_by_measure.csv", g6e),
]:
    t.to_csv(OUT / nm, index=False, encoding="utf-8-sig")
P("## G6. 傳統 beta、periodic beta、aperiodic 沿嚴重度的差別\n")
P("(a) 各指標在兩群（K = 2）與門檻上下的效果量\n")
P(g6a.to_markdown(index=False))
P("")
P(f"(b) 同一批人身上，傳統 beta 的 g 減去 periodic beta 的 g（bootstrap B = {B}）\n")
P(g6b.to_markdown(index=False))
P("")
P(
    "(c) 全樣本（含 HC）中，各指標的形狀。AIC 欄是相對於只有共變項的模型，越負越好；b_step9／b_step18 是兩個落差同時放進模型的係數（z 分數）；"
    "selected_threshold 是資料選出的單一落差位置，perm_p 已校正挑切點\n"
)
P(g6c.to_markdown(index=False))
P("")
P("(d) 把 aperiodic 的 exponent 與 offset 放進共變項後，兩群對 HC 的差距（HC 標準差）\n")
P(g6d.to_markdown(index=False))
P("")
P("數字全部由本腳本即時計算；尚未寫回 SSOT。")
(OUT / "07_summary.md").write_text("\n".join(MD), encoding="utf-8")
print(f"\n輸出 → {OUT}")
