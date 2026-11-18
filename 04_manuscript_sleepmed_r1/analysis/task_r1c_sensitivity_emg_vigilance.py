"""用既有受試者層級資料檢查肌電、清醒度與 aperiodic 的影響。

ISI ≥ 18 來自資料選點，這裡只作描述性分組。
逐 epoch 的清醒度另見 task_r1s；本支不讀原始 EEG。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from _r1_common import AN, COV, COV_COLS, OUTCOME, fmt_p, hedges_g, load_all, mwu_p

Y = OUTCOME
THR = 18
OUT = AN / "S_R1c_sensitivity"
OUT.mkdir(parents=True, exist_ok=True)
MD = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    MD.append(s)


df = load_all()
fq = pd.read_csv(AN / "S10_fooof_full" / "02_fit_quality.csv")[
    ["subject", "r2_mean", "fit_error_mean"]
]
df = df.merge(fq, on="subject", how="left")
df["logGamma_abs"] = np.log10(df["Gamma_abs_BRAIN"])
df["logHighBeta_abs"] = np.log10(df["HighBeta_abs_BRAIN"])
hc = df[df["group"] == "HC"]
ins = df[df["group"] != "HC"]

CONTRASTS = {
    "MEW": ins[ins["group"] == "MEW"],
    "MDI": ins[ins["group"] == "MDI"],
    "SCA": ins[ins["group"] == "SCA"],
    "K2 Lower": ins[ins["k2"] == "Lower"],
    "K2 Higher": ins[ins["k2"] == "Higher"],
    f"ISI 9-{THR - 1}": ins[ins["ISI_T1"] < THR],
    f"ISI >={THR}": ins[ins["ISI_T1"] >= THR],
    "All insomnia": ins,
}


def contrast_table(features):
    rows = []
    for f in features:
        try:
            H, pk = stats.kruskal(
                *[df.loc[df["group"] == g, f].dropna() for g in ["HC", "MEW", "MDI", "SCA"]]
            )
        except ValueError:  # 全部受試者數值相同（例如固定的 epoch 數）
            H, pk = np.nan, np.nan
        for name, d in CONTRASTS.items():
            g, lo, hi, n1, n2 = hedges_g(d[f], hc[f])
            rows.append(
                dict(
                    feature=f,
                    contrast=f"{name} vs HC",
                    n_a=n1,
                    n_HC=n2,
                    mean_a=round(d[f].mean(), 4),
                    mean_HC=round(hc[f].mean(), 4),
                    g=round(g, 3),
                    ci_low=round(lo, 3),
                    ci_high=round(hi, 3),
                    MWU_p=round(mwu_p(d[f], hc[f]), 4),
                    KW4_H=round(H, 2),
                    KW4_p=round(pk, 4),
                )
            )
    return pd.DataFrame(rows)


P("# S_R1c — 高嚴重度敏感度、肌電、清醒度、aperiodic\n")
P(f"N = {len(df)}；HC = {len(hc)}；失眠 = {len(ins)}\n")

# 1. 高嚴重度敏感度
rows = []
for f in ["beta_pw_BRAIN", "beta_pw_CENTRAL", "Beta_rel_BRAIN"]:
    for label, d in [
        ("ISI >= 15", ins[ins["ISI_T1"] >= 15]),
        ("ISI 9-14", ins[ins["ISI_T1"] < 15]),
        ("ISI >= 9 (all insomnia)", ins),
    ]:
        g, lo, hi, n1, n2 = hedges_g(d[f], hc[f])
        vc = d["group"].value_counts().to_dict()
        rows.append(
            dict(
                feature=f,
                subset=label,
                n=n1,
                n_HC=n2,
                g_vs_HC=round(g, 3),
                ci_low=round(lo, 3),
                ci_high=round(hi, 3),
                MWU_p=round(mwu_p(d[f], hc[f]), 4),
                n_MEW=vc.get("MEW", 0),
                n_MDI=vc.get("MDI", 0),
                n_SCA=vc.get("SCA", 0),
            )
        )
hs = pd.DataFrame(rows)
hs.to_csv(OUT / "01_high_severity.csv", index=False, encoding="utf-8-sig")
P("## 1. 高嚴重度敏感度（ISI >= 15 為慣用的臨床失眠切點）\n")
P(hs.to_markdown(index=False))
P("")

# 2. 肌電指標
EMG = ["Gamma_rel_BRAIN", "logGamma_abs", "HighBeta_rel_BRAIN", "logHighBeta_abs"]
emg = contrast_table(EMG)
cor = []
for f in EMG:
    for label, d in [("full sample", df), ("insomnia only", ins), ("HC only", hc)]:
        dd = d.dropna(subset=[Y, f])
        r, p = stats.spearmanr(dd[Y], dd[f])
        cor.append(
            dict(
                feature=f,
                contrast=f"Spearman with periodic beta ({label})",
                n_a=len(dd),
                g=round(r, 3),
                MWU_p=round(p, 4),
            )
        )
emg = pd.concat([emg, pd.DataFrame(cor)], ignore_index=True)
emg.to_csv(OUT / "02_emg_indices.csv", index=False, encoding="utf-8-sig")
P("## 2. 高頻功率（肌電的代理指標）：Gamma = 30–45 Hz，HighBeta = 30–40 Hz\n")
P("相關列的 g 欄位是 Spearman r，MWU_p 欄位是相關的 p。\n")
P(emg.to_markdown(index=False))
P("")

# 3. 肌電共變項敏感度
d = df.dropna(subset=[Y] + COV_COLS).copy()
d["k3"] = pd.Categorical(d["group"], categories=["HC", "MEW", "MDI", "SCA"])
d["k2c"] = pd.Categorical(d["k2"], categories=["HC", "Lower", "Higher"])
d["sev3"] = np.where(d["group"] == "HC", "HC", np.where(d["ISI_T1"] >= THR, "High", "Low"))
d["sev3"] = pd.Categorical(d["sev3"], categories=["HC", "Low", "High"])
SETS = {
    "none (age, sex, education only)": "",
    "+ 30–45 Hz relative power": " + Gamma_rel_BRAIN",
    "+ log10 30–45 Hz absolute power": " + logGamma_abs",
    "+ 30–40 Hz relative power": " + HighBeta_rel_BRAIN",
    "+ aperiodic exponent (whole-brain)": " + exponent_BRAIN",
    "+ log10 30–45 Hz absolute power + exponent": " + logGamma_abs + exponent_BRAIN",
}
TARGETS = [
    ("k3", "C(k3)[T.SCA]", "SCA vs HC"),
    ("k3", "C(k3)[T.MEW]", "MEW vs HC"),
    ("k3", "C(k3)[T.MDI]", "MDI vs HC"),
    ("k2c", "C(k2c)[T.Higher]", "K2 Higher vs HC"),
    ("k2c", "C(k2c)[T.Lower]", "K2 Lower vs HC"),
    ("sev3", "C(sev3)[T.High]", f"ISI >={THR} vs HC"),
    ("sev3", "C(sev3)[T.Low]", f"ISI 9-{THR - 1} vs HC"),
]
sd_hc = d.loc[d["group"] == "HC", Y].std(ddof=1)


def cov_table(sets):
    rows = []
    for sname, extra in sets.items():
        for fac, term, label in TARGETS:
            m = smf.ols(f"{Y} ~ {COV}{extra} + C({fac})", d).fit()
            ci = m.conf_int().loc[term]
            rows.append(
                dict(
                    added_covariate=sname,
                    contrast=label,
                    n=int(m.nobs),
                    b=round(m.params[term], 4),
                    ci_low=round(ci[0], 4),
                    ci_high=round(ci[1], 4),
                    b_in_HC_SD=round(m.params[term] / sd_hc, 3),
                    p=round(m.pvalues[term], 4),
                )
            )
    return pd.DataFrame(rows)


ec = cov_table(SETS)
ec.to_csv(OUT / "03_emg_covariate.csv", index=False, encoding="utf-8-sig")
P(
    "## 3. 加入高頻功率當共變項後的主要對比（b 為 periodic beta 的組間差；b_in_HC_SD＝以 HC 標準差為單位）\n"
)
P(
    ec[
        ec["contrast"].isin(["SCA vs HC", "K2 Higher vs HC", f"ISI >={THR} vs HC", "MEW vs HC"])
    ].to_markdown(index=False)
)
P("")

# 4. 腦區分布
REG = [
    f"beta_pw_{r}"
    for r in ["BRAIN", "FRONTAL", "CENTRAL", "PARIETAL", "TEMPORAL", "OCCIPITAL", "MIDLINE"]
]
reg = contrast_table(REG)
reg.to_csv(OUT / "04_regional.csv", index=False, encoding="utf-8-sig")
P("## 4. periodic beta 的腦區分布（MIDLINE = Fz、Cz、Pz）\n")
P(
    reg[
        reg["contrast"].isin(
            ["SCA vs HC", "K2 Higher vs HC", f"ISI >={THR} vs HC", "MEW vs HC", "K2 Lower vs HC"]
        )
    ][["feature", "contrast", "n_a", "n_HC", "g", "ci_low", "ci_high", "MWU_p"]].to_markdown(
        index=False
    )
)
P("")

# 5. 清醒度
VIG = [
    "Theta_rel_BRAIN",
    "Alpha_rel_BRAIN",
    "Delta_rel_BRAIN",
    "TAR_BRAIN",
    "theta_pw_BRAIN",
    "alpha_pw_BRAIN",
    "alpha_cf_BRAIN",
]
vig = contrast_table(VIG)
VSETS = {
    "none (age, sex, education only)": "",
    "+ theta/alpha ratio": " + TAR_BRAIN",
    "+ theta relative power": " + Theta_rel_BRAIN",
    "+ alpha peak frequency": " + alpha_cf_BRAIN",
    "+ theta/alpha ratio + alpha peak frequency": " + TAR_BRAIN + alpha_cf_BRAIN",
}
vc = cov_table(VSETS)
vig.to_csv(OUT / "05_vigilance.csv", index=False, encoding="utf-8-sig")
vc.to_csv(OUT / "05b_vigilance_covariate.csv", index=False, encoding="utf-8-sig")
P(
    "## 5. 清醒度指標（嗜睡的典型變化：theta 上升、alpha 下降、theta/alpha 上升、alpha 峰頻率變慢）\n"
)
P(
    vig[
        vig["contrast"].isin(
            ["SCA vs HC", "K2 Higher vs HC", f"ISI >={THR} vs HC", "MEW vs HC", "MDI vs HC"]
        )
    ][
        [
            "feature",
            "contrast",
            "n_a",
            "n_HC",
            "mean_a",
            "mean_HC",
            "g",
            "ci_low",
            "ci_high",
            "MWU_p",
            "KW4_p",
        ]
    ].to_markdown(index=False)
)
P("")
P("### 加入清醒度指標當共變項後的主要對比\n")
P(
    vc[vc["contrast"].isin(["SCA vs HC", "K2 Higher vs HC", f"ISI >={THR} vs HC"])].to_markdown(
        index=False
    )
)
P("")

# 6. aperiodic
AP = [
    "exponent_BRAIN",
    "offset_BRAIN",
    "exponent_CENTRAL",
    "offset_CENTRAL",
    "exponent_FRONTAL",
    "offset_FRONTAL",
    "exponent_PARIETAL",
    "offset_PARIETAL",
    "exponent_TEMPORAL",
    "offset_TEMPORAL",
    "exponent_OCCIPITAL",
    "offset_OCCIPITAL",
]
ap = contrast_table(AP)
ap.to_csv(OUT / "06_aperiodic.csv", index=False, encoding="utf-8-sig")
P("## 6. aperiodic exponent 與 offset\n")
P(
    ap[ap["contrast"].isin(["SCA vs HC", "K2 Higher vs HC", f"ISI >={THR} vs HC"])][
        ["feature", "contrast", "n_a", "n_HC", "g", "ci_low", "ci_high", "MWU_p", "KW4_p"]
    ].to_markdown(index=False)
)
P("")

# 7. 資料量與擬合品質
DQ = ["n_epochs", "r2_mean", "fit_error_mean"]
dq = contrast_table(DQ)
desc = df.groupby("group")[DQ].agg(["mean", "std", "min", "max"]).round(3)
dq.to_csv(OUT / "07_data_quality.csv", index=False, encoding="utf-8-sig")
P("## 7. 保留的 epoch 數與頻譜擬合品質\n")
P(desc.to_markdown())
P("")
P(
    dq[dq["contrast"].isin(["SCA vs HC", "MEW vs HC", "MDI vs HC", "K2 Higher vs HC"])][
        ["feature", "contrast", "n_a", "n_HC", "mean_a", "mean_HC", "g", "MWU_p", "KW4_p"]
    ].to_markdown(index=False)
)
P("")
dd = (
    d.merge(df[["subject", "r2_mean"]], on="subject", how="left")
    if "r2_mean" not in d.columns
    else d
)
m = smf.ols(f"{Y} ~ {COV} + n_epochs + r2_mean + C(k3)", dd).fit()
P(
    f"主要對比再加入 epoch 數與擬合 R² 當共變項：SCA vs HC b = {m.params['C(k3)[T.SCA]']:.4f}，"
    f"p = {fmt_p(m.pvalues['C(k3)[T.SCA]'])}（n = {int(m.nobs)}）\n"
)
P("數字全部由本腳本即時計算；CSV 為機讀來源。尚未寫回 SSOT。")

(OUT / "08_summary.md").write_text("\n".join(MD), encoding="utf-8")
print(f"\n輸出 → {OUT}")
