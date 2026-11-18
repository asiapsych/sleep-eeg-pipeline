"""讀 MATLAB 的 manifests 與 ICLabel decisions，彙整前處理品管。

beta_pct_change 包含去電源雜訊、ICA 剔除與改參考的整體變化，不能只歸因於 ICA。
逐人表供內部核對，補充資料只用組層級結果。"""

from __future__ import annotations
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from _r1_common import AN, COV, COV_COLS, OUTCOME, hedges_g, load_all, mwu_p

Y = OUTCOME
THR = 18
from _dynamic_n import local_path

PIPE = Path(
    local_path("eeg_pipeline_root", "EEG_PIPELINE_ROOT")
)
MAN = PIPE / "outputs" / "manifests"
DEC = PIPE / "outputs" / "decisions"
OUT = AN / "S_R1d_preproc_qc"
OUT.mkdir(parents=True, exist_ok=True)
MD = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    MD.append(s)


def flat(x):
    """MATLAB jsonencode 會把單一元素寫成純量、多個寫成巢狀 list；一律攤平成 list。"""
    if x is None:
        return []
    if isinstance(x, (list, tuple)):
        out = []
        for v in x:
            out += flat(v)
        return out
    return [x]


def read_subject(sid):
    r = dict(subject=sid)
    fm = MAN / f"{sid}.json"
    if not fm.exists():
        return r
    m = json.loads(fm.read_text())
    s1, s2, s2b, s3 = (m.get(k, {}) or {} for k in ("step1", "step2", "step2b", "step3"))
    r.update(
        n_bad_ch=s1.get("n_bad"),
        n_ic=s2.get("n_components"),
        n_ic_removed=s2b.get("n_removed"),
        beta_pct_change=s2b.get("beta_pct_change"),
        overcorr_flag=s2b.get("overcorr_flag"),
        n_epochs_total=s3.get("n_total"),
        n_after_reject=s3.get("n_after_layerC"),
        reject_rate=s3.get("reject_rate"),
        n_final=s3.get("n_final"),
    )
    removed = [int(v) for v in flat(s2b.get("final_remove"))]
    fi = DEC / f"{sid}_iclabel.json"
    if fi.exists():
        ic = json.loads(fi.read_text())
        classes = [str(c) for c in flat(ic.get("classes"))]
        sc = np.array(ic.get("scores"), dtype=float)
        if sc.ndim == 1:
            sc = sc.reshape(1, -1)
        if "Muscle" in classes and sc.shape[1] == len(classes):
            im, ie = classes.index("Muscle"), classes.index("Eye")
            keep = np.array([i + 1 not in removed for i in range(sc.shape[0])])
            top = sc.argmax(axis=1)
            r.update(
                n_ret_muscle_dom=int(np.sum(keep & (top == im))),
                max_ret_muscle_p=float(sc[keep, im].max()) if keep.any() else np.nan,
                n_ret_brain_dom=int(np.sum(keep & (top == classes.index("Brain")))),
            )
            rem = np.array([i for i in removed if 1 <= i <= sc.shape[0]], dtype=int) - 1
            r.update(
                n_removed_muscle=int(np.sum(top[rem] == im)) if len(rem) else 0,
                n_removed_eye=int(np.sum(top[rem] == ie)) if len(rem) else 0,
            )
        reasons = " ".join(str(v) for v in flat(ic.get("flag_reason")))
        r["n_flag_muscle"] = len(re.findall(r"Muscle=", reasons))
        r["n_flag_eye"] = len(re.findall(r"Eye=", reasons))
    return r


df = load_all()
qc = pd.DataFrame([read_subject(s) for s in df["subject"]])
for c in qc.columns:
    if c != "subject":
        qc[c] = pd.to_numeric(
            qc[c].astype(float) if qc[c].dtype == bool else qc[c], errors="coerce"
        )
df = df.merge(qc, on="subject", how="left")
df["any_ret_muscle"] = (
    (df["n_ret_muscle_dom"] > 0).astype(float).where(df["n_ret_muscle_dom"].notna())
)
df["abs_beta_pct_change"] = df["beta_pct_change"].abs()
keep_cols = ["subject", "group", "k2", "ISI_T1", Y] + [c for c in qc.columns if c != "subject"]
df[keep_cols].to_csv(OUT / "_subject_level_internal.csv", index=False, encoding="utf-8-sig")

hc = df[df["group"] == "HC"]
ins = df[df["group"] != "HC"]
P("# S_R1d — 前處理品管的組層級彙整\n")
P(
    f"N = {len(df)}；有 manifest 的人數 = {int(df['n_ic'].notna().sum())}；"
    f"有 ICLabel 紀錄的人數 = {int(df['n_ret_muscle_dom'].notna().sum())}\n"
)
P(
    f"pipeline：{PIPE.name}；ICLabel 自動剔除門檻 eye >= 0.70、muscle >= 0.80（config_thresholds.m）\n"
)

# 1. 各群描述
METRICS = [
    "n_bad_ch",
    "n_ic",
    "n_ic_removed",
    "n_removed_eye",
    "n_removed_muscle",
    "n_ret_muscle_dom",
    "any_ret_muscle",
    "max_ret_muscle_p",
    "n_ret_brain_dom",
    "beta_pct_change",
    "abs_beta_pct_change",
    "overcorr_flag",
    "n_epochs_total",
    "reject_rate",
]
METRICS = [m for m in METRICS if m in df.columns and df[m].notna().any()]
GROUPS = {
    "HC": hc,
    "MEW": ins[ins["group"] == "MEW"],
    "MDI": ins[ins["group"] == "MDI"],
    "SCA": ins[ins["group"] == "SCA"],
    "K2 Lower": ins[ins["k2"] == "Lower"],
    "K2 Higher": ins[ins["k2"] == "Higher"],
    f"ISI >={THR}": ins[ins["ISI_T1"] >= THR],
}
rows = []
for m in METRICS:
    try:
        H, pk = stats.kruskal(*[GROUPS[g][m].dropna() for g in ["HC", "MEW", "MDI", "SCA"]])
    except ValueError:
        H, pk = np.nan, np.nan
    r = dict(
        metric=m,
        KW4_H=round(H, 2) if H == H else np.nan,
        KW4_p=round(pk, 4) if pk == pk else np.nan,
    )
    for g, d in GROUPS.items():
        r[f"{g} mean"] = round(d[m].mean(), 3)
        r[f"{g} sd"] = round(d[m].std(ddof=1), 3)
    for g in ["SCA", "K2 Higher", f"ISI >={THR}"]:
        gg = hedges_g(GROUPS[g][m], hc[m])[0]
        r[f"g {g} vs HC"] = round(gg, 3) if gg == gg else np.nan
        r[f"p {g} vs HC"] = round(mwu_p(GROUPS[g][m], hc[m]), 4)
    rows.append(r)
gq = pd.DataFrame(rows)
gq.to_csv(OUT / "01_group_qc.csv", index=False, encoding="utf-8-sig")
show = gq[
    [
        "metric",
        "HC mean",
        "MEW mean",
        "MDI mean",
        "SCA mean",
        "K2 Higher mean",
        "KW4_p",
        "g SCA vs HC",
        "p SCA vs HC",
        "g K2 Higher vs HC",
        "p K2 Higher vs HC",
    ]
]
P("## 1. 各群的前處理指標（平均；any_ret_muscle、overcorr_flag 為比例）\n")
P(show.to_markdown(index=False))
P("")

cor = []
for m in METRICS:
    dd = df.dropna(subset=[Y, m])
    if dd[m].nunique() > 1:
        r_, p_ = stats.spearmanr(dd[Y], dd[m])
        cor.append(dict(metric=m, n=len(dd), spearman_r_with_beta=round(r_, 3), p=round(p_, 4)))
cor = pd.DataFrame(cor)
P("### 各指標與 periodic beta 的相關（全樣本）\n")
P(cor.to_markdown(index=False))
P("")

# 2. 共變項敏感度
d = df.dropna(subset=[Y] + COV_COLS).copy()
d["k3"] = pd.Categorical(d["group"], categories=["HC", "MEW", "MDI", "SCA"])
d["k2c"] = pd.Categorical(d["k2"], categories=["HC", "Lower", "Higher"])
d["sev3"] = pd.Categorical(
    np.where(d["group"] == "HC", "HC", np.where(d["ISI_T1"] >= THR, "High", "Low")),
    categories=["HC", "Low", "High"],
)
TARGETS = [
    ("k3", "C(k3)[T.SCA]", "SCA vs HC"),
    ("k3", "C(k3)[T.MEW]", "MEW vs HC"),
    ("k2c", "C(k2c)[T.Higher]", "K2 Higher vs HC"),
    ("k2c", "C(k2c)[T.Lower]", "K2 Lower vs HC"),
    ("sev3", "C(sev3)[T.High]", f"ISI >={THR} vs HC"),
]
SETS = {"none (age, sex, education only)": []}
for m in [
    "n_bad_ch",
    "n_ic_removed",
    "n_removed_muscle",
    "n_ret_muscle_dom",
    "max_ret_muscle_p",
    "beta_pct_change",
    "reject_rate",
]:
    if m in METRICS and d[m].nunique() > 1:
        SETS[f"+ {m}"] = [m]
allq = [
    m
    for m in [
        "n_bad_ch",
        "n_ic_removed",
        "n_ret_muscle_dom",
        "max_ret_muscle_p",
        "beta_pct_change",
        "reject_rate",
    ]
    if m in METRICS and d[m].nunique() > 1
]
SETS["+ all of the above QC metrics"] = allq
sd_hc = d.loc[d["group"] == "HC", Y].std(ddof=1)


def fit_rows(data, label, extra_cols):
    out = []
    dd = data.dropna(subset=extra_cols)
    extra = "".join(f" + {c}" for c in extra_cols)
    for fac, term, cname in TARGETS:
        m = smf.ols(f"{Y} ~ {COV}{extra} + C({fac})", dd).fit()
        if term not in m.params.index:
            continue
        ci = m.conf_int().loc[term]
        out.append(
            dict(
                setting=label,
                contrast=cname,
                n=int(m.nobs),
                b=round(m.params[term], 4),
                ci_low=round(ci[0], 4),
                ci_high=round(ci[1], 4),
                b_in_HC_SD=round(m.params[term] / sd_hc, 3),
                p=round(m.pvalues[term], 4),
            )
        )
    return out


cv = pd.DataFrame(sum((fit_rows(d, k, v) for k, v in SETS.items()), []))
cv.to_csv(OUT / "02_qc_covariate.csv", index=False, encoding="utf-8-sig")
P("## 2. 主要對比加入品管指標當共變項\n")
P(
    cv[cv["contrast"].isin(["SCA vs HC", "K2 Higher vs HC", f"ISI >={THR} vs HC"])].to_markdown(
        index=False
    )
)
P("")

# 3. 高品質子樣本
SUBS = {"full sample": d}
if "n_ret_muscle_dom" in d:
    SUBS["no retained muscle-dominant component"] = d[d["n_ret_muscle_dom"] == 0]
if "max_ret_muscle_p" in d:
    SUBS["max retained muscle probability < 0.50"] = d[d["max_ret_muscle_p"] < 0.50]
    SUBS["max retained muscle probability < 0.20"] = d[d["max_ret_muscle_p"] < 0.20]
if "overcorr_flag" in d:
    SUBS["no over-correction flag (|beta change| <= 20%)"] = d[d["overcorr_flag"] == 0]
if "n_bad_ch" in d:
    SUBS["no interpolated channel"] = d[d["n_bad_ch"] == 0]
    SUBS["<= 2 interpolated channels"] = d[d["n_bad_ch"] <= 2]
if "reject_rate" in d:
    SUBS["epoch rejection rate <= 10%"] = d[d["reject_rate"] <= 0.10]
rs = []
for label, dd in SUBS.items():
    h = dd.loc[dd["group"] == "HC", Y]
    for cname, mask in [
        ("SCA vs HC", dd["group"] == "SCA"),
        ("MEW vs HC", dd["group"] == "MEW"),
        ("K2 Higher vs HC", dd["k2"] == "Higher"),
        ("K2 Lower vs HC", dd["k2"] == "Lower"),
        (f"ISI >={THR} vs HC", (dd["group"] != "HC") & (dd["ISI_T1"] >= THR)),
    ]:
        g, lo, hi, n1, n2 = hedges_g(dd.loc[mask, Y], h)
        rs.append(
            dict(
                subsample=label,
                contrast=cname,
                n_a=n1,
                n_HC=n2,
                g=round(g, 3) if g == g else np.nan,
                ci_low=round(lo, 3) if lo == lo else np.nan,
                ci_high=round(hi, 3) if hi == hi else np.nan,
                MWU_p=round(mwu_p(dd.loc[mask, Y], h), 4),
            )
        )
rs = pd.DataFrame(rs)
rs.to_csv(OUT / "03_restricted_samples.csv", index=False, encoding="utf-8-sig")
P("## 3. 限定在高品質子樣本後的效果量\n")
P(
    rs[
        rs["contrast"].isin(["SCA vs HC", "K2 Higher vs HC", "K2 Lower vs HC", "MEW vs HC"])
    ].to_markdown(index=False)
)
P("")
P("數字全部由本腳本即時計算。")

(OUT / "04_summary.md").write_text("\n".join(MD), encoding="utf-8")
print(f"\n輸出 → {OUT}")
