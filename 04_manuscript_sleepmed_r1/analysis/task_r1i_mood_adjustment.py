"""比較加入 BDI-II、BAI 前後，各 EEG 指標的組間差異。

高嚴重度組另報 M2／M1 的保留比例及受試者層級 bootstrap 區間。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from _r1_common import AN, COV, COV_COLS, load_all

B = 2000
OUT = AN / "S_R1i_mood_adjustment"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261116)
MEAS = [
    ("beta_pw_BRAIN", "Periodic beta peak power, whole-scalp"),
    ("beta_pw_CENTRAL", "Periodic beta peak power, central"),
    ("Beta_rel_BRAIN", "Relative beta power (13–30 Hz), whole-scalp"),
    ("Beta_rel_CENTRAL", "Relative beta power, central"),
    ("logBeta_abs_BRAIN", "Absolute beta power (log), whole-scalp"),
    ("exponent_BRAIN", "Aperiodic exponent, whole-scalp"),
    ("exponent_CENTRAL", "Aperiodic exponent, central"),
    ("offset_BRAIN", "Aperiodic offset, whole-scalp"),
    ("offset_CENTRAL", "Aperiodic offset, central"),
]

df = load_all()
df["logBeta_abs_BRAIN"] = np.log10(df["Beta_abs_BRAIN"])
d0 = df.dropna(subset=COV_COLS + ["BDI_T1", "BAI_T1"]).copy()
d0["k2c"] = pd.Categorical(d0["k2"], categories=["HC", "Lower", "Higher"])


def coefs(d, f):
    sd = d.loc[d["k2"] == "HC", f].std(ddof=1)
    m1 = smf.ols(f"{f} ~ {COV} + C(k2c)", d).fit()
    m2 = smf.ols(f"{f} ~ {COV} + BDI_T1 + BAI_T1 + C(k2c)", d).fit()
    return m1, m2, sd


rows = []
for f, lab in MEAS:
    d = d0.dropna(subset=[f]).reset_index(drop=True)
    m1, m2, sd = coefs(d, f)
    for g in ["Higher", "Lower"]:
        t = f"C(k2c)[T.{g}]"
        c1, c2 = m1.conf_int().loc[t] / sd, m2.conf_int().loc[t] / sd
        r = dict(
            measure=lab,
            variable=f,
            group=f"{g} severity vs HC",
            n=int(m1.nobs),
            M1_diff_in_HC_SD=round(m1.params[t] / sd, 3),
            M1_ci=f"[{c1[0]:.2f}, {c1[1]:.2f}]",
            M1_p=round(m1.pvalues[t], 4),
            M2_diff_in_HC_SD=round(m2.params[t] / sd, 3),
            M2_ci=f"[{c2[0]:.2f}, {c2[1]:.2f}]",
            M2_p=round(m2.pvalues[t], 4),
        )
        if g == "Higher":
            ratio = m2.params[t] / m1.params[t]
            bs = []
            for _ in range(B):
                bd = d.iloc[rng.integers(0, len(d), len(d))]
                a, b_, _s = coefs(bd, f)
                if abs(a.params[t]) > 1e-9:
                    bs.append(b_.params[t] / a.params[t])
            r.update(
                retained=round(ratio, 2),
                retained_ci=f"[{np.percentile(bs, 2.5):.2f}, {np.percentile(bs, 97.5):.2f}]",
            )
        rows.append(r)
res = pd.DataFrame(rows)
res.to_csv(OUT / "01_mood_adjustment.csv", index=False, encoding="utf-8-sig")
md = [
    "# S_R1i — 控制憂鬱與焦慮後的組間差距\n",
    f"N = {len(d0)}（BDI-II、BAI 與共變項完整者）。M1 控制年齡、性別、教育；M2 再加 BDI-II、BAI。差距以健康組標準差為單位。",
    "retained＝M2 係數／M1 係數（只對高嚴重度群計算；bootstrap B = %d）。\n" % B,
    res.drop(columns="variable").to_markdown(index=False),
    "",
    "數字全部由本腳本即時計算；尚未寫回 SSOT。",
]
(OUT / "02_summary.md").write_text("\n".join(md), encoding="utf-8")
print("\n".join(md))
