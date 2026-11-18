"""
task_r1v_sleep_onset_check.py — Sleep Medicine R1：入睡困難與週期性 beta（探索性；用第二份問卷做一次事先講好的檢查）

背景：探索性分析看到，ISI 總分相近時，第 1 題（入睡困難）分數較高者的週期性 beta 較高
      （樣條總分下 b = 0.040，HC3 p = .025；七題校正後 q = .18）。第 1 題是看過七題的結果後才挑出來的，
      所以這個線索需要用另一個來源檢查一次。

分析計畫（2026-11-17 寫定，寫定時還沒有看過 PSQI 任何分項與腦波的關係；只核對過計分與分布）：
  指標    PSQI 第 2 分項（sleep latency，0–3 分；由「入睡所需時間」與「30 分鐘內無法入睡的頻率」兩題合成）。
          計分已核對：七個分項相加等於 PSQI 總分（268／268 人）。
  樣本    失眠症狀者（ISI ≥ 9）中有 PSQI 分項的人。
  主要檢定（只有這一個）
          週期性 beta ~ 年齡 + 性別 + 教育年數 + ISI 總分的自然三次樣條（3 個自由度）+ PSQI 第 2 分項（線性）
          看 PSQI 第 2 分項的係數；HC3 標準誤；雙尾 α = .05；預期方向為正。
  一併報告（不論主要檢定的結果）
          (a) 再加 BDI-II 與 BAI；(b) ISI 總分改用線性項；(c) 其餘六個 PSQI 分項各自放進同一個模型，
          用來看「入睡」是不是特別，而不是任何睡眠抱怨都有關（七個分項一起做 BH 校正）。
  解讀    主要檢定 p < .05 且方向為正：兩份問卷的入睡題都指向同一方向，可在補充資料與討論提一句，仍屬探索性。
          否則：這個線索沒有得到第二份問卷的支持，不寫進稿件。
  限制    PSQI 與 ISI 由同一批人在同一時間填寫，兩者的入睡題高度相關，這只是一致性的檢查，不是獨立樣本的重現。

輸出 outputs/補充資料/analysis/S_R1v_sleep_onset/
  01_primary.csv、02_all_components.csv、03_summary.md
執行：python task_r1v_sleep_onset_check.py
"""

from __future__ import annotations
import sqlite3
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from _dynamic_n import DB_PATH, ISI_ITEMS
from _r1_common import AN, COV, OUTCOME, load_all, fmt_p

OUT = AN / "S_R1v_sleep_onset"
OUT.mkdir(parents=True, exist_ok=True)
Y = OUTCOME
COMP = {
    f"PSQI_T1_C{i}": n
    for i, n in enumerate(
        [
            "subjective sleep quality",
            "sleep latency",
            "sleep duration",
            "habitual sleep efficiency",
            "sleep disturbances",
            "use of sleep medication",
            "daytime dysfunction",
        ],
        start=1,
    )
}
SPL = "cr(ISI_T1, df=3)"

df = load_all()
with sqlite3.connect(DB_PATH) as con:
    ps = pd.read_sql(
        "SELECT NEW_NUMBER AS subject, "
        + ", ".join(COMP)
        + " FROM isi_raw_data WHERE NEW_NUMBER IS NOT NULL",
        con,
    )
df = df.merge(ps, on="subject", how="left")
chk = df.dropna(subset=list(COMP) + ["PSQI_T1"])
assert (chk[list(COMP)].sum(axis=1) == chk["PSQI_T1"]).all(), "PSQI 分項相加不等於總分"
ins = df[df["ISI_T1"] >= 9].dropna(subset=["PSQI_T1_C2"]).copy()


def coef(formula, data, term, label):
    m = smf.ols(formula, data).fit()
    r = m.get_robustcov_results(cov_type="HC3", use_t=True)
    i = list(m.params.index).index(term)
    lo, hi = r.conf_int()[i]
    sd = data[Y].std(ddof=1)
    return dict(
        model=label,
        term=term,
        n=int(m.nobs),
        b=m.params[term],
        b_in_sd_of_beta=m.params[term] / sd,
        ci_low_hc3=lo,
        ci_high_hc3=hi,
        p_hc3=r.pvalues[i],
        p_ols=m.pvalues[term],
    )


m2 = ins.dropna(subset=["BDI_T1", "BAI_T1"])
prim = [
    coef(f"{Y} ~ {COV} + {SPL} + PSQI_T1_C2", ins, "PSQI_T1_C2", "PRIMARY: spline ISI total"),
    coef(
        f"{Y} ~ {COV} + BDI_T1 + BAI_T1 + {SPL} + PSQI_T1_C2", m2, "PSQI_T1_C2", "(a) + BDI-II, BAI"
    ),
    coef(f"{Y} ~ {COV} + ISI_T1 + PSQI_T1_C2", ins, "PSQI_T1_C2", "(b) linear ISI total"),
    coef(
        f"{Y} ~ {COV} + {SPL} + {ISI_ITEMS[0]}",
        ins,
        ISI_ITEMS[0],
        "reference: ISI item 1, same participants",
    ),
]
t1 = pd.DataFrame(prim)
t1.to_csv(OUT / "01_primary.csv", index=False, encoding="utf-8-sig")

rows = []
for c, name in COMP.items():
    d = ins.dropna(subset=[c])
    rows.append(
        dict(component=name, **coef(f"{Y} ~ {COV} + {SPL} + {c}", d, c, "spline ISI total"))
    )
t2 = pd.DataFrame(rows)
t2["q_bh_7"] = multipletests(t2["p_hc3"], method="fdr_bh")[1]
t2.to_csv(OUT / "02_all_components.csv", index=False, encoding="utf-8-sig")

rho = ins[[ISI_ITEMS[0], "PSQI_T1_C2"]].corr(method="spearman").iloc[0, 1]
L = [
    "# 入睡困難與週期性 beta：用 PSQI 入睡時間分項做的一次檢查（task_r1v；探索性）",
    "",
    f"樣本：失眠症狀者中有 PSQI 分項的 {len(ins)} 人。PSQI 入睡時間分項與 ISI 第 1 題的 Spearman ρ = {rho:.2f}。",
    "",
    "| 模型 | n | b | 以 beta 的標準差計 | HC3 95% CI | HC3 p |",
    "|---|---|---|---|---|---|",
]
for _, r in t1.iterrows():
    L.append(
        f"| {r['model']} | {r['n']} | {r['b']:+.4f} | {r['b_in_sd_of_beta']:+.2f} | [{r['ci_low_hc3']:+.4f}, {r['ci_high_hc3']:+.4f}] | {fmt_p(r['p_hc3'])} |"
    )
L += [
    "",
    "七個 PSQI 分項各自放進同一個模型（樣條總分下）：",
    "",
    "| 分項 | n | b | HC3 p | BH q |",
    "|---|---|---|---|---|",
]
for _, r in t2.iterrows():
    L.append(
        f"| {r['component']} | {r['n']} | {r['b']:+.4f} | {fmt_p(r['p_hc3'])} | {fmt_p(r['q_bh_7'])} |"
    )
(OUT / "03_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))
