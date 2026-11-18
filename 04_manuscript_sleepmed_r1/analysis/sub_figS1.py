"""Figure S1：用資料即時計算各階段的樣本人數。

前處理前的 308 筆是流程紀錄表的數字；若重跑前處理，這個常數要一併核對。"""

import sqlite3

from _dynamic_n import DB_PATH, EXCLUDED, EXCLUDED_ELIGIBILITY, ISI_ITEMS, S17_RELPOWER_CSV
from _r1_common import load_all
from _r1_submission import (
    AXIS,
    INK2,
    OUT,
    SURFACE,
    pd,
    plt,
    save_fig,
)

N_PREPROCESSED = 308  # 前處理流程 v2.3 的紀錄數（CHANGELOG_v2.3.md；outputs/logs）
eeg = set(pd.read_csv(S17_RELPOWER_CSV)["subject"])
with sqlite3.connect(DB_PATH) as con:
    db = pd.read_sql(
        "SELECT NEW_NUMBER AS subject, "
        + ", ".join(ISI_ITEMS)
        + " FROM isi_raw_data WHERE NEW_NUMBER IS NOT NULL",
        con,
    )
isi_ok = set(db.dropna(subset=ISI_ITEMS)["subject"])
both = eeg & isi_ok
n_orig = len(both - set(EXCLUDED))
assert n_orig == len(both), "初投前另有排除時，流程圖要多一步"
n_final = len(both - set(EXCLUDED) - set(EXCLUDED_ELIGIBILITY))
df = load_all()
assert len(df) == n_final, (len(df), n_final)
n_hc = int((df["ISI_T1"] < 9).sum())
n_ins = int((df["ISI_T1"] >= 9).sum())
steps = [
    (
        f"EEG recordings preprocessed\nn = {N_PREPROCESSED}",
        f"Fewer than 60 artefact-free epochs\nn = {N_PREPROCESSED - len(eeg)}",
    ),
    (
        f"Usable resting EEG (60 epochs)\nn = {len(eeg)}",
        f"Incomplete ISI items\nn = {len(eeg) - len(both)}",
    ),
    (
        f"EEG and complete ISI\nn = {n_orig}",
        f"Met exclusion criteria on audit\nof recruitment records\nn = {n_orig - n_final}",
    ),
    (f"Analytic sample\nn = {n_final}", None),
]
pd.DataFrame(
    [dict(step=a.replace("\n", ": "), excluded=(b or "").replace("\n", ": ")) for a, b in steps]
).to_csv(OUT / "FigureS1_counts.csv", index=False, encoding="utf-8-sig")

fig, ax = plt.subplots(figsize=(8.6, 6.6))
ax.axis("off")
ax.set_xlim(0, 10)
ax.set_ylim(2.6, 12.4)
box = dict(boxstyle="round,pad=0.5", fc=SURFACE, ec=INK2, lw=1.1)
y = 11.5
for i, (main, exc) in enumerate(steps):
    ax.text(3.2, y, main, ha="center", va="center", fontsize=10, bbox=box)
    if exc:
        ax.annotate(
            "",
            xy=(3.2, y - 1.45),
            xytext=(3.2, y - 0.55),
            arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.1),
        )
        ax.annotate(
            "",
            xy=(5.6, y - 1.0),
            xytext=(3.2, y - 1.0),
            arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.1),
        )
        ax.text(
            5.75,
            y - 1.0,
            exc,
            ha="left",
            va="center",
            fontsize=9,
            color=INK2,
            bbox=dict(boxstyle="round,pad=0.4", fc=SURFACE, ec=AXIS, lw=0.9),
        )
    y -= 2.0
y += 2.0
for x, lab in [
    (1.5, f"Healthy controls (ISI < 9)\nn = {n_hc}"),
    (5.0, f"Insomnia symptoms (ISI ≥ 9)\nn = {n_ins}"),
]:
    ax.annotate(
        "",
        xy=(x, y - 1.45),
        xytext=(3.2, y - 0.55),
        arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.1),
    )
    ax.text(x, y - 2.0, lab, ha="center", va="center", fontsize=10, bbox=box)
save_fig(fig, "FigureS1")
print(pd.read_csv(OUT / "FigureS1_counts.csv").to_string(index=False))
