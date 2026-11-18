"""送件圖表共用的樣式、數字格式與輸出工具。

圖存 PNG（300 dpi）及 PDF，表存 Markdown 及 CSV，都放在 S_R1_submission。
這裡設定的 rcParams 會套用到所有 sub_*.py。"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

from _r1_common import AN

OUT = AN / "S_R1_submission"
OUT.mkdir(parents=True, exist_ok=True)

SURFACE, INK, INK2, MUTED, GRID, AXIS = (
    "#ffffff",
    "#0b0b0b",
    "#52514e",
    "#898781",
    "#e1e0d9",
    "#c3c2b7",
)
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
DIV = LinearSegmentedColormap.from_list(
    "div", ["#1c5cab", "#86b6ef", "#f0efec", "#f3a583", "#c4481a"]
)
plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "Helvetica",
        "font.size": 10,
        "text.color": INK,
        "axes.labelcolor": INK2,
        "axes.edgecolor": AXIS,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "axes.axisbelow": True,
        "legend.frameon": False,
        "legend.fontsize": 9,
    }
)


def save_fig(fig, name):
    for ext, kw in (("png", dict(dpi=300)), ("pdf", {})):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)
    print("saved", OUT / f"{name}.png")


def panel_label(ax, s, x=-0.02, y=1.07):
    ax.text(
        x, y, s, transform=ax.transAxes, fontsize=13, fontweight="bold", ha="right", va="bottom"
    )


# 數字格式
def nl(x, d=2):
    """去前導 0（相關、p、比例這類上限為 1 的數）。"""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    s = f"{x:.{d}f}"
    return s.replace("-0.", "−.").replace("0.", ".", 1) if abs(x) < 1 else s.replace("-", "−")


def fp(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "–"
    return "< .001" if p < 0.001 else nl(p, 3)


def sg(x, d=2):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    return f"{0:.{d}f}" if round(x, d) == 0 else f"{x:+.{d}f}".replace("-", "−")


def uz(x, d=2):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    return f"{0:.{d}f}" if round(x, d) == 0 else f"{x:.{d}f}".replace("-", "−")


def ci(lo, hi, d=2):
    return f"[{sg(lo, d)}, {sg(hi, d)}]"


def gci(g, lo, hi, d=2):
    return f"{sg(g, d)} {ci(lo, hi, d)}"


def ci_txt(s, d=2):
    """分析輸出裡已經寫成 "[a, b]" 的信賴區間，重排成統一格式。"""
    lo, hi = (float(x) for x in str(s).strip("[] ").split(","))
    return ci(lo, hi, d)


def msd(s, d=2):
    s = pd.to_numeric(pd.Series(s), errors="coerce").dropna()
    return f"{s.mean():.{d}f} ({s.std(ddof=1):.{d}f})"


def kw3(a, b, c):
    a, b, c = (pd.to_numeric(pd.Series(v), errors="coerce").dropna() for v in (a, b, c))
    H, p = stats.kruskal(a, b, c)
    n = len(a) + len(b) + len(c)
    return H, p, max((H - 2) / (n - 3), 0)


def minus(s):
    """已經格式化好的字串：把連字號負號換成減號（只動數字前的）。"""
    import re

    return re.sub(r"(?<![\w.])-(?=\d|\.)", "−", str(s))


# 表
def md_table(headers, rows):
    esc = lambda c: str(c).replace("|", "\\|").replace("\n", " ")
    out = ["| " + " | ".join(esc(h) for h in headers) + " |", "|" + "---|" * len(headers)]
    out += ["| " + " | ".join(esc(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def write_table(name, label, title, parts, note=""):
    """parts：[(小標或 None, 欄名, 列), …]。寫出 <name>.md 與每一部分的 CSV。"""
    L = [f"**{label}.** {title}", ""]
    for i, (sub, headers, rows) in enumerate(parts):
        if sub:
            L += [f"*{sub}*", ""]
        L += [md_table(headers, rows), ""]
        suffix = "" if len(parts) == 1 else "_" + "abcdefghijklmnop"[i]
        pd.DataFrame(rows, columns=headers).to_csv(
            OUT / f"{name}{suffix}.csv", index=False, encoding="utf-8-sig"
        )
    if note:
        L += [f"*Note.* {note}", ""]
    (OUT / f"{name}.md").write_text("\n".join(L), encoding="utf-8")
    print("saved", OUT / f"{name}.md")


def csv_part(path, sub=None, drop_section_rows=False):
    """已經排好版的報告表（S_R1_tables/R*.csv）直接當成一個部分。"""
    t = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    rows = t.values.tolist()
    if drop_section_rows:
        rows = [r for r in rows if any(c != "" for c in r[1:])]
    return (sub, _headers(t), rows)


def _headers(t):
    """pandas 會把重複的欄名加上 .1、.2，這裡拿掉；上標記號 ^{a} 改成一般寫法。"""
    import re

    return [re.sub(r"\.\d+$", "", str(c)) for c in t.columns]


def csv_sections(path):
    """報告表裡「整列只有第一格有字」的列是小標。回傳 (欄名, {小標: 列})；小標之前的列放在 None。"""
    t = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    out, key = {None: []}, None
    for r in t.values.tolist():
        if all(c == "" for c in r[1:]):
            key = r[0]
            out[key] = []
        else:
            out[key].append(r)
    return _headers(t), out
