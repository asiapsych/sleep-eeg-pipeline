"""Word 三線表：Times New Roman、標題、表註及跨頁處理。

section／subheader 是表內的小標記號，寫檔時才轉成儲存格格式。"""

import re
from pathlib import Path
from typing import Optional, Sequence

from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


FONT_NAME = "Times New Roman"
FONT_SIZE = 10  # pt
TOP_BORDER_SZ = 8  # 1.0 pt × 8
MID_BORDER_SZ = 4  # 0.5 pt × 8
BOT_BORDER_SZ = 8  # 1.0 pt × 8

# 用標記區分跨欄小標與一般資料列；小標的格式比照主文 Table 4。
SECTION_PREFIX = "‖SEC‖"


def section(label: str):
    """Return a row that renders as a bold-italic section header spanning all columns."""
    return [SECTION_PREFIX + str(label)]


def _is_section(row) -> bool:
    return len(row) >= 1 and str(row[0]).startswith(SECTION_PREFIX)


# 同一張表各部分的欄名只加粗，不合併儲存格。
SUBHEADER_PREFIX = "‖SUB‖"


def subheader(cells):
    """Return a row whose cells render bold (per-panel column labels)."""
    cells = [str(c) for c in cells]
    cells[0] = SUBHEADER_PREFIX + cells[0]
    return cells


def _is_subheader(row) -> bool:
    return len(row) >= 1 and str(row[0]).startswith(SUBHEADER_PREFIX)


def _set_cell_border(cell, side: str, sz: int, color: str = "000000", val: str = "single"):
    """Set border for a single side of a cell."""
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_borders = tc_pr.find(qn("w:tcBorders"))
    if tc_borders is None:
        tc_borders = OxmlElement("w:tcBorders")
        tc_pr.append(tc_borders)
    edge = tc_borders.find(qn(f"w:{side}"))
    if edge is None:
        edge = OxmlElement(f"w:{side}")
        tc_borders.append(edge)
    edge.set(qn("w:val"), val)
    edge.set(qn("w:sz"), str(sz))
    edge.set(qn("w:color"), color)


def _clear_cell_borders(cell):
    """Remove all borders from cell (default = no borders)."""
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_borders = OxmlElement("w:tcBorders")
    for side in ("top", "left", "bottom", "right", "insideH", "insideV"):
        edge = OxmlElement(f"w:{side}")
        edge.set(qn("w:val"), "none")
        tc_borders.append(edge)
    existing = tc_pr.find(qn("w:tcBorders"))
    if existing is not None:
        tc_pr.remove(existing)
    tc_pr.append(tc_borders)


def _set_run_font(run, *, bold=False, italic=False):
    run.font.name = FONT_NAME
    run.font.size = Pt(FONT_SIZE)
    run.font.color.rgb = RGBColor(0, 0, 0)
    run.font.bold = bold
    run.font.italic = italic
    rPr = run._r.get_or_add_rPr()
    rFonts = rPr.find(qn("w:rFonts"))
    if rFonts is None:
        rFonts = OxmlElement("w:rFonts")
        rPr.append(rFonts)
    rFonts.set(qn("w:ascii"), FONT_NAME)
    rFonts.set(qn("w:hAnsi"), FONT_NAME)
    rFonts.set(qn("w:cs"), FONT_NAME)
    rFonts.set(qn("w:eastAsia"), FONT_NAME)


def _tighten(p):
    """1.15 line spacing, no space before/after — keeps table rows compact but readable."""
    pf = p.paragraph_format
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)
    pf.line_spacing = 1.15


def _add_paragraph(
    doc,
    text: str = "",
    *,
    bold=False,
    italic=False,
    align=WD_ALIGN_PARAGRAPH.LEFT,
    space_before=0,
    space_after=0,
):
    p = doc.add_paragraph()
    p.alignment = align
    p.paragraph_format.space_before = Pt(space_before)
    p.paragraph_format.space_after = Pt(space_after)
    if text:
        run = p.add_run(text)
        _set_run_font(run, bold=bold, italic=italic)
    return p


def _add_parsed_run(p, text: str):
    """Parse bold, italic, and ^{superscript} runs."""
    parts = re.split(r"(\^\{.*?\}|\*\*[^*\s].*?\*\*|\*[^*\s].*?\*)", text)
    for part in parts:
        if part.startswith("^{") and part.endswith("}"):
            r = p.add_run(part[2:-1])
            _set_run_font(r)
            r.font.superscript = True
        elif part.startswith("**") and part.endswith("**"):
            r = p.add_run(part[2:-2])
            _set_run_font(r, bold=True)
        elif part.startswith("*") and part.endswith("*"):
            r = p.add_run(part[1:-1])
            _set_run_font(r, italic=True)
        elif part:
            r = p.add_run(part)
            _set_run_font(r)


def write_apa_table(
    out_path: Path,
    table_number: str,  # e.g., "S20"
    title: str,  # italic 標題
    headers: Sequence[str],  # 表頭
    rows: Sequence[Sequence[str]],  # 表內容
    note: Optional[str] = None,
    col_widths: Optional[Sequence[float]] = None,  # 各欄寬 (inch)
    landscape: bool = False,
    row_group_size: Optional[int] = None,
):
    """
    Generate a docx file with one APA 7 style three-line table.

    Format:
      Table S20                  ← bold
      <Title in italic>          ← italic
      ┌──────────────────────┐
      │ Header1   Header2     │  ← top border 1.5pt
      ├──────────────────────┤
      │ row1col1  row1col2    │  ← middle border 1pt
      │ row2col1  row2col2    │
      └──────────────────────┘  ← bottom border 1.5pt
      Note. ...                  ← italic "Note." + text
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    doc = Document()
    if landscape:
        sec = doc.sections[0]
        sec.orientation = WD_ORIENT.LANDSCAPE
        sec.page_width, sec.page_height = sec.page_height, sec.page_width
        sec.left_margin = sec.right_margin = Inches(0.55)
        sec.top_margin = sec.bottom_margin = Inches(0.55)

    cap_p = doc.add_paragraph()
    cap_p.paragraph_format.space_after = Pt(0)
    r = cap_p.add_run(f"Supplementary Table {table_number}")
    _set_run_font(r, bold=True)

    title_p = doc.add_paragraph()
    title_p.paragraph_format.space_after = Pt(6)
    rt = title_p.add_run(title.rstrip().rstrip("."))
    _set_run_font(rt, italic=True)

    n_rows = len(rows) + 1
    n_cols = len(headers)
    table = doc.add_table(rows=n_rows, cols=n_cols)
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    # 同一筆資料不要拆頁；跨頁時重複欄名。
    for row in table.rows:
        tr_pr = row._tr.get_or_add_trPr()
        cant_split = OxmlElement("w:cantSplit")
        tr_pr.append(cant_split)
    header_pr = table.rows[0]._tr.get_or_add_trPr()
    repeat = OxmlElement("w:tblHeader")
    repeat.set(qn("w:val"), "true")
    header_pr.append(repeat)

    # APA 表格：數字置中，文字靠左。
    def _numlike(c):
        c = str(c).strip()
        if not c or c == "—":
            return None
        if re.search(r"[一-鿿]", c):
            return False
        if not re.search(r"\d", c):
            return False
        if len(re.findall(r"[A-Za-z]", c)) > 2:  # 多字母＝標籤/清單，視為文字
            return False
        return len(c) <= 18

    data_only = [r for r in rows if not _is_section(r) and not _is_subheader(r)]
    col_align = []
    for ci in range(n_cols):
        flags = [_numlike(row[ci]) for row in data_only if ci < len(row)]
        flags = [f for f in flags if f is not None]
        is_num = bool(flags) and (sum(flags) / len(flags) >= 0.6)
        col_align.append(WD_ALIGN_PARAGRAPH.CENTER if is_num else WD_ALIGN_PARAGRAPH.LEFT)

    for ci, h in enumerate(headers):
        cell = table.rows[0].cells[ci]
        _clear_cell_borders(cell)
        _set_cell_border(cell, "top", TOP_BORDER_SZ)
        _set_cell_border(cell, "bottom", MID_BORDER_SZ)
        for p in cell.paragraphs:
            p.text = ""
        p = cell.paragraphs[0]
        p.alignment = col_align[ci]
        _tighten(p)
        _add_parsed_run(p, str(h))  # APA: header not bold; statistical symbols italic

    section_rows = set()
    for ri, row in enumerate(rows, start=1):
        is_last = ri == n_rows - 1
        if _is_section(row):
            section_rows.add(ri)
            merged = table.rows[ri].cells[0]
            for ci in range(1, n_cols):
                merged = merged.merge(table.rows[ri].cells[ci])
            _clear_cell_borders(merged)
            if is_last:
                _set_cell_border(merged, "bottom", BOT_BORDER_SZ)
            for p in merged.paragraphs:
                p.text = ""
            p = merged.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            _tighten(p)
            r = p.add_run(str(row[0])[len(SECTION_PREFIX) :])
            _set_run_font(r, bold=True, italic=True)
            continue
        if _is_subheader(row):
            for ci, val in enumerate(row):
                cell = table.rows[ri].cells[ci]
                _clear_cell_borders(cell)
                if is_last:
                    _set_cell_border(cell, "bottom", BOT_BORDER_SZ)
                for p in cell.paragraphs:
                    p.text = ""
                p = cell.paragraphs[0]
                p.alignment = col_align[ci] if ci < len(col_align) else WD_ALIGN_PARAGRAPH.LEFT
                _tighten(p)
                txt = str(val)[len(SUBHEADER_PREFIX) :] if ci == 0 else str(val)
                r = p.add_run(txt)
                _set_run_font(r, bold=True)
            continue
        for ci, val in enumerate(row):
            cell = table.rows[ri].cells[ci]
            _clear_cell_borders(cell)
            if is_last:
                _set_cell_border(cell, "bottom", BOT_BORDER_SZ)
            for p in cell.paragraphs:
                p.text = ""
            p = cell.paragraphs[0]
            p.alignment = col_align[ci] if ci < len(col_align) else WD_ALIGN_PARAGRAPH.LEFT
            _tighten(p)
            _add_parsed_run(p, str(val))

    # 合併儲存格的小標列不套用個別欄寬。
    # 2026-08：未提供時以內容長度自動配寬（避免整表無明確欄寬）。
    if col_widths is None:
        usable = 9.0 if landscape else 6.5
        weights = []
        for ci in range(n_cols):
            lens = [len(str(headers[ci]))] if ci < len(headers) else [8]
            for row in rows:
                if _is_section(row):  # section row（合併列）不計權重
                    continue
                if ci < len(row):
                    lens.append(len(str(row[ci])))
            weights.append(min(max(max(lens), 5), 45))
        tot = sum(weights)
        col_widths = [usable * w / tot for w in weights]
    # 固定 table layout，避免 Word／LibreOffice 自動改掉欄寬。
    if col_widths is not None:
        table.autofit = False
        tblPr = table._tbl.tblPr
        layout = tblPr.find(qn("w:tblLayout"))
        if layout is None:
            layout = OxmlElement("w:tblLayout")
            tblPr.append(layout)
        layout.set(qn("w:type"), "fixed")
        total_dxa = int(sum(col_widths) * 1440)
        tblW = tblPr.find(qn("w:tblW"))
        if tblW is None:
            tblW = OxmlElement("w:tblW")
            tblPr.append(tblW)
        tblW.set(qn("w:type"), "dxa")
        tblW.set(qn("w:w"), str(total_dxa))
        grid = table._tbl.find(qn("w:tblGrid"))
        if grid is not None:
            for i, gc in enumerate(grid.findall(qn("w:gridCol"))):
                if i < len(col_widths):
                    gc.set(qn("w:w"), str(int(col_widths[i] * 1440)))
        for ci, w in enumerate(col_widths):
            for ri, row in enumerate(table.rows):
                if ri in section_rows:
                    continue
                row.cells[ci].width = Inches(w)

    # 同一指標的 Unadjusted／M1／M2 放在同一頁，方便對照。
    if row_group_size and row_group_size > 1:
        for group_start in range(1, n_rows, row_group_size):
            group_end = min(group_start + row_group_size, n_rows)
            for ri in range(group_start, group_end - 1):
                for cell in table.rows[ri].cells:
                    for p in cell.paragraphs:
                        p.paragraph_format.keep_with_next = True

    if note:
        note_p = doc.add_paragraph()
        note_p.paragraph_format.space_before = Pt(4)
        nr = note_p.add_run("Note. ")
        _set_run_font(nr, italic=True)
        _add_parsed_run(note_p, note)

    doc.save(str(out_path))


def write_csv_and_md(
    out_dir: Path,
    table_number: str,
    title: str,
    headers: Sequence[str],
    rows: Sequence[Sequence[str]],
    note: Optional[str] = None,
):
    """Also write csv + md alongside the docx for easy review."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n_cols = len(headers)

    def _csv_row(r):
        if _is_section(r):
            return [str(r[0])[len(SECTION_PREFIX) :]] + [""] * (n_cols - 1)
        if _is_subheader(r):
            return [str(r[0])[len(SUBHEADER_PREFIX) :]] + [str(c) for c in r[1:]]
        return [re.sub(r"\^\{(.*?)\}", r"\1", str(c)) for c in r]

    import csv as csvmod

    with (out_dir / f"{table_number}.csv").open("w", encoding="utf-8", newline="") as f:
        w = csvmod.writer(f)
        w.writerow(headers)
        for r in rows:
            w.writerow(_csv_row(r))
    lines = [
        f"# Supplementary Table {table_number}",
        "",
        f"*{title}*",
        "",
        "| " + " | ".join(headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for r in rows:
        if _is_section(r):
            label = str(r[0])[len(SECTION_PREFIX) :]
            cells = [f"**{label}**"] + [""] * (n_cols - 1)
            lines.append("| " + " | ".join(cells) + " |")
        elif _is_subheader(r):
            cells = [f"**{str(r[0])[len(SUBHEADER_PREFIX) :]}**"] + [f"**{c}**" for c in r[1:]]
            lines.append("| " + " | ".join(cells) + " |")
        else:
            cells = [re.sub(r"\^\{(.*?)\}", r"<sup>\1</sup>", str(c)) for c in r]
            lines.append("| " + " | ".join(cells) + " |")
    if note:
        lines.append("")
        lines.append(f"*Note.* {note}")
    (out_dir / f"{table_number}.md").write_text("\n".join(lines), encoding="utf-8")


def write_all(
    out_dir,
    table_number,
    title,
    headers,
    rows,
    note=None,
    col_widths=None,
    landscape=False,
    row_group_size=None,
):
    """Convenience: write .docx + .csv + .md in one call."""
    out_dir = Path(out_dir)
    write_apa_table(
        out_dir / f"{table_number}.docx",
        table_number,
        title,
        headers,
        rows,
        note=note,
        col_widths=col_widths,
        landscape=landscape,
        row_group_size=row_group_size,
    )
    write_csv_and_md(out_dir, table_number, title, headers, rows, note=note)
