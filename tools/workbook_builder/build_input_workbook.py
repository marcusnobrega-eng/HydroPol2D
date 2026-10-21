"""Build the HydroPol2D input workbook from input_catalog.yaml.

Values (and the map and gauge lists) are carried over from an existing workbook, so the
same script rebuilds the main workbook and each example. Every input sits in one row of
one block (setting, key, value, unit, default, used when) with a hover note on its value
cell; a block starts with its switch, and its rows grey out while the switch is 0.

The class tables and time series get one tab each, copied cell for cell from the old
separate workbooks next to SOURCE (or from SOURCE's own tabs when it is already a single
case workbook), so the model's readtable calls see the same tables. A series the case's
switches do not read keeps only its header rows, as a template.

Usage: python build_input_workbook.py input_catalog.yaml SOURCE.xlsx OUTPUT.xlsx
"""
import datetime as dt
import sys
import textwrap
from pathlib import Path

import openpyxl
import yaml
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
from openpyxl.comments import Comment
from openpyxl.formatting.rule import FormulaRule
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation

FONT = "Helvetica"
INK, SOFT, FAINT, RULE = "2C2C2A", "5F5E5A", "B4B2A9", "E3E1D9"
TEAL_DARK, TEAL, TEAL_LIGHT = "085041", "0F6E56", "E1F5EE"
GREEN_TEXT, GREEN_EDGE, GREEN_FILL = "27500A", "3B6D11", "EAF3DE"
HEADER_FILL, OFF_FILL, RED = "F1EFE8", "F7F6F2", "A32D2D"
# Block header colours: (fill, text, left bar)
BLOCK_STYLES = {
    "normal": ("E1F5EE", "085041", "0F6E56"),
    "advanced": ("F1EFE8", "444441", "888780"),
    "paused": ("FAEEDA", "633806", "BA7517"),
}


def font(size=11, bold=False, color=INK, italic=False):
    return Font(name=FONT, size=size, bold=bold, color=color, italic=italic)


def fill(color):
    return PatternFill("solid", start_color=color, end_color=color)


def side(color, style="thin"):
    return Side(style=style, color=color)


BOX = Border(left=side(FAINT), right=side(FAINT), top=side(FAINT), bottom=side(FAINT))
SWITCH_BOX = Border(left=side(GREEN_EDGE), right=side(GREEN_EDGE), top=side(GREEN_EDGE), bottom=side(GREEN_EDGE))
ROW_RULE = Border(bottom=side(RULE))
COLUMNS = ["Setting", "Key (model name)", "Value", "Unit", "Default", "Used when"]


# --------------------------------------------------------------------- source values
def source_values(path):
    """Cached cell values (data_only) of every key/flag and the four lists."""
    wb = openpyxl.load_workbook(path, data_only=True)
    values, lists = {}, {}
    for ws in wb.worksheets:
        for row in ws.iter_rows(max_row=400, max_col=60):
            for c in row:
                if isinstance(c.value, str) and c.value.strip() and c.column < ws.max_column:
                    key = c.value.strip()
                    if key not in values:
                        values[key] = ws.cell(c.row, c.column + 1).value

    def read_list(ws, anchor, columns, header_is_anchor=False, optional=()):
        # Optional columns may sit anywhere in the header row (or be absent: their default).
        for row in ws.iter_rows(max_row=50, max_col=60):
            for c in row:
                if c.value in (anchor, anchor.replace("Satellite", "Sattelite")):   # old spelling
                    header_row = c.row if header_is_anchor else c.row + 1
                    # The list's own columns: first match at or right of its title (lists
                    # side by side share a header row).
                    cols = []
                    for name in columns:
                        hit = [x.column for x in ws[header_row] if x.value == name and x.column >= c.column] or \
                              [x.column for x in ws[header_row] if x.value == name]
                        cols.append(hit[0] if hit else None)
                    if None in cols:
                        continue
                    extra = []
                    for spec in optional:
                        hit = [x.column for x in ws[header_row] if x.value == spec["name"] and x.column >= c.column] or \
                              [x.column for x in ws[header_row] if x.value == spec["name"]]
                        extra.append((hit[0] if hit else None, spec["default"]))
                    out, r = [], header_row + 1
                    while ws.cell(r, cols[0]).value not in (None, ""):
                        record = [ws.cell(r, k).value for k in cols]
                        for k, default in extra:
                            value = ws.cell(r, k).value if k else None
                            record.append(default if value in (None, "") else value)
                        out.append(record)
                        r += 1
                    return out
        return []

    gd = wb["General_Data"] if "General_Data" in wb.sheetnames else wb.worksheets[0]
    return values, lambda spec: read_list(gd, spec["anchor"], spec["columns"], spec.get("header_is_anchor", False),
                                          spec.get("optional_columns", ()))


# --------------------------------------------------------------------------- notes
def note_text(item):
    lines = [item["note"]]
    extra = [f"Unit: {item['unit']}"] if item.get("unit") else []
    if item["default"] == "required":
        extra.append("Required")
    elif item["default"] != "":
        extra.append(f"Default: {item['default']}")
    extra.append(f"Used when: {item['used_when']}")
    if item.get("options"):
        extra.append("Options: " + ", ".join(str(o) for o in item["options"]))
    lines.append(" · ".join(extra))
    return "\n\n".join(lines)


def add_note(cell, text):
    wrapped = sum(len(textwrap.wrap(p, 56)) or 1 for p in text.split("\n"))
    comment = Comment(text, "HydroPol2D")
    comment.width, comment.height = 380, max(90, 18 * wrapped + 24)
    cell.comment = comment


# ------------------------------------------------------------------------- layout
def title_rows(ws, title, subtitle, last_col="G"):
    ws.sheet_view.showGridLines = False
    ws.column_dimensions["A"].width = 2.5
    ws.row_dimensions[1].height = 36
    ws.merge_cells(f"B1:{last_col}1")
    c = ws["B1"]
    c.value, c.font = title, font(18, True, "FFFFFF")
    c.alignment = Alignment(vertical="center", indent=1)
    for col in range(1, openpyxl.utils.column_index_from_string(last_col) + 1):
        ws.cell(1, col).fill = fill(TEAL_DARK)
    ws.merge_cells(f"B2:{last_col}2")
    ws["B2"].value = subtitle
    ws["B2"].font = font(10, color=SOFT)
    ws["B2"].alignment = Alignment(vertical="center", indent=1)
    ws.row_dimensions[2].height = 22


def column_header(ws, row, names):
    for i, name in enumerate(names):
        c = ws.cell(row, 2 + i, name)
        c.font, c.fill = font(10, True, "444441"), fill(HEADER_FILL)
        c.border = Border(bottom=side("888780"))
        c.alignment = Alignment(vertical="center", indent=1 if i in (0, 1, 5) else 0,
                                horizontal="center" if i in (2, 3, 4) else None)
    ws.row_dimensions[row].height = 22


def block_header(ws, row, title, key_text=None, style="normal"):
    band, text, bar = BLOCK_STYLES[style]
    for col in range(2, 8):
        c = ws.cell(row, col)
        c.fill = fill(band)
    ws.cell(row, 2).border = Border(left=side(bar, "thick"))
    c = ws.cell(row, 2, title)
    c.font = font(12, True, text)
    c.alignment = Alignment(vertical="center", indent=1)
    if key_text:
        k = ws.cell(row, 3, key_text)
        k.font = font(10, color=TEAL)
        k.alignment = Alignment(vertical="center", indent=1)
    ws.row_dimensions[row].height = 24


def label_cell(cell, label, required):
    if required:
        cell.value = CellRichText([TextBlock(InlineFont(rFont=FONT, sz=11, color=INK), label),
                                   TextBlock(InlineFont(rFont=FONT, sz=11, color=RED, b=True), " *")])
    else:
        cell.value = label
        cell.font = font(11)
    cell.alignment = Alignment(vertical="center", indent=1)


def as_value(item, raw):
    if raw is None or (isinstance(raw, str) and raw.strip().lower() in ("", "none", "nan")):
        return None
    if item["kind"] in ("switch", "choice", "number"):
        try:
            number = float(raw)
            return int(number) if number.is_integer() and item["kind"] != "number" else number
        except (TypeError, ValueError):
            return raw
    return raw


def input_row(ws, row, item, value, validations):
    label_cell(ws.cell(row, 2), item["label"], item["default"] == "required" and item["kind"] != "switch")
    key = ws.cell(row, 3, item["key"])
    key.font, key.alignment = font(10, color=SOFT), Alignment(vertical="center", indent=1)
    v = ws.cell(row, 4, value)
    kind = item["kind"]
    if kind in ("switch", "choice"):
        v.font, v.fill, v.border = font(11, True, GREEN_TEXT), fill(GREEN_FILL), SWITCH_BOX
        v.alignment = Alignment(horizontal="center", vertical="center")
        options = item.get("options", [0, 1])
        validations.setdefault(",".join(str(o) for o in options), []).append(v.coordinate)
    else:
        v.font, v.fill, v.border = font(11), fill("FFFFFF"), BOX
        v.alignment = Alignment(vertical="center", horizontal="left" if kind in ("path", "text") else "right",
                                indent=1 if kind in ("path", "text") else 0)
        if kind == "date" and isinstance(value, (dt.datetime, dt.date)):
            v.number_format = "yyyy-mm-dd hh:mm"
    add_note(v, note_text(item))
    unit = ws.cell(row, 5, item.get("unit") or "")
    default = item["default"]
    dflt = ws.cell(row, 6, "required" if default == "required" else ("" if default == "" else default))
    used = ws.cell(row, 7, item["used_when"])
    for c, f in ((unit, font(10, color=SOFT)), (dflt, font(10, color=SOFT)), (used, font(10, color=SOFT, italic=True))):
        c.font = f
        c.alignment = Alignment(vertical="center", horizontal="center" if c in (unit, dflt) else None,
                                indent=1 if c is used else 0)
    for col in (2, 3, 5, 6, 7):
        ws.cell(row, col).border = ROW_RULE
    ws.row_dimensions[row].height = 20


def grey_out(ws, first, last, switch_cell):
    """Rows first..last turn grey while the block's switch is 0."""
    ws.conditional_formatting.add(
        f"B{first}:G{last}",
        FormulaRule(formula=[f"${switch_cell[0]}${switch_cell[1:]}=0"], font=Font(color=FAINT), fill=fill(OFF_FILL)),
    )


def build_tab(wb, tab, values, read_list, developer):
    ws = wb.create_sheet(tab["name"])
    title_rows(ws, f"HydroPol2D · {tab['name']}", f"{tab['purpose']}   Developed by {developer}.")
    widths = {"B": 36, "C": 40, "D": 64 if tab.get("wide_value") else 18, "E": 12, "F": 15, "G": 34}
    for col, width in widths.items():
        ws.column_dimensions[col].width = width
    column_header(ws, 4, COLUMNS)
    ws.freeze_panes = "A5"
    validations, switches, row = {}, {}, 6
    for block in tab["blocks"]:
        block_header(ws, row, block["title"], style=block.get("style", "normal"))
        row += 1
        first = row
        for item in block["inputs"]:
            value = as_value(item, values.get(item["key"]))
            if value is None and item["kind"] in ("switch", "choice") and item["default"] not in ("", "required"):
                value = item["default"]  # a switch always shows its state
            input_row(ws, row, item, value, validations)
            if block.get("switch") == item["key"]:
                switches[item["key"]] = f"D{row}"
            row += 1
        control = block.get("switch") or block.get("controlled_by")
        if control and control in switches:
            start = first + 1 if block.get("switch") else first
            if start <= row - 1:
                grey_out(ws, start, row - 1, switches[control])
        row += 1
    for spec in tab.get("lists", []):
        anchor_in_header = spec.get("header_is_anchor", False)
        block_header(ws, row, spec["title"], None if anchor_in_header else spec["anchor"])
        add_note(ws.cell(row, 2), spec["note"])
        row += 1
        optional = spec.get("optional_columns", [])
        names = spec["columns"] + [o["name"] for o in optional]
        for i, name in enumerate(names):
            c = ws.cell(row, 2 + i, name)
            c.font, c.fill, c.border = font(10, True, "444441"), fill(HEADER_FILL), Border(bottom=side("888780"))
            c.alignment = Alignment(vertical="center", indent=1, wrap_text=True)
        if any(len(n) + 2 > widths.get(get_column_letter(2 + i), 12) for i, n in enumerate(names)):
            ws.row_dimensions[row].height = 30   # a long header wraps onto two lines
        for i, o in enumerate(optional):
            add_note(ws.cell(row, 2 + len(spec["columns"]) + i), o["note"])
        row += 1
        first_record = row
        for record in read_list(spec):
            for i, value in enumerate(record):
                c = ws.cell(row, 2 + i, value)
                c.font, c.border = font(10), ROW_RULE
                c.alignment = Alignment(vertical="center", indent=1)
            row += 1
        for i, o in enumerate(optional):
            if o.get("options") and row > first_record:
                letter = get_column_letter(2 + len(spec["columns"]) + i)
                validations.setdefault(o["options"], []).append(f"{letter}{first_record}:{letter}{row - 1}")
        row += 1
    for options, cells in validations.items():
        dv = DataValidation(type="list", formula1=f'"{options}"', allow_blank=False,
                            showErrorMessage=True, errorTitle="HydroPol2D", error=f"Choose one of: {options}")
        for ref in cells:
            dv.add(ref)
        ws.add_data_validation(dv)
    return ws


# ------------------------------------------------------------ class tables and series
REPO_SHEETS = Path(__file__).resolve().parents[2] / "Input_Data_Sheets"
TAB_COLORS = {"settings": TEAL, "class": GREEN_EDGE, "series": "5F5E5A"}


def as_switch(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def table_plan(spec, source_path, values):
    """Where a table tab comes from and how much of it to copy."""
    source_path = Path(source_path)
    origin, where = None, None
    try:
        names = openpyxl.load_workbook(source_path, read_only=True).sheetnames
    except Exception:
        names = []
    if spec["name"] in names:
        origin, where = "case tab", (source_path, spec["name"])
    else:
        for rel in spec["files"]:
            if (source_path.parent / rel).is_file():
                origin, where = "case file", (source_path.parent / rel, None)
                break
        else:
            for rel in spec["files"]:
                if (REPO_SHEETS / rel).is_file():
                    origin, where = "default", (REPO_SHEETS / rel, None)
                    break
    rule = spec.get("read_when")
    read = rule is None or (all(as_switch(values.get(f)) == 1 for f in rule["on"]) and
                            not any(as_switch(values.get(f)) == 1 for f in rule["off"]))
    # A default series belongs to another case: keep only its headers.
    full = origin is not None and (origin != "default" or "template_rows" not in spec) and \
        (spec["kind"] != "series" or read)
    if origin is None:
        status = "no source found: empty tab"
    elif full:
        status = {"case tab": "carried over", "case file": f"from {where[0].name}",
                  "default": "repository defaults"}[origin]
    else:
        status = "template (header rows only): " + ("its switches are off" if not read else "no case data")
    return {"origin": origin, "where": where, "full": full, "read": read, "status": status}


def build_table_tab(wb, spec, plan, values):
    ws = wb.create_sheet(spec["name"])
    ws.sheet_properties.tabColor = TAB_COLORS[spec["kind"]]
    ws.sheet_view.showGridLines = False
    if plan["origin"] is None:
        return ws
    path, sheet = plan["where"]
    big = path.stat().st_size > 2e6
    src_wb = openpyxl.load_workbook(path, data_only=True, read_only=big)
    src = src_wb[sheet] if sheet else src_wb.worksheets[0]
    last_row = None if plan["full"] else spec.get("template_rows", 3)
    header_rows = set(spec.get("header_rows", []))
    notes = spec.get("column_notes", {})
    max_col, title_cell = 0, None
    for row in src.iter_rows(max_row=last_row):
        for c in row:
            if c.value is None:
                continue
            d = ws.cell(c.row, c.column, c.value)   # same row and column: readtable sees the same table
            max_col = max(max_col, c.column)
            if c.number_format and c.number_format != "General":
                d.number_format = c.number_format
            if c.row == 1 and title_cell is None:
                title_cell = d
                d.font = font(16, True, "FFFFFF")
            elif c.row in header_rows:
                d.font = font(10, True, "444441")
                d.alignment = Alignment(vertical="center", wrap_text=True)
            else:
                d.font = font(10)
            comment = None if big else c.comment
            if comment is not None:
                add_note(d, comment.text.strip())
            elif c.row in header_rows and isinstance(c.value, str) and c.value.strip() in notes:
                add_note(d, notes[c.value.strip()])
    if not big:
        for rng in src.merged_cells.ranges:
            if last_row is None or rng.max_row <= last_row:
                ws.merge_cells(str(rng))
    # Mesh from the repository defaults: take the case's values from General_Data.
    if spec["name"] == "Mesh" and plan["origin"] == "default":
        for r in range(1, ws.max_row + 1):
            label = ws.cell(r, 2).value
            if isinstance(label, str) and label.strip() in values and values[label.strip()] is not None:
                ws.cell(r, 3).value = values[label.strip()]
    max_col = max(max_col, 1)
    for col in range(1, max_col + 1):   # bands inside the table's own columns only
        ws.cell(1, col).fill = fill(TEAL_DARK)
        for r in header_rows:
            if last_row is None or r <= last_row:
                ws.cell(r, col).fill = fill(HEADER_FILL)
                ws.cell(r, col).border = Border(bottom=side("888780"))
    ws.row_dimensions[1].height = 30
    if title_cell is not None:
        title_cell.alignment = Alignment(vertical="center")
        note = spec["note"] + "\n\nThis tab: " + plan["status"] + "."
        add_note(title_cell, note)
    for col in range(1, min(max_col, 60) + 1):
        letter = get_column_letter(col)
        width = src.column_dimensions[letter].width if not big and letter in src.column_dimensions else None
        head = max((len(str(ws.cell(r, col).value or "")) for r in header_rows), default=0)
        ws.column_dimensions[letter].width = min(max(width or 0, head + 2, 10), 42)
    if header_rows:
        ws.freeze_panes = ws.cell(max(header_rows) + 1, 1)
    return ws


def start_tab(wb, catalog, built):
    ws = wb.active
    ws.title = "Start here"
    dev = catalog["developer"]
    title_rows(ws, "HydroPol2D input workbook", f"Developed by {dev}", last_col="E")
    for col, width in {"B": 30, "C": 62, "D": 18, "E": 18}.items():
        ws.column_dimensions[col].width = width
    row = 4
    ws.cell(row, 2, "Developed by").font = font(10, True, SOFT)
    ws.cell(row, 3, dev).font = font(14, True, TEAL_DARK)
    row += 2

    def heading(text):
        nonlocal row
        block_header(ws, row, text)
        for col in (6, 7):
            ws.cell(row, col).fill = PatternFill(None)
        row += 1

    heading("How to use this workbook")
    steps = [
        "Go through the tabs from left to right; each groups one part of the model.",
        "Type only in the white boxes (values) and pick 0 or 1 in the green boxes (switches).",
        "Hover over any value to read what it means, its unit, its default and when it is used.",
        "A block greys out while its switch is 0: its values are then ignored.",
        "Settings marked * are required whenever their block is on.",
    ]
    for i, text in enumerate(steps, 1):
        ws.cell(row, 2, f"{i}.").font = font(11, True, TEAL)
        ws.cell(row, 2).alignment = Alignment(horizontal="right")
        c = ws.cell(row, 3, text)
        c.font, c.alignment = font(11), Alignment(wrap_text=True, vertical="center")
        row += 1
    row += 1
    heading("Colour legend")
    samples = [
        ("White box", "A value you type", fill("FFFFFF"), BOX, font(11)),
        ("Green box", "A switch: 0 = off, 1 = on (some take 1, 2, 3)", fill(GREEN_FILL), SWITCH_BOX, font(11, True, GREEN_TEXT)),
        ("Grey row", "Its block is switched off", fill(OFF_FILL), None, font(11, color=FAINT)),
        ("Red *", "Required while the block is on", None, None, font(11, True, RED)),
        ("Advanced", "Grey block: legacy settings, keep the defaults unless you know why", fill(BLOCK_STYLES["advanced"][0]), None, font(11, True, BLOCK_STYLES["advanced"][1])),
        ("Paused", "Amber block: a paused feature; keep it at 0", fill(BLOCK_STYLES["paused"][0]), None, font(11, True, BLOCK_STYLES["paused"][1])),
    ]
    for name, meaning, f, b, ft in samples:
        c = ws.cell(row, 2, name)
        c.font = ft
        if f is not None:
            c.fill = f
        if b is not None:
            c.border = b
        c.alignment = Alignment(indent=1, vertical="center")
        ws.cell(row, 3, meaning).font = font(11)
        ws.row_dimensions[row].height = 20
        row += 1
    row += 1
    heading("Tabs")
    for tab in catalog["tabs"]:
        c = ws.cell(row, 2, tab["name"])
        c.hyperlink = f"#'{tab['name']}'!A1"
        c.font = Font(name=FONT, size=11, bold=True, color=TEAL, underline="single")
        c.alignment = Alignment(indent=1)
        p = ws.cell(row, 3, tab["purpose"])
        p.font, p.alignment = font(11), Alignment(wrap_text=True, vertical="center")
        row += 1
    row += 1
    heading("Tables and series")
    for spec in catalog.get("tables", []):
        c = ws.cell(row, 2, spec["name"])
        c.hyperlink = f"#'{spec['name']}'!A1"
        c.font = Font(name=FONT, size=11, bold=True, color=TEAL, underline="single")
        c.alignment = Alignment(indent=1)
        status = (built or {}).get(spec["name"], {}).get("status", "")
        p = ws.cell(row, 3, spec["purpose"] + (f"  ({status})" if status else ""))
        p.font, p.alignment = font(11), Alignment(wrap_text=True, vertical="center")
        row += 1
    row += 1
    ws.cell(row, 2, f"Built {dt.date.today().isoformat()} from input_catalog.yaml "
                    f"({sum(len(b['inputs']) for t in catalog['tabs'] for b in t['blocks'])} inputs).").font = font(9, color=FAINT)


def main(catalog_path, source_path, out_path):
    catalog = yaml.safe_load(open(catalog_path))
    values, read_list = source_values(source_path)
    plans = {spec["name"]: table_plan(spec, source_path, values) for spec in catalog.get("tables", [])}
    wb = openpyxl.Workbook()
    start_tab(wb, catalog, plans)
    for tab in catalog["tabs"]:
        build_tab(wb, tab, values, read_list, catalog["developer"]).sheet_properties.tabColor = TAB_COLORS["settings"]
    for spec in catalog.get("tables", []):
        build_table_tab(wb, spec, plans[spec["name"]], values)
        print(f"  {spec['name']}: {plans[spec['name']]['status']}")
    wb.properties.creator = catalog["developer"]
    wb.properties.title = catalog["title"]
    wb.properties.subject = "HydroPol2D model inputs"
    wb.save(out_path)
    print(f"wrote {out_path} ({len(wb.sheetnames)} tabs)")


if __name__ == "__main__":
    main(*sys.argv[1:4])
