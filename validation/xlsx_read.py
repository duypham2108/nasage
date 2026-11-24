"""Minimal .xlsx reader: zipfile + XML, no openpyxl.

openpyxl is not installed in this environment and this project does not install packages
into the user's environment to read one file. An .xlsx is a zip of XML parts, so the sheet
grid can be recovered directly.

Handles: shared strings, inline strings, numbers, booleans, and gaps in the cell grid.
Does not handle: dates as dates (they come back as Excel serial numbers), formulas
(cached values are returned), styles.
"""
import re
import zipfile
from xml.etree import ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RELNS = "{http://schemas.openxmlformats.org/package/2006/relationships}"


def _col_to_idx(ref):
    m = re.match(r"([A-Z]+)", ref or "")
    if not m:
        return None
    n = 0
    for ch in m.group(1):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def sheet_names(path):
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        return [s.get("name") for s in wb.iter(f"{NS}sheet")]


def read_sheet(path, sheet=0):
    """Return a list of rows (list of cell values as str/float/None)."""
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in ss.iter(f"{NS}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{NS}t")))

        wb = ET.fromstring(z.read("xl/workbook.xml"))
        sheets = list(wb.iter(f"{NS}sheet"))
        names = [s.get("name") for s in sheets]
        idx = names.index(sheet) if isinstance(sheet, str) else sheet
        rid = sheets[idx].get(
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = None
        for r in rels.iter(f"{RELNS}Relationship"):
            if r.get("Id") == rid:
                target = r.get("Target")
        target = (target or f"worksheets/sheet{idx + 1}.xml").lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target

        root = ET.fromstring(z.read(target))
        rows = []
        for row in root.iter(f"{NS}row"):
            cells = {}
            for c in row.iter(f"{NS}c"):
                i = _col_to_idx(c.get("r"))
                t = c.get("t")
                if t == "inlineStr":
                    v = "".join(x.text or "" for x in c.iter(f"{NS}t"))
                else:
                    ve = c.find(f"{NS}v")
                    v = ve.text if ve is not None else None
                    if v is not None and t == "s":
                        v = shared[int(v)]
                    elif v is not None and t not in ("str", "e"):
                        try:
                            v = float(v)
                        except ValueError:
                            pass
                if i is not None:
                    cells[i] = v
            rows.append([cells.get(i) for i in range(max(cells) + 1)] if cells else [])
        return rows


def to_frame(path, sheet=0, header_row=0):
    import pandas as pd
    rows = read_sheet(path, sheet)
    if not rows:
        return pd.DataFrame()
    width = max(len(r) for r in rows)
    rows = [r + [None] * (width - len(r)) for r in rows]
    hdr = [str(h) if h is not None else f"col{i}" for i, h in enumerate(rows[header_row])]
    return pd.DataFrame(rows[header_row + 1:], columns=hdr)


def iter_rows(path, sheet=0, max_rows=None):
    """Stream rows from a large sheet without holding the XML in memory.

    Needed for the NOVA normalized matrix, whose sheet1.xml is 460 MB uncompressed.
    """
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in ss.iter(f"{NS}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{NS}t")))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        sheets = list(wb.iter(f"{NS}sheet"))
        names = [s.get("name") for s in sheets]
        idx = names.index(sheet) if isinstance(sheet, str) else sheet
        rid = sheets[idx].get(
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = None
        for r in rels.iter(f"{RELNS}Relationship"):
            if r.get("Id") == rid:
                target = r.get("Target")
        target = (target or f"worksheets/sheet{idx + 1}.xml").lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        n = 0
        with z.open(target) as fh:
            for event, el in ET.iterparse(fh, events=("end",)):
                if el.tag != f"{NS}row":
                    continue
                cells = {}
                for c in el.iter(f"{NS}c"):
                    i = _col_to_idx(c.get("r"))
                    t = c.get("t")
                    if t == "inlineStr":
                        v = "".join(x.text or "" for x in c.iter(f"{NS}t"))
                    else:
                        ve = c.find(f"{NS}v")
                        v = ve.text if ve is not None else None
                        if v is not None and t == "s":
                            v = shared[int(v)]
                        elif v is not None and t not in ("str", "e"):
                            try:
                                v = float(v)
                            except ValueError:
                                pass
                    if i is not None:
                        cells[i] = v
                yield [cells.get(i) for i in range(max(cells) + 1)] if cells else []
                el.clear()
                n += 1
                if max_rows and n >= max_rows:
                    return


if __name__ == "__main__":
    import sys
    p = sys.argv[1]
    print("sheets:", sheet_names(p))
    for s in sheet_names(p):
        rows = read_sheet(p, s)
        print(f"--- {s!r}: {len(rows)} rows x "
              f"{max((len(r) for r in rows), default=0)} cols")
        for r in rows[:4]:
            print("   ", [str(x)[:20] if x is not None else "" for x in r[:14]])
