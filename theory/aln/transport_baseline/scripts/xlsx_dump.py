"""Minimal dependency-free .xlsx reader (zip + XML) for the Phonon Olympics spreadsheets.

read_sheets(path) -> {sheet_name: {row_number: {col_index: value}}}
Cached values are returned (formulas are not evaluated); numeric strings are
converted to float.  Command line: python xlsx_dump.py file.xlsx [...] prints rows.
"""
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def _col_index(ref):
    letters = re.match(r"([A-Z]+)", ref).group(1)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def read_sheets(path):
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    relmap = {r.get("Id"): r.get("Target") for r in rels}
    out = {}
    for sh in wb.find("m:sheets", NS):
        name = sh.get("name")
        target = relmap[sh.get("{%s}id" % NS["r"])].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        root = ET.fromstring(z.read(target))
        rows = {}
        for row in root.iter("{%s}row" % NS["m"]):
            cells = {}
            for c in row.findall("m:c", NS):
                t = c.get("t")
                v = c.find("m:v", NS)
                isel = c.find("m:is", NS)
                if t == "s" and v is not None:
                    val = shared[int(v.text)]
                elif t == "inlineStr" and isel is not None:
                    val = "".join(x.text or "" for x in isel.iter("{%s}t" % NS["m"]))
                else:
                    val = _num(v.text) if v is not None else None
                if val is not None:
                    cells[_col_index(c.get("r"))] = val
            if cells:
                rows[int(row.get("r"))] = cells
        out[name] = rows
    return out


if __name__ == "__main__":
    for p in sys.argv[1:]:
        print("=" * 20, p)
        for sheet, rows in read_sheets(p).items():
            print("--- sheet", sheet)
            for r, cells in rows.items():
                mx = max(cells)
                print(r, "|", " | ".join(str(cells.get(i, "")) for i in range(mx + 1)))
