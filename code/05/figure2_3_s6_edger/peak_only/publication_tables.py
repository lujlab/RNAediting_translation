#!/usr/bin/env python3
"""Create manuscript Table S6 and editable, auditable source tables."""

import argparse, csv, json, hashlib, math
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib.utils import simpleSplit
from docx import Document
from docx.shared import Inches, Pt
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, Alignment, Border, Side
from publication_figures import load, read, write_tsv, ORDER, LABELS

TITLE = "Table S6. Enrichment of edited codons among edgeR-defined ribosome-occupancy peaks."
HEAD = [
    "Library",
    "Unedited outside peak",
    "Unedited inside peak",
    "Unedited within peaks (%)",
    "Edited outside peak",
    "Edited inside peak",
    "Edited within peaks (%)",
    "P value (one-sided Fisher)",
]


def label(key):
    a, b = LABELS[key]
    if a == "Adar5G1":
        return "Dmel (Adar5G1) M"
    if a == "iso-1":
        return "Dmel (iso-1) F+M"
    if a == "WT":
        return "Dmel (WT) " + b
    return a + " " + b


def numeric_rows(data):
    return [
        [
            label(k),
            data[k]["unedited_total"] - data[k]["unedited_peak"],
            data[k]["unedited_peak"],
            data[k]["unedited_percent"],
            data[k]["edited_total"] - data[k]["edited_peak"],
            data[k]["edited_peak"],
            data[k]["edited_percent"],
            data[k]["p_value"],
        ]
        for k in ORDER
    ]


def formatted(rows):
    return [
        [
            r[0],
            f"{r[1]:,}",
            f"{r[2]:,}",
            f"{r[3]:.4f}",
            f"{r[4]:,}",
            f"{r[5]:,}",
            f"{r[6]:.4f}",
            f"{r[7]:.2E}",
        ]
        for r in rows
    ]


def make_pdf(out, rows):
    # Landscape A4, original three-line table structure with grouped column headings.
    c = canvas.Canvas(
        str(out / "TableS6.pdf"), pagesize=(841.89, 595.28), pageCompression=1
    )
    c.setTitle(TITLE)
    c.setFont("Times-Bold", 11)
    c.drawString(35, 555, TITLE)
    edges = [35, 181, 272, 351, 448, 529, 601, 697, 807]

    def center(col, y, s, size=9, bold=False):
        c.setFont("Times-Bold" if bold else "Times-Roman", size)
        c.drawCentredString((edges[col] + edges[col + 1]) / 2, y, s)

    c.setLineWidth(0.8)
    c.line(35, 540, 807, 540)
    center(0, 524, "Library", 10, True)
    for a, b, s in [(181, 448, "Unedited codon"), (448, 697, "Edited codon")]:
        c.setFont("Times-Bold", 10)
        c.drawCentredString((a + b) / 2, 524, s)
    center(7, 524, "P value", 10, True)
    c.setLineWidth(0.4)
    c.line(184, 517, 445, 517)
    c.line(451, 517, 694, 517)
    for col, s in [
        (1, ["Outside peak", "number"]),
        (2, ["Inside peak", "number"]),
        (3, ["Percentage within", "peaks (%)"]),
        (4, ["Outside peak", "number"]),
        (5, ["Inside peak", "number"]),
        (6, ["Percentage within", "peaks (%)"]),
    ]:
        for j, t in enumerate(s):
            center(col, 504 - j * 11, t)
    c.setLineWidth(0.5)
    c.line(35, 486, 807, 486)
    for i, row in enumerate(rows):
        for col, s in enumerate(row):
            center(col, 468 - i * 25, s, 9.3)
    bottom = 468 - (len(rows) - 1) * 25 - 12
    c.setLineWidth(0.8)
    c.line(35, bottom, 807, bottom)
    c.setFont("Times-Roman", 9)
    c.showPage()
    c.save()


def make_docx(out, rows):
    d = Document()
    sec = d.sections[0]
    sec.page_width = Inches(11.69)
    sec.page_height = Inches(8.27)
    sec.top_margin = sec.bottom_margin = Inches(0.55)
    sec.left_margin = sec.right_margin = Inches(0.48)
    style = d.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(9)
    style.paragraph_format.space_after = Pt(0)
    p = d.add_paragraph()
    p.add_run(TITLE).bold = True
    t = d.add_table(rows=2, cols=8)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.cell(0, 0).merge(t.cell(1, 0)).text = "Library"
    t.cell(0, 1).merge(t.cell(0, 3)).text = "Unedited codon"
    t.cell(0, 4).merge(t.cell(0, 6)).text = "Edited codon"
    t.cell(0, 7).merge(t.cell(1, 7)).text = "P value"
    for offset in [1, 4]:
        for j, s in enumerate(
            [
                "Outside peak\nnumber",
                "Inside peak\nnumber",
                "Percentage within\npeaks (%)",
            ]
        ):
            t.cell(1, offset + j).text = s
    for row in rows:
        cells = t.add_row().cells
        for cell, s in zip(cells, row):
            cell.text = s
    for i, row in enumerate(t.rows):
        for cell in row.cells:
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                for run in p.runs:
                    run.bold = i < 2
            pr = cell._tc.get_or_add_tcPr()
            border = OxmlElement("w:tcBorders")
            for side in ["top", "left", "bottom", "right"]:
                node = OxmlElement("w:" + side)
                node.set(
                    qn("w:val"),
                    (
                        "single"
                        if (side == "top" and i == 0)
                        or (side == "bottom" and i in [1, len(t.rows) - 1])
                        else "nil"
                    ),
                )
                node.set(qn("w:sz"), "6")
                border.append(node)
            pr.append(border)
    d.save(out / "TableS6.docx")


def add_sheet(wb, name, head, rows):
    ws = wb.create_sheet(name)
    ws.append(head)
    for row in rows:
        ws.append(row)
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = ws.dimensions
    for c in ws[1]:
        c.font = Font(bold=True, name="Times New Roman")
        c.alignment = Alignment(wrap_text=True)
        c.border = Border(bottom=Side(style="thin"))
    for col in ws.columns:
        ws.column_dimensions[col[0].column_letter].width = (
            27 if col[0].column == 1 else 23
        )
    return ws


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edgeR-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    src = (
        args.edgeR_dir
        / "edited_unedited_summary_codonA/edited_unedited_peaks_by_library.tsv"
    )
    hairpin = (
        args.edgeR_dir
        / "edited_unedited_hairpin_summary/edited_unedited_hairpin_in_peaks_by_library.long.tsv"
    )
    data = load(src)
    hd = load(hairpin, True)
    assert len(data) == len(hd) == 14
    nums = numeric_rows(data)
    fmt = formatted(nums)
    for file, rows in [
        ("TableS6.tsv", nums),
        ("TableS6.formatted.tsv", fmt),
        ("F3A_hairpin_enrichment.tsv", numeric_rows(hd)),
    ]:
        with (out / file).open("w") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(HEAD)
            w.writerows(rows)
    make_pdf(out, fmt)
    make_docx(out, fmt)
    wb = Workbook()
    wb.remove(wb.active)
    ws = add_sheet(wb, "TableS6", HEAD, nums)
    for row in ws.iter_rows(min_row=2):
        for i in [1, 2, 4, 5]:
            row[i].number_format = "#,##0"
        for i in [3, 6]:
            row[i].number_format = "0.0000"
        row[7].number_format = "0.00E+00"
    add_sheet(wb, "F3A_hairpin", HEAD, numeric_rows(hd))
    for name, path in [
        ("All14_peak_source", src),
        ("All14_hairpin_in", hairpin),
        (
            "All14_hairpin_out",
            hairpin.with_name("edited_unedited_hairpin_out_peaks_by_library.long.tsv"),
        ),
    ]:
        rr = read(path)
        add_sheet(
            wb,
            name,
            list(rr[0]),
            [
                [
                    (
                        float(v)
                        if v
                        and v not in ["NA", "NaN"]
                        and all(ch in "-+.eE0123456789" for ch in v)
                        else v
                    )
                    for v in r.values()
                ]
                for r in rr
            ],
        )
    wb.save(out / "TableS6_and_hairpin.xlsx")
    reread = load_workbook(out / "TableS6_and_hairpin.xlsx", data_only=True)
    assert reread["TableS6"].max_row == 11
    for a, b in zip(nums, reread["TableS6"].iter_rows(min_row=2, values_only=True)):
        for x, y in zip(a, b):
            assert (
                x == y
                if isinstance(x, str)
                else math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-14)
            )
    pass
    print("Table S6 and hairpin tables complete:", out, flush=True)


if __name__ == "__main__":
    main()
