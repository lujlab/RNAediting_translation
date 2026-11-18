#!/usr/bin/env python3
"""Draw manuscript panels as native vector PDF objects (no raster, masks or clips)."""

import argparse, csv, hashlib, json, math, subprocess, textwrap
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from fontTools.ttLib import TTFont as FontToolsFont
from scipy.stats import fisher_exact

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[3]
FONTS = CODE.parents[1] / "figure_s6_s8_pausepred/reference_style_fonts"
ORDER = [
    "dmel_Adar5G1_M",
    "dmel_ISO1_FM",
    "dmel_WBg_F",
    "dmel_WBg_M",
    "dsim_w501_F",
    "dsim_w501_M",
    "dsec_Dsec_F",
    "dsec_Dsec_M",
    "dere_Dere_F",
    "dere_Dere_M",
]
MALE = [ORDER[i] for i in [0, 1, 3, 5, 7, 9]]
FEMALE = [ORDER[i] for i in [2, 4, 6, 8]]
LABELS = dict(
    zip(
        ORDER,
        [
            ("Adar5G1", "M"),
            ("iso-1", "M+F"),
            ("WT", "F"),
            ("WT", "M"),
            ("Dsim", "F"),
            ("Dsim", "M"),
            ("Dsec", "F"),
            ("Dsec", "M"),
            ("Dere", "F"),
            ("Dere", "M"),
        ],
    )
)
EDITED = (0.8, 0.078, 0.498)
UNEDIT = (0.29, 0.631, 0.710)
INK = (0.08, 0.08, 0.08)
CMAPS = {}


def setup_fonts():
    for suffix in ["Regular", "Bold", "Italic"]:
        p = FONTS / ("Arial-" + suffix + ".ttf")
        font = FontToolsFont(str(p))
        cm = font.getBestCmap()
        CMAPS[suffix] = {
            n
            for n, v in cm.items()
            if n == 32 or getattr(font["glyf"][v], "numberOfContours", 0)
        }
        pdfmetrics.registerFont(TTFont("Arial-" + suffix, str(p)))
        fallback = subprocess.check_output(
            ["fc-match", "-f", "%{file}", "Liberation Sans:style=" + suffix],
            universal_newlines=True,
        ).strip()
        pdfmetrics.registerFont(TTFont("Fallback-" + suffix, fallback))


def font_for(s, bold=False, italic=False):
    suffix = "Bold" if bold else "Italic" if italic else "Regular"
    return (
        "Arial-" if all(ord(ch) in CMAPS[suffix] for ch in s) else "Fallback-"
    ) + suffix


def text(c, x, y, s, size=7, color=INK, align="left", bold=False, italic=False):
    c.setFillColorRGB(*color)
    c.setFont(font_for(s, bold, italic), size)
    {"left": c.drawString, "center": c.drawCentredString, "right": c.drawRightString}[
        align
    ](x, y, s)


def line(c, x1, y1, x2, y2, width=0.45, color=INK):
    c.setStrokeColorRGB(*color)
    c.setLineWidth(width)
    c.line(x1, y1, x2, y2)


def rect(c, x, y, w, h, color):
    c.setFillColorRGB(*color)
    c.rect(x, y, w, h, fill=1, stroke=0)


def arrow(c, x1, y1, x2, y2, color=INK):
    line(c, x1, y1, x2, y2, 0.45, color)
    angle = math.atan2(y2 - y1, x2 - x1)
    for d in [-0.55, 0.55]:
        line(
            c,
            x2,
            y2,
            x2 - 2.7 * math.cos(angle + d),
            y2 - 2.7 * math.sin(angle + d),
            0.45,
            color,
        )


def mark(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def load(path, hairpin=False):
    result = {}
    for r in read(path):
        assert r["event"] == "stalling"
        if hairpin:
            key = r["library"]
            e, et, u, ut = [
                int(r[k])
                for k in [
                    "edited_inside_event",
                    "edited_total_region_tested",
                    "unedited_inside_event",
                    "unedited_total_region_tested",
                ]
            ]
            p = float(r["event_fisher_p_greater"])
        else:
            key = r["group_id"]
            e, et, u, ut = [
                int(r[k])
                for k in [
                    "edited_event",
                    "edited_total",
                    "unedited_event",
                    "unedited_total",
                ]
            ]
            p = float(r["fisher_p"])
        assert min(et, ut) > 0 and 0 <= e <= et and 0 <= u <= ut
        check = fisher_exact([[e, et - e], [u, ut - u]], alternative="greater")[1]
        assert math.isclose(check, p, rel_tol=1e-7, abs_tol=1e-12), (key, p, check)
        result[key] = dict(
            group_id=key,
            edited_peak=e,
            edited_total=et,
            edited_percent=e / et * 100,
            unedited_peak=u,
            unedited_total=ut,
            unedited_percent=u / ut * 100,
            p_value=p,
            significance=mark(p),
        )
    return result


def group_label(c, x, y, key, size):
    a, b = LABELS[key]
    sex_size = size
    if a == "Adar5G1":
        size *= 0.86
        base = font_for("Adar", italic=True)
        small = font_for("5G1", italic=True)
        w = pdfmetrics.stringWidth("Adar", base, size) + pdfmetrics.stringWidth(
            "5G1", small, size * 0.63
        )
        text(c, x - w / 2, y, "Adar", size, italic=True)
        text(
            c,
            x - w / 2 + pdfmetrics.stringWidth("Adar", base, size),
            y + size * 0.35,
            "5G1",
            size * 0.63,
            italic=True,
        )
    else:
        text(c, x, y, a, size, align="center")
    text(c, x, y - 8.3, b, sex_size, align="center")


def bar_chart(c, rows, order, x, y, w, h, ymax, label_size=7):
    ymax = max(
        ymax,
        5
        * math.ceil(
            (
                max(
                    max(rows[g]["edited_percent"], rows[g]["unedited_percent"])
                    for g in order
                )
                + 3
            )
            / 5
        ),
    )
    line(c, x, y, x + w, y, 0.55)
    line(c, x, y, x, y + h, 0.55)
    for v in range(0, int(ymax) + 1, 10):
        yy = y + h * v / ymax
        line(c, x - 1.5, yy, x, yy, 0.4)
        text(c, x - 3, yy - 2.2, str(v), 6.8, align="right")
    c.saveState()
    c.translate(x - 20, y + h / 2)
    c.rotate(90)
    text(c, 0, 0, "Codons in edgeR peak (%)", 7.6, align="center")
    c.restoreState()
    n = len(order)
    step = w / n
    bw = step * 0.30
    gap = step * 0.045
    for i, key in enumerate(order):
        row = rows[key]
        cx = x + (i + 0.5) * step
        e = row["edited_percent"]
        u = row["unedited_percent"]
        rect(c, cx - gap / 2 - bw, y, bw, h * e / ymax, EDITED)
        rect(c, cx + gap / 2, y, bw, h * u / ymax, UNEDIT)
        text(
            c,
            cx,
            y + h * max(e, u) / ymax + 3.2,
            row["significance"],
            6.7,
            EDITED,
            "center",
            bold=True,
        )
        group_label(c, cx, y - 10, key, label_size)
    # Match manuscript's two compact legend keys above the plotting area.
    lx = x + w - 69
    ly = y + h + 3.8
    rect(c, lx, ly, 4.5, 4.5, EDITED)
    text(c, lx + 6.4, ly - 0.3, "Edited", 6.8)
    rect(c, lx + 34, ly, 4.5, 4.5, UNEDIT)
    text(c, lx + 40.4, ly - 0.3, "Unedit", 6.8)
    return ymax


def schematic(c, x, y):
    # Values in this schematic are illustrative; all quantitative panels use fitted results.
    heads = [
        ("Extract matched", "coverage"),
        ("Model focal versus", "flanking counts"),
        ("Ribosome", "occupancy peak"),
    ]
    for cx, (a, b) in zip([x + 34, x + 102, x + 171], heads):
        text(c, cx, y + 107, a, 6.7, align="center")
        text(c, cx, y + 100, b, 6.7, align="center")
    for bx in [x + 75, x + 143]:
        c.setFillColorRGB(0.95, 0.95, 0.95)
        c.roundRect(bx, y + 6, 57, 87, 3, stroke=0, fill=1)
    # Matched mRNA and A-site Ribo-seq tracks.
    for assay, base, col in [
        ("mRNA-seq", y + 57, (0.57, 0.57, 0.57)),
        ("A-site Ribo-seq", y + 10, UNEDIT),
    ]:
        text(c, x + 2, base + 29, assay, 6.3)
        for j in range(37):
            ht = (5 + (j * 7 % 9)) * 0.78
            if assay.startswith("A-site") and j == 18:
                ht = 23
            rect(c, x + 3 + j * 1.52, base, 1.05, ht, EDITED if j == 18 else col)
        line(c, x + 2, base, x + 60, base, 0.4)
        yy = base + 21
        for aa, bb in [(x + 4, x + 27), (x + 34, x + 59)]:
            line(c, aa, yy, bb, yy, 0.4, (0.1, 0.48, 0.72))
            line(c, aa, yy - 1.7, aa, yy + 1.7, 0.4, (0.1, 0.48, 0.72))
            line(c, bb, yy - 1.7, bb, yy + 1.7, 0.4, (0.1, 0.48, 0.72))
        prefix = "m" if assay == "mRNA-seq" else "R"
        for cx, s in [
            (x + 15, prefix + "f"),
            (x + 30, prefix + "c"),
            (x + 47, prefix + "f"),
        ]:
            text(c, cx, yy + 3, s, 6.0, (0.1, 0.48, 0.72), "center", italic=True)
    # Display the four observations rather than imply a test on calculated ratios.
    text(c, x + 103, y + 82, "edgeR QL model", 6.1, align="center")
    text(c, x + 104, y + 68, "Focal", 5.7, align="center")
    text(c, x + 126, y + 68, "Flank", 5.7, align="center")
    for yy, lab, a, b in [(y + 53, "RNA", "mc", "mf"), (y + 37, "RPF", "Rc", "Rf")]:
        text(c, x + 78, yy, lab, 5.7)
        text(c, x + 104, yy, a, 6.5, align="center", italic=True)
        text(c, x + 126, yy, b, 6.5, align="center", italic=True)
    text(c, x + 103, y + 16, "Biological replicates", 5.5, align="center")
    text(c, x + 171, y + 82, "Focal mRC > flank", 5.9, align="center")
    text(c, x + 171, y + 73, "FDR ≤ 0.05", 6.1, align="center")
    vals = [5, 7, 9, 7, 11, 22, 12, 9, 8, 6, 5]
    for j, ht in enumerate(vals):
        rect(c, x + 148 + j * 4.4, y + 33, 2.8, ht, EDITED if j == 5 else UNEDIT)
    line(c, x + 146, y + 33, x + 197, y + 33, 0.4)
    arrow(c, x + 171, y + 27, x + 148, y + 27, (0.1, 0.48, 0.72))
    arrow(c, x + 173, y + 27, x + 195, y + 27, (0.1, 0.48, 0.72))
    text(c, x + 171, y + 17, "±30 codons", 5.8, (0.1, 0.48, 0.72), "center")
    arrow(c, x + 64, y + 56, x + 72, y + 56)
    arrow(c, x + 135, y + 56, x + 141, y + 56)


def save_panel(
    out, name, data, order, size, plotbox, letter, ymax, with_schematic=False
):
    p = out / (name + ".pdf")
    c = canvas.Canvas(str(p), pagesize=size, pageCompression=1)
    c.setTitle(name + " — revised peak-only edgeR analysis")
    c.setAuthor("RNA editing analysis")
    text(c, 2, size[1] - 11, letter, 10, bold=True)
    if with_schematic:
        schematic(c, 12, 5)
    top = bar_chart(c, data, order, *plotbox, ymax, 6.8 if len(order) > 6 else 7)
    c.showPage()
    c.save()
    rows = [dict(panel=name, **data[k]) for k in order]
    write_tsv(out / (name + ".plot_data.tsv"), rows)
    return dict(
        panel=name,
        pdf=p.name,
        groups=order,
        page_points=size,
        plot_points=plotbox,
        ymax=top,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edgeR-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    setup_fonts()
    src = (
        args.edgeR_dir
        / "edited_unedited_summary_codonA/edited_unedited_peaks_by_library.tsv"
    )
    hs = (
        args.edgeR_dir
        / "edited_unedited_hairpin_summary/edited_unedited_hairpin_in_peaks_by_library.long.tsv"
    )
    data = load(src)
    hairpin = load(hs, True)
    assert len(data) == len(hairpin) == 14
    outside = load(
        hs.with_name("edited_unedited_hairpin_out_peaks_by_library.long.tsv"), True
    )
    for key in data:
        for metric in [
            "edited_peak",
            "edited_total",
            "unedited_peak",
            "unedited_total",
        ]:
            assert hairpin[key][metric] + outside[key][metric] == data[key][metric], (
                key,
                metric,
                "hairpin partition mismatch",
            )
    rec = []
    rec.append(
        save_panel(
            out, "F2B", data, MALE, (420, 137), (254, 25, 160, 94), "B", 30, True
        )
    )
    rec.append(
        save_panel(out, "F2B_chart", data, MALE, (199, 137), (32, 25, 160, 94), "B", 30)
    )
    rec.append(
        save_panel(out, "S6A", data, FEMALE, (208, 137), (32, 25, 169, 94), "A", 25)
    )
    rec.append(
        save_panel(out, "F3A", hairpin, ORDER, (238, 137), (32, 25, 200, 94), "A", 30)
    )
    report = {
        "panels": rec,
        "reference_manuscript": "RNA_editing_MC_manuscript-0909-LUJ-0920.pdf pp.48,50",
        "reference_supplement": "RNA_editing_MC_SOM-0909-LUJ-0918.pdf p.8",
        "colors": {"edited": EDITED, "unedited": UNEDIT},
        "rendering": "Direct ReportLab PDF text, paths and rectangles; no images, clipping paths, masks or transparency",
        "fonts": "Embedded Arial subsets from reference figures; unsupported glyphs use metrically compatible Liberation Sans",
        "hairpin_partition_checks": "All 14 groups: in + out equals all tested codons and peaks, separately for edited and unedited controls",
        "statistics": "Original one-sided Fisher greater P, independently recomputed from each 2x2 count table; no multiplicity correction for enrichment stars",
        "input_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [src, hs]
        },
    }
    pass
    print("Publication panels complete:", out, flush=True)


if __name__ == "__main__":
    main()
