#!/usr/bin/env python3
"""Native vector PDF rendering. No images, clipping paths, transparency or masks.

Run with /usr/bin/python3 and runtime/python36 on PYTHONPATH.
Statistical results are read from R outputs; only display geometry is calculated here.
"""

import csv, math, random, statistics, json, argparse, os
from pathlib import Path
from collections import defaultdict
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor, Color

FIG = Path(os.environ.get("FIGURE_ROOT", Path(__file__).resolve().parents[1])).resolve()
GLM = Path(os.environ.get("GLM_RESULT_DIR", FIG / "GLM")).resolve()
OUT = FIG / "figures"
TABLES = FIG / "tables"
OUT.mkdir(exist_ok=True)
TABLES.mkdir(exist_ok=True)
for name, file in [("Sans", "Regular"), ("SansBold", "Bold"), ("SansItalic", "Italic")]:
    pdfmetrics.registerFont(
        TTFont(name, "/usr/share/fonts/liberation/LiberationSans-" + file + ".ttf")
    )
BLACK = "#231815"
EDIT = "#BE1681"
CONTROL = "#279D9F"
MALE = "#4C78A8"
FEMALE = "#E78AC3"


def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(name, rows):
    with (TABLES / name).open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def num(x):
    return float(x)


def new(name, w, h):
    c = canvas.Canvas(
        str(OUT / (name + ".pdf")),
        pagesize=(w, h),
        pageCompression=1,
        initialFontName="Sans",
    )
    c.setTitle(name + " | recalculated gmRC")
    c.setAuthor("RNA editing analysis")
    return c


def text(c, x, y, s, size=7, align="left", font="Sans", color=BLACK, angle=0):
    c.saveState()
    c.setFillColor(HexColor(color))
    c.setFont(font, size)
    c.translate(x, y)
    c.rotate(angle)
    {"left": c.drawString, "center": c.drawCentredString, "right": c.drawRightString}[
        align
    ](0, 0, str(s))
    c.restoreState()


def line(c, x1, y1, x2, y2, color=BLACK, width=0.55, dash=None):
    c.saveState()
    c.setStrokeColor(HexColor(color))
    c.setLineWidth(width)
    if dash:
        c.setDash(*dash)
    c.line(x1, y1, x2, y2)
    c.restoreState()


def poly(c, points, fill=None, stroke=None, width=0.5, close=False):
    p = c.beginPath()
    p.moveTo(*points[0])
    for xy in points[1:]:
        p.lineTo(*xy)
    if close:
        p.close()
    c.saveState()
    c.setLineWidth(width)
    if fill:
        c.setFillColor(HexColor(fill))
    if stroke:
        c.setStrokeColor(HexColor(stroke))
    c.drawPath(p, fill=int(fill is not None), stroke=int(stroke is not None))
    c.restoreState()


def rect(c, x, y, w, h, fill, stroke=None, width=0.5):
    poly(c, [(x, y), (x + w, y), (x + w, y + h), (x, y + h)], fill, stroke, width, True)


def marker(c, x, y, color, shape="Female", size=1.6):
    if shape == "Male":
        poly(
            c,
            [(x, y + size), (x - size, y - size), (x + size, y - size)],
            color,
            None,
            close=True,
        )
    else:
        c.saveState()
        c.setFillColor(HexColor(color))
        c.circle(x, y, size, stroke=0, fill=1)
        c.restoreState()


def axis(c, x, y, w, h, xlim, ylim, xticks, yticks, xlabel=None, ylabel=None, size=6):
    X = lambda a: x + (a - xlim[0]) / (xlim[1] - xlim[0]) * w
    Y = lambda a: y + (a - ylim[0]) / (ylim[1] - ylim[0]) * h
    line(c, x, y, x + w, y)
    line(c, x, y, x, y + h)
    for t in xticks:
        line(c, X(t), y, X(t), y - 2.5)
        text(c, X(t), y - 10, ("%g" % t), size, "center")
    for t in yticks:
        line(c, x - 2.5, Y(t), x, Y(t))
        text(c, x - 5, Y(t) - 2, ("%g" % t), size, "right")
    if xlabel:
        text(c, x + w / 2, y - 24, xlabel, 7, "center")
    if ylabel:
        text(c, x - 29, y + h / 2, ylabel, 7, "center", angle=90)
    return X, Y


def quant(v, p):
    v = sorted(v)
    idx = (len(v) - 1) * p
    a = int(idx)
    b = min(a + 1, len(v) - 1)
    return v[a] + (idx - a) * (v[b] - v[a])


def fmtp(p):
    return "<1e-300" if p < 1e-300 else ("%.2g" % p)


def metagene():
    rows = read(FIG / "metagene/metagene_gmRC_summary.tsv")
    specs = [
        (
            "F2C",
            [
                ("Dmel_Adar5G1__male", "Dmel (Adar5G1) M"),
                ("Dmel_ISO1_WBg__male", "Dmel (ISO1 + WT) M"),
                ("Dsim__male", "Dsim M"),
            ],
            3,
        ),
        (
            "FS7A",
            [
                ("Dmel_ISO1_WBg__female", "Dmel (ISO1 + WT) F"),
                ("Dsec__male", "Dsec M"),
                ("Dere__male", "Dere M"),
                ("Dsim__female", "Dsim F"),
                ("Dsec__female", "Dsec F"),
                ("Dere__female", "Dere F"),
            ],
            3,
        ),
    ]
    for name, units, ncol in specs:
        width = 510
        height = 148 if name == "F2C" else 260
        c = new(name, width, height)
        text(c, 6, height - 13, "C" if name == "F2C" else "A", 11, font="SansBold")
        line(c, 370, height - 10, 385, height - 10, EDIT, 1.4)
        text(c, 389, height - 12, "Edited", 6.5)
        line(c, 433, height - 10, 448, height - 10, CONTROL, 1.4)
        text(c, 452, height - 12, "Unedited", 6.5)
        selected = []
        pw = (width - 51) / ncol
        ph = 76
        for i, (unit, title) in enumerate(units):
            d = sorted(
                [r for r in rows if r["unit_id"] == unit],
                key=lambda r: num(r["Position"]),
            )
            selected.extend(d)
            x = 38 + (i % ncol) * pw
            y = 44 + (len(units) // ncol - 1 - i // ncol) * 112
            ymax = max(
                2.2,
                math.ceil(max(num(r[k]) for r in d for k in ["edited", "control"]) * 5)
                / 5,
            )
            X, Y = axis(
                c,
                x,
                y,
                pw - 15,
                ph,
                (-100, 100),
                (0.2, ymax),
                [-100, -50, 0, 50, 100],
                [0.2, 0.6, 1, 1.4, 1.8, 2.2] if ymax == 2.2 else [0.2, 1, 2, ymax],
                ylabel="gmRC" if i % ncol == 0 else None,
                size=5.7,
            )
            line(c, X(0), Y(0.2), X(0), Y(ymax), "#888888", 0.5, [2, 2])
            for key, color in [("control", CONTROL), ("edited", EDIT)]:
                poly(
                    c,
                    [(X(num(r["Position"])), Y(num(r[key]))) for r in d],
                    stroke=color,
                    width=0.65,
                )
            text(c, x + (pw - 15) / 2, y + ph + 8, title, 7, "center")
            center = next(r for r in d if num(r["Position"]) == 0)
            text(
                c,
                x + pw - 18,
                y + ph - 10,
                "n = " + center["edited_n_before_trim"],
                6,
                "right",
            )
        text(c, width / 2, 12, "Distance from focal codon (codons)", 7.5, "center")
        c.save()
        write(name + "_profile_data.tsv", selected)


def f3c():
    # F3C is a four-species WT comparison; exclude Adar5G1 from plot and export.
    d = [r for r in read(FIG / "F3C/F3C_plot_data.tsv") if r["plot_group"] != "adar5G1"]
    aa = sorted(
        set(r["AA"] for r in d) - {"edited", "AA_mean"},
        key=lambda a: -statistics.median(num(r["gmRC"]) for r in d if r["AA"] == a),
    )
    order = ["edited", "AA_mean"] + aa
    col = {"Dmel": "#F2B37B", "Dsim": "#D6456C", "Dsec": "#7BBCD5", "Dere": "#6E7CB8"}
    c = new("F3C", 510, 210)
    text(c, 6, 195, "C", 11, font="SansBold")
    values = [math.log(num(r["gmRC"]), 2) for r in d]
    lo = math.floor(min(values) * 2) / 2 - 0.1
    hi = math.ceil(max(values) * 2) / 2 + 0.28
    X, Y = axis(
        c,
        39,
        52,
        454,
        130,
        (0.5, len(order) + 0.5),
        (lo, hi),
        [],
        [v / 2 for v in range(math.ceil(lo * 2), math.floor(hi * 2) + 1)],
        ylabel="log2(mean gmRC)",
    )
    rng = random.Random(18907675)
    coords = []
    for j, a in enumerate(order, 1):
        vals = [
            math.log(num(r["gmRC"]), 2)
            for r in d
            if r["AA"] == a and r["plot_group"] != "adar5G1"
        ]
        q1, med, q3 = quant(vals, 0.25), quant(vals, 0.5), quant(vals, 0.75)
        iqr = q3 - q1
        low = min(v for v in vals if v >= q1 - 1.5 * iqr)
        high = max(v for v in vals if v <= q3 + 1.5 * iqr)
        line(c, X(j), Y(low), X(j), Y(high))
        line(c, X(j - 0.14), Y(low), X(j + 0.14), Y(low))
        line(c, X(j - 0.14), Y(high), X(j + 0.14), Y(high))
        rect(
            c,
            X(j - 0.32),
            Y(q1),
            X(j + 0.32) - X(j - 0.32),
            Y(q3) - Y(q1),
            "#D9D9D9",
            BLACK,
        )
        line(c, X(j - 0.32), Y(med), X(j + 0.32), Y(med), BLACK, 0.8)
        for r in [r for r in d if r["AA"] == a]:
            jitter = rng.uniform(-0.12, 0.12)
            val = math.log(num(r["gmRC"]), 2)
            marker(c, X(j + jitter), Y(val), col[r["plot_group"]], r["sex"], 1.5)
            coords.append(dict(r, category_order=j, log2_mean_gmRC=val, jitter=jitter))
        text(
            c,
            X(j),
            44,
            {"edited": "Edited", "AA_mean": "AA mean"}.get(a, a),
            6,
            "right",
            angle=45 if a == "AA_mean" else 0,
        )
    p = num(read(FIG / "F3C/AA_gmRC_v10_edited_vs_AAmean_pairedWilcox.tsv")[0]["p_adj"])
    by = hi - 0.10
    line(c, X(1), Y(by) - 3, X(1), Y(by))
    line(c, X(1), Y(by), X(2), Y(by))
    line(c, X(2), Y(by), X(2), Y(by) - 3)
    text(c, X(1) - 3, Y(by) + 4, "P = " + fmtp(p), 6)
    for i, (s, color) in enumerate(col.items()):
        marker(c, 264 + i * 45, 199, color, "Female", 1.7)
        text(c, 268 + i * 45, 197, "Adar5G1" if s == "adar5G1" else s, 6)
    text(c, 350, 187, "Female: circle   Male: triangle", 6)
    c.save()
    write("F3C_plot_data.tsv", coords)


def blend(a, b, t):
    a = HexColor(a)
    b = HexColor(b)
    return (
        Color(
            a.red * (1 - t) + b.red * t,
            a.green * (1 - t) + b.green * t,
            a.blue * (1 - t) + b.blue * t,
        )
        .hexval()
        .replace("0x", "#")
    )


def heatcolor(v):
    v = max(0, min(2.5, v))
    return (
        blend("#2166AC", "#FFFFFF", v / 1.25)
        if v <= 1.25
        else blend("#FFFFFF", "#B2182B", (v - 1.25) / 1.25)
    )


def f3d():
    d = read(FIG / "AA/AA_three_group_four_species_pool_summary.tsv")
    stats = read(FIG / "AA/AA_three_group_codon_pooled_simple_group_effect_by_AA.tsv")
    order = [
        a
        for a in [
            "Ala",
            "Arg",
            "Asn",
            "Asp",
            "Cys",
            "Gln",
            "Glu",
            "Gly",
            "His",
            "Ile",
            "Leu",
            "Lys",
            "Met",
            "Phe",
            "Pro",
            "Ser",
            "Thr",
            "Trp",
            "Tyr",
            "Val",
        ]
        if any(r["AA"] == a for r in d)
    ]
    c = new("F3D", 310, 194)
    text(c, 6, 179, "D", 11, font="SansBold")
    x = 62
    y = 54
    w = 232 / len(order)
    h = 32
    groups = ["Unedited_preferred", "Unedited_unpreferred", "Edited"]
    selected = []
    for i, a in enumerate(order):
        s = next(r for r in stats if r["AA"] == a)
        positive = s["edited_higher_than_any_unedited"] == "TRUE"
        star = s["significance"]
        for j, g in enumerate(groups):
            r = next(r for r in d if r["AA"] == a and r["group"] == g)
            v = num(r["mean_gmRC"])
            rect(c, x + i * w, y + (2 - j) * h, w, h, heatcolor(v), "#FFFFFF", 0.5)
            text(
                c,
                x + (i + 0.5) * w,
                y + (2 - j) * h + 14,
                "%.2f" % v,
                5.8,
                "center",
                color="#FFFFFF" if v >= 1.75 or v <= 0.25 else BLACK,
            )
            selected.append(
                dict(
                    r,
                    p_adj_BH=s["p_adj_BH"],
                    significance=star,
                    edited_higher_than_any_unedited=positive,
                )
            )
        text(
            c,
            x + (i + 0.5) * w,
            y - 7,
            a,
            6,
            "right",
            color="#BE2235" if positive else BLACK,
            angle=45,
        )
        if star != "ns":
            text(
                c,
                x + (i + 0.5) * w,
                y - 31,
                star,
                6,
                "center",
                color="#BE2235" if positive else BLACK,
            )
    for j, labels in enumerate(
        [["Unedited", "preferred"], ["Unedited", "unpreferred"], ["Edited"]]
    ):
        for k, label in enumerate(labels):
            text(c, x - 5, y + (2 - j) * h + 18 - k * 8, label, 6, "right")
    for i in range(80):
        rect(c, 185 + i * 1.35, 171, 1.4, 5, heatcolor(i / 79 * 2.5))
    text(c, 177, 170, "gmRC", 6, "right")
    for v in [0, 0.5, 1, 1.25, 1.5, 2, 2.5]:
        text(c, 185 + v / 2.5 * 108, 180, "%g" % v, 5, "center")
    c.save()
    write("F3D_heatmap_data.tsv", selected)


def glm():
    d = read(
        GLM
        / "WT_sex_BH_plots/latest_all_variables_BH/all_WT_fixed_effects_BH_all_variables.tsv"
    )
    terms = [
        ("edit_statusedited", "Editing status"),
        ("editing_level_within_editSD", "Editing level (+1 SD; edited)"),
        ("category_beforeoptimal", "Focal codon preference"),
        ("phyloP_z", "Codon phyloP"),
        ("secondary_structureinside", "Stable hairpin"),
        ("logTE_z", "Translation efficiency"),
        ("transcript_length_z", "Transcript length"),
        ("upstream_categoryoptimal", "Upstream codon preference"),
        ("downstream_categoryoptimal", "Downstream codon preference"),
    ]
    for name, model in [("F3E", "GroupMean_Full"), ("FS11", "GroupMean_Matched")]:
        shown = (
            terms
            if name == "F3E"
            else [x for x in terms if x[0] != "editing_level_within_editSD"]
        )
        rows = [r for r in d if r["model_type"] == model]
        c = new(name, 440, 232)
        text(c, 6, 218, "E" if name == "F3E" else "", 11, font="SansBold")
        estimates = [
            num(r["estimate"])
            for r in rows
            if r["term"] in dict(shown) or r["term"] == "edit_status_matchededited"
        ]
        xmin = min(-0.05, math.floor(min(estimates) / 0.05) * 0.05)
        xmax = max(0.20, math.ceil((max(estimates) + 0.014) / 0.05) * 0.05)
        X, Y = axis(
            c,
            183,
            42,
            235,
            166,
            (xmin, xmax),
            (-0.5, len(shown) - 0.5),
            [i * 0.05 for i in range(round(xmin / 0.05), round(xmax / 0.05) + 1)],
            [],
            xlabel="Fixed-effect coefficient for ln(1 + gmRC)",
        )
        line(c, X(0), 42, X(0), 208, "#B5B5B5", 0.45)
        plotted = []
        for j, (term, label) in enumerate(shown):
            yc = Y(len(shown) - 1 - j)
            text(c, 177, yc - 2, label, 7, "right")
            for sex, color, off in [("Male", MALE, 3), ("Female", FEMALE, -3)]:
                match = [
                    r
                    for r in rows
                    if r["sex_label"] == sex
                    and r["term"]
                    in (
                        [term, "edit_status_matchededited"]
                        if term == "edit_statusedited"
                        else [term]
                    )
                ]
                assert len(match) == 1, (name, sex, term)
                r = match[0]
                v = num(r["estimate"])
                rect(c, min(X(0), X(v)), yc + off - 2.3, abs(X(v) - X(0)), 4.6, color)
                if num(r["p_adjust_bh_all"]) < 0.05:
                    text(
                        c,
                        X(v) + (3 if v >= 0 else -3),
                        yc + off - 2,
                        "*",
                        8,
                        "left" if v >= 0 else "right",
                    )
                plotted.append(dict(r, display_variable=label))
        for i, (s, color) in enumerate([("Male", MALE), ("Female", FEMALE)]):
            rect(c, 342, 80 - i * 11, 6, 6, color)
            text(c, 351, 80 - i * 11, s, 6)
        c.save()
        write(name + "_fixed_effects.tsv", plotted)




if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "stage", choices=["metagene", "f3c", "f3d", "glm", "all"]
    )
    a = p.parse_args()
    for name, fn in [
        ("metagene", metagene),
        ("f3c", f3c),
        ("f3d", f3d),
        ("glm", glm),
    ]:
        if a.stage in [name, "all"]:
            fn()
            print("Rendered", name)
