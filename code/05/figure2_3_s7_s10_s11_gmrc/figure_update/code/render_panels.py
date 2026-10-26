#!/usr/bin/env python3
"""Draw replacement panels at their final physical size.

The supplied manuscript PDFs use Arial at 6 pt for ticks, 7 pt for legends,
8 pt for strip titles and 9--10 pt for panel letters.  This renderer keeps
those sizes and draws directly in the final panel coordinate system.  The
composite builder therefore only translates panels; it never scales text.
"""

import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

import render_figures as rf
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

FIG = Path(__file__).resolve().parents[1]
OUT = FIG / "panels"
OUT.mkdir(exist_ok=True)
rf.OUT = OUT

PAGE_W = 467.717
FONT_DIR = FIG / "reference_fonts"
for prefix in ("FS7", "F3", "FS11"):
    for style in ("Regular", "Bold", "Italic"):
        path = FONT_DIR / (prefix + "-Arial-" + style + ".ttf")
        if path.exists():
            pdfmetrics.registerFont(TTFont(prefix + "Arial" + style, str(path)))
ACTIVE_REG, ACTIVE_BOLD, ACTIVE_ITALIC = "Sans", "SansBold", "SansItalic"


def use_font_family(prefix=None):
    global ACTIVE_REG, ACTIVE_BOLD, ACTIVE_ITALIC
    if prefix is None:
        ACTIVE_REG, ACTIVE_BOLD, ACTIVE_ITALIC = "Sans", "SansBold", "SansItalic"
    else:
        ACTIVE_REG = prefix + "ArialRegular"
        ACTIVE_BOLD = prefix + "ArialBold"
        ACTIVE_ITALIC = (
            prefix + "ArialItalic"
            if (FONT_DIR / (prefix + "-Arial-Italic.ttf")).exists()
            else ACTIVE_REG
        )


def label(c, x, y, value, size=7, align="left", font=None, color=rf.BLACK, angle=0):
    if font is None or font == "Sans":
        font = ACTIVE_REG
    elif font == "SansBold":
        font = ACTIVE_BOLD
    elif font == "SansItalic":
        font = ACTIVE_ITALIC
    rf.text(c, x, y, value, size, align, font, color, angle)


def axes(
    c,
    x,
    y,
    w,
    h,
    xlim,
    ylim,
    xticks,
    yticks,
    ylabel=None,
    tick_size=6,
    label_size=8,
    x_tick_offset=-10.0,
    ylabel_gap=16.5,
):
    X = lambda value: x + (value - xlim[0]) / (xlim[1] - xlim[0]) * w
    Y = lambda value: y + (value - ylim[0]) / (ylim[1] - ylim[0]) * h
    rf.line(c, x, y, x + w, y, width=0.55)
    rf.line(c, x, y, x, y + h, width=0.55)
    for value in xticks:
        rf.line(c, X(value), y, X(value), y - 2.2, width=0.5)
        label(c, X(value), y + x_tick_offset, "%g" % value, tick_size, "center")
    for value in yticks:
        rf.line(c, x - 2.2, Y(value), x, Y(value), width=0.5)
        label(c, x - 3.8, Y(value) - 2.0, "%g" % value, tick_size, "right")
    if ylabel:
        label(c, x - ylabel_gap, y + h / 2.0, ylabel, label_size, "center", angle=90)
    return X, Y


def panel_letter(c, letter, h):
    label(c, 5.77, h - 8.2, letter, 10, font="SansBold")


def line_legend(c, h, centered=False):
    if centered:
        rf.line(c, 195, h - 8, 208, h - 8, rf.EDIT, 1.3)
        label(c, 211.786, h - 10.679, "Edited", 7)
        rf.line(c, 235.696, h - 8.309, 244.3, h - 8.309, rf.CONTROL, 1.3)
        label(c, 249.036, h - 10.727, "Unedit", 7)
    else:
        rf.line(c, 367, h - 8.5, 381, h - 8.5, rf.EDIT, 1.3)
        label(c, 384.9, h - 10.5, "Edited", 7)
        rf.line(c, 411, h - 8.5, 425, h - 8.5, rf.CONTROL, 1.3)
        label(c, 428.5, h - 10.5, "Unedited", 7)


def render_metagene(name, letter, units, h, two_rows=False):
    rows = rf.read(FIG / "metagene/metagene_gmRC_summary.tsv")
    c = rf.new(name, PAGE_W, h)
    if two_rows:
        label(c, 5.996, h - 10.679, letter, 10, font="SansBold")
    else:
        panel_letter(c, letter, h)
    line_legend(c, h, two_rows)
    x_positions = [24.9249, 174.8078, 326.1921] if two_rows else [26.5, 179.0, 329.7]
    if two_rows:
        y_positions = [99.5, 25.5]
        plot_h = 48.33
    else:
        y_positions = [24.0]
        plot_h = 52.6
    plot_w = 137.232 if two_rows else 130.0
    for index, (unit, title) in enumerate(units):
        data = sorted(
            (r for r in rows if r["unit_id"] == unit),
            key=lambda r: rf.num(r["Position"]),
        )
        col, row = index % 3, index // 3
        ymax = max(
            2.2,
            math.ceil(
                max(rf.num(r[key]) for r in data for key in ("edited", "control")) * 5
            )
            / 5.0,
        )
        yticks = [0.2, 0.6, 1, 1.4, 1.8, 2.2] if ymax <= 2.2 else [0.2, 1, 2, ymax]
        X, Y = axes(
            c,
            x_positions[col],
            y_positions[row],
            plot_w,
            plot_h,
            (-100, 100),
            (0.2, ymax),
            [-100, -50, 0, 50, 100],
            yticks,
            "gmRC" if col == 0 else None,
            tick_size=6.5 if two_rows else 6,
            ylabel_gap=13.1 if two_rows else 16.5,
        )
        rf.line(c, X(0), Y(0.2), X(0), Y(ymax), "#888888", 0.5, [2, 2])
        for key, color in (("control", rf.CONTROL), ("edited", rf.EDIT)):
            rf.poly(
                c,
                [(X(rf.num(r["Position"])), Y(rf.num(r[key]))) for r in data],
                stroke=color,
                width=0.65,
            )
        title_y = [151.754, 77.046][row] if two_rows else y_positions[row] + plot_h + 8
        label(c, x_positions[col] + plot_w / 2, title_y, title, 8, "center")
        center = next(r for r in data if rf.num(r["Position"]) == 0)
        n_y = y_positions[row] + plot_h - (5.5 if two_rows else 9)
        label(
            c,
            x_positions[col] + plot_w - 4,
            n_y,
            "n = " + center["edited_n_before_trim"],
            7,
            "right",
        )
    label(
        c,
        233.858 if two_rows else 250.0,
        4.533 if two_rows else 3.1,
        "Distance from focal codon (codons)",
        8,
        "center",
    )
    c.save()


def track_groups(path):
    grouped = defaultdict(list)
    for row in rf.read(path):
        grouped[row["group_id"]].append(row)
    return grouped


def draw_track(c, rows, x, y, w, h, ylabel=False, title_y=None):
    rows = sorted(rows, key=lambda r: int(r["codonpos"]))
    lo, hi = min(int(r["codonpos"]) for r in rows), max(
        int(r["codonpos"]) for r in rows
    )
    # All jef windows are 560--760; retaining explicit endpoints keeps old axes fixed.
    lo, hi = (560, 760) if lo >= 500 else (164, 364)
    xticks = [600, 650, 700, 750] if lo == 560 else [200, 250, 300, 350]
    X, Y = axes(
        c,
        x,
        y,
        w,
        h,
        (lo, hi),
        (0, 10),
        xticks,
        [0, 2, 4, 6, 8, 10],
        "gmRC" if ylabel else None,
        tick_size=6.5,
        x_tick_offset=-8.0,
        ylabel_gap=13.1,
    )
    bw = max(0.42, min(0.72, w / (hi - lo) * 0.86))
    focal = None
    for row in rows:
        pos = int(row["codonpos"])
        value = rf.num(row["gmRC"])
        edited = row["is_focal_edited_codon"] == "1"
        if edited:
            focal = (pos, value)
        if value > 0:
            rf.rect(
                c,
                X(pos) - bw / 2,
                Y(0),
                bw,
                Y(min(10, value)) - Y(0),
                rf.EDIT if edited else rf.CONTROL,
            )
    if focal:
        pos, value = focal
        anchor = Y(min(10, value))
        lx = min(x + w - 3, X(pos) + 30)
        ly = min(y + h - 5, max(y + h * 0.63, anchor + 8))
        rf.line(c, X(pos) + bw / 2, anchor, lx - 3, ly - 2, rf.EDIT, 0.55)
        label(c, lx, ly, "%.2f" % value, 7, "right", color=rf.EDIT)
    ty = y + h + 8 if title_y is None else title_y
    title = rows[0]["title"]
    if ACTIVE_REG.startswith("FS7Arial") and title.endswith("(jef)"):
        base = title[:-5] + "("
        parts = [(base, ACTIVE_REG), ("jef", ACTIVE_ITALIC), (")", ACTIVE_REG)]
        widths = [pdfmetrics.stringWidth(value, font, 8) for value, font in parts]
        cursor = x + w / 2 - sum(widths) / 2
        for (value, font), width in zip(parts, widths):
            label(c, cursor, ty, value, 8, font=font)
            cursor += width
    else:
        label(c, x + w / 2, ty, title, 8, "center")


def bar_legend(c, h, centered=False):
    if centered:
        rf.line(c, 198.447, h - 7.449, 207.05, h - 7.449, rf.EDIT, 1.3)
        label(c, 211.787, h - 8.868, "Edited", 7)
        rf.line(c, 235.697, h - 7.498, 244.30, h - 7.498, rf.CONTROL, 1.3)
        label(c, 249.037, h - 8.987, "Unedit", 7)
    else:
        rf.rect(c, 383, h - 11.8, 7, 7, rf.EDIT)
        label(c, 395.2, h - 10.8, "Edited", 7)
        rf.rect(c, 417, h - 11.8, 7, 7, rf.CONTROL)
        label(c, 428.9, h - 10.8, "Unedited", 7)


def render_tracks(name, letter, table, order, h, two_rows=False):
    grouped = track_groups(FIG / "tables" / table)
    c = rf.new(name, PAGE_W, h)
    (
        label(c, 5.996, h - 10.098, letter, 10, font="SansBold")
        if two_rows
        else label(c, 5.77, h - 7.3, letter, 10, font="SansBold")
    )
    bar_legend(c, h, two_rows)
    xs = [26.0, 178.4, 328.9] if not two_rows else [24.9249, 174.8078, 326.1921]
    ys = [35.2] if not two_rows else [110.214, 23.356]
    ph = 47.6 if not two_rows else 63.372
    pw = 128.5 if not two_rows else 137.232
    fs7_titles = {
        "Dmel__female": "Dmel (iso-1 + WT) F (jef)",
        "Dsec__male": "Dsec M (jef)",
        "Dere__male": "Dere M (jef)",
        "Dsim__female": "Dsim F (jef)",
        "Dsec__female": "Dsec F (jef)",
        "Dere__female": "Dere F (jef)",
    }
    for i, group in enumerate(order):
        current = [
            dict(row, title=fs7_titles.get(group, row["title"]))
            for row in grouped[group]
        ]
        draw_track(
            c,
            current,
            xs[i % 3],
            ys[i // 3],
            pw,
            ph,
            i % 3 == 0,
            89.34 if not two_rows else [177.915, 89.752][i // 3],
        )
    label(
        c,
        233.858,
        15.7 if not two_rows else 3.939,
        "Transcript codon position",
        8,
        "center",
    )
    c.save()


def render_f3c():
    h = 172.12
    data = [
        r
        for r in rf.read(FIG / "tables/F3C_plot_data.tsv")
        if r["plot_group"] != "adar5G1"
    ]
    order = [a for _, a in sorted({(int(r["category_order"]), r["AA"]) for r in data})]
    colors = {
        "Dmel": "#F2B37B",
        "Dsim": "#D6456C",
        "Dsec": "#7BBCD5",
        "Dere": "#6E7CB8",
    }
    c = rf.new("F3C", PAGE_W, h)
    panel_letter(c, "C", h)
    values = [rf.num(r["log2_mean_gmRC"]) for r in data]
    lo = math.floor(min(values) * 2) / 2 - 0.1
    hi = math.ceil(max(values) * 2) / 2 + 0.25
    X, Y = axes(
        c,
        42.0,
        31.0,
        400.0,
        108.0,
        (0.5, len(order) + 0.5),
        (lo, hi),
        [],
        [v / 2.0 for v in range(math.ceil(lo * 2), math.floor(hi * 2) + 1)],
        "log2(mean gmRC)",
        tick_size=7,
        label_size=8,
    )
    rng = random.Random(18907675)
    for j, aa in enumerate(order, 1):
        vals = [rf.num(r["log2_mean_gmRC"]) for r in data if r["AA"] == aa]
        q1, med, q3 = rf.quant(vals, 0.25), rf.quant(vals, 0.5), rf.quant(vals, 0.75)
        iqr = q3 - q1
        low = min(v for v in vals if v >= q1 - 1.5 * iqr)
        high = max(v for v in vals if v <= q3 + 1.5 * iqr)
        rf.line(c, X(j), Y(low), X(j), Y(high))
        rf.line(c, X(j - 0.14), Y(low), X(j + 0.14), Y(low))
        rf.line(c, X(j - 0.14), Y(high), X(j + 0.14), Y(high))
        rf.rect(
            c,
            X(j - 0.31),
            Y(q1),
            X(j + 0.31) - X(j - 0.31),
            Y(q3) - Y(q1),
            "#D9D9D9",
            rf.BLACK,
        )
        rf.line(c, X(j - 0.31), Y(med), X(j + 0.31), Y(med), width=0.8)
        for row in (r for r in data if r["AA"] == aa):
            rf.marker(
                c,
                X(j + rng.uniform(-0.12, 0.12)),
                Y(rf.num(row["log2_mean_gmRC"])),
                colors[row["plot_group"]],
                row["sex"],
                1.55,
            )
        label(
            c,
            X(j),
            22.0,
            {"edited": "Edited", "AA_mean": "AA_mean"}.get(aa, aa),
            7,
            "right",
            angle=45 if aa == "AA_mean" else 0,
        )
    p = rf.num(
        rf.read(FIG / "F3C/AA_gmRC_v10_edited_vs_AAmean_pairedWilcox.tsv")[0]["p_adj"]
    )
    by = hi - 0.10
    rf.line(c, X(1), Y(by) - 3, X(1), Y(by))
    rf.line(c, X(1), Y(by), X(2), Y(by))
    rf.line(c, X(2), Y(by), X(2), Y(by) - 3)
    label(c, X(1) - 2, Y(by) + 4, "P = " + rf.fmtp(p), 7)
    label(c, 389, h - 18, "Sex", 7)
    rf.marker(c, 394, h - 27, colors["Dmel"], "Female", 1.7)
    label(c, 400, h - 29, "Female", 7)
    rf.marker(c, 394, h - 38, colors["Dmel"], "Male", 1.7)
    label(c, 400, h - 40, "Male", 7)
    label(c, 423, h - 18, "Species", 7)
    for i, sp in enumerate(("Dmel", "Dsim", "Dsec", "Dere")):
        rf.marker(c, 426, h - 29 - i * 9, colors[sp], "Female", 1.7)
        label(c, 432, h - 31 - i * 9, sp, 7)
    c.save()


def heatcolor(value):
    value = max(0.0, min(2.5, value))
    if value <= 1.25:
        return rf.blend("#2166AC", "#FFFFFF", value / 1.25)
    return rf.blend("#FFFFFF", "#B2182B", (value - 1.25) / 1.25)


def render_f3d():
    h, wpage = 240.009, 245.87
    data = rf.read(FIG / "AA/AA_three_group_four_species_pool_summary.tsv")
    stats = rf.read(
        FIG / "AA/AA_three_group_codon_pooled_simple_group_effect_by_AA.tsv"
    )
    order = [
        aa
        for aa in (
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
        )
        if any(r["AA"] == aa for r in data)
    ]
    c = rf.new("F3D", wpage, h)
    label(c, 7.96, h - 14.3, "D", 10, font="SansBold")
    x, y, total_w, row_h = 50.0, 115.0, 182.0, 34.0
    cell_w = total_w / len(order)
    groups = ["Unedited_preferred", "Unedited_unpreferred", "Edited"]
    for i, aa in enumerate(order):
        stat = next(r for r in stats if r["AA"] == aa)
        positive = stat["edited_higher_than_any_unedited"] == "TRUE"
        for j, group in enumerate(groups):
            row = next(r for r in data if r["AA"] == aa and r["group"] == group)
            value = rf.num(row["mean_gmRC"])
            rf.rect(
                c,
                x + i * cell_w,
                y + (2 - j) * row_h,
                cell_w,
                row_h,
                heatcolor(value),
                "#FFFFFF",
                0.45,
            )
            label(
                c,
                x + (i + 0.5) * cell_w,
                y + (2 - j) * row_h + 13,
                "%.2f" % value,
                5.7,
                "center",
                color="#FFFFFF" if value <= 0.25 or value >= 2.15 else rf.BLACK,
            )
        label(
            c,
            x + (i + 0.5) * cell_w,
            y - 8,
            aa,
            6,
            "right",
            color="#BE2235" if positive else rf.BLACK,
            angle=45,
        )
        if stat["significance"] != "ns":
            label(
                c,
                x + (i + 0.5) * cell_w,
                y - 31,
                stat["significance"],
                6,
                "center",
                color="#BE2235" if positive else rf.BLACK,
            )
    labels = [("Unedited\npreferred", 2), ("Unedited\nunpreferred", 1), ("Edited", 0)]
    for value, row in labels:
        for k, part in enumerate(value.split("\n")):
            label(c, x - 5, y + row * row_h + 18 - k * 8, part, 6, "right")
    bx, by, bw = 141.0, 221.0, 90.0
    label(c, bx - 6, by - 1, "gmRC", 7, "right")
    for i in range(100):
        rf.rect(c, bx + i * bw / 100, by, bw / 100 + 0.05, 5, heatcolor(i / 99.0 * 2.5))
    for value in (0, 0.5, 1, 1.25, 1.5, 2, 2.5):
        label(c, bx + value / 2.5 * bw, by + 9, "%g" % value, 5.5, "center")
    c.save()


TERMS = [
    ("Editing status", "edit_statusedited"),
    ("Editing level (+1 SD; edited)", "editing_level_within_editSD"),
    ("Focal codon preference", "category_beforeoptimal"),
    ("Codon phyloP", "phyloP_z"),
    ("Stable hairpin", "secondary_structureinside"),
    ("Translation efficiency", "logTE_z"),
    ("Transcript length", "transcript_length_z"),
    ("Upstream codon preference", "upstream_categoryoptimal"),
    ("Downstream codon preference", "downstream_categoryoptimal"),
]


def glm_panel(name, w, h, compact=False):
    """Reproduce the supplied Arial GLM layout; replace only bar values/stars."""
    rows = rf.read(FIG / "tables" / (name + "_fixed_effects.tsv"))
    shown = (
        TERMS
        if name == "F3E"
        else [x for x in TERMS if x[1] != "editing_level_within_editSD"]
    )
    c = rf.new(name, w, h)
    if name == "F3E":
        label(c, 0.001, 229.453, "E", 10, font="SansBold")
        label_x = 98.208
        label_y = [
            221.6846,
            207.8359,
            193.9902,
            180.1445,
            166.2949,
            152.4463,
            138.5996,
            124.7549,
            110.9053,
        ]
        axis_x0, axis_y, axis_x1, axis_top = 102.328, 105.8581, 207.344, 230.4882
        zero_x, step = 123.3312, 21.0033
        ticks = (-0.05, 0, 0.05, 0.10, 0.15, 0.20)
        tick_size, bar_h = 6, 4.506
        male_bottom = lambda ly: ly + 2.442
        female_bottom = lambda ly: ly - 3.191
        star_y = {"Male": lambda ly: ly + 1.005, "Female": lambda ly: ly - 4.628}
        star_size = 8
    else:
        label_x = 105.63
        label_y = [
            170.4064,
            151.4626,
            132.5202,
            113.5769,
            94.6345,
            75.6911,
            56.7487,
            37.8054,
        ]
        axis_x0, axis_y, axis_x1, axis_top = 111.7822, 29.9278, 446.7172, 182.3228
        zero_x, step = 182.1809, 70.3986
        ticks = (-0.05, 0, 0.05, 0.10, 0.15)
        tick_size, bar_h = 6.5, 7.16
        male_bottom = lambda ly: ly + 2.203
        female_bottom = lambda ly: ly - 5.426
        star_y = {"Male": lambda ly: ly + 2.528, "Female": lambda ly: ly - 6.101}
        star_size = 8
    X = lambda value: zero_x + value / 0.05 * step
    rf.line(c, axis_x0, axis_y, axis_x1, axis_y, width=0.5)
    rf.line(c, axis_x0, axis_y, axis_x0, axis_top, width=0.5)
    for tick in ticks:
        rf.line(
            c,
            X(tick),
            axis_y,
            X(tick),
            axis_y - (1.878 if name == "F3E" else 2.26),
            width=0.4,
        )
        label(
            c,
            X(tick),
            axis_y - (8.2614 if name == "F3E" else 8.9193),
            "%.2f" % tick,
            tick_size,
            "center",
        )
    for index, (display, term) in enumerate(shown):
        ly = label_y[index]
        # Added editing-level text needs glyphs absent from legacy PDF subsets.
        if term == "editing_level_within_editSD":
            rf.text(c, label_x, ly, display, 7, "right", font="Sans")
        else:
            label(c, label_x, ly, display, 7, "right")
        rf.line(
            c,
            axis_x0 - (1.878 if name == "F3E" else 4.093),
            ly + (1.878 if name == "F3E" else 2.393),
            axis_x0,
            ly + (1.878 if name == "F3E" else 2.393),
            width=0.45,
        )
        for sex, color in (("Male", rf.MALE), ("Female", rf.FEMALE)):
            match = [
                r
                for r in rows
                if r["sex_label"] == sex and r["display_variable"] == display
            ]
            if not match:
                match = [r for r in rows if r["sex_label"] == sex and r["term"] == term]
            row = match[0]
            value = rf.num(row["estimate"])
            yb = male_bottom(ly) if sex == "Male" else female_bottom(ly)
            rf.rect(c, min(X(0), X(value)), yb, abs(X(value) - X(0)), bar_h, color)
            if rf.num(row["p_adjust_bh_all"]) < 0.05:
                label(
                    c,
                    X(value) + (3 if value >= 0 else -3),
                    star_y[sex](ly),
                    "*",
                    star_size,
                    "left" if value >= 0 else "right",
                    font="SansBold",
                )
    if name == "F3E":
        label(c, 76.92, 84.7061, "Fixed-effect coefficient for ln(1 + gmRC)", 8)
        rf.rect(c, 179.171, 126.426, 7.593, 5.695, rf.MALE)
        label(c, 191.034, 127.375, "Male", 7)
        rf.rect(c, 179.171, 115.512, 7.593, 5.695, rf.FEMALE)
        label(c, 191.034, 116.4609, "Female", 7)
    else:
        label(c, 200.1934, 8.0651, "Fixed-effect coefficient for ln(1 + gmRC)", 8)
        label(c, 408.2896, 70.3376, "Sex", 7)
        rf.rect(c, 408.681, 59.815, 9.547, 7.16, rf.MALE)
        label(c, 421.5562, 61.3151, "Male", 7)
        rf.rect(c, 408.681, 50.316, 9.547, 7.16, rf.FEMALE)
        label(c, 421.5562, 51.5886, "Female", 7)
    c.save()




if __name__ == "__main__":
    use_font_family(None)
    render_metagene(
        "F2C",
        "C",
        [
            ("Dmel_Adar5G1__male", "Dmel (Adar5G1) M"),
            ("Dmel_ISO1_WBg__male", "Dmel (ISO1 + WT) M"),
            ("Dsim__male", "Dsim M"),
        ],
        104.43,
    )
    render_tracks(
        "F2D",
        "D",
        "F2D_profile_data.tsv",
        ["Dmel_Adar5G1__male", "Dmel__male", "Dsim__male"],
        112.499,
    )
    use_font_family("FS7")
    render_metagene(
        "FS7A",
        "A",
        [
            ("Dmel_ISO1_WBg__female", "Dmel (iso-1 + WT) F"),
            ("Dsec__male", "Dsec M"),
            ("Dere__male", "Dere M"),
            ("Dsim__female", "Dsim F"),
            ("Dsec__female", "Dsec F"),
            ("Dere__female", "Dere F"),
        ],
        174.23,
        True,
    )
    render_tracks(
        "FS7B",
        "B",
        "FS7B_profile_data.tsv",
        [
            "Dmel__female",
            "Dsec__male",
            "Dere__male",
            "Dsim__female",
            "Dsec__female",
            "Dere__female",
        ],
        202.77,
        True,
    )
    use_font_family("F3")
    render_f3c()
    render_f3d()
    glm_panel("F3E", 221.847, 240.009, True)
    use_font_family("FS11")
    glm_panel("FS11", PAGE_W, 190.0, False)
    use_font_family(None)
    print("Rendered individual panels in", OUT)
