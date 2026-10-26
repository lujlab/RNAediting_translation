#!/usr/bin/env python3
"""Render the new nine-row FS11 panel at native size with the FS11 Arial fonts."""

import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGURE_ROOT = HERE.parent
if not (FIGURE_ROOT / "reference_fonts").exists():
    FIGURE_ROOT = HERE.parents[1]
RESULT = FIGURE_ROOT / "FS11_with_editing_level"
TABLE = RESULT / "FS11_with_editing_level_fixed_effects.tsv"

sys.path.insert(0, str(FIGURE_ROOT / "code"))
import render_panels as native  # noqa: E402
import render_figures as rf  # noqa: E402

TERMS = [
    ("Editing status", "edit_status_matchededited"),
    ("Editing level (+1 SD; edited)", "editing_level_within_editSD"),
    ("Focal codon preference", "category_beforeoptimal"),
    ("Codon phyloP", "phyloP_z"),
    ("Stable hairpin", "secondary_structureinside"),
    ("Translation efficiency", "logTE_z"),
    ("Transcript length", "transcript_length_z"),
    ("Upstream codon preference", "upstream_categoryoptimal"),
    ("Downstream codon preference", "downstream_categoryoptimal"),
]


def render() -> Path:
    rows = rf.read(TABLE)
    if len(rows) != 18:
        raise RuntimeError(f"Expected 18 displayed rows, found {len(rows)}")

    output_name = "FS11_with_editing_level"
    rf.OUT = RESULT
    native.use_font_family("FS11")

    page_w, page_h = 467.717, 209.0
    c = rf.new(output_name, page_w, page_h)

    label_x = 105.63
    label_y = [189.35 - i * 18.943 for i in range(9)]
    axis_x0, axis_y, axis_x1, axis_top = 111.7822, 29.9278, 446.7172, 201.33
    zero_x, step = 182.1809, 70.3986
    ticks = (-0.05, 0, 0.05, 0.10, 0.15)
    bar_h = 7.16
    xcoord = lambda value: zero_x + value / 0.05 * step

    rf.line(c, axis_x0, axis_y, axis_x1, axis_y, width=0.5)
    rf.line(c, axis_x0, axis_y, axis_x0, axis_top, width=0.5)
    for tick in ticks:
        rf.line(c, xcoord(tick), axis_y, xcoord(tick), axis_y - 2.26, width=0.4)
        native.label(c, xcoord(tick), axis_y - 8.9193, f"{tick:.2f}", 6.5, "center")

    for index, (display, term) in enumerate(TERMS):
        ly = label_y[index]
        # The extracted legacy Arial subset lacks glyphs in this added label.
        if term == "editing_level_within_editSD":
            rf.text(c, label_x, ly, display, 7, "right", font="Sans")
        else:
            native.label(c, label_x, ly, display, 7, "right")
        rf.line(c, axis_x0 - 4.093, ly + 2.393, axis_x0, ly + 2.393, width=0.45)
        for sex, color in (("Male", rf.MALE), ("Female", rf.FEMALE)):
            matches = [
                row
                for row in rows
                if row["sex_label"] == sex
                and (row["display_variable"] == display or row["term"] == term)
            ]
            if len(matches) != 1:
                raise RuntimeError(
                    f"Expected one row for {sex}/{display}; found {len(matches)}"
                )
            row = matches[0]
            value = rf.num(row["estimate"])
            yb = ly + 2.203 if sex == "Male" else ly - 5.426
            rf.rect(
                c,
                min(xcoord(0), xcoord(value)),
                yb,
                abs(xcoord(value) - xcoord(0)),
                bar_h,
                color,
            )
            if rf.num(row["p_adjust_bh_all"]) < 0.05:
                star_y = ly + 2.528 if sex == "Male" else ly - 6.101
                native.label(
                    c,
                    xcoord(value) + (3 if value >= 0 else -3),
                    star_y,
                    "*",
                    8,
                    "left" if value >= 0 else "right",
                    font="SansBold",
                )

    native.label(c, 200.1934, 8.0651, "Fixed-effect coefficient for ln(1 + gmRC)", 8)
    native.label(c, 408.2896, 78.84, "Sex", 7)
    rf.rect(c, 408.681, 68.32, 9.547, 7.16, rf.MALE)
    native.label(c, 421.5562, 69.82, "Male", 7)
    rf.rect(c, 408.681, 58.82, 9.547, 7.16, rf.FEMALE)
    native.label(c, 421.5562, 60.09, "Female", 7)
    c.save()

    source = RESULT / f"{output_name}.pdf"
    destinations = [
        FIGURE_ROOT / "figures" / f"{output_name}.pdf",
        FIGURE_ROOT / "panels" / f"{output_name}.pdf",
        FIGURE_ROOT / "figures" / "FS11.pdf",
        FIGURE_ROOT / "panels" / "FS11.pdf",
    ]
    for destination in destinations:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    return source


if __name__ == "__main__":
    print(render())
