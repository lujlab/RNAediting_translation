#!/usr/bin/env python3
"""Fresh editable vector pages; no template overlays, images or clipping masks."""

import os
import json
from pathlib import Path
import fitz

OUT = Path(
    os.environ.get(
        "FS14_OUTPUT_DIR",
        str(Path(__file__).resolve().parents[3] / "results/06/FS14_F4B"),
    )
)
D = json.loads((OUT / "figure_data.json").read_text())
W, H = 467.7170104980469, 566.9290161132812
FONT_PATH = Path(
    os.environ.get(
        "FIGURE_FONT_FILE", "/usr/share/fonts/liberation/LiberationSans-Regular.ttf"
    )
)
if "FIGURE_FONT_FILE" in os.environ and not FONT_PATH.is_file():
    raise FileNotFoundError(str(FONT_PATH))
F = fitz.Font(fontfile=str(FONT_PATH)) if FONT_PATH.is_file() else fitz.Font("helv")
FONT = F.buffer
C = [
    (0.792, 0.800, 0.800),
    (0.976, 0.647, 0.169),
    (0.396, 0.784, 0.816),
    (0.604, 0.478, 0.733),
]
LABELS = ["AAG", "IAG", "AIG", "IIG"]
STATES = ["AA", "IA", "AI", "II"]


def text(p, x, y, s, size=6, align="left", rotate=0):
    p.insert_font(fontname="ArialFigure", fontbuffer=FONT)
    if align == "center":
        x -= F.text_length(s, fontsize=size) / 2
    if align == "right":
        x -= F.text_length(s, fontsize=size)
    p.insert_text((x, y), s, fontname="ArialFigure", fontsize=size, rotate=rotate)


def line(p, x, y, x2, y2, width=0.5):
    p.draw_line((x, y), (x2, y2), color=(0, 0, 0), width=width)


def legend(p):
    for i, (label, color) in enumerate(zip(LABELS, C)):
        x = 184 + i * 39
        p.draw_rect(fitz.Rect(x, 6, x + 9, 13), color=None, fill=color)
        text(p, x + 12, 12.98, label, 7)


def panel(p, key, col, row, s14=False):
    sp, sex = key.split("_")
    dx = 219.32 * col
    top = [39.35, 171.84, 307.37][row]
    tx = {
        "Dmel": "FBtr0334522",
        "Dsim": "FBtr0361051",
        "Dsec": "FBtr0200880",
        "Dere": "FBtr0143136",
    }[sp]
    title = (
        ("Dmel (WT)" if sp == "Dmel" else sp) + " " + sex + " (" + tx + " 137th codon)"
    )
    text(p, 143 + dx, top - 9, title, 8, "center")
    x0, x1 = 42.444 + dx, 131.504 + dx
    y0, y1 = top, top + (78.04 if s14 else 91.357)
    vals = D[key]["mRC"]
    upper, step = {
        "Dmel_F": (0.35, 0.1),
        "Dmel_M": (0.55, 0.1),
        "Dsim_F": (9, 2),
        "Dsim_M": (6, 1),
        "Dsec_M": (0.275, 0.05),
        "Dsec_F": (0.44, 0.1),
        "Dere_M": (0.33, 0.1),
        "Dere_F": (0.44, 0.1),
    }[key]
    line(p, x0, y0, x0, y1)
    line(p, x0, y1, x1, y1)
    for i in range(int(upper / step) + 1):
        v = i * step
        y = y1 - v / upper * (y1 - y0)
        line(p, x0 - 1.89, y, x0, y)
        label = (
            ("{:.2f}".format(v) if step == 0.05 else "{:.1f}".format(v))
            if step < 1
            else str(int(v))
        )
        text(p, x0 - 3.4, y + 2.4, label, 6.5, "right")
    text(p, x0 - 20, (y0 + y1) / 2 + 10, "mRC", 8, rotate=90)
    for i, (v, color, label) in enumerate(zip(vals, C, LABELS)):
        center = x0 + (x1 - x0) * (1 + 2 * i) / 8
        bw = (x1 - x0) / 8
        if v is None:
            text(p, center, y1 - 2.05, "NA", 6, "center")
        else:
            yy = y1 - v / upper * (y1 - y0)
            if v > 0:
                p.draw_rect(
                    fitz.Rect(center - bw / 2, yy, center + bw / 2, y1),
                    color=(0, 0, 0),
                    fill=color,
                    width=0.75,
                )
            else:
                line(p, center - bw / 2, y1, center + bw / 2, y1, 0.75)
            text(p, center, yy - 2.05, "{:.3f}".format(v), 6, "center")
        line(p, center, y1, center, y1 + 1.89)
        text(p, center, y1 + 9, label, 7, "center")
    x0, x1 = 162.227 + dx, 227.246 + dx
    y1 = top + 78.04
    for i in range(6):
        v = i / 5
        y = y1 - v * 78.04
        line(p, x0 - 1.89, y, x0, y)
        text(p, x0 - 3.4, y + 2.4, "{:.1f}".format(v), 6.5, "right")
    text(p, x0 - 18, top + 80.54, "Editing-state proportion", 8, rotate=90)
    for assay, x, label in zip(
        ["RNA", "RPF", "MS"],
        [166.493 + dx, 186.811 + dx, 207.130 + dx],
        ["mRNA-seq", "Ribo-seq", "MS"],
    ):
        yy = top
        for state, color in zip(STATES, C):
            value = D[key]["row"][assay + "_" + state + "_percent"] / 100
            if value > 0:
                p.draw_rect(
                    fitz.Rect(x, yy, x + 15.847, yy + value * 78.04),
                    color=None,
                    fill=color,
                )
            yy += value * 78.04
        center = x + 15.847 / 2
        line(p, center, y1, center, y1 + 1.89)
        width = F.text_length(label, fontsize=6.5)
        if s14 and assay == "MS":
            text(p, center, y1 + 10, label, 6.5, "center")
        else:
            cosine, sine, angle = (0.8660254, 0.5, 30) if s14 else (0.7071, 0.7071, 45)
            origin = fitz.Point(center - width * cosine, y1 + 6 + width * sine)
            p.insert_text(
                origin,
                label,
                fontname="ArialFigure",
                fontsize=6.5,
                morph=(origin, fitz.Matrix(angle)),
            )
    line(p, x0, top, x0, y1)
    line(p, x0, y1, x1, y1)


for name, keys, height in [
    (
        "FS14_DIA_WBg_only.pdf",
        ["Dmel_F", "Dsim_F", "Dsec_M", "Dsec_F", "Dere_M", "Dere_F"],
        H,
    ),
    ("F4B_DIA_WBg_only.pdf", ["Dmel_M", "Dsim_M"], 148),
]:
    doc = fitz.open()
    p = doc.new_page(width=W, height=height)
    legend(p)
    for i, key in enumerate(keys):
        panel(p, key, i % 2, i // 2, s14=name.startswith("FS14"))
    doc.set_metadata(
        {
            "title": name[:-4]
            + "; Dmel WT=WBg only; DIA-NN 100-aa, no duplicate unlisted rows"
        }
    )
    doc.save(str(OUT / name), garbage=4, deflate=True)
    assert not p.get_images()
    doc.close()
print("Created two fresh vector PDFs.")
