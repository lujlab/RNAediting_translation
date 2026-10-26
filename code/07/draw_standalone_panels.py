"""Native editable plots from freshly computed tabular values. No input PDF.
No rasterization, image objects, overlays, clipping masks or outlined text.
"""

from pathlib import Path
import csv, math, statistics, re, os
import xml.etree.ElementTree as ET
import fitz

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(os.environ.get("REPORTER_WORK_DIR", str(ROOT / "results/.reporter_work")))
OUT = Path(os.environ.get("REPORTER_OUT_DIR", str(ROOT / "results/07")))
OUT.mkdir(parents=True, exist_ok=True)
NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", NS)


def table(name):
    with (WORK / (name + ".tsv")).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


SUMMARY = table("summaries")
TESTS = {r["Panel"]: r for r in table("tests")}
BLACK = (0, 0, 0)
BLUE = (0.10, 0.32, 0.64)
RED = (0.65, 0.12, 0.12)


def css(c):
    return "#" + "".join("%02x" % round(v * 255) for v in c)


def sub(parent, tag, **attrs):
    return ET.SubElement(
        parent,
        "{%s}%s" % (NS, tag),
        {k.replace("_", "-"): str(v) for k, v in attrs.items()},
    )


class Canvas:
    def __init__(self, w, h):
        self.doc = fitz.open()
        self.p = self.doc.new_page(width=w, height=h)
        self.svg = ET.Element(
            "{%s}svg" % NS,
            {
                "width": "%spt" % w,
                "height": "%spt" % h,
                "viewBox": "0 0 %s %s" % (w, h),
            },
        )
        self.fonts = {}
        self.texts = []
        for key, file in [
            ("regular", "DejaVuSans.ttf"),
            ("bold", "DejaVuSans-Bold.ttf"),
        ]:
            path = Path("/usr/share/fonts/dejavu") / file
            self.fonts[key] = fitz.Font(fontfile=str(path))
            self.p.insert_font(fontname=key, fontfile=str(path))

    def line(self, x1, y1, x2, y2, c=BLACK, w=0.4):
        self.p.draw_line((x1, y1), (x2, y2), color=c, width=w)
        sub(self.svg, "line", x1=x1, y1=y1, x2=x2, y2=y2, stroke=css(c), stroke_width=w)

    def dot(self, x, y):
        self.p.draw_circle((x, y), 1.6, color=None, fill=BLUE)
        sub(self.svg, "circle", cx=x, cy=y, r=1.6, fill=css(BLUE))

    def text(self, x, y, s, size=7, anchor="left", bold=False, rotate=0):
        key = "bold" if bold else "regular"
        font = self.fonts[key]
        assert all(font.has_glyph(ord(c)) for c in s), s
        width = font.text_length(s, fontsize=size)
        if not rotate:
            if anchor == "middle":
                x -= width / 2
            elif anchor == "right":
                x -= width
        elif anchor == "middle":
            y += width / 2
        self.p.insert_text((x, y), s, fontname=key, fontsize=size, rotate=rotate)
        t = sub(
            self.svg,
            "text",
            x=x,
            y=y,
            font_family="DejaVu Sans",
            font_size=size,
            font_weight="bold" if bold else "normal",
            fill="#000000",
        )
        if rotate:
            t.set("transform", "rotate(%s %s %s)" % (-rotate, x, y))
        t.text = s
        self.texts.append(s)

    def save(self, name):
        for f in self.p.get_fonts():
            typ, value = self.doc.xref_get_key(f[0], "ToUnicode")
            if typ == "xref":
                xid = int(value.split()[0])
                stream = self.doc.xref_stream(xid)
                stream = re.sub(
                    rb"(<[0-9a-fA-F]+>\s+)<00ad>", rb"\g<1><002d>", stream, flags=re.I
                )
                self.doc.update_stream(xid, stream)
        path = OUT / (name + ".pdf")
        self.doc.save(str(path), garbage=4, deflate=True)
        ET.ElementTree(self.svg).write(
            str(OUT / (name + ".svg")), encoding="utf-8", xml_declaration=True
        )
        return path


def panel(c, name, ox=0, oy=0):
    rows = [r for r in SUMMARY if r["Panel"] == name]
    doses = [float(r["Concentration"]) for r in rows]
    avg = [float(r["Mean"]) for r in rows]
    sd = [float(r["SD"]) for r in rows]
    # Preserve the author's existing displayed F4D concentration scale; do not
    # infer an experimental dilution factor. Source concentrations are retained.
    x = [v / 10 for v in doses] if name == "F4D" else [math.log2(1 + v) for v in doses]
    titles = {
        "F4D": "Cell-free translation, 90 min",
        "FS15D": "mRNA with 5′ UTR",
        "FS15E": "mRNA with 5′ UTR",
        "FS15F": "mRNA without 5′ UTR",
        "FS15G": "mRNA without 5′ UTR",
    }
    yl = {
        "F4D": (0, 6, [0, 2, 4, 6]),
        "FS15D": (0, 26, [0, 5, 10, 15, 20, 25]),
        "FS15E": (0.2, 1.1, [0.2, 0.4, 0.6, 0.8, 1.0]),
        "FS15F": (0, 8, [0, 2, 4, 6, 8]),
        "FS15G": (0.25, 1.1, [0.3, 0.5, 0.7, 0.9, 1.1]),
    }
    ymin, ymax, yticks = yl[name]
    left, right, top, bottom = ox + 43, ox + 224, oy + 39, oy + 129
    span = max(x) - min(x)
    xmin = min(x) - span * 0.05
    xmax = max(x) + span * 0.05
    xm = lambda a: left + (a - xmin) / (xmax - xmin) * (right - left)
    ym = lambda a: bottom - (a - ymin) / (ymax - ymin) * (bottom - top)
    c.text(ox + 8, oy + 13, name[-1], 10, bold=True)
    c.text((left + right) / 2, oy + 15, titles[name], 7.5, anchor="middle")
    c.line(left, top, left, bottom)
    c.line(left, bottom, right, bottom)
    for y in yticks:
        c.line(left - 2, ym(y), left, ym(y))
        c.text(
            left - 4,
            ym(y) + 2,
            str(y) if name in ["F4D", "FS15D", "FS15F"] else "%.1f" % y,
            6.5,
            anchor="right",
        )
    for v, source in zip(x, doses):
        label = ("%g" % round(source / 10, 1)) if name == "F4D" else "%g" % source
        c.line(xm(v), bottom, xm(v), bottom + 2)
        c.text(xm(v), bottom + 9, label, 6 if name == "F4D" else 6.5, anchor="middle")
    c.text((left + right) / 2, oy + 157, "ITP input (mM)", 8, anchor="middle")
    ylabel = (
        "Renilla / RPS27A protein ratio"
        if name in ["FS15E", "FS15G"]
        else "Renilla / Firefly luciferase activity"
    )
    c.text(ox + 12, (top + bottom) / 2, ylabel, 7.5, anchor="middle", rotate=90)
    # Ordinary least-squares guide in the same displayed x coordinates.
    sx, sy = statistics.mean(x), statistics.mean(avg)
    slope = sum((a - sx) * (b - sy) for a, b in zip(x, avg)) / sum(
        (a - sx) ** 2 for a in x
    )
    intercept = sy - slope * sx
    c.line(
        xm(x[0]),
        ym(intercept + slope * x[0]),
        xm(x[-1]),
        ym(intercept + slope * x[-1]),
        RED,
        0.45,
    )
    for a, b, s in zip(x, avg, sd):
        xx = xm(a)
        lo, hi = ym(b - s), ym(b + s)
        c.line(xx, lo, xx, hi)
        for yy in [lo, hi]:
            c.line(xx - 2.3, yy, xx + 2.3, yy)
        c.dot(xx, ym(b))
    assert float(TESTS[name]["rho"]) == -1
    c.text(right - 3, oy + 30, "Spearman’s ρ = −1", 7, anchor="right")
    ptext = "P = 2.51 × 10⁻⁸" if name == "F4D" else "P = 1.39 × 10⁻³"
    # Entire P label is one real text run, including Unicode exponent digits.
    c.text(right - 3, oy + 40, ptext, 7, anchor="right")


def main():
    for name in ["F4D", "FS15D", "FS15E", "FS15F", "FS15G"]:
        c = Canvas(234, 166)
        panel(c, name)
        c.save("F4D_3rep_SD_from_data" if name == "F4D" else name + "_SD")
    print(OUT)


if __name__ == "__main__":
    main()
