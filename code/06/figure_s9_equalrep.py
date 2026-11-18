#!/usr/bin/env python3
"""Final FS9: coverage-qualified pairs, bilateral any-edited support,
equal-total rescaled counts pooled without pseudocounts. Python stdlib only.
The publication layout is fixed to the approved 84 codons; a changed set fails
explicitly rather than silently omitting/adding rows.
"""

import argparse
import base64
import csv
import json
import math
import os
import re
from collections import defaultdict
from pathlib import Path

COHORTS = ["Dmel_WT", "Dsim", "Dsec", "Dere"]
SPECIES = {"mel": "Dmel_WT", "sim": "Dsim", "sec": "Dsec", "ere": "Dere"}
STATES = ["AA", "IA", "AI", "II"]
STEM = "Supplementary_Figure_9_equalrep_both_edited_EDITABLE"


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def table(path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def shade(v):
    if v is None:
        return (0.906, 0.906, 0.906)
    low, mid, high = (0.145, 0.380, 0.627), (0.969, 0.965, 0.965), (0.851, 0.302, 0.157)
    t = max(-1, min(1, v / 4))
    a, b, f = (low, mid, t + 1) if t < 0 else (mid, high, t)
    return tuple(x + (y - x) * f for x, y in zip(a, b))


def draw_pdf(path, layout, fonts, cells):
    """Write a plain PDF with independent text/rectangles, no private AI data."""
    W, H = layout["size"]
    ops = []

    def text(x, y, value, size, font, color=(34 / 255, 23 / 255, 20 / 255)):
        escaped = value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
        ops.append(
            "BT /%s %.6f Tf %.6f %.6f %.6f rg 1 0 0 1 %.6f %.6f Tm (%s) Tj ET"
            % (font, size, *color, x, H - y, escaped)
        )

    def centered(x, y, value, size, font="TT1"):
        width = (
            sum(fonts["metrics"][font].get(str(ord(ch)), 0) for ch in value)
            * size
            / 1000
        )
        text(x - width / 2, y, value, size, font)

    def rectangle(rect, fill, stroke=None, sw=0.16):
        x0, y0, x1, y1 = rect
        ops.append("q %.6f %.6f %.6f rg" % tuple(fill))
        if stroke:
            ops.append("%.6f %.6f %.6f RG %.6f w" % (*stroke, sw))
        ops.append(
            "%.6f %.6f %.6f %.6f re %s Q"
            % (x0, H - y1, x1 - x0, y1 - y0, "B" if stroke else "f")
        )

    labels = []
    for slot in layout["cells"]:
        d = cells[(slot["cohort"], slot["target_handle"], slot["allele"])]
        v = None if d["status"] == "both_zero" else d["value"]
        r = slot["rect"]
        x = (r[0] + r[2]) / 2
        y = (r[1] + r[3]) / 2 + 2.04
        rectangle(r, shade(v), (0.827, 0.863, 0.890))
        label = (
            "NA"
            if v is None
            else (
                "-inf"
                if v == -math.inf
                else (
                    "+inf"
                    if v == math.inf
                    else "{:.2f}".format(0 if abs(v) < 0.005 else v)
                )
            )
        )
        if v is not None and math.isinf(v):
            centered(x - 2.3, y, "-" if v < 0 else "+", 6.35)
            text(x + 1.1 - 0.713 * 6.35 / 2, y, "\xa5", 6.35, "symb")
        else:
            centered(x, y, label, 5.8 if v is None else 6.35)
        labels.append(
            dict(
                cohort=slot["cohort"],
                target_handle=slot["target_handle"],
                allele=slot["allele"],
                label=label,
            )
        )
    for d in layout["legend"]:
        rectangle(d["rect"], d["fill"], d["color"], d["width"] or 0)
    for s in layout["static_text"]:
        f = s["font"]
        alias = "TT0" if "Italic" in f else "TT2" if "Bold" in f else "TT1"
        color = tuple(((s["color"] >> shift) & 255) / 255 for shift in [16, 8, 0])
        text(*s["origin"], s["text"], s["size"], alias, color)
    objects = {}
    old_ids = sorted(int(k) for k in fonts["objects"])
    mapping = {old: 5 + i for i, old in enumerate(old_ids)}
    for old, new in mapping.items():
        item = fonts["objects"][str(old)]
        obj = re.sub(
            r"\b(\d+) 0 R\b",
            lambda m: "%d 0 R" % mapping[int(m.group(1))],
            item["object"],
        ).encode("latin1")
        if item["stream"] is not None:
            obj += b"\nstream\n" + base64.b64decode(item["stream"]) + b"\nendstream"
        objects[new] = obj
    resources = (
        "<< /Font << "
        + " ".join(
            "/%s %d 0 R" % (name, mapping[old]) for name, old in fonts["roots"].items()
        )
        + " >> >>"
    )
    objects[1] = b"<< /Type /Catalog /Pages 2 0 R >>"
    objects[2] = b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>"
    objects[3] = (
        "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.6f %.6f] /Resources %s /Contents 4 0 R >>"
        % (W, H, resources)
    ).encode()
    stream = ("\n".join(ops) + "\n").encode("latin1")
    objects[4] = (
        ("<< /Length %d >>\nstream\n" % len(stream)).encode() + stream + b"endstream"
    )
    result = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = [0]
    for i in range(1, len(objects) + 1):
        offsets.append(len(result))
        result.extend(("%d 0 obj\n" % i).encode() + objects[i] + b"\nendobj\n")
    start = len(result)
    result.extend(("xref\n0 %d\n0000000000 65535 f \n" % len(offsets)).encode())
    for n in offsets[1:]:
        result.extend(("%010d 00000 n \n" % n).encode())
    result.extend(
        (
            "trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n"
            % (len(offsets), start)
        ).encode()
    )
    tmp = path.with_suffix(".pdf.tmp")
    tmp.write_bytes(result)
    os.replace(str(tmp), str(path))
    return labels


def main():
    script = Path(__file__).resolve()
    root = script.parents[2]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog", type=Path, required=True)
    p.add_argument("--counts", type=Path, required=True)
    p.add_argument(
        "--outdir", type=Path, default=root / "FigureS9_equalrep_both_edited_20261005"
    )
    args = p.parse_args()
    assets = script.parent / "figure_s9_assets"
    layout = json.loads((assets / "layout.json").read_text())
    fonts = json.loads((assets / "fonts.json").read_text())
    catalog = read(args.catalog)
    raw = read(args.counts)
    byrep = defaultdict(dict)
    for r in raw:
        c = "Dmel_WT" if r["genotype"] in ("ISO1", "WBg") else r["genotype"]
        if c not in COHORTS or int(r["edit_site_count"]) != 2:
            continue
        k = (c, r["target_handle"], r["genotype"], r["sex"], r["rep"], r["sample"])
        assert r["allele"] not in byrep[k], ("duplicate", k)
        byrep[k][r["allele"]] = r
    qualified = defaultdict(list)
    rep_rows = []
    for k, a in sorted(byrep.items()):
        assert set(a) == set(STATES), k
        rt = sum(float(a[s]["RNA_count"]) for s in STATES)
        pt = sum(float(a[s]["RPF_count"]) for s in STATES)
        assert all(
            float(r["RNA_total"]) == rt and float(r["RPF_total"]) == pt
            for r in a.values()
        ), k
        if rt < 5 or pt < 5:
            continue
        state = {}
        for s in STATES:
            rn, pn = float(a[s]["RNA_count"]), float(a[s]["RPF_count"])
            assert rn >= 0 and pn >= 0
            state[s] = (rn / rt, pn / pt, rn, pn)
            rep_rows.append(
                dict(
                    cohort=k[0],
                    target_handle=k[1],
                    genotype=k[2],
                    sex=k[3],
                    rep=k[4],
                    sample=k[5],
                    allele=s,
                    RNA_count=rn,
                    RPF_count=pn,
                    RNA_total=rt,
                    RPF_total=pt,
                    RNA_fraction=rn / rt,
                    RPF_fraction=pn / pt,
                )
            )
        qualified[k[:2]].append(state)
    disposition = []
    cells = {}
    retained = set()
    for m in catalog:
        c, h = SPECIES[m["species"]], m["target_handle"]
        reps = qualified.get((c, h), [])
        n = len(reps)
        er = sum(r[s][2] for r in reps for s in STATES[1:])
        ep = sum(r[s][3] for r in reps for s in STATES[1:])
        keep = n > 0 and er > 0 and ep > 0
        site = dict(
            cohort=c,
            target_handle=h,
            gene=m["gene"],
            target_transcript=m["target_transcript"],
            aa_position_1based=m["aa_position_1based"],
            n_pairs=n,
            edited_RNA_sum=er,
            edited_RPF_sum=ep,
            retained=int(keep),
        )
        disposition.append(site)
        if not keep:
            continue
        retained.add((c, h))
        for s in STATES:
            # Pool counts rescaled to one total per pair, then divide by n.
            rp = sum(r[s][0] for r in reps) / n
            pp = sum(r[s][1] for r in reps) / n
            status = (
                "both_zero"
                if rp == pp == 0
                else "RNA_zero" if rp == 0 else "RPF_zero" if pp == 0 else "finite"
            )
            v = (
                "NA"
                if status == "both_zero"
                else (
                    math.inf
                    if rp == 0
                    else -math.inf if pp == 0 else math.log2(pp / rp)
                )
            )
            cells[c, h, s] = dict(
                **site,
                allele=s,
                RNA_mean_fraction=rp,
                RPF_mean_fraction=pp,
                value=v,
                status=status,
                RNA_count_sum=sum(r[s][2] for r in reps),
                RPF_count_sum=sum(r[s][3] for r in reps)
            )
    expected = {(r["cohort"], r["target_handle"]) for r in layout["cells"]}
    if retained != expected:
        raise ValueError(
            "Codon set differs from publication layout. Missing=%r; added=%r"
            % (expected - retained, retained - expected)
        )
    summary = []
    for c in COHORTS:
        ss = [r for r in disposition if r["cohort"] == c and r["retained"]]
        ii = [r for r in cells.values() if r["cohort"] == c and r["allele"] == "II"]
        summary.append(
            dict(
                cohort=c,
                codons=len(ss),
                one_pair=sum(r["n_pairs"] == 1 for r in ss),
                II_positive_finite=sum(
                    r["status"] == "finite" and r["value"] > 0 for r in ii
                ),
                II_negative_finite=sum(
                    r["status"] == "finite" and r["value"] < 0 for r in ii
                ),
                II_minus_inf=sum(r["status"] == "RPF_zero" for r in ii),
                II_plus_inf=sum(r["status"] == "RNA_zero" for r in ii),
                II_NA=sum(r["status"] == "both_zero" for r in ii),
            )
        )
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    labels = draw_pdf(out / "figures" / (STEM + ".pdf"), layout, fonts, cells)
    for name, rows in [
        ("summary", summary),
        ("all_candidate_disposition", disposition),
        ("retained_codons", [r for r in disposition if r["retained"]]),
        ("heatmap_data", list(cells.values())),
        (
            "qualified_replicate_fractions",
            [r for r in rep_rows if (r["cohort"], r["target_handle"]) in retained],
        ),
        ("displayed_cell_labels", labels),
    ]:
        table(out / (name + ".tsv"), rows)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
