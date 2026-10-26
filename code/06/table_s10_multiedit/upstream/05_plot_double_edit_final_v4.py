#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Final vector plotting for the rebuilt multi-edit codon pipeline v4.

Double-edit main figure (one two-page PDF per species):
  page 1: AA / IA / AI / II mean log2(RPF/RNA) enrichment + Lm
  page 2: AA / IA / AI / II mean MS intensity proportion + Lp

Main RNA/RPF cohorts:
  mel -> Dmel_WT
  sim -> Dsim
  sec -> Dsec
  ere -> Dere

MS panels:
  mel -> WBg_ISO1, AS, AG
  sim -> Dsim
  sec -> Dsec
  ere -> Dere

The RNA/RPF heatmap contains only rows with RNA_RPF_heatmap_eligible == yes.
MS rows are displayed when at least one within-site composition statistic is
available. Raw MS intensity is never used as the color scale.

A separate two-page PDF is generated for the triple-edit codon(s), when present.
All output is vector PDF via ReportLab.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path
from typing import Any, Iterable

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas

SPECIES_ORDER = ("mel", "sim", "sec", "ere")
MAIN_COHORT = {"mel": "Dmel_WT", "sim": "Dsim", "sec": "Dsec", "ere": "Dere"}
MS_GROUPS = {
    "mel": ("WBg_ISO1", "AS", "AG"),
    "sim": ("Dsim",),
    "sec": ("Dsec",),
    "ere": ("Dere",),
}
DOUBLE_ALLELES = ("AA", "IA", "AI", "II")
TRIPLE_ALLELES = ("AAA", "IAA", "AIA", "AAI", "IIA", "IAI", "AII", "III")


def clean(x: Any) -> str:
    return "" if x is None else str(x).strip().strip("\ufeff")


def optional_float(x: Any) -> float | None:
    s = clean(x)
    if not s or s.upper() in {"NA", "NAN", "NULL"}:
        return None
    try:
        value = float(s)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def natural_key(value: str) -> list[Any]:
    import re

    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", value)]


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 1.0
    pos = (len(ordered) - 1) * fraction
    lo, hi = math.floor(pos), math.ceil(pos)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - pos) + ordered[hi] * (pos - lo)


def register_fonts() -> tuple[str, str]:
    candidates = [
        (
            Path("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"),
            Path("/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf"),
        ),
        (
            Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"),
            Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"),
        ),
        (
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        ),
    ]
    for regular, bold in candidates:
        if regular.is_file() and bold.is_file():
            pdfmetrics.registerFont(TTFont("EditRegular", str(regular)))
            pdfmetrics.registerFont(TTFont("EditBold", str(bold)))
            return "EditRegular", "EditBold"
    return "Helvetica", "Helvetica-Bold"


def mix(low: tuple[float, float, float], high: tuple[float, float, float], t: float):
    t = max(0.0, min(1.0, t))
    return tuple(low[i] + (high[i] - low[i]) * t for i in range(3))


def rpf_color(value: float, limit: float) -> tuple[float, float, float]:
    # blue -> nearly white -> red
    low = (0.13, 0.40, 0.67)
    mid = (0.97, 0.97, 0.97)
    high = (0.90, 0.30, 0.15)
    z = max(-1.0, min(1.0, value / limit if limit else 0.0))
    if z < 0:
        return mix(low, mid, z + 1.0)
    return mix(mid, high, z)


def ms_color(value: float) -> tuple[float, float, float]:
    # fixed 0..1 scale, because values are within-site proportions.
    return mix((0.97, 0.98, 1.00), (0.03, 0.19, 0.42), value**0.75)


def level_color(value: float) -> tuple[float, float, float]:
    # fixed 0..1 scale for Lm/Lp.
    return mix((1.00, 0.98, 0.88), (0.70, 0.36, 0.02), value**0.65)


def luminance(color: tuple[float, float, float]) -> float:
    return 0.2126 * color[0] + 0.7152 * color[1] + 0.0722 * color[2]


def set_text_for_fill(pdf: canvas.Canvas, color: tuple[float, float, float]) -> None:
    if luminance(color) < 0.48:
        pdf.setFillColorRGB(1, 1, 1)
    else:
        pdf.setFillColorRGB(0.05, 0.12, 0.22)


def format_value(value: float | None, digits: int = 2) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def row_label(row: dict[str, str], max_chars: int = 58) -> str:
    gene = clean(row.get("gene", ""))
    handle = clean(row.get("target_handle", ""))
    genomic = clean(row.get("genomic_sites", ""))
    label = " | ".join(x for x in (gene, handle, genomic) if x)
    if len(label) <= max_chars:
        return label
    return label[: max_chars - 1] + "…"


def draw_legend_diverging(
    pdf: canvas.Canvas,
    x: float,
    y: float,
    width: float,
    limit: float,
    regular: str,
) -> None:
    steps = 120
    for i in range(steps):
        value = -limit + 2 * limit * i / max(1, steps - 1)
        pdf.setFillColorRGB(*rpf_color(value, limit))
        pdf.rect(x + width * i / steps, y, width / steps + 0.4, 9, fill=1, stroke=0)
    pdf.setFillColorRGB(0.12, 0.18, 0.27)
    pdf.setFont(regular, 6.8)
    pdf.drawString(x, y - 10, f"-{limit:.2g}")
    pdf.drawCentredString(x + width / 2, y - 10, "0")
    pdf.drawRightString(x + width, y - 10, f"+{limit:.2g}")
    pdf.drawString(x + width + 12, y + 1, "mean log2(RPF/RNA)")


def draw_legend_01(
    pdf: canvas.Canvas,
    x: float,
    y: float,
    width: float,
    regular: str,
    label: str,
    color_fun,
) -> None:
    steps = 100
    for i in range(steps):
        v = i / max(1, steps - 1)
        pdf.setFillColorRGB(*color_fun(v))
        pdf.rect(x + width * i / steps, y, width / steps + 0.4, 9, fill=1, stroke=0)
    pdf.setFillColorRGB(0.12, 0.18, 0.27)
    pdf.setFont(regular, 6.8)
    pdf.drawString(x, y - 10, "0")
    pdf.drawCentredString(x + width / 2, y - 10, "0.5")
    pdf.drawRightString(x + width, y - 10, "1")
    pdf.drawString(x + width + 12, y + 1, label)


def draw_rna_rpf_page(
    pdf: canvas.Canvas,
    species: str,
    cohort: str,
    rows: list[dict[str, str]],
    alleles: tuple[str, ...],
    rpf_limit: float,
    regular: str,
    bold: str,
    title_prefix: str = "Double-edit",
) -> None:
    shown = [
        r for r in rows if clean(r.get("RNA_RPF_heatmap_eligible", "")).lower() == "yes"
    ]
    shown.sort(
        key=lambda r: (
            -(optional_float(r.get("RNA_editing_level_mean")) or -1.0),
            natural_key(clean(r.get("target_handle", ""))),
        )
    )

    row_h = 15.0
    left = 330.0
    cell_w = 82.0 if len(alleles) <= 4 else 62.0
    level_w = 68.0
    level_x = left + len(alleles) * cell_w + 9
    width = level_x + level_w + 42
    top = 112.0
    bottom = 72.0
    n_display = max(1, len(shown))
    height = max(620.0, top + bottom + n_display * row_h)
    pdf.setPageSize((width, height))

    pdf.setFillColorRGB(0.05, 0.16, 0.30)
    pdf.setFont(bold, 14)
    pdf.drawString(34, height - 34, f"{species} {title_prefix} codons: RNA/RPF")
    pdf.setFont(regular, 8.2)
    pdf.setFillColorRGB(0.30, 0.38, 0.47)
    pdf.drawString(
        34,
        height - 52,
        f"{cohort}; allele-wise mean log2(RPF/RNA); only RNA/RPF-eligible codons are shown.",
    )

    header_y = height - 86
    pdf.setFillColorRGB(0.05, 0.16, 0.30)
    pdf.setFont(bold, 8.5)
    pdf.drawRightString(left - 9, header_y, "Gene | transcript_codon | genomic sites")
    for i, allele in enumerate(alleles):
        pdf.drawCentredString(left + i * cell_w + cell_w / 2, header_y, allele)
    pdf.drawCentredString(level_x + level_w / 2, header_y, "Lm")

    if not shown:
        pdf.setFont(regular, 9)
        pdf.setFillColorRGB(0.35, 0.42, 0.50)
        pdf.drawString(
            left, header_y - 35, "No codons passed the RNA/RPF inclusion criteria."
        )
    else:
        for j, row in enumerate(shown):
            y = header_y - 14 - (j + 1) * row_h
            pdf.setFillColorRGB(0.08, 0.16, 0.27)
            pdf.setFont(regular, 6.3)
            pdf.drawRightString(left - 9, y + 4.0, row_label(row))
            for i, allele in enumerate(alleles):
                x = left + i * cell_w
                value = optional_float(
                    row.get(f"RNA_RPF_{allele}_mean_log2_enrichment")
                )
                color = (
                    (0.93, 0.93, 0.93) if value is None else rpf_color(value, rpf_limit)
                )
                pdf.setFillColorRGB(*color)
                pdf.rect(x, y, cell_w, row_h, fill=1, stroke=0)
                pdf.setStrokeColorRGB(0.83, 0.87, 0.90)
                pdf.setLineWidth(0.25)
                pdf.rect(x, y, cell_w, row_h, fill=0, stroke=1)
                set_text_for_fill(pdf, color)
                pdf.setFont(regular, 6.7)
                pdf.drawCentredString(x + cell_w / 2, y + 4.0, format_value(value))

            lm = optional_float(row.get("RNA_editing_level_mean"))
            lcolor = (
                (0.93, 0.93, 0.93)
                if lm is None
                else level_color(max(0.0, min(1.0, lm)))
            )
            pdf.setFillColorRGB(*lcolor)
            pdf.rect(level_x, y, level_w, row_h, fill=1, stroke=0)
            pdf.setStrokeColorRGB(0.83, 0.87, 0.90)
            pdf.rect(level_x, y, level_w, row_h, fill=0, stroke=1)
            set_text_for_fill(pdf, lcolor)
            pdf.setFont(regular, 6.7)
            pdf.drawCentredString(level_x + level_w / 2, y + 4.0, format_value(lm))

    legend_y = 37
    draw_legend_diverging(pdf, left, legend_y, 210, rpf_limit, regular)
    draw_legend_01(pdf, level_x, legend_y, level_w, regular, "Lm", level_color)
    pdf.setFillColorRGB(0.24, 0.31, 0.39)
    pdf.setFont(regular, 6.8)
    pdf.drawString(
        34,
        17,
        f"n = {len(shown)} codons. Lm = mean RNA edited-allele fraction across eligible RNA samples.",
    )


def has_ms_signal(row: dict[str, str], alleles: tuple[str, ...]) -> bool:
    if optional_float(row.get("Lp_mean")) is not None:
        return True
    return any(
        optional_float(row.get(f"{a}_mean_intensity_proportion")) is not None
        for a in alleles
    )


def draw_ms_page(
    pdf: canvas.Canvas,
    species: str,
    ms_rows: list[dict[str, str]],
    main_order_rows: list[dict[str, str]],
    groups: tuple[str, ...],
    alleles: tuple[str, ...],
    regular: str,
    bold: str,
    title_prefix: str = "Double-edit",
) -> None:
    lm_by_handle = {
        clean(r.get("target_handle")): optional_float(r.get("RNA_editing_level_mean"))
        for r in main_order_rows
    }

    group_rows: dict[str, list[dict[str, str]]] = {}
    for group in groups:
        selected = [
            r
            for r in ms_rows
            if clean(r.get("MS_group")) == group and has_ms_signal(r, alleles)
        ]
        selected.sort(
            key=lambda r: (
                -(
                    lm_by_handle.get(clean(r.get("target_handle")))
                    if lm_by_handle.get(clean(r.get("target_handle"))) is not None
                    else -1.0
                ),
                natural_key(clean(r.get("target_handle", ""))),
            )
        )
        group_rows[group] = selected

    row_h = 14.0
    left = 330.0
    cell_w = 82.0 if len(alleles) <= 4 else 62.0
    level_w = 68.0
    level_x = left + len(alleles) * cell_w + 9
    width = level_x + level_w + 42
    rows_needed = sum(max(1, len(group_rows[g])) for g in groups)
    group_overhead = 51 * len(groups)
    height = max(650.0, 120 + rows_needed * row_h + group_overhead + 75)
    pdf.setPageSize((width, height))

    pdf.setFillColorRGB(0.05, 0.16, 0.30)
    pdf.setFont(bold, 14)
    pdf.drawString(
        34, height - 34, f"{species} {title_prefix} codons: MS allele composition"
    )
    pdf.setFont(regular, 8.2)
    pdf.setFillColorRGB(0.30, 0.38, 0.47)
    pdf.drawString(
        34,
        height - 52,
        "Cells are within-site allele intensity proportions; raw unnormalized intensity is not used as a display scale.",
    )

    cursor = height - 82
    for group in groups:
        selected = group_rows[group]
        pdf.setFillColorRGB(0.05, 0.16, 0.30)
        pdf.setFont(bold, 9.5)
        pdf.drawString(
            34, cursor, f"{group}  (n={len(selected)} codons with MS composition)"
        )
        cursor -= 19

        pdf.setFont(bold, 8.2)
        pdf.drawRightString(left - 9, cursor, "Gene | transcript_codon | genomic sites")
        for i, allele in enumerate(alleles):
            pdf.drawCentredString(left + i * cell_w + cell_w / 2, cursor, allele)
        pdf.drawCentredString(level_x + level_w / 2, cursor, "Lp")
        cursor -= 13

        if not selected:
            pdf.setFont(regular, 8)
            pdf.setFillColorRGB(0.35, 0.42, 0.50)
            pdf.drawString(left, cursor - 4, "No MS signal available.")
            cursor -= row_h
        else:
            for row in selected:
                y = cursor - row_h
                pdf.setFillColorRGB(0.08, 0.16, 0.27)
                pdf.setFont(regular, 6.2)
                pdf.drawRightString(left - 9, y + 3.6, row_label(row))
                for i, allele in enumerate(alleles):
                    x = left + i * cell_w
                    value = optional_float(
                        row.get(f"{allele}_mean_intensity_proportion")
                    )
                    color = (
                        (0.93, 0.93, 0.93)
                        if value is None
                        else ms_color(max(0.0, min(1.0, value)))
                    )
                    pdf.setFillColorRGB(*color)
                    pdf.rect(x, y, cell_w, row_h, fill=1, stroke=0)
                    pdf.setStrokeColorRGB(0.83, 0.87, 0.90)
                    pdf.setLineWidth(0.25)
                    pdf.rect(x, y, cell_w, row_h, fill=0, stroke=1)
                    set_text_for_fill(pdf, color)
                    pdf.setFont(regular, 6.6)
                    pdf.drawCentredString(x + cell_w / 2, y + 3.6, format_value(value))

                lp = optional_float(row.get("Lp_mean"))
                lcolor = (
                    (0.93, 0.93, 0.93)
                    if lp is None
                    else level_color(max(0.0, min(1.0, lp)))
                )
                pdf.setFillColorRGB(*lcolor)
                pdf.rect(level_x, y, level_w, row_h, fill=1, stroke=0)
                pdf.setStrokeColorRGB(0.83, 0.87, 0.90)
                pdf.rect(level_x, y, level_w, row_h, fill=0, stroke=1)
                set_text_for_fill(pdf, lcolor)
                pdf.setFont(regular, 6.6)
                pdf.drawCentredString(level_x + level_w / 2, y + 3.6, format_value(lp))
                cursor = y
        cursor -= 25

    legend_y = 37
    draw_legend_01(
        pdf, left, legend_y, 210, regular, "MS allele intensity proportion", ms_color
    )
    draw_legend_01(pdf, level_x, legend_y, level_w, regular, "Lp", level_color)
    pdf.setFillColorRGB(0.24, 0.31, 0.39)
    pdf.setFont(regular, 6.8)
    pdf.drawString(
        34,
        17,
        "Lp = edited protein-allele intensity / total within-site intensity; samples with total intensity = 0 are NA.",
    )


def create_species_pdf(
    species: str,
    integrated_rows: list[dict[str, str]],
    ms_rows: list[dict[str, str]],
    out: Path,
    rpf_limit: float,
    regular: str,
    bold: str,
) -> int:
    cohort = MAIN_COHORT[species]
    main_rows = [
        r
        for r in integrated_rows
        if clean(r.get("species")) == species and clean(r.get("cohort")) == cohort
    ]
    if not main_rows and not any(clean(r.get("species")) == species for r in ms_rows):
        return 0
    out.parent.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(str(out), pageCompression=1)
    draw_rna_rpf_page(
        pdf, species, cohort, main_rows, DOUBLE_ALLELES, rpf_limit, regular, bold
    )
    pdf.showPage()
    draw_ms_page(
        pdf,
        species,
        [r for r in ms_rows if clean(r.get("species")) == species],
        main_rows,
        MS_GROUPS[species],
        DOUBLE_ALLELES,
        regular,
        bold,
    )
    pdf.save()
    return 2


def create_all_species_pdf(
    integrated_rows: list[dict[str, str]],
    ms_rows: list[dict[str, str]],
    out: Path,
    rpf_limit: float,
    regular: str,
    bold: str,
) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(str(out), pageCompression=1)
    page_count = 0
    for species in SPECIES_ORDER:
        cohort = MAIN_COHORT[species]
        main_rows = [
            r
            for r in integrated_rows
            if clean(r.get("species")) == species and clean(r.get("cohort")) == cohort
        ]
        species_ms = [r for r in ms_rows if clean(r.get("species")) == species]
        if not main_rows and not species_ms:
            continue
        if page_count:
            pdf.showPage()
        draw_rna_rpf_page(
            pdf, species, cohort, main_rows, DOUBLE_ALLELES, rpf_limit, regular, bold
        )
        page_count += 1
        pdf.showPage()
        draw_ms_page(
            pdf,
            species,
            species_ms,
            main_rows,
            MS_GROUPS[species],
            DOUBLE_ALLELES,
            regular,
            bold,
        )
        page_count += 1
    if page_count == 0:
        pdf.setFont(regular, 10)
        pdf.drawString(36, 700, "No double-edit plotting data were found.")
        page_count = 1
    pdf.save()
    return page_count


def create_triple_pdf(
    integrated_rows: list[dict[str, str]],
    ms_rows: list[dict[str, str]],
    out: Path,
    rpf_limit: float,
    regular: str,
    bold: str,
) -> int:
    if not integrated_rows and not ms_rows:
        return 0
    out.parent.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(str(out), pageCompression=1)

    # Triple is deliberately treated separately. If several cohorts exist, use all
    # cohort rows on page 1, grouped implicitly by the cohort in the subtitle/label.
    # Since the current data contain only one triple-edit codon, this remains compact.
    rows = sorted(
        integrated_rows,
        key=lambda r: (
            clean(r.get("cohort")),
            natural_key(clean(r.get("target_handle", ""))),
        ),
    )
    # Reuse the standard page routine cohort by cohort; usually only one WT/species row is relevant.
    cohort = clean(rows[0].get("cohort")) if rows else "triple-edit"
    species = (
        clean(rows[0].get("species"))
        if rows
        else (clean(ms_rows[0].get("species")) if ms_rows else "")
    )
    draw_rna_rpf_page(
        pdf,
        species,
        cohort,
        rows,
        TRIPLE_ALLELES,
        rpf_limit,
        regular,
        bold,
        title_prefix="Triple-edit",
    )
    pdf.showPage()
    groups = tuple(
        sorted({clean(r.get("MS_group")) for r in ms_rows if clean(r.get("MS_group"))})
    ) or ("MS",)
    draw_ms_page(
        pdf,
        species,
        ms_rows,
        rows,
        groups,
        TRIPLE_ALLELES,
        regular,
        bold,
        title_prefix="Triple-edit",
    )
    pdf.save()
    return 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--integrated-dir",
        type=Path,
        default=Path("double_edit_rebuild_v4/04_integrated"),
    )
    parser.add_argument(
        "--ms-dir", type=Path, default=Path("double_edit_rebuild_v4/03_ms")
    )
    parser.add_argument(
        "--outdir", type=Path, default=Path("double_edit_rebuild_v4/05_plots")
    )
    parser.add_argument(
        "--rpf-color-limit",
        type=float,
        default=0.0,
        help="Symmetric +/- limit for RNA/RPF heatmaps. <=0 means auto 95th percentile with minimum 1.0.",
    )
    args = parser.parse_args()

    integrated_dir = args.integrated_dir.resolve()
    ms_dir = args.ms_dir.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    double_integrated = read_tsv(
        integrated_dir / "double_edit_RNA_RPF_MS.display_wide.tsv"
    )
    triple_integrated = read_tsv(
        integrated_dir / "triple_edit_RNA_RPF_MS.display_wide.tsv"
    )
    double_ms = read_tsv(ms_dir / "double_edit_MS_display_wide.tsv")
    triple_ms = read_tsv(ms_dir / "triple_edit_MS_display_wide.tsv")

    score_values: list[float] = []
    for row in double_integrated:
        if clean(row.get("cohort")) != MAIN_COHORT.get(clean(row.get("species")), ""):
            continue
        if clean(row.get("RNA_RPF_heatmap_eligible", "")).lower() != "yes":
            continue
        for allele in DOUBLE_ALLELES:
            v = optional_float(row.get(f"RNA_RPF_{allele}_mean_log2_enrichment"))
            if v is not None:
                score_values.append(abs(v))
    auto_limit = max(1.0, percentile(score_values, 0.95))
    rpf_limit = args.rpf_color_limit if args.rpf_color_limit > 0 else auto_limit

    regular, bold = register_fonts()
    page_counts: dict[str, int] = {}

    for species in SPECIES_ORDER:
        output = outdir / "double_edit" / f"{species}.double_edit_RNA_RPF_MS.vector.pdf"
        n = create_species_pdf(
            species,
            double_integrated,
            double_ms,
            output,
            rpf_limit,
            regular,
            bold,
        )
        if n:
            page_counts[str(output)] = n

    combined = outdir / "double_edit" / "all_species.double_edit_RNA_RPF_MS.vector.pdf"
    page_counts[str(combined)] = create_all_species_pdf(
        double_integrated, double_ms, combined, rpf_limit, regular, bold
    )

    triple_scores: list[float] = []
    for row in triple_integrated:
        for allele in TRIPLE_ALLELES:
            v = optional_float(row.get(f"RNA_RPF_{allele}_mean_log2_enrichment"))
            if v is not None:
                triple_scores.append(abs(v))
    triple_limit = (
        max(1.0, percentile(triple_scores, 0.95)) if triple_scores else rpf_limit
    )
    triple_out = outdir / "triple_edit" / "triple_edit_RNA_RPF_MS.vector.pdf"
    n_triple = create_triple_pdf(
        triple_integrated, triple_ms, triple_out, triple_limit, regular, bold
    )
    if n_triple:
        page_counts[str(triple_out)] = n_triple

    summary_path = outdir / "plot_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["file", "pages", "RNA_RPF_color_limit"])
        for path, pages in page_counts.items():
            writer.writerow([path, pages, f"{rpf_limit:.10g}"])

    print(f"[plot] RNA/RPF symmetric color limit: +/-{rpf_limit:.4g}")
    for path, pages in page_counts.items():
        print(f"[plot] {pages} pages: {path}")
    print(f"[plot] summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
