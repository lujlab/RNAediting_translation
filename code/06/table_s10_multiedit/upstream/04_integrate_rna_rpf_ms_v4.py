#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Join RNA/RPF joint-allele enrichment with MS composition metrics.

v4 MS display metrics:
- per-allele mean intensity proportion = intensity / total_intensity
- site-level Lp = total edited protein intensity / total_intensity

Raw unnormalized MS intensity is intentionally excluded from the final
integrated display tables.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any, Iterable


def clean(x: Any) -> str:
    return "" if x is None else str(x).strip().strip("\ufeff")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({f: row.get(f, "") for f in fields})
            n += 1
    return n


def ms_group_for_cohort(cohort: str) -> str:
    return "WBg_ISO1" if cohort == "Dmel_WT" else cohort


def allele_names_from_rna_row(row: dict[str, str]) -> list[str]:
    suffix = "_mean_log2_RPF_RNA"
    return [key[: -len(suffix)] for key in row if key.endswith(suffix)]


def process_long(rna_path: Path, ms_allele_path: Path, out: Path) -> None:
    rna_rows = read_tsv(rna_path)
    ms_rows = read_tsv(ms_allele_path)
    ms_idx = {(r["MS_group"], r["target_handle"], r["allele"]): r for r in ms_rows}

    rows: list[dict[str, Any]] = []
    for rna in rna_rows:
        cohort = rna["cohort"]
        group = ms_group_for_cohort(cohort)
        site_count = int(float(rna["edit_site_count"]))
        for allele in allele_names_from_rna_row(rna):
            ms = ms_idx.get((group, rna["target_handle"], allele), {})
            rows.append(
                {
                    "species": rna.get("species", ""),
                    "cohort": cohort,
                    "MS_group": group,
                    "target_handle": rna["target_handle"],
                    "edit_site_count": site_count,
                    "edit_cds_positions": rna.get("edit_cds_positions", ""),
                    "edit_codon_positions": rna.get("edit_codon_positions", ""),
                    "gene": rna.get("gene", ""),
                    "genomic_sites": rna.get("genomic_sites", ""),
                    "allele": allele,
                    "RNA_RPF_heatmap_eligible": rna.get("heatmap_eligible", ""),
                    "valid_RNA_RPF_sample_count": rna.get("valid_sample_count", ""),
                    # This is the RNA-level editing fraction and is useful beside Lp.
                    "RNA_editing_level_mean": rna.get("RNA_editing_level_mean", ""),
                    "RNA_editing_level_sd": rna.get("RNA_editing_level_sd", ""),
                    "RNA_editing_level_sem": rna.get("RNA_editing_level_sem", ""),
                    "mean_log2_RPF_RNA": rna.get(f"{allele}_mean_log2_RPF_RNA", ""),
                    "sd_log2_RPF_RNA": rna.get(f"{allele}_sd_log2_RPF_RNA", ""),
                    "sem_log2_RPF_RNA": rna.get(f"{allele}_sem_log2_RPF_RNA", ""),
                    "RNA_RPF_n": rna.get(f"{allele}_n", ""),
                    "MS_available": "yes" if ms else "no",
                    "MS_n_samples_with_signal": ms.get("n_MS_samples_with_signal", ""),
                    "MS_detected_sample_count": ms.get("detected_sample_count", ""),
                    "MS_mean_intensity_proportion": ms.get(
                        "mean_intensity_proportion", ""
                    ),
                    "MS_sd_intensity_proportion": ms.get("sd_intensity_proportion", ""),
                    "MS_sem_intensity_proportion": ms.get(
                        "sem_intensity_proportion", ""
                    ),
                    "Lp_mean": ms.get("Lp_mean", ""),
                    "Lp_sd": ms.get("Lp_sd", ""),
                    "Lp_sem": ms.get("Lp_sem", ""),
                    "Lp_n": ms.get("Lp_n", ""),
                }
            )

    fields = [
        "species",
        "cohort",
        "MS_group",
        "target_handle",
        "edit_site_count",
        "edit_cds_positions",
        "edit_codon_positions",
        "gene",
        "genomic_sites",
        "allele",
        "RNA_RPF_heatmap_eligible",
        "valid_RNA_RPF_sample_count",
        "RNA_editing_level_mean",
        "RNA_editing_level_sd",
        "RNA_editing_level_sem",
        "mean_log2_RPF_RNA",
        "sd_log2_RPF_RNA",
        "sem_log2_RPF_RNA",
        "RNA_RPF_n",
        "MS_available",
        "MS_n_samples_with_signal",
        "MS_detected_sample_count",
        "MS_mean_intensity_proportion",
        "MS_sd_intensity_proportion",
        "MS_sem_intensity_proportion",
        "Lp_mean",
        "Lp_sd",
        "Lp_sem",
        "Lp_n",
    ]
    write_tsv(out, rows, fields)


def process_wide(
    rna_path: Path,
    ms_wide_path: Path,
    out: Path,
    site_count: int,
) -> None:
    rna_rows = read_tsv(rna_path)
    ms_rows = read_tsv(ms_wide_path)
    ms_idx = {(r["MS_group"], r["target_handle"]): r for r in ms_rows}

    rows: list[dict[str, Any]] = []
    allele_order = (
        ["AA", "IA", "AI", "II"]
        if site_count == 2
        else ["AAA", "IAA", "AIA", "AAI", "IIA", "IAI", "AII", "III"]
    )

    for rna in rna_rows:
        if int(float(rna["edit_site_count"])) != site_count:
            continue
        cohort = rna["cohort"]
        group = ms_group_for_cohort(cohort)
        ms = ms_idx.get((group, rna["target_handle"]), {})
        row: dict[str, Any] = {
            "species": rna.get("species", ""),
            "cohort": cohort,
            "MS_group": group,
            "target_handle": rna["target_handle"],
            "edit_site_count": site_count,
            "edit_cds_positions": rna.get("edit_cds_positions", ""),
            "edit_codon_positions": rna.get("edit_codon_positions", ""),
            "gene": rna.get("gene", ""),
            "genomic_sites": rna.get("genomic_sites", ""),
            "RNA_RPF_heatmap_eligible": rna.get("heatmap_eligible", ""),
            "RNA_editing_level_mean": rna.get("RNA_editing_level_mean", ""),
            "Lp_mean": ms.get("Lp_mean", ""),
            "Lp_sd": ms.get("Lp_sd", ""),
            "Lp_sem": ms.get("Lp_sem", ""),
            "Lp_n": ms.get("Lp_n", ""),
        }
        for allele in allele_order:
            row[f"RNA_RPF_{allele}_mean_log2_enrichment"] = rna.get(
                f"{allele}_mean_log2_RPF_RNA", ""
            )
            row[f"MS_{allele}_mean_intensity_proportion"] = ms.get(
                f"{allele}_mean_intensity_proportion", ""
            )
        rows.append(row)

    fields = [
        "species",
        "cohort",
        "MS_group",
        "target_handle",
        "edit_site_count",
        "edit_cds_positions",
        "edit_codon_positions",
        "gene",
        "genomic_sites",
        "RNA_RPF_heatmap_eligible",
        "RNA_editing_level_mean",
    ]
    for allele in allele_order:
        fields.append(f"RNA_RPF_{allele}_mean_log2_enrichment")
    fields += ["Lp_mean", "Lp_sd", "Lp_sem", "Lp_n"]
    for allele in allele_order:
        fields.append(f"MS_{allele}_mean_intensity_proportion")

    write_tsv(out, rows, fields)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--normalized-dir",
        type=Path,
        default=Path("double_edit_rebuild_v4/02_normalized"),
    )
    parser.add_argument(
        "--ms-dir", type=Path, default=Path("double_edit_rebuild_v4/03_ms")
    )
    parser.add_argument(
        "--outdir", type=Path, default=Path("double_edit_rebuild_v4/04_integrated")
    )
    args = parser.parse_args()

    normalized_dir = args.normalized_dir.resolve()
    ms_dir = args.ms_dir.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    process_long(
        normalized_dir / "double_edit_scores.cohort_mean.tsv",
        ms_dir / "double_edit_MS_allele_summary.tsv",
        outdir / "double_edit_RNA_RPF_MS.integrated.long.tsv",
    )
    process_long(
        normalized_dir / "triple_edit_scores.cohort_mean.tsv",
        ms_dir / "triple_edit_MS_allele_summary.tsv",
        outdir / "triple_edit_RNA_RPF_MS.integrated.long.tsv",
    )
    process_wide(
        normalized_dir / "double_edit_scores.cohort_mean.tsv",
        ms_dir / "double_edit_MS_display_wide.tsv",
        outdir / "double_edit_RNA_RPF_MS.display_wide.tsv",
        2,
    )
    process_wide(
        normalized_dir / "triple_edit_scores.cohort_mean.tsv",
        ms_dir / "triple_edit_MS_display_wide.tsv",
        outdir / "triple_edit_RNA_RPF_MS.display_wide.tsv",
        3,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
