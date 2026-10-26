#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Normalize joint double-/triple-edit RNA/RPF haplotype counts per biological
sample and summarize them without depth weighting.

For allele a in one sample:
    pRNA(a) = (RNA_a + alpha) / (RNA_total + K*alpha)
    pRPF(a) = (RPF_a + alpha) / (RPF_total + K*alpha)
    score(a) = log2[pRPF(a) / pRNA(a)]

K=4 for double-edit codons and K=8 for triple-edit codons.

Coverage filter is applied per real RNA/RPF sample. Main WT/species cohorts:
    ISO1 + WBg -> Dmel_WT
    Dsim       -> Dsim
    Dsec       -> Dsec
    Dere       -> Dere
AG/AS/Adar5G1 are also summarized separately when present but are not part of
`main_double_edit_heatmap.tsv` unless the code is changed explicitly.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


def clean(x: Any) -> str:
    return "" if x is None else str(x).strip().strip("\ufeff")


def as_float(x: Any) -> float:
    s = clean(x)
    if not s or s.upper() in {"NA", "NAN", "NULL"}:
        return 0.0
    return float(s)


def as_int(x: Any) -> int:
    return int(float(clean(x)))


def expected_patterns(site_count: int) -> list[str]:
    patterns: list[str] = []
    for edited_count in range(site_count + 1):
        for edited_positions in itertools.combinations(range(site_count), edited_count):
            chars = ["A"] * site_count
            for idx in edited_positions:
                chars[idx] = "I"
            patterns.append("".join(chars))
    return patterns


def cohort_for_genotype(genotype: str) -> str:
    g = clean(genotype)
    if g in {"ISO1", "WBg"}:
        return "Dmel_WT"
    if g == "AG":
        return "AG"
    if g == "AS":
        return "AS"
    if g == "Adar5G1":
        return "Adar5G1"
    if g in {"Dsim", "Dsec", "Dere"}:
        return g
    return g


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


def mean_sd_sem(values: list[float]) -> tuple[str, str, str]:
    if not values:
        return "", "", ""
    m = statistics.mean(values)
    if len(values) == 1:
        return f"{m:.10g}", "", ""
    sd = statistics.stdev(values)
    sem = sd / math.sqrt(len(values))
    return f"{m:.10g}", f"{sd:.10g}", f"{sem:.10g}"


def process_table(
    rows: list[dict[str, str]],
    site_count: int,
    alpha: float,
    min_rna_total: float,
    min_rpf_total: float,
    min_valid_samples: int,
):
    patterns = expected_patterns(site_count)

    # One row per sample x codon x allele is expected from step 01.
    by_sample: dict[tuple[str, str, str, int, str], dict[str, dict[str, str]]] = (
        defaultdict(dict)
    )
    meta_by_handle: dict[str, dict[str, str]] = {}
    for row in rows:
        if as_int(row["edit_site_count"]) != site_count:
            continue
        key = (
            row["species"],
            row["genotype"],
            row["sex"],
            as_int(row["rep"]),
            row["target_handle"],
        )
        by_sample[key][row["allele"]] = row
        meta_by_handle.setdefault(row["target_handle"], row)

    sample_score_rows: list[dict[str, Any]] = []
    sample_codon_rows: list[dict[str, Any]] = []

    # Cohort/codon buckets.
    score_bucket: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    editing_bucket: dict[tuple[str, str], list[float]] = defaultdict(list)
    valid_sample_bucket: dict[tuple[str, str], set[tuple[str, str, int]]] = defaultdict(
        set
    )
    rna_detect_bucket: dict[tuple[str, str], bool] = defaultdict(bool)
    rpf_detect_bucket: dict[tuple[str, str], bool] = defaultdict(bool)
    matched_sample_bucket: dict[tuple[str, str], set[tuple[str, str, int]]] = (
        defaultdict(set)
    )

    for (species, genotype, sex, rep, handle), allele_map in sorted(by_sample.items()):
        meta = meta_by_handle[handle]
        cohort = cohort_for_genotype(genotype)

        rna_counts = {
            a: as_float(allele_map.get(a, {}).get("RNA_count", 0)) for a in patterns
        }
        rpf_counts = {
            a: as_float(allele_map.get(a, {}).get("RPF_count", 0)) for a in patterns
        }
        rna_total = sum(rna_counts.values())
        rpf_total = sum(rpf_counts.values())
        pass_cov = rna_total >= min_rna_total and rpf_total >= min_rpf_total
        unedited = "A" * site_count
        edited_rna = rna_total - rna_counts[unedited]
        edited_rpf = rpf_total - rpf_counts[unedited]
        rna_edit_level = (
            (edited_rna / rna_total)
            if rna_total >= min_rna_total and rna_total > 0
            else None
        )

        sample_key = (genotype, sex, rep)
        matched_sample_bucket[(cohort, handle)].add(sample_key)
        if pass_cov:
            valid_sample_bucket[(cohort, handle)].add(sample_key)
            rna_detect_bucket[(cohort, handle)] |= edited_rna > 0
            rpf_detect_bucket[(cohort, handle)] |= edited_rpf > 0
        if rna_edit_level is not None:
            editing_bucket[(cohort, handle)].append(rna_edit_level)

        sample_codon_rows.append(
            {
                "species": species,
                "cohort": cohort,
                "genotype": genotype,
                "sex": sex,
                "rep": rep,
                "sample": (
                    next(iter(allele_map.values())).get("sample", "")
                    if allele_map
                    else ""
                ),
                "target_handle": handle,
                "gene": meta.get("gene", ""),
                "genomic_sites": meta.get("genomic_sites", ""),
                "edit_site_count": site_count,
                "RNA_total": f"{rna_total:.10g}",
                "RPF_total": f"{rpf_total:.10g}",
                "edited_RNA_count": f"{edited_rna:.10g}",
                "edited_RPF_count": f"{edited_rpf:.10g}",
                "pair_pass_coverage": "yes" if pass_cov else "no",
                "RNA_editing_level_exact": (
                    f"{rna_edit_level:.10g}" if rna_edit_level is not None else ""
                ),
            }
        )

        k = len(patterns)
        for allele in patterns:
            p_rna = (rna_counts[allele] + alpha) / (rna_total + k * alpha)
            p_rpf = (rpf_counts[allele] + alpha) / (rpf_total + k * alpha)
            score = math.log2(p_rpf / p_rna)
            if pass_cov:
                score_bucket[(cohort, handle, allele)].append(score)
            sample_score_rows.append(
                {
                    "species": species,
                    "cohort": cohort,
                    "genotype": genotype,
                    "sex": sex,
                    "rep": rep,
                    "sample": (
                        next(iter(allele_map.values())).get("sample", "")
                        if allele_map
                        else ""
                    ),
                    "target_handle": handle,
                    "gene": meta.get("gene", ""),
                    "genomic_sites": meta.get("genomic_sites", ""),
                    "edit_site_count": site_count,
                    "allele": allele,
                    "RNA_count": f"{rna_counts[allele]:.10g}",
                    "RPF_count": f"{rpf_counts[allele]:.10g}",
                    "RNA_total": f"{rna_total:.10g}",
                    "RPF_total": f"{rpf_total:.10g}",
                    "pair_pass_coverage": "yes" if pass_cov else "no",
                    "RNA_fraction_smoothed": f"{p_rna:.10g}",
                    "RPF_fraction_smoothed": f"{p_rpf:.10g}",
                    "log2_RPF_RNA_enrichment": f"{score:.10g}",
                }
            )

    cohort_rows: list[dict[str, Any]] = []
    cohorts_handles = sorted(
        {(c, h) for c, h, _ in score_bucket} | set(matched_sample_bucket)
    )
    for cohort, handle in cohorts_handles:
        meta = meta_by_handle[handle]
        valid_n = len(valid_sample_bucket[(cohort, handle)])
        coverage = valid_n >= min_valid_samples
        rna_detected = rna_detect_bucket[(cohort, handle)]
        rpf_detected = rpf_detect_bucket[(cohort, handle)]
        eligible = coverage and rna_detected and rpf_detected
        edit_vals = editing_bucket[(cohort, handle)]
        edit_mean, edit_sd, edit_sem = mean_sd_sem(edit_vals)

        base: dict[str, Any] = {
            "species": meta.get("species", ""),
            "cohort": cohort,
            "target_handle": handle,
            "target_transcript": meta.get("target_transcript", ""),
            "codon_start_1based": meta.get("codon_start_1based", ""),
            "aa_position_1based": meta.get("aa_position_1based", ""),
            "edit_site_count": site_count,
            "edit_cds_positions": meta.get("edit_cds_positions", ""),
            "edit_codon_positions": meta.get("edit_codon_positions", ""),
            "genomic_sites": meta.get("genomic_sites", ""),
            "gene": meta.get("gene", ""),
            "ms_site_aliases": meta.get("ms_site_aliases", ""),
            "matched_sample_count": len(matched_sample_bucket[(cohort, handle)]),
            "valid_sample_count": valid_n,
            "coverage_qualified": "yes" if coverage else "no",
            "RNA_edited_detected_in_valid_samples": "yes" if rna_detected else "no",
            "RPF_edited_detected_in_valid_samples": "yes" if rpf_detected else "no",
            "heatmap_eligible": "yes" if eligible else "no",
            "RNA_editing_level_mean": edit_mean,
            "RNA_editing_level_sd": edit_sd,
            "RNA_editing_level_sem": edit_sem,
            "RNA_editing_level_n": len(edit_vals),
        }
        for allele in patterns:
            m, sd, sem = mean_sd_sem(score_bucket.get((cohort, handle, allele), []))
            base[f"{allele}_mean_log2_RPF_RNA"] = m
            base[f"{allele}_sd_log2_RPF_RNA"] = sd
            base[f"{allele}_sem_log2_RPF_RNA"] = sem
            base[f"{allele}_n"] = len(score_bucket.get((cohort, handle, allele), []))
        cohort_rows.append(base)

    return sample_score_rows, sample_codon_rows, cohort_rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--counts",
        type=Path,
        default=Path("double_edit_rebuild_v4/01_counts/allele_counts.by_sample.tsv"),
    )
    p.add_argument(
        "--outdir", type=Path, default=Path("double_edit_rebuild_v4/02_normalized")
    )
    p.add_argument("--alpha", type=float, default=0.5)
    p.add_argument("--min-rna-total", type=float, default=5)
    p.add_argument("--min-rpf-total", type=float, default=5)
    p.add_argument("--min-valid-samples", type=int, default=1)
    args = p.parse_args()

    rows = read_tsv(args.counts.resolve())
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    summary = {}
    for n_sites, prefix in ((2, "double_edit"), (3, "triple_edit")):
        sample_scores, sample_codons, cohort_rows = process_table(
            rows,
            site_count=n_sites,
            alpha=args.alpha,
            min_rna_total=args.min_rna_total,
            min_rpf_total=args.min_rpf_total,
            min_valid_samples=args.min_valid_samples,
        )
        patterns = expected_patterns(n_sites)
        sample_score_fields = [
            "species",
            "cohort",
            "genotype",
            "sex",
            "rep",
            "sample",
            "target_handle",
            "gene",
            "genomic_sites",
            "edit_site_count",
            "allele",
            "RNA_count",
            "RPF_count",
            "RNA_total",
            "RPF_total",
            "pair_pass_coverage",
            "RNA_fraction_smoothed",
            "RPF_fraction_smoothed",
            "log2_RPF_RNA_enrichment",
        ]
        sample_codon_fields = [
            "species",
            "cohort",
            "genotype",
            "sex",
            "rep",
            "sample",
            "target_handle",
            "gene",
            "genomic_sites",
            "edit_site_count",
            "RNA_total",
            "RPF_total",
            "edited_RNA_count",
            "edited_RPF_count",
            "pair_pass_coverage",
            "RNA_editing_level_exact",
        ]
        cohort_fields = [
            "species",
            "cohort",
            "target_handle",
            "target_transcript",
            "codon_start_1based",
            "aa_position_1based",
            "edit_site_count",
            "edit_cds_positions",
            "edit_codon_positions",
            "genomic_sites",
            "gene",
            "ms_site_aliases",
            "matched_sample_count",
            "valid_sample_count",
            "coverage_qualified",
            "RNA_edited_detected_in_valid_samples",
            "RPF_edited_detected_in_valid_samples",
            "heatmap_eligible",
            "RNA_editing_level_mean",
            "RNA_editing_level_sd",
            "RNA_editing_level_sem",
            "RNA_editing_level_n",
        ]
        for allele in patterns:
            cohort_fields += [
                f"{allele}_mean_log2_RPF_RNA",
                f"{allele}_sd_log2_RPF_RNA",
                f"{allele}_sem_log2_RPF_RNA",
                f"{allele}_n",
            ]

        write_tsv(
            outdir / f"{prefix}_scores.by_sample.tsv",
            sample_scores,
            sample_score_fields,
        )
        write_tsv(
            outdir / f"{prefix}_codon_qc.by_sample.tsv",
            sample_codons,
            sample_codon_fields,
        )
        write_tsv(
            outdir / f"{prefix}_scores.cohort_mean.tsv", cohort_rows, cohort_fields
        )

        if n_sites == 2:
            main_rows = [
                r
                for r in cohort_rows
                if r["cohort"] in {"Dmel_WT", "Dsim", "Dsec", "Dere"}
                and r["heatmap_eligible"] == "yes"
            ]
            write_tsv(outdir / "main_double_edit_heatmap.tsv", main_rows, cohort_fields)

        summary[prefix] = {
            "sample_score_rows": len(sample_scores),
            "sample_codon_rows": len(sample_codons),
            "cohort_rows": len(cohort_rows),
            "eligible_cohort_rows": sum(
                r["heatmap_eligible"] == "yes" for r in cohort_rows
            ),
        }

    run = {
        "counts": str(args.counts.resolve()),
        "outdir": str(outdir),
        "alpha": args.alpha,
        "min_rna_total": args.min_rna_total,
        "min_rpf_total": args.min_rpf_total,
        "min_valid_samples": args.min_valid_samples,
        **summary,
    }
    (outdir / "analysis_summary.json").write_text(
        json.dumps(run, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(run, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
