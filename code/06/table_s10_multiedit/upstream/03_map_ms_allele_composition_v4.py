#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map precomputed allele-level MS intensities to double-/triple-edit codons and
summarize only within-site MS composition.

Input columns:
    site genotype sex rep sample allele intensity total_intensity proportion

Scientific rule used in v4
--------------------------
Raw `intensity` is NOT used as a cross-library display quantity because the
input intensities have not been library-normalized. Raw intensity is retained
only for mapping/QC and for reconstructing within-site fractions.

For one mapped site/sample:
    allele_proportion(a) = intensity(a) / total_intensity

Protein editing level (Lp):
    double-edit: Lp = (IA + AI + II) intensity / total_intensity = 1 - AA fraction
    triple-edit: Lp = sum(all non-AAA allele intensities) / total_intensity

If total_intensity == 0, Lp and all allele proportions are NA. A missing allele
row is treated as proportion 0 only when the same site/sample has positive
`total_intensity`.

The supplied input `proportion` is retained for QC, while the v4 display value
is recomputed from `intensity / total_intensity` so its definition is explicit.
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

DEFAULT_MS = Path(
    "/data/dongxn/RNAediting/codex/03_MS/DIA_NN_analysis/"
    "site_genotype_rep_allele_intensity.tsv"
)


def clean(x: Any) -> str:
    return "" if x is None else str(x).strip().strip("\ufeff")


def optional_float(x: Any) -> float | None:
    s = clean(x)
    if not s or s.upper() in {"NA", "NAN", "NULL"}:
        return None
    return float(s)


def as_float(x: Any) -> float:
    value = optional_float(x)
    return 0.0 if value is None else value


def as_int(x: Any) -> int:
    return int(float(clean(x)))


def fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.10g}"


def normalize_genotype(value: str) -> str:
    x = clean(value)
    u = x.upper()
    if u in {"ISO", "ISO1", "DML_ISO1", "DMEL_ISO1"}:
        return "ISO1"
    if u in {"WBG", "WBGISO", "WBG_ISO1"}:
        return "WBg"
    if u in {"AG", "ADARG"}:
        return "AG"
    if u in {"AS", "ADARS"}:
        return "AS"
    if u in {"DSIM", "SIM"}:
        return "Dsim"
    if u in {"DSEC", "SEC"}:
        return "Dsec"
    if u in {"DERE", "ERE"}:
        return "Dere"
    return x


def normalize_sex(value: str) -> str:
    x = clean(value).upper()
    if x in {"F", "FEMALE"}:
        return "F"
    if x in {"M", "MALE"}:
        return "M"
    return x


def ms_group(genotype: str) -> str:
    g = normalize_genotype(genotype)
    if g in {"ISO1", "WBg"}:
        return "WBg_ISO1"
    return g


def expected_patterns(site_count: int) -> list[str]:
    patterns: list[str] = []
    for edited_count in range(site_count + 1):
        for edited_positions in itertools.combinations(range(site_count), edited_count):
            chars = ["A"] * site_count
            for idx in edited_positions:
                chars[idx] = "I"
            patterns.append("".join(chars))
    return patterns


def unedited_pattern(site_count: int) -> str:
    return "A" * site_count


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
            writer.writerow({field: row.get(field, "") for field in fields})
            n += 1
    return n


def mean_sd_sem(values: list[float]) -> tuple[str, str, str]:
    if not values:
        return "", "", ""
    mean = statistics.mean(values)
    if len(values) == 1:
        return f"{mean:.10g}", "", ""
    sd = statistics.stdev(values)
    sem = sd / math.sqrt(len(values))
    return f"{mean:.10g}", f"{sd:.10g}", f"{sem:.10g}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path("double_edit_rebuild_v4/01_counts/multiedit_codons.tsv"),
    )
    parser.add_argument("--ms", type=Path, default=DEFAULT_MS)
    parser.add_argument(
        "--outdir", type=Path, default=Path("double_edit_rebuild_v4/03_ms")
    )
    parser.add_argument(
        "--proportion-tolerance",
        type=float,
        default=1e-6,
        help="QC tolerance for supplied proportion vs intensity/total_intensity.",
    )
    parser.add_argument(
        "--total-relative-tolerance",
        type=float,
        default=1e-6,
        help="QC tolerance for sum(mapped allele intensities) vs supplied total_intensity.",
    )
    args = parser.parse_args()

    catalog_path = args.catalog.resolve()
    ms_path = args.ms.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    catalog = read_tsv(catalog_path)
    ms_rows = read_tsv(ms_path)

    required = {
        "site",
        "genotype",
        "sex",
        "rep",
        "sample",
        "allele",
        "intensity",
        "total_intensity",
        "proportion",
    }
    if ms_rows:
        missing = required - set(ms_rows[0])
        if missing:
            raise ValueError(f"MS table missing columns: {sorted(missing)}")

    alias_to_handles: dict[str, set[str]] = defaultdict(set)
    meta_by_handle: dict[str, dict[str, str]] = {}
    for row in catalog:
        handle = row["target_handle"]
        meta_by_handle[handle] = row
        for alias in clean(row.get("ms_site_aliases")).split(";"):
            alias = alias.strip()
            if alias:
                alias_to_handles[alias].add(handle)

    alias_conflicts = [
        {
            "site_alias": alias,
            "target_handles": ";".join(sorted(handles)),
            "n_targets": len(handles),
        }
        for alias, handles in alias_to_handles.items()
        if len(handles) > 1
    ]

    unmatched: list[dict[str, Any]] = []
    incompatible: list[dict[str, Any]] = []
    duplicate_rows: list[dict[str, Any]] = []

    # One retained raw record per target/sample/allele.
    values: dict[tuple[str, str, str, int, str, str], dict[str, Any]] = {}
    site_present_handles: set[str] = set()
    sample_keys_by_handle: dict[str, set[tuple[str, str, int, str]]] = defaultdict(set)

    for row in ms_rows:
        site = clean(row["site"])
        targets = alias_to_handles.get(site, set())
        if len(targets) != 1:
            unmatched.append(
                {
                    **row,
                    "mapping_status": "unmatched" if not targets else "ambiguous_alias",
                    "candidate_target_handles": ";".join(sorted(targets)),
                }
            )
            continue

        handle = next(iter(targets))
        meta = meta_by_handle[handle]
        n_sites = as_int(meta["edit_site_count"])
        allele = clean(row["allele"]).upper()
        expected = set(expected_patterns(n_sites))
        if allele not in expected:
            incompatible.append(
                {
                    **row,
                    "target_handle": handle,
                    "edit_site_count": n_sites,
                    "expected_alleles": ";".join(expected_patterns(n_sites)),
                    "reason": "allele label incompatible with candidate site count",
                }
            )
            continue

        genotype = normalize_genotype(row["genotype"])
        sex = normalize_sex(row["sex"])
        rep = as_int(row["rep"])
        sample = clean(row["sample"])
        key = (handle, genotype, sex, rep, sample, allele)
        value = {
            "source_site": site,
            "intensity": as_float(row["intensity"]),
            "total_intensity": as_float(row["total_intensity"]),
            "input_proportion": optional_float(row["proportion"]),
        }
        if key in values:
            old = values[key]
            duplicate_rows.append(
                {
                    "target_handle": handle,
                    "source_site": site,
                    "genotype": genotype,
                    "sex": sex,
                    "rep": rep,
                    "sample": sample,
                    "allele": allele,
                    "old_intensity": old["intensity"],
                    "new_intensity": value["intensity"],
                    "action": "keep maximum intensity",
                }
            )
            if value["intensity"] > old["intensity"]:
                values[key] = value
        else:
            values[key] = value

        site_present_handles.add(handle)
        sample_keys_by_handle[handle].add((genotype, sex, rep, sample))

    # Index retained rows by site/sample for efficient reconstruction.
    values_by_sample: dict[
        tuple[str, str, str, int, str], dict[str, dict[str, Any]]
    ] = defaultdict(dict)
    for (handle, genotype, sex, rep, sample, allele), value in values.items():
        values_by_sample[(handle, genotype, sex, rep, sample)][allele] = value

    mapped_rep: list[dict[str, Any]] = []
    lp_rep: list[dict[str, Any]] = []
    proportion_qc: list[dict[str, Any]] = []
    total_qc: list[dict[str, Any]] = []

    for handle in sorted(site_present_handles):
        meta = meta_by_handle[handle]
        n_sites = as_int(meta["edit_site_count"])
        patterns = expected_patterns(n_sites)
        unedited = unedited_pattern(n_sites)

        for genotype, sex, rep, sample in sorted(sample_keys_by_handle[handle]):
            sample_key = (handle, genotype, sex, rep, sample)
            allele_values = values_by_sample[sample_key]
            observed_values = list(allele_values.values())
            reported_totals = sorted({v["total_intensity"] for v in observed_values})
            total_intensity = max(reported_totals, default=0.0)
            sum_expected_intensity = sum(
                allele_values.get(allele, {}).get("intensity", 0.0)
                for allele in patterns
            )

            if total_intensity > 0:
                relative_total_diff = (
                    abs(sum_expected_intensity - total_intensity) / total_intensity
                )
                total_status = (
                    "PASS"
                    if relative_total_diff <= args.total_relative_tolerance
                    else "CHECK"
                )
            elif sum_expected_intensity == 0:
                relative_total_diff = 0.0
                total_status = "NO_SIGNAL"
            else:
                relative_total_diff = None
                total_status = "CHECK_TOTAL_ZERO_WITH_POSITIVE_ALLELES"

            total_qc.append(
                {
                    "species": meta["species"],
                    "MS_group": ms_group(genotype),
                    "genotype": genotype,
                    "sex": sex,
                    "rep": rep,
                    "sample": sample,
                    "target_handle": handle,
                    "edit_site_count": n_sites,
                    "reported_total_min": fmt(
                        min(reported_totals) if reported_totals else None
                    ),
                    "reported_total_max": fmt(
                        max(reported_totals) if reported_totals else None
                    ),
                    "n_distinct_reported_totals": len(reported_totals),
                    "sum_expected_allele_intensity": fmt(sum_expected_intensity),
                    "relative_difference_sum_vs_total": fmt(relative_total_diff),
                    "status": total_status,
                }
            )

            proportions: dict[str, float | None] = {}
            edited_intensity = 0.0
            for allele in patterns:
                raw = allele_values.get(allele)
                if raw is None:
                    intensity = 0.0
                    input_prop = None
                    source_site = ""
                    observed_flag = "no"
                else:
                    intensity = raw["intensity"]
                    input_prop = raw["input_proportion"]
                    source_site = raw["source_site"]
                    observed_flag = "yes"

                computed_prop = (
                    intensity / total_intensity if total_intensity > 0 else None
                )
                proportions[allele] = computed_prop
                if allele != unedited:
                    edited_intensity += intensity

                prop_diff = (
                    computed_prop - input_prop
                    if computed_prop is not None and input_prop is not None
                    else None
                )
                prop_qc_status = (
                    "PASS"
                    if prop_diff is not None
                    and abs(prop_diff) <= args.proportion_tolerance
                    else "CHECK" if prop_diff is not None else "NA"
                )
                proportion_qc.append(
                    {
                        "species": meta["species"],
                        "MS_group": ms_group(genotype),
                        "genotype": genotype,
                        "sex": sex,
                        "rep": rep,
                        "sample": sample,
                        "target_handle": handle,
                        "allele": allele,
                        "input_proportion": fmt(input_prop),
                        "computed_intensity_proportion": fmt(computed_prop),
                        "difference_computed_minus_input": fmt(prop_diff),
                        "status": prop_qc_status,
                    }
                )

                mapped_rep.append(
                    {
                        "species": meta["species"],
                        "MS_group": ms_group(genotype),
                        "genotype": genotype,
                        "sex": sex,
                        "rep": rep,
                        "sample": sample,
                        "target_handle": handle,
                        "source_site": source_site,
                        "edit_site_count": n_sites,
                        "gene": meta.get("gene", ""),
                        "genomic_sites": meta.get("genomic_sites", ""),
                        "allele": allele,
                        "allele_row_present_in_input": observed_flag,
                        "raw_intensity_for_QC": fmt(intensity),
                        "raw_total_intensity_for_QC": fmt(total_intensity),
                        "input_proportion_for_QC": fmt(input_prop),
                        "intensity_proportion": fmt(computed_prop),
                        "MS_signal_available": "yes" if total_intensity > 0 else "no",
                    }
                )

            lp = edited_intensity / total_intensity if total_intensity > 0 else None
            # Mathematically equivalent check against 1 - unedited proportion.
            lp_from_unedited = (
                1.0 - proportions[unedited]
                if proportions.get(unedited) is not None
                else None
            )
            lp_diff = (
                lp - lp_from_unedited
                if lp is not None and lp_from_unedited is not None
                else None
            )
            lp_rep.append(
                {
                    "species": meta["species"],
                    "MS_group": ms_group(genotype),
                    "genotype": genotype,
                    "sex": sex,
                    "rep": rep,
                    "sample": sample,
                    "target_handle": handle,
                    "edit_site_count": n_sites,
                    "gene": meta.get("gene", ""),
                    "genomic_sites": meta.get("genomic_sites", ""),
                    "unedited_allele": unedited,
                    "raw_total_intensity_for_QC": fmt(total_intensity),
                    "raw_edited_intensity_for_QC": fmt(edited_intensity),
                    "Lp": fmt(lp),
                    "Lp_from_1_minus_unedited_proportion": fmt(lp_from_unedited),
                    "Lp_internal_difference": fmt(lp_diff),
                    "MS_signal_available": "yes" if total_intensity > 0 else "no",
                }
            )

    # Allele proportion summary: only samples with positive total intensity contribute.
    allele_bucket: dict[tuple[str, str, str], dict[str, list[float]]] = defaultdict(
        lambda: {"proportion": [], "detected": []}
    )
    for row in mapped_rep:
        prop = optional_float(row["intensity_proportion"])
        if prop is None:
            continue
        key = (row["MS_group"], row["target_handle"], row["allele"])
        allele_bucket[key]["proportion"].append(prop)
        allele_bucket[key]["detected"].append(
            1.0 if as_float(row["raw_intensity_for_QC"]) > 0 else 0.0
        )

    # Lp summary: only samples with positive total intensity contribute.
    lp_bucket: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in lp_rep:
        lp = optional_float(row["Lp"])
        if lp is not None:
            lp_bucket[(row["MS_group"], row["target_handle"])].append(lp)

    allele_summary: list[dict[str, Any]] = []
    for (group, handle, allele), bucket in sorted(allele_bucket.items()):
        meta = meta_by_handle[handle]
        mean_prop, sd_prop, sem_prop = mean_sd_sem(bucket["proportion"])
        lp_mean, lp_sd, lp_sem = mean_sd_sem(lp_bucket[(group, handle)])
        allele_summary.append(
            {
                "species": meta["species"],
                "MS_group": group,
                "target_handle": handle,
                "edit_site_count": meta["edit_site_count"],
                "gene": meta.get("gene", ""),
                "genomic_sites": meta.get("genomic_sites", ""),
                "allele": allele,
                "n_MS_samples_with_signal": len(bucket["proportion"]),
                "detected_sample_count": int(sum(bucket["detected"])),
                "mean_intensity_proportion": mean_prop,
                "sd_intensity_proportion": sd_prop,
                "sem_intensity_proportion": sem_prop,
                "Lp_mean": lp_mean,
                "Lp_sd": lp_sd,
                "Lp_sem": lp_sem,
                "Lp_n": len(lp_bucket[(group, handle)]),
            }
        )

    site_summary: list[dict[str, Any]] = []
    groups_handles = sorted({(r["MS_group"], r["target_handle"]) for r in mapped_rep})
    allele_summary_idx = {
        (r["MS_group"], r["target_handle"], r["allele"]): r for r in allele_summary
    }
    for group, handle in groups_handles:
        meta = meta_by_handle[handle]
        n_sites = as_int(meta["edit_site_count"])
        patterns = expected_patterns(n_sites)
        lp_mean, lp_sd, lp_sem = mean_sd_sem(lp_bucket[(group, handle)])
        row: dict[str, Any] = {
            "species": meta["species"],
            "MS_group": group,
            "target_handle": handle,
            "edit_site_count": n_sites,
            "gene": meta.get("gene", ""),
            "genomic_sites": meta.get("genomic_sites", ""),
            "Lp_mean": lp_mean,
            "Lp_sd": lp_sd,
            "Lp_sem": lp_sem,
            "Lp_n": len(lp_bucket[(group, handle)]),
        }
        for allele in patterns:
            summary = allele_summary_idx.get((group, handle, allele), {})
            row[f"{allele}_mean_intensity_proportion"] = summary.get(
                "mean_intensity_proportion", ""
            )
            row[f"{allele}_sd_intensity_proportion"] = summary.get(
                "sd_intensity_proportion", ""
            )
            row[f"{allele}_sem_intensity_proportion"] = summary.get(
                "sem_intensity_proportion", ""
            )
            row[f"{allele}_detected_sample_count"] = summary.get(
                "detected_sample_count", ""
            )
        site_summary.append(row)

    rep_fields = [
        "species",
        "MS_group",
        "genotype",
        "sex",
        "rep",
        "sample",
        "target_handle",
        "source_site",
        "edit_site_count",
        "gene",
        "genomic_sites",
        "allele",
        "allele_row_present_in_input",
        "raw_intensity_for_QC",
        "raw_total_intensity_for_QC",
        "input_proportion_for_QC",
        "intensity_proportion",
        "MS_signal_available",
    ]
    lp_rep_fields = [
        "species",
        "MS_group",
        "genotype",
        "sex",
        "rep",
        "sample",
        "target_handle",
        "edit_site_count",
        "gene",
        "genomic_sites",
        "unedited_allele",
        "raw_total_intensity_for_QC",
        "raw_edited_intensity_for_QC",
        "Lp",
        "Lp_from_1_minus_unedited_proportion",
        "Lp_internal_difference",
        "MS_signal_available",
    ]
    allele_summary_fields = [
        "species",
        "MS_group",
        "target_handle",
        "edit_site_count",
        "gene",
        "genomic_sites",
        "allele",
        "n_MS_samples_with_signal",
        "detected_sample_count",
        "mean_intensity_proportion",
        "sd_intensity_proportion",
        "sem_intensity_proportion",
        "Lp_mean",
        "Lp_sd",
        "Lp_sem",
        "Lp_n",
    ]

    write_tsv(outdir / "MS_allele_composition.by_sample.tsv", mapped_rep, rep_fields)
    write_tsv(outdir / "MS_Lp.by_sample.tsv", lp_rep, lp_rep_fields)
    write_tsv(
        outdir / "double_edit_MS_allele_summary.tsv",
        [r for r in allele_summary if as_int(r["edit_site_count"]) == 2],
        allele_summary_fields,
    )
    write_tsv(
        outdir / "triple_edit_MS_allele_summary.tsv",
        [r for r in allele_summary if as_int(r["edit_site_count"]) == 3],
        allele_summary_fields,
    )
    write_tsv(
        outdir / "MS_allele_summary.all.tsv", allele_summary, allele_summary_fields
    )

    # Wide plotting-ready MS tables: allele proportions + one Lp column.
    for n_sites, prefix in ((2, "double_edit"), (3, "triple_edit")):
        patterns = expected_patterns(n_sites)
        rows = [r for r in site_summary if as_int(r["edit_site_count"]) == n_sites]
        fields = [
            "species",
            "MS_group",
            "target_handle",
            "edit_site_count",
            "gene",
            "genomic_sites",
            "Lp_mean",
            "Lp_sd",
            "Lp_sem",
            "Lp_n",
        ]
        for allele in patterns:
            fields.extend(
                [
                    f"{allele}_mean_intensity_proportion",
                    f"{allele}_sd_intensity_proportion",
                    f"{allele}_sem_intensity_proportion",
                    f"{allele}_detected_sample_count",
                ]
            )
        write_tsv(outdir / f"{prefix}_MS_display_wide.tsv", rows, fields)

    unmatched_fields = (
        list(ms_rows[0].keys()) + ["mapping_status", "candidate_target_handles"]
        if ms_rows
        else [
            "site",
            "genotype",
            "sex",
            "rep",
            "sample",
            "allele",
            "intensity",
            "total_intensity",
            "proportion",
            "mapping_status",
            "candidate_target_handles",
        ]
    )
    write_tsv(outdir / "MS_unmatched_rows.tsv", unmatched, unmatched_fields)

    incompatible_fields = (
        list(ms_rows[0].keys())
        + ["target_handle", "edit_site_count", "expected_alleles", "reason"]
        if ms_rows
        else []
    )
    write_tsv(
        outdir / "MS_incompatible_allele_rows.tsv", incompatible, incompatible_fields
    )
    write_tsv(
        outdir / "MS_alias_conflicts.tsv",
        alias_conflicts,
        ["site_alias", "target_handles", "n_targets"],
    )
    write_tsv(
        outdir / "MS_duplicate_rows.tsv",
        duplicate_rows,
        [
            "target_handle",
            "source_site",
            "genotype",
            "sex",
            "rep",
            "sample",
            "allele",
            "old_intensity",
            "new_intensity",
            "action",
        ],
    )
    write_tsv(
        outdir / "MS_proportion_QC.tsv",
        proportion_qc,
        [
            "species",
            "MS_group",
            "genotype",
            "sex",
            "rep",
            "sample",
            "target_handle",
            "allele",
            "input_proportion",
            "computed_intensity_proportion",
            "difference_computed_minus_input",
            "status",
        ],
    )
    write_tsv(
        outdir / "MS_total_intensity_QC.tsv",
        total_qc,
        [
            "species",
            "MS_group",
            "genotype",
            "sex",
            "rep",
            "sample",
            "target_handle",
            "edit_site_count",
            "reported_total_min",
            "reported_total_max",
            "n_distinct_reported_totals",
            "sum_expected_allele_intensity",
            "relative_difference_sum_vs_total",
            "status",
        ],
    )

    catalog_status = [
        {
            "species": row["species"],
            "target_handle": row["target_handle"],
            "edit_site_count": row["edit_site_count"],
            "ms_site_aliases": row.get("ms_site_aliases", ""),
            "present_in_MS_input": (
                "yes" if row["target_handle"] in site_present_handles else "no"
            ),
        }
        for row in catalog
    ]
    write_tsv(
        outdir / "MS_candidate_coverage.tsv",
        catalog_status,
        [
            "species",
            "target_handle",
            "edit_site_count",
            "ms_site_aliases",
            "present_in_MS_input",
        ],
    )

    run = {
        "catalog": str(catalog_path),
        "ms": str(ms_path),
        "outdir": str(outdir),
        "catalog_codons": len(catalog),
        "candidate_codons_present_in_MS_input": len(site_present_handles),
        "mapped_MS_sample_allele_rows": len(mapped_rep),
        "Lp_sample_rows": len(lp_rep),
        "unmatched_input_rows": len(unmatched),
        "incompatible_allele_rows": len(incompatible),
        "alias_conflicts": len(alias_conflicts),
        "duplicate_input_keys": len(duplicate_rows),
        "proportion_QC_check_rows": sum(r["status"] == "CHECK" for r in proportion_qc),
        "total_intensity_QC_check_rows": sum(
            r["status"].startswith("CHECK") for r in total_qc
        ),
        "MS_display_definition": "allele intensity / within-site total intensity",
        "Lp_definition_double": "(IA+AI+II)/total_intensity = 1-AA proportion",
        "Lp_definition_triple": "sum(non-AAA allele intensities)/total_intensity = 1-AAA proportion",
        "note": "raw intensity is retained only for QC; it is not used as a cross-library display metric",
    }
    (outdir / "analysis_summary.json").write_text(
        json.dumps(run, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(run, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
