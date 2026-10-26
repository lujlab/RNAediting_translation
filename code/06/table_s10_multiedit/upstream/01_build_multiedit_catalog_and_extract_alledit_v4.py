#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build double-/triple-edit codon candidates from longest-projection candidate tables
and extract joint RNA/RPF haplotype counts from *.Asite_df_alledit.

Main design
-----------
1) Candidate tables ONLY define which A-to-I sites belong to the same codon.
2) All four species use the same longest-projection coordinate columns:
   target_region + new_transcript + new_cds_pos_1based.
3) transcript_id + hgvs_c are retained only as source/original-coordinate aliases
   for downstream MS site matching; they do NOT define codon membership.
4) For each target codon, editing sites are ordered by CDS coordinate.
   - 2 sites: AA, IA, AI, II
   - 3 sites: AAA, IAA, AIA, AAI, IIA, IAI, AII, III
5) RNA/RPF counts come from Asite_df_alledit, not from raw *.num files.
   Each exact observed codon can occur in several overlapping transition rows;
   therefore, for the same exact codon and assay, the maximum reported count is
   retained before converting codons to A/I haplotype labels. This prevents the
   same haplotype from being counted several times.
6) The one/few triple-edit codons are written to separate output tables.

The script writes extensive QC so count conflicts and coordinate mismatches can
be checked before normalization.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

SPECIES = ("mel", "sim", "sec", "ere")

DEFAULT_DML = Path(
    "/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/"
    "01_rebuild/result/sites_candidate_rescue_3R_18907675_all_test2/"
    "raw_rows_by_group/ISO1.candidate.AtoG.annotations.long.tsv.new"
)
DEFAULT_PROJECTION_ROOT = Path(
    "/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/"
    "03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed"
)
DEFAULT_ALLEDIT_DIR = Path(
    "/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/"
    "04_alledit/result/alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed_newdata/results"
)


def clean(x: Any) -> str:
    return "" if x is None else str(x).strip().strip("\ufeff")


def natural_key(value: str) -> list[Any]:
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", value)]


def as_int(value: Any) -> int | None:
    s = clean(value)
    if not s:
        return None
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return None


def as_float(value: Any) -> float:
    s = clean(value)
    if not s or s.upper() in {"NA", "NAN", "NULL"}:
        return 0.0
    try:
        return float(s)
    except ValueError:
        return 0.0


def boolish(value: Any) -> bool:
    return clean(value).lower() in {"1", "true", "t", "yes", "y"}


def codon_start_1based(cds_pos: int) -> int:
    return ((cds_pos - 1) // 3) * 3 + 1


def aa_pos_1based(cds_pos: int) -> int:
    return ((cds_pos - 1) // 3) + 1


def parse_coding_hgvs(value: str) -> int | None:
    """Accept coding A>G coordinates such as c.313A>G; reject UTR c.* / c.-."""
    m = re.search(r"(?:^|\|)c\.(\d+)A>G(?:$|\|)", clean(value), flags=re.IGNORECASE)
    if not m:
        m = re.fullmatch(r"c\.(\d+)A>G", clean(value), flags=re.IGNORECASE)
    return int(m.group(1)) if m else None


def expected_patterns(site_count: int) -> list[str]:
    patterns: list[str] = []
    for edited_count in range(site_count + 1):
        for edited_positions in itertools.combinations(range(site_count), edited_count):
            chars = ["A"] * site_count
            for idx in edited_positions:
                chars[idx] = "I"
            patterns.append("".join(chars))
    return patterns


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


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        return fields, list(reader)


def candidate_row_to_site(species: str, row: dict[str, str]) -> dict[str, Any] | None:
    """Convert one longest-projection candidate row to one target CDS editing site.

    Codon membership is defined identically for mel/sim/sec/ere by the projected
    longest-transcript coordinate:
        target_region == CDS
        new_transcript
        new_cds_pos_1based

    transcript_id + hgvs_c are kept only as source/original-coordinate aliases
    for downstream MS matching and auditing.
    """
    if clean(row.get("change2")) and clean(row.get("change2")).upper() != "A>G":
        return None

    if clean(row.get("target_region")).upper() != "CDS":
        return None

    target_tx = clean(row.get("new_transcript"))
    target_cds = as_int(row.get("new_cds_pos_1based"))
    if not target_tx or target_cds is None or target_cds < 1:
        return None

    # Source/original annotation is useful for MS aliases, but it does not define
    # which sites are grouped into a codon in this analysis.
    source_tx = clean(row.get("transcript_id")) or target_tx
    source_cds = parse_coding_hgvs(clean(row.get("hgvs_c")))
    if source_cds is None and source_tx == target_tx:
        source_cds = target_cds

    genomic = clean(row.get("species_coord")) or clean(row.get("handle"))
    gene = clean(row.get("gene"))
    ann = clean(row.get("ann"))
    lib = clean(row.get("lib"))

    return {
        "species": species,
        "target_transcript": target_tx,
        "target_cds_pos_1based": target_cds,
        "target_codon_start_1based": codon_start_1based(target_cds),
        "target_handle": f"{target_tx}_{codon_start_1based(target_cds)}",
        "source_transcript": source_tx,
        "source_cds_pos_1based": source_cds if source_cds is not None else "",
        "source_ms_site": (
            f"{source_tx}:{aa_pos_1based(source_cds)}"
            if source_tx and source_cds is not None
            else ""
        ),
        "target_ms_site": f"{target_tx}:{aa_pos_1based(target_cds)}",
        "genomic_site": genomic,
        "gene": gene,
        "ann": ann,
        "lib_info": lib,
    }


def build_candidate_catalog(candidate_paths: dict[str, Path], outdir: Path):
    site_rows_by_species: dict[str, list[dict[str, Any]]] = {}
    candidate_qc: list[dict[str, Any]] = []

    for species in SPECIES:
        path = candidate_paths[species]
        if not path.is_file():
            raise FileNotFoundError(path)
        fields, rows = read_tsv(path)

        converted: list[dict[str, Any]] = []
        for row in rows:
            site = candidate_row_to_site(species, row)
            if site is not None:
                converted.append(site)

        # Deduplicate repeated biological sites across libraries/rows while preserving aliases.
        grouped: dict[tuple[str, int], dict[str, Any]] = {}
        for site in converted:
            key = (site["target_transcript"], int(site["target_cds_pos_1based"]))
            rec = grouped.setdefault(
                key,
                {
                    "species": species,
                    "target_transcript": key[0],
                    "target_cds_pos_1based": key[1],
                    "target_codon_start_1based": codon_start_1based(key[1]),
                    "target_handle": f"{key[0]}_{codon_start_1based(key[1])}",
                    "genomic_sites": set(),
                    "genes": set(),
                    "source_ms_sites": set(),
                    "target_ms_sites": set(),
                    "source_coordinates": set(),
                    "source_row_count": 0,
                    "candidate_lib_info": set(),
                },
            )
            if site["genomic_site"]:
                rec["genomic_sites"].add(site["genomic_site"])
            if site["gene"]:
                rec["genes"].add(site["gene"])
            if site["source_ms_site"]:
                rec["source_ms_sites"].add(site["source_ms_site"])
            rec["target_ms_sites"].add(site["target_ms_site"])
            if site["source_transcript"] and site["source_cds_pos_1based"] != "":
                rec["source_coordinates"].add(
                    f"{site['source_transcript']}_{site['source_cds_pos_1based']}"
                )
            if site["lib_info"]:
                rec["candidate_lib_info"].add(site["lib_info"])
            rec["source_row_count"] += 1

        site_rows = []
        for rec in grouped.values():
            site_rows.append(
                {
                    "species": rec["species"],
                    "target_transcript": rec["target_transcript"],
                    "target_cds_pos_1based": rec["target_cds_pos_1based"],
                    "target_codon_start_1based": rec["target_codon_start_1based"],
                    "target_handle": rec["target_handle"],
                    "genomic_sites": ";".join(
                        sorted(rec["genomic_sites"], key=natural_key)
                    ),
                    "genes": ";".join(sorted(rec["genes"], key=natural_key)),
                    "source_ms_sites": ";".join(
                        sorted(rec["source_ms_sites"], key=natural_key)
                    ),
                    "target_ms_sites": ";".join(
                        sorted(rec["target_ms_sites"], key=natural_key)
                    ),
                    "source_coordinates": ";".join(
                        sorted(rec["source_coordinates"], key=natural_key)
                    ),
                    "source_row_count": rec["source_row_count"],
                    "candidate_lib_info": ";".join(
                        sorted(rec["candidate_lib_info"], key=natural_key)
                    ),
                }
            )
        site_rows.sort(
            key=lambda r: (
                natural_key(r["target_transcript"]),
                int(r["target_cds_pos_1based"]),
            )
        )
        site_rows_by_species[species] = site_rows
        candidate_qc.append(
            {
                "species": species,
                "candidate_file": str(path),
                "input_rows": len(rows),
                "coding_AtoG_rows_parsed": len(converted),
                "unique_target_CDS_sites": len(site_rows),
                "input_columns": ",".join(fields),
            }
        )

    # Group unique sites by target transcript + codon start.
    codons: list[dict[str, Any]] = []
    site_output: list[dict[str, Any]] = []
    for species in SPECIES:
        by_codon: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
        for site in site_rows_by_species[species]:
            by_codon[
                (site["target_transcript"], int(site["target_codon_start_1based"]))
            ].append(site)

        for (tx, start), sites in by_codon.items():
            sites.sort(key=lambda r: int(r["target_cds_pos_1based"]))
            if len(sites) < 2:
                continue
            if len(sites) > 3:
                raise ValueError(
                    f"{species} {tx}_{start}: found {len(sites)} distinct editing sites in one codon"
                )

            offsets = [int(s["target_cds_pos_1based"]) - start + 1 for s in sites]
            if any(x not in (1, 2, 3) for x in offsets):
                raise ValueError(
                    f"Bad within-codon offsets for {species} {tx}_{start}: {offsets}"
                )

            # MS site aliases: the target/longest-transcript alias always exists.
            # For source/original-transcript aliases, require an alias to be shared by
            # ALL editing sites in the codon. This avoids mapping a peptide-site alias
            # that represents only one of the two editing nucleotides.
            target_alias = f"{tx}:{aa_pos_1based(start)}"
            source_sets = [
                {x for x in s["source_ms_sites"].split(";") if x} for s in sites
            ]
            nonempty_source_sets = [x for x in source_sets if x]
            common_source_aliases = (
                set.intersection(*nonempty_source_sets)
                if len(nonempty_source_sets) == len(sites) and nonempty_source_sets
                else set()
            )
            aliases: set[str] = {target_alias} | common_source_aliases
            genes: set[str] = set()
            genomic: set[str] = set()
            source_coords: set[str] = set()
            for s in sites:
                genes.update(x for x in s["genes"].split(";") if x)
                genomic.update(x for x in s["genomic_sites"].split(";") if x)
                source_coords.update(x for x in s["source_coordinates"].split(";") if x)

            n_sites = len(sites)
            rec = {
                "species": species,
                "target_handle": f"{tx}_{start}",
                "target_transcript": tx,
                "codon_start_1based": start,
                "aa_position_1based": aa_pos_1based(start),
                "edit_site_count": n_sites,
                "edit_cds_positions": ";".join(
                    str(s["target_cds_pos_1based"]) for s in sites
                ),
                "edit_codon_positions": ";".join(str(x) for x in offsets),
                "pos1_cds": sites[0]["target_cds_pos_1based"],
                "pos2_cds": sites[1]["target_cds_pos_1based"],
                "pos3_cds": sites[2]["target_cds_pos_1based"] if n_sites == 3 else "",
                "genomic_sites": ";".join(sorted(genomic, key=natural_key)),
                "gene": ";".join(sorted(genes, key=natural_key)),
                "source_coordinates": ";".join(sorted(source_coords, key=natural_key)),
                "ms_common_source_aliases": ";".join(
                    sorted(common_source_aliases, key=natural_key)
                ),
                "ms_source_alias_status": (
                    "common_source_alias"
                    if common_source_aliases
                    else "target_alias_only"
                ),
                "ms_site_aliases": ";".join(sorted(aliases, key=natural_key)),
                "expected_alleles": ";".join(expected_patterns(n_sites)),
            }
            codons.append(rec)

            for rank, site in enumerate(sites, 1):
                site_output.append(
                    {
                        **site,
                        "multi_edit_target_handle": rec["target_handle"],
                        "multi_edit_site_count": n_sites,
                        "site_rank_in_codon": rank,
                        "codon_position_1based": offsets[rank - 1],
                    }
                )

    codons.sort(
        key=lambda r: (
            SPECIES.index(r["species"]),
            natural_key(r["target_transcript"]),
            int(r["codon_start_1based"]),
        )
    )
    site_output.sort(
        key=lambda r: (
            SPECIES.index(r["species"]),
            natural_key(r["target_transcript"]),
            int(r["target_cds_pos_1based"]),
        )
    )

    candidate_qc_by_species = {r["species"]: r for r in candidate_qc}
    for species in SPECIES:
        sc = [r for r in codons if r["species"] == species]
        candidate_qc_by_species[species]["multi_edit_codons"] = len(sc)
        candidate_qc_by_species[species]["double_edit_codons"] = sum(
            int(r["edit_site_count"]) == 2 for r in sc
        )
        candidate_qc_by_species[species]["triple_edit_codons"] = sum(
            int(r["edit_site_count"]) == 3 for r in sc
        )

    site_fields = [
        "species",
        "target_transcript",
        "target_cds_pos_1based",
        "target_codon_start_1based",
        "target_handle",
        "genomic_sites",
        "genes",
        "source_ms_sites",
        "target_ms_sites",
        "source_coordinates",
        "source_row_count",
        "candidate_lib_info",
        "multi_edit_target_handle",
        "multi_edit_site_count",
        "site_rank_in_codon",
        "codon_position_1based",
    ]
    codon_fields = [
        "species",
        "target_handle",
        "target_transcript",
        "codon_start_1based",
        "aa_position_1based",
        "edit_site_count",
        "edit_cds_positions",
        "edit_codon_positions",
        "pos1_cds",
        "pos2_cds",
        "pos3_cds",
        "genomic_sites",
        "gene",
        "source_coordinates",
        "ms_common_source_aliases",
        "ms_source_alias_status",
        "ms_site_aliases",
        "expected_alleles",
    ]
    qc_fields = [
        "species",
        "candidate_file",
        "input_rows",
        "coding_AtoG_rows_parsed",
        "unique_target_CDS_sites",
        "multi_edit_codons",
        "double_edit_codons",
        "triple_edit_codons",
        "input_columns",
    ]

    write_tsv(
        outdir / "candidate_sites_in_multiedit_codons.tsv", site_output, site_fields
    )
    write_tsv(outdir / "multiedit_codons.tsv", codons, codon_fields)
    write_tsv(
        outdir / "double_edit_codons.tsv",
        [r for r in codons if int(r["edit_site_count"]) == 2],
        codon_fields,
    )
    write_tsv(
        outdir / "triple_edit_codons.tsv",
        [r for r in codons if int(r["edit_site_count"]) == 3],
        codon_fields,
    )
    write_tsv(outdir / "candidate_QC.tsv", candidate_qc, qc_fields)
    return codons


def parse_alledit_sample(path: Path) -> dict[str, Any] | None:
    name = path.name
    suffix = ".Asite_df_alledit"
    if not name.endswith(suffix):
        return None
    stem = name[: -len(suffix)]

    # ISO1 has one real RNA/RPF sample per sex.
    m = re.fullmatch(r"Dmel_ISO1_head_(female|male)_ribo", stem, flags=re.IGNORECASE)
    if m:
        sex = "F" if m.group(1).lower() == "female" else "M"
        return {
            "species": "mel",
            "genotype": "ISO1",
            "sex": sex,
            "rep": 1,
            "sample": stem,
        }

    patterns = [
        (r"WBg_([FM])(\d+)", "mel", "WBg"),
        (r"AdarG_([FM])(\d+)", "mel", "AG"),
        (r"AdarS_([FM])(\d+)", "mel", "AS"),
        (r"Dmel_Adar5G1_([FM])(\d+)(?:_ribo)?", "mel", "Adar5G1"),
        (r"Dsim_([FM])(\d+)", "sim", "Dsim"),
        (r"Dsec_([FM])(\d+)", "sec", "Dsec"),
        (r"Dere_([FM])(\d+)", "ere", "Dere"),
    ]
    for pattern, species, genotype in patterns:
        m = re.fullmatch(pattern, stem, flags=re.IGNORECASE)
        if m:
            return {
                "species": species,
                "genotype": genotype,
                "sex": m.group(1).upper(),
                "rep": int(m.group(2)),
                "sample": stem,
            }

    # Optional head samples are real additional RNA/RPF samples; treat as rep4.
    m = re.fullmatch(
        r"D(sim|sec|ere)_(female|male)_head_ribo", stem, flags=re.IGNORECASE
    )
    if m:
        sp = m.group(1).lower()
        species = {"sim": "sim", "sec": "sec", "ere": "ere"}[sp]
        genotype = {"sim": "Dsim", "sec": "Dsec", "ere": "Dere"}[sp]
        sex = "F" if m.group(2).lower() == "female" else "M"
        return {
            "species": species,
            "genotype": genotype,
            "sex": sex,
            "rep": 4,
            "sample": stem,
        }

    return None


def classify_codon_haplotype(codon: str, codon_positions: list[int]) -> str | None:
    c = clean(codon).upper().replace("U", "T")
    if len(c) != 3:
        return None
    states: list[str] = []
    for pos in codon_positions:
        base = c[pos - 1]
        if base == "A":
            states.append("A")
        elif base == "G":
            states.append("I")
        else:
            return None
    return "".join(states)


def extract_counts(codons: list[dict[str, Any]], alledit_dir: Path, outdir: Path):
    if not alledit_dir.is_dir():
        raise NotADirectoryError(alledit_dir)

    target_by_species_handle: dict[str, dict[str, dict[str, Any]]] = {
        sp: {} for sp in SPECIES
    }
    for codon in codons:
        target_by_species_handle[codon["species"]][codon["target_handle"]] = codon

    files = sorted(
        alledit_dir.glob("*.Asite_df_alledit"), key=lambda p: natural_key(p.name)
    )
    if not files:
        raise FileNotFoundError(f"No *.Asite_df_alledit in {alledit_dir}")

    manifest: list[dict[str, Any]] = []
    output_rows: list[dict[str, Any]] = []
    conflict_rows: list[dict[str, Any]] = []
    unexpected_rows: list[dict[str, Any]] = []

    required = {
        "handle",
        "ref_codon",
        "alt_codon",
        "ref_number_rna",
        "alt_number_rna",
        "ref_number_ribo",
        "alt_number_ribo",
        "edit_pos",
        "FBt_cdspos",
    }

    for path in files:
        sm = parse_alledit_sample(path)
        if sm is None:
            manifest.append(
                {
                    "file": str(path),
                    "used": "no",
                    "reason": "unrecognized_sample_name",
                    "species": "",
                    "genotype": "",
                    "sex": "",
                    "rep": "",
                    "rows_total": "",
                    "rows_target_handle": "",
                    "target_codons_with_rows": "",
                }
            )
            continue

        targets = target_by_species_handle[sm["species"]]
        # exact_codon_observations[handle][codon][assay] -> list[count]
        exact_obs: dict[str, dict[str, dict[str, list[float]]]] = defaultdict(
            lambda: defaultdict(lambda: {"RNA": [], "RPF": []})
        )
        rows_total = 0
        rows_target = 0

        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = set(reader.fieldnames or [])
            missing = required - fields
            if missing:
                raise ValueError(f"{path} missing columns: {sorted(missing)}")

            for row in reader:
                rows_total += 1
                h = clean(row.get("handle"))
                if h not in targets:
                    continue
                rows_target += 1

                ref_codon = clean(row.get("ref_codon")).upper().replace("U", "T")
                alt_codon = clean(row.get("alt_codon")).upper().replace("U", "T")
                if ref_codon:
                    exact_obs[h][ref_codon]["RNA"].append(
                        as_float(row.get("ref_number_rna"))
                    )
                    exact_obs[h][ref_codon]["RPF"].append(
                        as_float(row.get("ref_number_ribo"))
                    )
                if alt_codon:
                    exact_obs[h][alt_codon]["RNA"].append(
                        as_float(row.get("alt_number_rna"))
                    )
                    exact_obs[h][alt_codon]["RPF"].append(
                        as_float(row.get("alt_number_ribo"))
                    )

        # Count conflicts are expected when the same haplotype appears in overlapping
        # transition rows. Keep max, but expose all disagreements in QC.
        exact_counts: dict[str, dict[str, dict[str, float]]] = defaultdict(dict)
        for h, codon_map in exact_obs.items():
            exact_counts[h] = {}
            for raw_codon, assay_map in codon_map.items():
                exact_counts[h][raw_codon] = {}
                for assay in ("RNA", "RPF"):
                    vals = assay_map[assay]
                    chosen = max(vals) if vals else 0.0
                    exact_counts[h][raw_codon][assay] = chosen
                    uniq = sorted(set(vals))
                    if len(uniq) > 1:
                        conflict_rows.append(
                            {
                                **sm,
                                "file": path.name,
                                "target_handle": h,
                                "observed_codon": raw_codon,
                                "assay": assay,
                                "n_observations": len(vals),
                                "unique_observed_counts": ";".join(
                                    f"{x:.10g}" for x in uniq
                                ),
                                "chosen_max_count": f"{chosen:.10g}",
                            }
                        )

        for h, meta in sorted(targets.items(), key=lambda kv: natural_key(kv[0])):
            n_sites = int(meta["edit_site_count"])
            patterns = expected_patterns(n_sites)
            codon_positions = [
                int(x) for x in str(meta["edit_codon_positions"]).split(";")
            ]
            allele_counts = {
                allele: {"RNA": 0.0, "RPF": 0.0, "codons": []} for allele in patterns
            }

            for raw_codon, assay_counts in exact_counts.get(h, {}).items():
                allele = classify_codon_haplotype(raw_codon, codon_positions)
                if allele is None or allele not in allele_counts:
                    unexpected_rows.append(
                        {
                            **sm,
                            "file": path.name,
                            "target_handle": h,
                            "edit_codon_positions": meta["edit_codon_positions"],
                            "observed_codon": raw_codon,
                            "RNA_count_max": f"{assay_counts['RNA']:.10g}",
                            "RPF_count_max": f"{assay_counts['RPF']:.10g}",
                            "reason": "editing positions are not all A/G",
                        }
                    )
                    continue
                allele_counts[allele]["RNA"] += assay_counts["RNA"]
                allele_counts[allele]["RPF"] += assay_counts["RPF"]
                allele_counts[allele]["codons"].append(raw_codon)

            rna_total = sum(v["RNA"] for v in allele_counts.values())
            rpf_total = sum(v["RPF"] for v in allele_counts.values())
            for allele in patterns:
                v = allele_counts[allele]
                output_rows.append(
                    {
                        **sm,
                        "target_handle": h,
                        "target_transcript": meta["target_transcript"],
                        "codon_start_1based": meta["codon_start_1based"],
                        "aa_position_1based": meta["aa_position_1based"],
                        "edit_site_count": n_sites,
                        "edit_cds_positions": meta["edit_cds_positions"],
                        "edit_codon_positions": meta["edit_codon_positions"],
                        "genomic_sites": meta["genomic_sites"],
                        "gene": meta["gene"],
                        "ms_site_aliases": meta["ms_site_aliases"],
                        "allele": allele,
                        "observed_codons_for_allele": ";".join(sorted(v["codons"])),
                        "RNA_count": f"{v['RNA']:.10g}",
                        "RPF_count": f"{v['RPF']:.10g}",
                        "RNA_total": f"{rna_total:.10g}",
                        "RPF_total": f"{rpf_total:.10g}",
                        "RNA_proportion": f"{(v['RNA']/rna_total) if rna_total > 0 else 0:.10g}",
                        "RPF_proportion": f"{(v['RPF']/rpf_total) if rpf_total > 0 else 0:.10g}",
                    }
                )

        manifest.append(
            {
                "file": str(path),
                "used": "yes",
                "reason": "included",
                **sm,
                "rows_total": rows_total,
                "rows_target_handle": rows_target,
                "target_codons_with_rows": len(exact_obs),
            }
        )
        print(
            f"[alledit] {path.name}: {rows_target:,} target transition rows; "
            f"{len(exact_obs)} target codons with rows"
        )

    output_rows.sort(
        key=lambda r: (
            SPECIES.index(r["species"]),
            natural_key(r["target_handle"]),
            natural_key(r["genotype"]),
            r["sex"],
            int(r["rep"]),
            expected_patterns(int(r["edit_site_count"])).index(r["allele"]),
        )
    )

    fields = [
        "species",
        "genotype",
        "sex",
        "rep",
        "sample",
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
        "allele",
        "observed_codons_for_allele",
        "RNA_count",
        "RPF_count",
        "RNA_total",
        "RPF_total",
        "RNA_proportion",
        "RPF_proportion",
    ]
    write_tsv(outdir / "allele_counts.by_sample.tsv", output_rows, fields)
    write_tsv(
        outdir / "double_edit_allele_counts.by_sample.tsv",
        [r for r in output_rows if int(r["edit_site_count"]) == 2],
        fields,
    )
    write_tsv(
        outdir / "triple_edit_allele_counts.by_sample.tsv",
        [r for r in output_rows if int(r["edit_site_count"]) == 3],
        fields,
    )

    manifest_fields = [
        "file",
        "used",
        "reason",
        "species",
        "genotype",
        "sex",
        "rep",
        "sample",
        "rows_total",
        "rows_target_handle",
        "target_codons_with_rows",
    ]
    write_tsv(outdir / "alledit_manifest.tsv", manifest, manifest_fields)

    conflict_fields = [
        "species",
        "genotype",
        "sex",
        "rep",
        "sample",
        "file",
        "target_handle",
        "observed_codon",
        "assay",
        "n_observations",
        "unique_observed_counts",
        "chosen_max_count",
    ]
    write_tsv(outdir / "alledit_count_conflicts.tsv", conflict_rows, conflict_fields)

    unexpected_fields = [
        "species",
        "genotype",
        "sex",
        "rep",
        "sample",
        "file",
        "target_handle",
        "edit_codon_positions",
        "observed_codon",
        "RNA_count_max",
        "RPF_count_max",
        "reason",
    ]
    write_tsv(
        outdir / "alledit_unexpected_codons.tsv", unexpected_rows, unexpected_fields
    )

    return {
        "alledit_files_found": len(files),
        "alledit_files_used": sum(r["used"] == "yes" for r in manifest),
        "allele_rows": len(output_rows),
        "count_conflict_rows": len(conflict_rows),
        "unexpected_codon_rows": len(unexpected_rows),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-dmel", type=Path, default=DEFAULT_DML)
    p.add_argument(
        "--candidate-sim",
        type=Path,
        default=DEFAULT_PROJECTION_ROOT
        / "sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
    )
    p.add_argument(
        "--candidate-sec",
        type=Path,
        default=DEFAULT_PROJECTION_ROOT
        / "sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
    )
    p.add_argument(
        "--candidate-ere",
        type=Path,
        default=DEFAULT_PROJECTION_ROOT
        / "ere.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
    )
    p.add_argument("--alledit-dir", type=Path, default=DEFAULT_ALLEDIT_DIR)
    p.add_argument(
        "--outdir", type=Path, default=Path("double_edit_rebuild_v4/01_counts")
    )
    args = p.parse_args()

    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    candidate_paths = {
        "mel": args.candidate_dmel.resolve(),
        "sim": args.candidate_sim.resolve(),
        "sec": args.candidate_sec.resolve(),
        "ere": args.candidate_ere.resolve(),
    }

    codons = build_candidate_catalog(candidate_paths, outdir)
    stats = extract_counts(codons, args.alledit_dir.resolve(), outdir)

    summary = {
        "candidate_files": {k: str(v) for k, v in candidate_paths.items()},
        "alledit_dir": str(args.alledit_dir.resolve()),
        "outdir": str(outdir),
        "multi_edit_codons": len(codons),
        "double_edit_codons": sum(int(r["edit_site_count"]) == 2 for r in codons),
        "triple_edit_codons": sum(int(r["edit_site_count"]) == 3 for r in codons),
        **stats,
    }
    (outdir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
