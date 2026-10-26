#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""PausePred peak enrichment inside CDS-coordinate hairpins.

PausePred positions are transcript-coordinate, 1-based nucleotide positions.
Hairpin intervals and the explicit editing-site coordinates are CDS-coordinate,
1-based inclusive. Every PausePred and unedited position is therefore mapped
from transcript position to cumulative CDS position before hairpin membership
is tested.

The two reported scopes reproduce the definitions of the existing new
PausePred analyses:

* all_edit_cds: all editing sites retained in the PausePred universe are
  compared with all non-editing universe sites.
* unedit_A_base: all retained editing sites are kept, while the unedited
  comparison background is restricted to reference A bases in CDS.

For each scope, sample, Pause_score cutoff, and transcript-coordinate window,
the script reports counts and a one-sided Fisher exact enrichment test for
editing sites in a peak. Peak targets are made from the complete PausePred
universe, then tested for membership of the hairpin subset; this preserves the
definition used by the existing enrichment scripts, including window_nt.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

try:
    from scipy.stats import fisher_exact
except Exception as exc:  # pragma: no cover
    fisher_exact = None
    SCIPY_IMPORT_ERROR = exc


SAMPLE_TO_SPECIES = {
    "AdarG_F": "Dmel",
    "AdarG_M": "Dmel",
    "AdarS_F": "Dmel",
    "AdarS_M": "Dmel",
    "Dere_F": "Dere",
    "Dere_M": "Dere",
    "Dmel_Adar5G1_M": "Dmel",
    "Dmel_ISO1_head_female_ribo": "Dmel",
    "Dmel_ISO1_head_male_ribo": "Dmel",
    "Dsec_F": "Dsec",
    "Dsec_M": "Dsec",
    "Dsim_F": "Dsim",
    "Dsim_M": "Dsim",
    "WBg_F": "Dmel",
    "WBg_M": "Dmel",
}

SPECIES_CONFIG = {
    "Dmel": {
        "hairpin": "/data/dongxn/RNAediting/codex/02_stalling/3flank/hairpin/Dmel.fa.T25Z3.2.merged.tsv",
        "cds": "/path/to/RNA_editing_DATAprocessing_new/ref_genome/dmel-r6.12_GTF.CDS_position.tsv",
    },
    "Dsec": {
        "hairpin": "/data/dongxn/RNAediting/codex/02_stalling/3flank/hairpin/Dsec.fa.T25Z3.2.merged.tsv",
        "cds": "/path/to/Proj1_editing/code/02_stalling/2pausepred/dsec_all_r1.3_with_UTR_CDSstop_UTR_position_CDS",
    },
    "Dere": {
        "hairpin": "/data/dongxn/RNAediting/codex/02_stalling/3flank/hairpin/Dere.fa.T25Z3.2.merged.tsv",
        "cds": "/path/to/RNA_editing_DATAprocessing_new/ref_genome/dere-all-r1.05_GTF.CDS_position.tsv",
    },
    "Dsim": {
        "hairpin": "/data/dongxn/RNAediting/codex/02_stalling/3flank/hairpin/Dsim.fa.T25Z3.2.merged.tsv",
        "cds": "/path/to/RNA_editing_DATAprocessing_new/ref_genome/dsim-all-r2.02_GTF.CDS_position.tsv",
    },
}

SCOPES = ("all_edit_cds", "unedit_A_base")


def as_int(value: object) -> Optional[int]:
    try:
        text = "" if value is None else str(value).strip()
        if not text or text.lower() in {"nan", "na", "none", "null", "."}:
            return None
        number = float(text)
        if not math.isfinite(number):
            return None
        return int(number)
    except (TypeError, ValueError):
        return None


def as_float(value: object) -> Optional[float]:
    try:
        text = "" if value is None else str(value).strip()
        if not text or text.lower() in {"nan", "na", "none", "null", "."}:
            return None
        number = float(text)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def normalize_tx(value: object) -> Optional[str]:
    text = "" if value is None else str(value).strip()
    if not text or text.lower() in {"nan", "na", "none", "null", "."}:
        return None
    return re.sub(r"\.[0-9]+$", "", text)


def parse_site_key(value: object) -> Optional[Tuple[str, int]]:
    text = "" if value is None else str(value).strip()
    if not text or text.lower() in {"nan", "na", "none", "null", "."}:
        return None
    if "_" not in text:
        return None
    tx, pos = text.rsplit("_", 1)
    tx_norm = normalize_tx(tx)
    pos_int = as_int(pos)
    if tx_norm is None or pos_int is None or pos_int <= 0:
        return None
    return tx_norm, pos_int


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


class CdsCoordinateMap:
    """Map transcript coordinates to cumulative CDS coordinates.

    This follows the coordinate construction used by the existing hairpin
    codon-level analysis: CDS segments are sorted by transcript coordinate and
    cumulative CDS positions are assigned across segments.
    """

    def __init__(self, path: Path):
        raw: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        with path.open("r", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if len(row) < 3:
                    continue
                tx = normalize_tx(row[0])
                start = as_int(row[1])
                end = as_int(row[2])
                feature = str(row[3]).strip().upper() if len(row) >= 4 else "CDS"
                if tx is None or start is None or end is None or end < start:
                    continue
                if feature not in {"CDS", "", ".", "NA"}:
                    continue
                raw[tx].append((start, end))

        self.segments: Dict[str, List[Tuple[int, int, int, int]]] = {}
        for tx, values in raw.items():
            values.sort()
            cumulative = 1
            built = []
            for tx_start, tx_end in values:
                cds_start = cumulative
                cds_end = cumulative + tx_end - tx_start
                built.append((tx_start, tx_end, cds_start, cds_end))
                cumulative = cds_end + 1
            self.segments[tx] = built
        self.starts = {
            tx: [item[0] for item in items] for tx, items in self.segments.items()
        }

    def tx_to_cds(self, tx: str, tx_pos: int) -> Optional[int]:
        segments = self.segments.get(tx)
        if not segments:
            return None
        index = bisect_right(self.starts[tx], tx_pos) - 1
        if index < 0:
            return None
        for i in range(max(0, index - 1), min(len(segments), index + 2)):
            tx_start, tx_end, cds_start, _ = segments[i]
            if tx_start <= tx_pos <= tx_end:
                return cds_start + tx_pos - tx_start
        return None

    def cds_length(self, tx: str) -> Optional[int]:
        segments = self.segments.get(tx)
        if not segments:
            return None
        return max(item[3] for item in segments)


class HairpinIndex:
    """Merged CDS-coordinate hairpin intervals with fast membership lookup."""

    def __init__(self, path: Path):
        raw: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        with path.open("r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"transcript", "start", "end"}
            if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
                raise RuntimeError(
                    f"Hairpin file must contain transcript/start/end columns: {path}"
                )
            for row in reader:
                tx = normalize_tx(row.get("transcript"))
                start = as_int(row.get("start"))
                end = as_int(row.get("end"))
                if tx is None or start is None or end is None or end < start:
                    continue
                raw[tx].append((start, end))

        self.intervals: Dict[str, List[Tuple[int, int, str]]] = {}
        for tx, values in raw.items():
            values.sort()
            merged: List[Tuple[int, int, str]] = []
            start, end = values[0]
            ids = [f"{tx}:{start}-{end}"]
            for next_start, next_end in values[1:]:
                if next_start <= end + 1:
                    end = max(end, next_end)
                    ids.append(f"{tx}:{next_start}-{next_end}")
                else:
                    merged.append((start, end, ";".join(ids)))
                    start, end = next_start, next_end
                    ids = [f"{tx}:{start}-{end}"]
            merged.append((start, end, ";".join(ids)))
            self.intervals[tx] = merged

        self.starts = {
            tx: [item[0] for item in items] for tx, items in self.intervals.items()
        }
        self.interval_count = sum(len(items) for items in self.intervals.values())
        self.union_nt = sum(
            end - start + 1
            for items in self.intervals.values()
            for start, end, _ in items
        )

    def find(self, tx: str, cds_pos: int) -> Optional[str]:
        items = self.intervals.get(tx)
        if not items:
            return None
        index = bisect_right(self.starts[tx], cds_pos) - 1
        if index < 0:
            return None
        start, end, identifier = items[index]
        if start <= cds_pos <= end:
            return identifier
        return None


def read_used_editing_keys(path: Path) -> Set[Tuple[str, int]]:
    keys: Set[Tuple[str, int]] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "site_key" not in reader.fieldnames:
            raise RuntimeError(f"Editing-used file has no site_key column: {path}")
        for row in reader:
            parsed = parse_site_key(row.get("site_key"))
            if parsed is not None:
                keys.add(parsed)
    return keys


def first_value(row: dict, names: Sequence[str]) -> Optional[str]:
    for name in names:
        if name in row:
            value = row.get(name)
            if value is not None and str(value).strip() != "":
                return str(value).strip()
    return None


def read_editing_coordinates(
    path: Path,
    used_keys: Set[Tuple[str, int]],
    cds_map: CdsCoordinateMap,
) -> Tuple[List[Tuple[str, int, int, str]], dict]:
    """Return used editing sites with CDS coordinates and coordinate QC."""

    coord_by_key: Dict[Tuple[str, int], Tuple[int, str]] = {}
    qc = {
        "editing_used_keys": len(used_keys),
        "editing_extracted_rows": 0,
        "editing_keys_found_in_extracted": 0,
        "editing_keys_with_explicit_cds_coordinate": 0,
        "editing_keys_with_cds_map_coordinate": 0,
        "editing_keys_with_invalid_explicit_coordinate_fallback": 0,
        "editing_keys_with_conflicting_coordinates": 0,
        "editing_keys_unmapped_to_cds": 0,
    }

    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"Empty editing extracted file: {path}")
        for row in reader:
            qc["editing_extracted_rows"] += 1
            parsed = parse_site_key(row.get("site_key"))
            if parsed is None:
                tx = normalize_tx(
                    first_value(row, ("new_transcript", "transcript_id", "transcript"))
                )
                tx_pos = as_int(
                    first_value(
                        row, ("new_tx_pos_1based", "transcript_pos", "tx_pos_1based")
                    )
                )
                if tx is not None and tx_pos is not None and tx_pos > 0:
                    parsed = (tx, tx_pos)
            if parsed is None or parsed not in used_keys:
                continue

            tx, tx_pos = parsed
            explicit_text = first_value(
                row,
                (
                    "new_cds_pos_1based",
                    "alledit_cds_pos",
                    "cds_pos_1based",
                    "cds_position",
                ),
            )
            explicit = as_int(explicit_text)
            source = ""
            cds_pos = None
            if explicit is not None and explicit > 0:
                cds_len = cds_map.cds_length(tx)
                if cds_len is None or explicit <= cds_len:
                    cds_pos = explicit
                    source = "explicit_cds_coordinate"
                    qc["editing_keys_with_explicit_cds_coordinate"] += 1
                else:
                    qc["editing_keys_with_invalid_explicit_coordinate_fallback"] += 1

            if cds_pos is None:
                cds_pos = cds_map.tx_to_cds(tx, tx_pos)
                if cds_pos is not None:
                    source = "transcript_to_cds_map"
                    qc["editing_keys_with_cds_map_coordinate"] += 1

            if cds_pos is None:
                continue

            if parsed in coord_by_key:
                qc["editing_keys_with_conflicting_coordinates"] += int(
                    coord_by_key[parsed][0] != cds_pos
                )
                # Prefer an explicitly supplied coordinate over a fallback.
                if coord_by_key[parsed][1] == "explicit_cds_coordinate":
                    continue
                if source != "explicit_cds_coordinate":
                    continue
            else:
                qc["editing_keys_found_in_extracted"] += 1
            coord_by_key[parsed] = (cds_pos, source)

    # A small number of used keys can be absent from the extracted file. The
    # transcript/CDS reference is an explicit fallback, while QC exposes use.
    for key in sorted(used_keys):
        if key in coord_by_key:
            continue
        tx, tx_pos = key
        cds_pos = cds_map.tx_to_cds(tx, tx_pos)
        if cds_pos is not None:
            coord_by_key[key] = (cds_pos, "transcript_to_cds_map_missing_extracted")
            qc["editing_keys_with_cds_map_coordinate"] += 1

    qc["editing_keys_unmapped_to_cds"] = len(used_keys - set(coord_by_key))
    records = [
        (tx, tx_pos, cds_pos, source)
        for (tx, tx_pos), (cds_pos, source) in sorted(coord_by_key.items())
    ]
    return records, qc


def read_a_base_keys(path: Path) -> Set[Tuple[str, int]]:
    keys: Set[Tuple[str, int]] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"transcript_id", "transcript_pos_1based", "reference_base"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise RuntimeError(f"A-base map is missing required columns: {path}")
        for row in reader:
            if str(row.get("reference_base", "")).strip().upper() != "A":
                continue
            tx = normalize_tx(row.get("transcript_id"))
            pos = as_int(row.get("transcript_pos_1based"))
            if tx is not None and pos is not None and pos > 0:
                keys.add((tx, pos))
    return keys


def read_universe_and_peak_index(
    path: Path,
    cutoffs: Sequence[float],
    editing_keys: Set[Tuple[str, int]],
    a_base_keys: Set[Tuple[str, int]],
) -> Tuple[Dict[float, Dict[str, List[int]]], dict]:
    """Read the already filtered universe and index peak positions."""

    peak_index: Dict[float, Dict[str, List[int]]] = {
        cutoff: defaultdict(list) for cutoff in cutoffs
    }
    qc = {
        "pausepred_universe_rows": 0,
        "pausepred_universe_editing_keys_seen": 0,
        "pausepred_universe_A_base_keys_seen": 0,
        "pausepred_universe_invalid_rows": 0,
    }
    editing_seen: Set[Tuple[str, int]] = set()
    a_seen: Set[Tuple[str, int]] = set()

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"transcript_id", "transcript_pos", "Pause_score"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise RuntimeError(
                f"PausePred universe is missing required columns: {path}"
            )
        for row in reader:
            tx = normalize_tx(row.get("transcript_id"))
            pos = as_int(row.get("transcript_pos"))
            score = as_float(row.get("Pause_score"))
            if tx is None or pos is None or pos <= 0 or score is None:
                qc["pausepred_universe_invalid_rows"] += 1
                continue
            qc["pausepred_universe_rows"] += 1
            key = (tx, pos)
            if key in editing_keys:
                editing_seen.add(key)
            if key in a_base_keys:
                a_seen.add(key)
            for cutoff in cutoffs:
                if score >= cutoff:
                    peak_index[cutoff][tx].append(pos)

    for cutoff in cutoffs:
        for positions in peak_index[cutoff].values():
            positions.sort()

    qc["pausepred_universe_editing_keys_seen"] = len(editing_seen)
    qc["pausepred_universe_A_base_keys_seen"] = len(a_seen)
    return peak_index, qc


def peak_hit(
    peak_index: Dict[float, Dict[str, List[int]]],
    cutoff: float,
    window_nt: int,
    tx: str,
    tx_pos: int,
) -> bool:
    positions = peak_index[cutoff].get(tx)
    if not positions:
        return False
    left = bisect_left(positions, tx_pos - window_nt)
    return left < len(positions) and positions[left] <= tx_pos + window_nt


def initialize_counts(cutoffs: Sequence[float], windows: Sequence[int]) -> dict:
    output = {}
    for scope in SCOPES:
        output[scope] = {}
        for cutoff in cutoffs:
            for window_nt in windows:
                output[scope][(cutoff, window_nt)] = {
                    "edited_total": 0,
                    "edited_peak": 0,
                    "unedited_total": 0,
                    "unedited_peak": 0,
                }
    return output


def add_site_to_counts(
    counts: dict,
    scope: str,
    tx: str,
    tx_pos: int,
    peak_index: Dict[float, Dict[str, List[int]]],
    cutoffs: Sequence[float],
    windows: Sequence[int],
    edited: bool,
) -> None:
    for cutoff in cutoffs:
        for window_nt in windows:
            record = counts[scope][(cutoff, window_nt)]
            field = "edited_total" if edited else "unedited_total"
            peak_field = "edited_peak" if edited else "unedited_peak"
            record[field] += 1
            if peak_hit(peak_index, cutoff, window_nt, tx, tx_pos):
                record[peak_field] += 1


def process_sample(
    sample: str,
    all_edit_dir: Path,
    a_base_dir: Path,
    hairpin_dir: Path,
    cutoffs: Sequence[float],
    windows: Sequence[int],
    progress_every: int,
) -> Tuple[List[dict], dict]:
    if sample not in SAMPLE_TO_SPECIES:
        raise RuntimeError(f"No species configuration for sample: {sample}")
    species = SAMPLE_TO_SPECIES[sample]
    config = SPECIES_CONFIG[species]
    hairpin_path = hairpin_dir / Path(config["hairpin"]).name
    cds_path = Path(config["cds"])

    prefix = f"{sample}.pausepred_enrichment"
    universe_path = all_edit_dir / f"{prefix}.pausepred_universe.tsv"
    used_path = all_edit_dir / f"{prefix}.editing_sites_used_in_universe.tsv"
    extracted_path = all_edit_dir / f"{prefix}.editing_sites_extracted.tsv"
    unedited_path = all_edit_dir / f"{prefix}.unedited_sites.csv"
    a_map_path = a_base_dir / "a_base_maps" / f"{species}.A_base_map.tsv.gz"

    for path in (
        hairpin_path,
        cds_path,
        universe_path,
        used_path,
        extracted_path,
        unedited_path,
    ):
        if not path.exists():
            raise FileNotFoundError(f"Required file is missing: {path}")
    if not a_map_path.exists():
        raise FileNotFoundError(f"Required A-base map is missing: {a_map_path}")

    print(f"[start] {sample} ({species})", flush=True)
    hairpins = HairpinIndex(hairpin_path)
    cds_map = CdsCoordinateMap(cds_path)
    used_keys = read_used_editing_keys(used_path)
    editing_records, edit_qc = read_editing_coordinates(
        extracted_path, used_keys, cds_map
    )
    a_base_keys = read_a_base_keys(a_map_path)
    peak_index, universe_qc = read_universe_and_peak_index(
        universe_path, cutoffs, used_keys, a_base_keys
    )

    counts = initialize_counts(cutoffs, windows)
    qc = {
        "sample": sample,
        "species": species,
        "hairpin_file": str(hairpin_path),
        "cds_coordinate_file": str(cds_path),
        "pausepred_universe_file": str(universe_path),
        "editing_sites_used_file": str(used_path),
        "editing_sites_extracted_file": str(extracted_path),
        "unedited_sites_file": str(unedited_path),
        "a_base_map_file": str(a_map_path),
        "hairpin_coordinate_system": "CDS_1based_inclusive",
        "pausepred_coordinate_system": "transcript_1based",
        "editing_coordinate_system": "CDS_1based",
        "hairpin_merged_interval_count": hairpins.interval_count,
        "hairpin_transcript_count": len(hairpins.intervals),
        "hairpin_union_nt": hairpins.union_nt,
        "cds_map_transcript_count": len(cds_map.segments),
        "hairpin_transcripts_missing_cds_map": len(
            set(hairpins.intervals) - set(cds_map.segments)
        ),
        "editing_records_with_cds_coordinate": len(editing_records),
        "editing_records_inside_hairpin": 0,
        "editing_records_outside_hairpin": 0,
        "editing_records_cds_unmapped": edit_qc["editing_keys_unmapped_to_cds"],
        "a_base_keys_total": len(a_base_keys),
        "unedited_input_rows": 0,
        "unedited_rows_with_cds_coordinate": 0,
        "unedited_rows_cds_unmapped": 0,
        "unedited_rows_on_transcripts_without_hairpin": 0,
        "unedited_rows_outside_hairpin": 0,
        "unedited_rows_inside_hairpin_all_edit": 0,
        "unedited_rows_inside_hairpin_A_base": 0,
        "unedited_A_base_rows_cds_unmapped": 0,
        "unedited_A_base_rows_not_in_hairpin": 0,
    }
    qc.update(edit_qc)
    qc.update(universe_qc)

    # Editing coordinates are explicitly CDS coordinates. Peak membership
    # still uses the original transcript key, as in the source enrichment.
    for tx, tx_pos, cds_pos, _source in editing_records:
        hp_id = hairpins.find(tx, cds_pos)
        if hp_id is None:
            qc["editing_records_outside_hairpin"] += 1
            continue
        qc["editing_records_inside_hairpin"] += 1
        for scope in SCOPES:
            add_site_to_counts(
                counts,
                scope,
                tx,
                tx_pos,
                peak_index,
                cutoffs,
                windows,
                edited=True,
            )

    # The unedited file is the complement of editing sites in the filtered
    # PausePred universe. It is streamed because it can contain millions of
    # rows. A-base filtering is applied only to this unedited background.
    with unedited_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"transcript_id", "transcript_pos", "site_key"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise RuntimeError(
                f"Unedited file is missing required columns: {unedited_path}"
            )
        for row in reader:
            qc["unedited_input_rows"] += 1
            if progress_every and qc["unedited_input_rows"] % progress_every == 0:
                print(
                    f"[progress] {sample}: unedited rows {qc['unedited_input_rows']:,}",
                    flush=True,
                )
            tx = normalize_tx(row.get("transcript_id"))
            tx_pos = as_int(row.get("transcript_pos"))
            if tx is None or tx_pos is None or tx_pos <= 0:
                qc["unedited_rows_cds_unmapped"] += 1
                continue
            key = (tx, tx_pos)
            cds_pos = cds_map.tx_to_cds(tx, tx_pos)
            if cds_pos is None:
                qc["unedited_rows_cds_unmapped"] += 1
                if key in a_base_keys:
                    qc["unedited_A_base_rows_cds_unmapped"] += 1
                continue
            qc["unedited_rows_with_cds_coordinate"] += 1
            hp_id = hairpins.find(tx, cds_pos)
            if hp_id is None:
                if tx not in hairpins.intervals:
                    qc["unedited_rows_on_transcripts_without_hairpin"] += 1
                else:
                    qc["unedited_rows_outside_hairpin"] += 1
                if key in a_base_keys:
                    qc["unedited_A_base_rows_not_in_hairpin"] += 1
                continue

            qc["unedited_rows_inside_hairpin_all_edit"] += 1
            add_site_to_counts(
                counts,
                "all_edit_cds",
                tx,
                tx_pos,
                peak_index,
                cutoffs,
                windows,
                edited=False,
            )

            if key in a_base_keys:
                qc["unedited_rows_inside_hairpin_A_base"] += 1
                add_site_to_counts(
                    counts,
                    "unedit_A_base",
                    tx,
                    tx_pos,
                    peak_index,
                    cutoffs,
                    windows,
                    edited=False,
                )

    output_rows = []
    for scope in SCOPES:
        for cutoff in cutoffs:
            for window_nt in windows:
                record = counts[scope][(cutoff, window_nt)]
                edited_total = record["edited_total"]
                edited_peak = record["edited_peak"]
                unedited_total = record["unedited_total"]
                unedited_peak = record["unedited_peak"]
                edited_nonpeak = edited_total - edited_peak
                unedited_nonpeak = unedited_total - unedited_peak
                if edited_total and unedited_total:
                    if fisher_exact is None:
                        raise RuntimeError(
                            f"scipy.stats.fisher_exact is required: {SCIPY_IMPORT_ERROR}"
                        )
                    odds_ratio, p_value = fisher_exact(
                        [
                            [edited_peak, edited_nonpeak],
                            [unedited_peak, unedited_nonpeak],
                        ],
                        alternative="greater",
                    )
                    odds_ratio = float(odds_ratio)
                    p_value = float(p_value)
                else:
                    odds_ratio = float("nan")
                    p_value = float("nan")
                output_rows.append(
                    {
                        "sample": sample,
                        "species": species,
                        "analysis_scope": scope,
                        "Pause_score_cutoff": cutoff,
                        "window_nt": window_nt,
                        "edited_total_in_hairpin": edited_total,
                        "edited_in_peak_in_hairpin": edited_peak,
                        "edited_not_in_peak_in_hairpin": edited_nonpeak,
                        "edited_peak_fraction": (
                            edited_peak / edited_total if edited_total else float("nan")
                        ),
                        "unedited_total_in_hairpin": unedited_total,
                        "unedited_in_peak_in_hairpin": unedited_peak,
                        "unedited_not_in_peak_in_hairpin": unedited_nonpeak,
                        "unedited_peak_fraction": (
                            unedited_peak / unedited_total
                            if unedited_total
                            else float("nan")
                        ),
                        "total_sites_in_hairpin": edited_total + unedited_total,
                        "total_peak_sites_in_hairpin": edited_peak + unedited_peak,
                        "odds_ratio_editing_in_peak": odds_ratio,
                        "p_value_one_sided_fisher_greater": p_value,
                        "p_adj_BH_within_scope_cutoff_window": float("nan"),
                        "peak_definition": (
                            "Pause_score>=cutoff in complete filtered universe; "
                            "window expanded on transcript coordinate"
                        ),
                        "hairpin_coordinate_system": "CDS_1based_inclusive",
                        "pausepred_coordinate_system": "transcript_1based_mapped_to_CDS",
                        "edited_background_filter": "all_editing_sites_used_in_universe",
                        "unedited_background_filter": (
                            "all_unedited_universe_sites"
                            if scope == "all_edit_cds"
                            else "reference_CDS_base_is_A"
                        ),
                    }
                )

    print(
        f"[done] {sample}: universe={qc['pausepred_universe_rows']:,}; "
        f"editing_in_hairpin={qc['editing_records_inside_hairpin']:,}; "
        f"unedited_in_hairpin(all/A)="
        f"{qc['unedited_rows_inside_hairpin_all_edit']:,}/"
        f"{qc['unedited_rows_inside_hairpin_A_base']:,}",
        flush=True,
    )
    return output_rows, qc


def process_sample_worker(payload):
    """Top-level worker entry point for ProcessPoolExecutor pickling."""
    return process_sample(*payload)


def bh_adjust(rows: List[dict]) -> None:
    groups: Dict[Tuple[str, float, int], List[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        p_value = row["p_value_one_sided_fisher_greater"]
        if p_value == p_value and math.isfinite(p_value):
            groups[
                (
                    row["analysis_scope"],
                    float(row["Pause_score_cutoff"]),
                    int(row["window_nt"]),
                )
            ].append(index)
    for indices in groups.values():
        indices.sort(key=lambda i: rows[i]["p_value_one_sided_fisher_greater"])
        n = len(indices)
        running = 1.0
        for rank_from_end, index in enumerate(reversed(indices), start=1):
            rank = n - rank_from_end + 1
            p_value = rows[index]["p_value_one_sided_fisher_greater"]
            adjusted = min(running, p_value * n / rank)
            running = adjusted
            rows[index]["p_adj_BH_within_scope_cutoff_window"] = adjusted


def write_tsv(path: Path, rows: Sequence[dict]) -> None:
    if not rows:
        path.write_text("\n")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_float_list(value: str) -> List[float]:
    output = []
    for token in value.split(","):
        number = as_float(token)
        if number is None:
            raise ValueError(f"Invalid numeric value: {token}")
        output.append(number)
    return output


def parse_int_list(value: str) -> List[int]:
    output = []
    for token in value.split(","):
        number = as_int(token)
        if number is None:
            raise ValueError(f"Invalid integer value: {token}")
        output.append(number)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--all-edit-dir",
        type=Path,
        default=Path(
            "/data/dongxn/RNAediting/codex/02_stalling/2pausepred/"
            "pausepred_results_new/enrichment_alledit_cds"
        ),
    )
    parser.add_argument(
        "--a-base-dir",
        type=Path,
        default=Path(
            "/data/dongxn/RNAediting/codex/02_stalling/2pausepred/"
            "pausepred_results_new/enrichment_unedit_A_base"
        ),
    )
    parser.add_argument(
        "--hairpin-dir",
        type=Path,
        default=Path("/data/dongxn/RNAediting/codex/02_stalling/3flank/hairpin"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(
            "/data/dongxn/RNAediting/codex/02_stalling/2pausepred/"
            "pausepred_results_new/hairpin_edit_unedit_enrichment"
        ),
    )
    parser.add_argument(
        "--samples",
        default=",".join(
            sample
            for sample in SAMPLE_TO_SPECIES
            if sample not in {"AdarG_F", "AdarG_M", "AdarS_F", "AdarS_M"}
        ),
    )
    parser.add_argument("--pause-cutoffs", default="5,10")
    parser.add_argument("--windows", default="0,2")
    parser.add_argument("--progress-every", type=int, default=1000000)
    parser.add_argument(
        "--jobs",
        "--threads",
        dest="jobs",
        type=int,
        default=int(os.environ.get("PAUSEPRED_JOBS", "1")),
        help=(
            "Number of samples processed in parallel. Process workers are "
            "used so CSV parsing is not limited by the Python GIL."
        ),
    )
    parser.add_argument(
        "--executor",
        choices=("process", "thread"),
        default="process",
        help="Parallel backend when --jobs > 1 (default: process).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cutoffs = parse_float_list(args.pause_cutoffs)
    windows = parse_int_list(args.windows)
    samples = [item.strip() for item in args.samples.split(",") if item.strip()]
    unknown = [sample for sample in samples if sample not in SAMPLE_TO_SPECIES]
    if unknown:
        raise SystemExit(f"Unknown samples: {', '.join(unknown)}")
    if args.jobs < 1:
        raise SystemExit("--jobs/--threads must be >= 1")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    all_rows: List[dict] = []
    qc_rows: List[dict] = []
    payloads = [
        (
            sample,
            args.all_edit_dir,
            args.a_base_dir,
            args.hairpin_dir,
            cutoffs,
            windows,
            args.progress_every,
        )
        for sample in samples
    ]

    def record_completed(sample, result):
        rows, qc = result
        all_rows.extend(rows)
        qc_rows.append(qc)
        # Checkpoints make an interrupted multi-sample run recoverable and
        # provide progress files before the combined table is finalized.
        write_tsv(
            args.out_dir / f"{sample}.hairpin_enrichment_summary.tsv",
            rows,
        )
        write_tsv(
            args.out_dir / f"{sample}.hairpin_enrichment_QC.tsv",
            [qc],
        )
        print(f"[checkpoint] {sample}", flush=True)

    if args.jobs == 1:
        for sample, payload in zip(samples, payloads):
            record_completed(sample, process_sample_worker(payload))
    else:
        executor_class = (
            ProcessPoolExecutor if args.executor == "process" else ThreadPoolExecutor
        )
        worker_count = min(args.jobs, len(payloads))
        print(
            f"[parallel] samples={len(payloads)} workers={worker_count} "
            f"executor={args.executor} windows={windows}",
            flush=True,
        )
        with executor_class(max_workers=worker_count) as executor:
            futures = {
                executor.submit(process_sample_worker, payload): sample
                for sample, payload in zip(samples, payloads)
            }
            for future in as_completed(futures):
                sample = futures[future]
                record_completed(sample, future.result())

    all_rows.sort(
        key=lambda row: (
            row["analysis_scope"],
            row["species"],
            row["sample"],
            float(row["Pause_score_cutoff"]),
            int(row["window_nt"]),
        )
    )
    qc_rows.sort(key=lambda row: (row["species"], row["sample"]))
    bh_adjust(all_rows)
    summary_path = args.out_dir / "hairpin_pausepred_edit_unedit_enrichment_summary.tsv"
    qc_path = args.out_dir / "hairpin_pausepred_edit_unedit_enrichment_QC.tsv"
    write_tsv(summary_path, all_rows)
    write_tsv(qc_path, qc_rows)

    readme = args.out_dir / "README.md"
    pass
    print(f"[write] {summary_path}", flush=True)
    print(f"[write] {qc_path}", flush=True)


if __name__ == "__main__":
    main()
