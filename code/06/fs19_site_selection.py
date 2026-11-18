#!/usr/bin/env python3
"""Site-selection helpers used only by the Figure S19 analysis."""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# Standard genetic code in T/C/A/G lexicographic codon order.
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: aa
    for aa, (a, b, c) in zip(
        _AAS,
        ((a, b, c) for a in _BASES for b in _BASES for c in _BASES),
    )
}

GENE_RE = re.compile(r"\bFBgn\d+\b")
TRANSCRIPT_RE = re.compile(r"\bFBtr\d+\b")
PARENT_RE = re.compile(r"(?:^|[;\s])parent=([^;]+)")


@dataclass
class CdsRecord:
    header: str
    sequence: str
    gene_id: str
    transcript_id: str
    nt_length: int
    trimmed_nt: int
    ambiguous_codon_count: int
    internal_stop_count: int
    terminal_stop_present: bool
    protein_sequence: str


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    parts: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                parts.append(re.sub(r"\s+", "", line).upper().replace("U", "T"))
    if header is not None:
        yield header, "".join(parts)


def extract_parent_ids(header: str) -> Tuple[Optional[str], Optional[str]]:
    parent_match = PARENT_RE.search(header)
    text = parent_match.group(1) if parent_match else header
    gene_match = GENE_RE.search(text)
    transcript_match = TRANSCRIPT_RE.search(text)
    return (
        gene_match.group(0) if gene_match else None,
        transcript_match.group(0) if transcript_match else None,
    )


def translate_cds(sequence: str) -> Tuple[str, int, int, int, bool]:
    clean = sequence.upper().replace("U", "T")
    trimmed_nt = len(clean) % 3
    usable = clean[: len(clean) - trimmed_nt] if trimmed_nt else clean
    aas: List[str] = []
    ambiguous = 0
    for i in range(0, len(usable), 3):
        aa = CODON_TABLE.get(usable[i : i + 3])
        if aa is None:
            aa = "X"
            ambiguous += 1
        aas.append(aa)
    terminal_stop = bool(aas and aas[-1] == "*")
    if terminal_stop:
        aas.pop()
    internal_stops = sum(aa == "*" for aa in aas)
    return (
        "".join("X" if aa == "*" else aa for aa in aas),
        trimmed_nt,
        ambiguous,
        internal_stops,
        terminal_stop,
    )


def build_cds_record(header: str, sequence: str) -> Optional[CdsRecord]:
    gene_id, transcript_id = extract_parent_ids(header)
    if gene_id is None or transcript_id is None:
        return None
    protein, trimmed_nt, ambiguous, internal_stops, terminal_stop = translate_cds(
        sequence
    )
    return CdsRecord(
        header=header,
        sequence=sequence,
        gene_id=gene_id,
        transcript_id=transcript_id,
        nt_length=len(sequence),
        trimmed_nt=trimmed_nt,
        ambiguous_codon_count=ambiguous,
        internal_stop_count=internal_stops,
        terminal_stop_present=terminal_stop,
        protein_sequence=protein,
    )


def choose_longest(records: Sequence[CdsRecord]) -> CdsRecord:
    candidates = list(records)
    max_length = max(r.nt_length for r in candidates)
    candidates = [r for r in candidates if r.nt_length == max_length]
    min_internal = min(r.internal_stop_count for r in candidates)
    candidates = [r for r in candidates if r.internal_stop_count == min_internal]
    min_ambiguous = min(r.ambiguous_codon_count for r in candidates)
    candidates = [r for r in candidates if r.ambiguous_codon_count == min_ambiguous]
    min_trimmed = min(r.trimmed_nt for r in candidates)
    candidates = [r for r in candidates if r.trimmed_nt == min_trimmed]
    return sorted(candidates, key=lambda r: r.transcript_id)[0]


def select_longest_cds(
    cds_fasta: Path,
) -> Tuple[Dict[str, CdsRecord], pd.DataFrame, pd.DataFrame]:
    by_gene: Dict[str, List[CdsRecord]] = defaultdict(list)
    missing_rows: List[dict] = []
    for header, sequence in read_fasta(cds_fasta):
        record = build_cds_record(header, sequence)
        if record is None:
            missing_rows.append(
                {
                    "header": header,
                    "nt_length": len(sequence),
                    "status": "missing_FBgn_or_FBtr",
                }
            )
            continue
        by_gene[record.gene_id].append(record)

    selected: Dict[str, CdsRecord] = {}
    selection_rows: List[dict] = []
    for gene_id, records in by_gene.items():
        winner = choose_longest(records)
        selected[gene_id] = winner
        for record in records:
            selection_rows.append(
                {
                    "gene_id": record.gene_id,
                    "transcript_id": record.transcript_id,
                    "nt_length": record.nt_length,
                    "protein_length": len(record.protein_sequence),
                    "trimmed_nt": record.trimmed_nt,
                    "ambiguous_codon_count": record.ambiguous_codon_count,
                    "internal_stop_count": record.internal_stop_count,
                    "selected": int(record.transcript_id == winner.transcript_id),
                }
            )
    return selected, pd.DataFrame(selection_rows), pd.DataFrame(missing_rows)


def read_table_auto(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep=None, engine="python")


def is_cds_region(value: object) -> bool:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    region = str(value).strip().upper()
    return region == "CDS" or region.startswith("CDS_") or "CODING" in region


def gene_ids_for_rows(
    edit_df: pd.DataFrame,
    gene_col: Optional[str],
    transcript_col: str,
    transcript_to_gene: Dict[str, str],
) -> pd.Series:
    values: List[str] = []
    for _, row in edit_df.iterrows():
        gene_id = ""
        if gene_col is not None:
            match = GENE_RE.search(str(row.get(gene_col, "")).strip())
            if match:
                gene_id = match.group(0)
        if not gene_id:
            match = TRANSCRIPT_RE.search(str(row.get(transcript_col, "")).strip())
            if match:
                gene_id = transcript_to_gene.get(match.group(0), "")
        values.append(gene_id)
    return pd.Series(values, index=edit_df.index, dtype="object")


def _position_to_nt0(value: object, position_type: str) -> Optional[int]:
    try:
        numeric = int(float(value))
    except (TypeError, ValueError):
        return None
    if position_type == "cds_nt_1based":
        return numeric - 1
    if position_type == "cds_nt_0based":
        return numeric
    if position_type == "aa_1based":
        return (numeric - 1) * 3
    if position_type == "aa_0based":
        return numeric * 3
    raise ValueError(position_type)


def prepare_editing_sites(
    edit_df: pd.DataFrame,
    selected_by_gene: Dict[str, CdsRecord],
    transcript_col: str,
    gene_col: Optional[str],
    position_col: str,
    editing_level_col: str,
    editing_sample_col: Optional[str],
    position_type: str,
    require_a: bool,
    mode: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Validate CDS sites and retain the highest mean-editing site per gene."""
    if mode != "one-per-gene":
        raise ValueError("Only one-per-gene site selection is supported")
    transcript_to_record = {r.transcript_id: r for r in selected_by_gene.values()}
    qc_rows: List[dict] = []
    records: List[dict] = []

    for row_index, row in edit_df.iterrows():
        tx_match = TRANSCRIPT_RE.search(str(row.get(transcript_col, "")).strip())
        transcript_id = tx_match.group(0) if tx_match else ""
        record = transcript_to_record.get(transcript_id)
        if record is None:
            qc_rows.append(
                {
                    "input_row": row_index,
                    "transcript_id": transcript_id,
                    "status": "transcript_not_selected_longest_CDS",
                }
            )
            continue
        if gene_col:
            supplied = GENE_RE.search(str(row.get(gene_col, "")))
            if supplied and supplied.group(0) != record.gene_id:
                qc_rows.append(
                    {
                        "input_row": row_index,
                        "transcript_id": transcript_id,
                        "status": "gene_transcript_mismatch",
                    }
                )
                continue
        nt0 = _position_to_nt0(row.get(position_col), position_type)
        if nt0 is None or nt0 < 0 or nt0 >= record.nt_length:
            qc_rows.append(
                {
                    "input_row": row_index,
                    "transcript_id": transcript_id,
                    "status": "invalid_or_out_of_range_position",
                }
            )
            continue
        base = record.sequence[nt0]
        if require_a and position_type.startswith("cds_nt") and base != "A":
            qc_rows.append(
                {
                    "input_row": row_index,
                    "transcript_id": transcript_id,
                    "status": "editing_position_not_A",
                }
            )
            continue
        level = pd.to_numeric(
            pd.Series([row.get(editing_level_col)]), errors="coerce"
        ).iloc[0]
        if not np.isfinite(level):
            qc_rows.append(
                {
                    "input_row": row_index,
                    "transcript_id": transcript_id,
                    "status": "invalid_editing_level",
                }
            )
            continue
        sample = (
            str(row.get(editing_sample_col, "")).strip() if editing_sample_col else ""
        )
        if not sample or sample.lower() == "nan":
            sample = f"__record_{row_index}"
        records.append(
            {
                "group": "cds_edited",
                "gene_id": record.gene_id,
                "transcript_id": record.transcript_id,
                "site_nt_0based": nt0,
                "site_nt_1based": nt0 + 1,
                "site_aa_index0": nt0 // 3,
                "site_aa_1based": nt0 // 3 + 1,
                "site_base": base,
                "editing_sample": sample,
                "editing_level": float(level),
                "input_row": row_index,
            }
        )

    records_df = pd.DataFrame(records)
    if records_df.empty:
        return records_df, pd.DataFrame(qc_rows), records_df
    site_keys = [
        "group",
        "gene_id",
        "transcript_id",
        "site_nt_0based",
        "site_nt_1based",
        "site_aa_index0",
        "site_aa_1based",
        "site_base",
    ]
    per_sample = records_df.groupby(
        site_keys + ["editing_sample"], as_index=False, dropna=False
    ).agg(
        editing_level_sample_mean=("editing_level", "mean"),
        editing_level_sample_n_records=("editing_level", "size"),
        first_input_row=("input_row", "min"),
    )
    all_sites = per_sample.groupby(site_keys, as_index=False, dropna=False).agg(
        editing_level_mean=("editing_level_sample_mean", "mean"),
        editing_level_sd_across_samples=("editing_level_sample_mean", "std"),
        editing_level_n_samples=("editing_sample", "nunique"),
        editing_level_n_records=("editing_level_sample_n_records", "sum"),
        input_row=("first_input_row", "min"),
    )
    all_sites["site_id"] = (
        "EDIT_"
        + all_sites["transcript_id"].astype(str)
        + "_"
        + all_sites["site_nt_1based"].astype(str)
    )
    ranked = all_sites.sort_values(
        [
            "gene_id",
            "editing_level_mean",
            "editing_level_n_samples",
            "transcript_id",
            "site_nt_0based",
        ],
        ascending=[True, False, False, True, True],
        kind="mergesort",
    )
    selected = ranked.groupby("gene_id", sort=False, as_index=False).head(1).copy()
    return (
        selected.sort_values(["gene_id", "site_nt_0based"]).reset_index(drop=True),
        pd.DataFrame(qc_rows),
        all_sites.sort_values(["gene_id", "site_nt_0based"]).reset_index(drop=True),
    )
