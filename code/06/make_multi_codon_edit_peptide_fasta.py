#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate peptide FASTA records for multiple A-to-G edits in one codon.

The ordinary peptide generator creates one edited peptide per editing event.
This script adds the missing joint state: when two or more entries in a
missense.pos file belong to the same transcript and codon, all those edits
are applied simultaneously and the resulting amino acid is used to make one
combined edited peptide.

The default FASTA header follows the existing peptide-database convention:

    FBpp_start_edit1_edit2_end

For two edits affecting the same amino-acid position this intentionally looks
like, for example, ``FBpp0071261_1008_1030_1030_1030``.  The distinct CDS
coordinates and site IDs are retained in the summary and QC tables.

The input missense.pos format does not contain the reference nucleotide
codon.  Therefore, unless --cds-fasta is supplied, the reference codon is
inferred by enumerating codons whose translation and individual A-to-G
mutations agree with refAA and each altAA in the input.  Ambiguous or
inconsistent groups are reported and skipped rather than guessed.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

AA3_TO_AA1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Sec": "U",
    "Pyl": "O",
    "Ter": "*",
    "Stop": "*",
}

CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}
CODONS = tuple(CODON_TABLE)
TRANSCRIPT_RE = re.compile(r"(?:^|[|;\s])(FBtr\d+)(?:$|[|;\s])")
AA_PATTERN = re.compile(r"^p\.([A-Za-z*]+?)([0-9]+)([A-Za-z*]+)$")
MISSENSE_TYPES = {
    "missense_variant",
    "missense_variant&splice_region_variant",
}


def aa_to_one(value: str) -> str:
    value = value.strip()
    if len(value) == 1:
        return value.upper()
    if value in AA3_TO_AA1:
        return AA3_TO_AA1[value]
    raise ValueError(f"Unknown amino-acid code: {value}")


def read_tsv(path: Path) -> Iterator[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


@dataclass(frozen=True)
class Edit:
    site_id: str
    fbtr: str
    cds_pos: int
    aa_pos: int
    ref: str
    alt: str


@dataclass(frozen=True)
class Peptide:
    sequence: str
    start: int
    end: int
    missed_cleavages: str
    enzyme: str


def load_fbtr_to_fbpp(path: Path) -> Dict[str, str]:
    mapping: Dict[str, Set[str]] = defaultdict(set)
    for row in read_tsv(path):
        fbtr = row.get("FBtr", "").strip()
        fbpp = row.get("FBpp", "").strip()
        if fbtr and fbpp:
            mapping[fbtr].add(fbpp)
    # The supplied complete FlyBase map is expected to be one-to-one.  Keep
    # one deterministic value if a future map contains a duplicate row.
    return {fbtr: sorted(fbpps)[0] for fbtr, fbpps in mapping.items()}


def load_edits(path: Path) -> List[Edit]:
    edits: List[Edit] = []
    seen: Set[Tuple[str, str, int, int, str, str]] = set()
    required = {"site_id", "FBtr", "cds_pos", "aa_pos", "refAA", "altAA"}
    for row in read_tsv(path):
        missing = required.difference(row)
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        edit = Edit(
            site_id=row["site_id"].strip(),
            fbtr=row["FBtr"].strip(),
            cds_pos=int(row["cds_pos"]),
            aa_pos=int(row["aa_pos"]),
            ref=aa_to_one(row["refAA"]),
            alt=aa_to_one(row["altAA"]),
        )
        key = (edit.site_id, edit.fbtr, edit.cds_pos, edit.aa_pos, edit.ref, edit.alt)
        if key not in seen:
            edits.append(edit)
            seen.add(key)
    return edits


def split_position(value: str) -> Optional[Tuple[str, str]]:
    parts = value.strip().split("/", 1)
    if len(parts) != 2 or not parts[0].isdigit() or not parts[1].isdigit():
        return None
    return parts[0], parts[1]


def parse_ann_record(
    ann: str,
    effect: str,
    expected_transcript: str,
) -> Optional[Dict[str, str]]:
    """Parse the same ANN fields used to make missense_pos_new."""
    candidates: List[Dict[str, str]] = []
    for raw_record in ann.split(","):
        fields = raw_record.split("|")
        if len(fields) < 14 or fields[1] != effect:
            continue
        transcript = fields[6]
        if not transcript.startswith("FBtr"):
            continue
        trans = split_position(fields[11])
        cds = split_position(fields[12])
        aa = split_position(fields[13])
        protein = AA_PATTERN.match(fields[10])
        if trans is None or cds is None or aa is None or protein is None:
            continue
        candidates.append(
            {
                "gene": fields[3],
                "FBgn": fields[4],
                "FBtr": transcript,
                "refAA": protein.group(1),
                "altAA": protein.group(3),
                "trans_pos": trans[0],
                "trans_len": trans[1],
                "cds_pos": cds[0],
                "cds_len": cds[1],
                "aa_pos": aa[0],
                "aa_len": aa[1],
            }
        )
    if not candidates:
        return None
    if expected_transcript:
        matching = [row for row in candidates if row["FBtr"] == expected_transcript]
        if matching:
            return matching[0]
    return candidates[0]


def load_annotation_source(path: Path) -> List[Edit]:
    """Read the raw long annotation table used to create missense_pos_new.

    This is deliberately separate from load_edits(): the raw table retains
    every original candidate row, while missense.pos is a derived one-row-per-
    site table.  Rows from different libraries are collapsed by handle using
    the same deterministic rule as convert_candidate_to_missense_pos.py.
    """
    records: Dict[str, Edit] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"handle", "type", "ann"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        transcript_column = next(
            (
                name
                for name in ("transcript", "transcript_id")
                if name in (reader.fieldnames or [])
            ),
            None,
        )
        if transcript_column is None:
            raise ValueError(f"{path}: missing transcript/transcript_id column")

        for row in reader:
            effect = row.get("type", "").strip()
            if effect not in MISSENSE_TYPES:
                continue
            # The source file is an A-to-G candidate table.  Keep this guard
            # so a mixed annotation table cannot silently introduce another
            # nucleotide substitution into the joint-edit logic.
            change2 = row.get("change2", "").strip()
            if change2 and change2 != "A>G":
                continue
            parsed = parse_ann_record(
                row.get("ann", ""),
                effect,
                row.get(transcript_column, "").strip(),
            )
            if parsed is None:
                continue
            try:
                edit = Edit(
                    site_id=row["handle"].strip(),
                    fbtr=parsed["FBtr"],
                    cds_pos=int(parsed["cds_pos"]),
                    aa_pos=int(parsed["aa_pos"]),
                    ref=aa_to_one(parsed["refAA"]),
                    alt=aa_to_one(parsed["altAA"]),
                )
            except (KeyError, ValueError) as exc:
                raise ValueError(
                    f"{path}: invalid row for handle {row.get('handle', '')}"
                ) from exc
            old = records.get(edit.site_id)
            if old is None:
                records[edit.site_id] = edit
            elif old != edit:
                # Match the conversion script's deterministic collapse of
                # multiple annotations for one handle.
                records[edit.site_id] = min(
                    (old, edit),
                    key=lambda value: (
                        value.site_id,
                        value.fbtr,
                        value.cds_pos,
                        value.aa_pos,
                        value.ref,
                        value.alt,
                    ),
                )
    return [records[key] for key in sorted(records)]


def load_peptides(path: Path, wanted_fbpps: Set[str]) -> Dict[str, List[Peptide]]:
    peptides: Dict[str, List[Peptide]] = defaultdict(list)
    for row in read_tsv(path):
        fbpp = row.get("source_id", "").strip()
        if fbpp not in wanted_fbpps:
            continue
        peptides[fbpp].append(
            Peptide(
                sequence=row["peptide_seq"].strip(),
                start=int(row["aa_start_1based"]),
                end=int(row["aa_end_1based"]),
                missed_cleavages=row.get("missed_cleavages", ""),
                enzyme=row.get("enzyme", ""),
            )
        )
    return peptides


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    sequence: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:].split()[0]
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def load_reference_cds(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    result: Dict[str, str] = {}
    for header, sequence in read_fasta(path):
        match = re.search(r"FBtr\d+", header)
        if match is None:
            continue
        sequence = sequence.upper().replace("U", "T")
        if re.search(r"[^ACGTN]", sequence):
            raise ValueError(
                f"{path} does not look like a nucleotide CDS FASTA: {header}"
            )
        result[match.group(0)] = sequence
    return result


def mutate_codon(codon: str, offsets: Iterable[int]) -> Optional[str]:
    chars = list(codon)
    for offset in offsets:
        if chars[offset] != "A":
            return None
        chars[offset] = "G"
    return "".join(chars)


def infer_combined_aa(
    events: Sequence[Edit],
    reference_cds: Dict[str, str],
) -> Tuple[Optional[str], str, str]:
    """Return combined AA, method, and detail/QC message."""
    refs = {event.ref for event in events}
    aa_positions = {event.aa_pos for event in events}
    offsets = [(event.cds_pos - 1) % 3 for event in events]
    if len(refs) != 1:
        return None, "", "different_refAA_within_codon"
    if len(aa_positions) != 1:
        return None, "", "different_aa_pos_within_codon"
    if len(offsets) != len(set(offsets)):
        return None, "", "duplicate_nucleotide_position"
    ref_aa = next(iter(refs))

    if events[0].fbtr in reference_cds:
        cds = reference_cds[events[0].fbtr]
        codon_index = (events[0].cds_pos - 1) // 3
        start = codon_index * 3
        codon = cds[start : start + 3]
        if len(codon) != 3 or re.search(r"[^ACGT]", codon):
            return None, "reference_cds", "reference_codon_missing_or_ambiguous"
        if CODON_TABLE.get(codon) != ref_aa:
            return (
                None,
                "reference_cds",
                f"reference_codon_translates_{CODON_TABLE.get(codon)}_not_{ref_aa}",
            )
        for event, offset in zip(events, offsets):
            single = mutate_codon(codon, [offset])
            if single is None or CODON_TABLE.get(single) != event.alt:
                return (
                    None,
                    "reference_cds",
                    f"single_edit_does_not_match_{event.site_id}",
                )
        combined_codon = mutate_codon(codon, offsets)
        if combined_codon is None:
            return None, "reference_cds", "reference_base_is_not_A"
        combined = CODON_TABLE.get(combined_codon)
        if combined is None:
            return None, "reference_cds", "combined_codon_unknown"
        return combined, "reference_cds", codon + ">" + combined_codon

    candidates: List[Tuple[str, str]] = []
    for codon in CODONS:
        if CODON_TABLE[codon] != ref_aa:
            continue
        combined_codon = mutate_codon(codon, offsets)
        if combined_codon is None:
            continue
        valid = True
        for event, offset in zip(events, offsets):
            single = mutate_codon(codon, [offset])
            if single is None or CODON_TABLE.get(single) != event.alt:
                valid = False
                break
        if valid:
            candidates.append((codon, combined_codon))
    combined_aas = {CODON_TABLE[codon] for _, codon in candidates}
    if not candidates:
        return (
            None,
            "inferred_from_single_edits",
            "no_reference_codon_matches_single_edits",
        )
    if len(combined_aas) != 1:
        details = ";".join(f"{ref}>{mut}" for ref, mut in candidates)
        return None, "inferred_from_single_edits", "ambiguous_combined_AA:" + details
    combined = next(iter(combined_aas))
    return (
        combined,
        "inferred_from_single_edits",
        ";".join(f"{ref}>{mut}" for ref, mut in candidates),
    )


def group_multi_edits(edits: Sequence[Edit]) -> Dict[Tuple[str, int], List[Edit]]:
    grouped: Dict[Tuple[str, int], List[Edit]] = defaultdict(list)
    for edit in edits:
        codon = (edit.cds_pos - 1) // 3 + 1
        grouped[(edit.fbtr, codon)].append(edit)
    return {
        key: sorted(value, key=lambda event: (event.cds_pos, event.site_id))
        for key, value in grouped.items()
        if len(value) >= 2
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate combined edited peptides for multiple edits in one codon."
    )
    parser.add_argument("--map", required=True, type=Path)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--edits",
        type=Path,
        help="Derived missense.pos table.",
    )
    input_group.add_argument(
        "--source-annotations",
        type=Path,
        help="Raw *.candidate.AtoG.annotations.long.tsv(.new) used to make missense_pos_new.",
    )
    parser.add_argument("--peptides", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument(
        "--cds-fasta",
        type=Path,
        default=None,
        help="Optional reference nucleotide CDS FASTA. Without it, infer codons from single-edit AA changes.",
    )
    parser.add_argument(
        "--header-id",
        choices=("fbpp", "fbtr"),
        default="fbpp",
        help="Identifier used as the first FASTA header field. Default: fbpp.",
    )
    parser.add_argument(
        "--header-position",
        choices=("aa", "cds"),
        default="aa",
        help="Use amino-acid or CDS positions for edit1/edit2 in the header. Default: aa.",
    )
    parser.add_argument(
        "--allow-stop",
        action="store_true",
        help="Allow a combined edit that translates to a stop codon. By default it is reported and skipped.",
    )
    parser.add_argument(
        "--no-strict-ref-check",
        action="store_true",
        help="Do not skip a peptide when its residue does not match the reference AA.",
    )
    return parser.parse_args()


def append_suffix(prefix: Path, suffix: str) -> Path:
    """Append a result suffix without replacing dots in the prefix name."""
    return prefix.parent / f"{prefix.name}{suffix}"


def main() -> None:
    args = parse_args()
    fbtr_to_fbpp = load_fbtr_to_fbpp(args.map)
    if args.source_annotations is not None:
        edits = load_annotation_source(args.source_annotations)
        input_description = f"raw annotation source: {args.source_annotations}"
    else:
        edits = load_edits(args.edits)
        input_description = f"missense.pos: {args.edits}"
    reference_cds = load_reference_cds(args.cds_fasta)
    groups = group_multi_edits(edits)

    wanted_fbpps = {fbtr_to_fbpp[fbtr] for fbtr, _ in groups if fbtr in fbtr_to_fbpp}
    peptides_by_fbpp = load_peptides(args.peptides, wanted_fbpps)

    edit_fasta = append_suffix(args.out_prefix, ".edit.fasta")
    summary_tsv = append_suffix(args.out_prefix, ".summary.tsv")
    skipped_tsv = append_suffix(args.out_prefix, ".skipped.tsv")
    groups_tsv = append_suffix(args.out_prefix, ".codon_groups.tsv")
    for path in (edit_fasta, summary_tsv, skipped_tsv, groups_tsv):
        path.parent.mkdir(parents=True, exist_ok=True)

    summary_fields = [
        "FBtr",
        "FBpp",
        "codon_index",
        "cds_positions",
        "aa_positions",
        "site_ids",
        "refAA",
        "individual_altAAs",
        "combinedAA",
        "pep_start",
        "pep_end",
        "local_pos_in_peptide",
        "header",
        "wt_peptide",
        "edit_peptide",
        "missed_cleavages",
        "enzyme",
        "inference_method",
        "inference_detail",
    ]
    skip_fields = [
        "FBtr",
        "FBpp",
        "codon_index",
        "cds_positions",
        "aa_positions",
        "site_ids",
        "refAA",
        "individual_altAAs",
        "reason",
        "details",
    ]
    group_fields = [
        "FBtr",
        "FBpp",
        "codon_index",
        "cds_positions",
        "aa_positions",
        "site_ids",
        "refAA",
        "individual_altAAs",
        "combinedAA",
        "status",
        "details",
    ]

    n_records = 0
    with edit_fasta.open("w", encoding="utf-8") as fasta_out, summary_tsv.open(
        "w", encoding="utf-8", newline=""
    ) as summary_out, skipped_tsv.open(
        "w", encoding="utf-8", newline=""
    ) as skipped_out, groups_tsv.open(
        "w", encoding="utf-8", newline=""
    ) as groups_out:
        summary_writer = csv.DictWriter(
            summary_out, fieldnames=summary_fields, delimiter="\t", lineterminator="\n"
        )
        skipped_writer = csv.DictWriter(
            skipped_out, fieldnames=skip_fields, delimiter="\t", lineterminator="\n"
        )
        group_writer = csv.DictWriter(
            groups_out, fieldnames=group_fields, delimiter="\t", lineterminator="\n"
        )
        summary_writer.writeheader()
        skipped_writer.writeheader()
        group_writer.writeheader()

        for (fbtr, codon_index), events in sorted(groups.items()):
            fbpp = fbtr_to_fbpp.get(fbtr, "")
            cds_positions = sorted({event.cds_pos for event in events})
            aa_positions = sorted({event.aa_pos for event in events})
            site_ids = sorted({event.site_id for event in events})
            ref_aas = sorted({event.ref for event in events})
            alt_aas = [event.alt for event in events]
            common = {
                "FBtr": fbtr,
                "FBpp": fbpp,
                "codon_index": codon_index,
                "cds_positions": ",".join(map(str, cds_positions)),
                "aa_positions": ",".join(map(str, aa_positions)),
                "site_ids": ",".join(site_ids),
                "refAA": ",".join(ref_aas),
                "individual_altAAs": ",".join(alt_aas),
            }
            if not fbpp:
                group_writer.writerow(
                    {
                        **common,
                        "combinedAA": "",
                        "status": "skipped",
                        "details": "FBtr_not_found_in_map",
                    }
                )
                skipped_writer.writerow(
                    {**common, "reason": "FBtr_not_found_in_map", "details": ""}
                )
                continue
            combined_aa, method, detail = infer_combined_aa(events, reference_cds)
            if combined_aa is None:
                group_writer.writerow(
                    {**common, "combinedAA": "", "status": "skipped", "details": detail}
                )
                skipped_writer.writerow(
                    {**common, "reason": "cannot_infer_combined_AA", "details": detail}
                )
                continue
            if combined_aa == "*" and not args.allow_stop:
                group_writer.writerow(
                    {
                        **common,
                        "combinedAA": combined_aa,
                        "status": "skipped",
                        "details": "combined_stop_codon",
                    }
                )
                skipped_writer.writerow(
                    {**common, "reason": "combined_stop_codon", "details": detail}
                )
                continue

            group_writer.writerow(
                {
                    **common,
                    "combinedAA": combined_aa,
                    "status": "candidate",
                    "details": detail,
                }
            )
            aa_pos = aa_positions[0]
            seen_records: Set[Tuple[str, int, int, str]] = set()
            for peptide in peptides_by_fbpp.get(fbpp, []):
                if not (peptide.start <= aa_pos <= peptide.end):
                    continue
                local_pos = aa_pos - peptide.start + 1
                observed_ref = peptide.sequence[local_pos - 1]
                ref_aa = events[0].ref
                if observed_ref != ref_aa and not args.no_strict_ref_check:
                    skipped_writer.writerow(
                        {
                            **common,
                            "reason": "refAA_mismatch",
                            "details": f"peptide_has={observed_ref};expected={ref_aa}",
                        }
                    )
                    continue
                edit_peptide = (
                    peptide.sequence[: local_pos - 1]
                    + combined_aa
                    + peptide.sequence[local_pos:]
                )
                key = (peptide.sequence, peptide.start, peptide.end, edit_peptide)
                if key in seen_records:
                    continue
                seen_records.add(key)
                labels = [
                    str(event.aa_pos if args.header_position == "aa" else event.cds_pos)
                    for event in events
                ]
                name_id = fbpp if args.header_id == "fbpp" else fbtr
                header = f"{name_id}_{peptide.start}_{'_'.join(labels)}_{peptide.end}"
                fasta_out.write(f">{header}\n{edit_peptide}\n")
                summary_writer.writerow(
                    {
                        **common,
                        "combinedAA": combined_aa,
                        "pep_start": peptide.start,
                        "pep_end": peptide.end,
                        "local_pos_in_peptide": local_pos,
                        "header": header,
                        "wt_peptide": peptide.sequence,
                        "edit_peptide": edit_peptide,
                        "missed_cleavages": peptide.missed_cleavages,
                        "enzyme": peptide.enzyme,
                        "inference_method": method,
                        "inference_detail": detail,
                    }
                )
                n_records += 1
            if not seen_records:
                skipped_writer.writerow(
                    {**common, "reason": "no_covering_peptide", "details": ""}
                )

    candidate_groups = sum(
        1 for row in read_tsv(groups_tsv) if row["status"] == "candidate"
    )
    print(f"Input: {input_description}")
    print(f"Input missense sites: {len(edits)}")
    print(f"Multi-edit codon groups: {len(groups)}")
    print(f"Candidate groups: {candidate_groups}")
    print(f"Combined edited peptide records: {n_records}")
    print(f"EDIT FASTA: {edit_fasta}")
    print(f"SUMMARY: {summary_tsv}")
    print(f"SKIPPED: {skipped_tsv}")
    print(f"CODON GROUPS: {groups_tsv}")


if __name__ == "__main__":
    main()
