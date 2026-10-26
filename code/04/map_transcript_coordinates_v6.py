#!/usr/bin/env python3
"""
Map transcript-coordinate BED intervals to genomic coordinates and project the
same genomic bases onto every other transcript of the same gene. If no other
transcript overlaps the query bases, the source transcript is emitted as an
automatic fallback so the query remains in the complete mapping output.

Input transcript BED uses standard BED coordinates:
    transcript_id  start  end  [name]
where start is 0-based and end is half-open. Coordinates include UTRs because
mapping is based on GTF exon records, not CDS records.

Outputs:
1) <prefix>.genome.bed
   One BED row per genomic exon block of each query interval.
2) <prefix>.source_coordinates.tsv
   Source transcript, transcript-coordinate, CDS-coordinate, and genomic blocks.
3) <prefix>.other_transcript_blocks.tsv
   A complete source-to-target mapping table. Each row contains the exact source
   transcript/CDS coordinates, genomic coordinates, and the projected target
   transcript/CDS coordinates for the same genomic segment.
4) <prefix>.other_transcript_summary.tsv
   One row per query/target-transcript pair, including FULL/PARTIAL/NONE status,
   target CDS length, and a summary of projected regions/CDS blocks. By default,
   the source transcript is added as a target only when no other transcript
   overlaps the query interval. Use --include-source to always include it.

All transcript and CDS coordinates in the output are 0-based and half-open.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass
class Exon:
    chrom: str
    start: int  # genomic, 0-based inclusive
    end: int  # genomic, 0-based exclusive
    strand: str
    tx_start: int = 0
    tx_end: int = 0

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: List[Exon] = field(default_factory=list)
    cds_genomic_intervals: List[Tuple[int, int]] = field(default_factory=list)
    cds_tx_intervals: List[Tuple[int, int]] = field(default_factory=list)
    length: int = 0
    cds_length: int = 0

    def finalize(self) -> None:
        # Remove exact duplicate exon records, which can occur in some generated GTFs.
        unique_exons = {}
        for exon in self.exons:
            unique_exons[(exon.chrom, exon.start, exon.end, exon.strand)] = exon
        self.exons = list(unique_exons.values())

        if self.strand == "+":
            self.exons.sort(key=lambda x: (x.start, x.end))
        elif self.strand == "-":
            self.exons.sort(key=lambda x: (x.start, x.end), reverse=True)
        else:
            raise ValueError(
                f"Unsupported strand {self.strand!r} for {self.transcript_id}"
            )

        cursor = 0
        for exon in self.exons:
            exon.tx_start = cursor
            exon.tx_end = cursor + exon.length
            cursor = exon.tx_end
        self.length = cursor

        # Project GTF CDS records to mature-transcript coordinates. The resulting
        # coordinates are oriented 5' to 3' for both positive and negative strands.
        cds_tx_parts: List[Tuple[int, int]] = []
        for cds_start, cds_end in sorted(set(self.cds_genomic_intervals)):
            for exon in self.exons:
                ov_start = max(cds_start, exon.start)
                ov_end = min(cds_end, exon.end)
                if ov_start >= ov_end:
                    continue
                if self.strand == "+":
                    t_start = exon.tx_start + (ov_start - exon.start)
                    t_end = exon.tx_start + (ov_end - exon.start)
                else:
                    t_start = exon.tx_start + (exon.end - ov_end)
                    t_end = exon.tx_start + (exon.end - ov_start)
                cds_tx_parts.append((t_start, t_end))

        self.cds_tx_intervals = merge_half_open(cds_tx_parts)
        self.cds_length = sum(end - start for start, end in self.cds_tx_intervals)


@dataclass
class Query:
    transcript_id: str
    start: int
    end: int
    name: str
    line_no: int

    @property
    def length(self) -> int:
        return self.end - self.start


ATTR_RE = re.compile(r'\s*([^\s;]+)\s+"([^"]*)"\s*')


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attr_text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        match = ATTR_RE.fullmatch(item)
        if match:
            attrs[match.group(1)] = match.group(2)
            continue
        # Also tolerate key=value GFF3-like attributes.
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
    return attrs


def normalize_id(value: str, strip_version: bool) -> str:
    value = value.strip()
    if strip_version:
        value = re.sub(r"\.\d+$", "", value)
    return value


def load_gtf(
    gtf_path: str,
    strip_version: bool = False,
) -> Tuple[Dict[str, Transcript], Dict[str, List[str]], Dict[str, List[Transcript]]]:
    """
    Load exon annotations from a GTF.

    A transcript_id can occasionally occur at more than one gene/chromosome/strand
    combination in malformed or merged annotation files. Such IDs are retained in
    ``ambiguous_transcripts`` and excluded from coordinate mapping rather than
    causing the entire run to abort or silently selecting one locus.
    """
    candidates: Dict[str, Dict[Tuple[str, str, str], Transcript]] = {}

    with open_text(gtf_path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {line_no} does not have 9 columns")

            chrom, _, feature, start1, end1, _, strand, _, attr_text = fields
            if feature not in {"exon", "CDS"}:
                continue
            if strand not in {"+", "-"}:
                continue

            attrs = parse_attributes(attr_text)
            raw_tx = (
                attrs.get("transcript_id")
                or attrs.get("transcript")
                or attrs.get("Parent")
            )
            raw_gene = (
                attrs.get("gene_id") or attrs.get("gene") or attrs.get("gene_name")
            )
            if raw_tx is None or raw_gene is None:
                continue

            tx_id = normalize_id(raw_tx, strip_version)
            gene_id = normalize_id(raw_gene, strip_version)
            start0 = int(start1) - 1
            end0 = int(end1)
            if end0 <= start0:
                raise ValueError(
                    f"Invalid exon coordinates at GTF line {line_no}: {chrom}:{start1}-{end1}"
                )

            locus_key = (gene_id, chrom, strand)
            tx_candidates = candidates.setdefault(tx_id, {})
            if locus_key not in tx_candidates:
                tx_candidates[locus_key] = Transcript(
                    transcript_id=tx_id,
                    gene_id=gene_id,
                    chrom=chrom,
                    strand=strand,
                )
            tx_candidate = tx_candidates[locus_key]
            if feature == "exon":
                tx_candidate.exons.append(Exon(chrom, start0, end0, strand))
            else:
                tx_candidate.cds_genomic_intervals.append((start0, end0))

    transcripts: Dict[str, Transcript] = {}
    ambiguous_transcripts: Dict[str, List[Transcript]] = {}

    for tx_id, by_locus in candidates.items():
        tx_list = list(by_locus.values())
        if len(tx_list) != 1:
            ambiguous_transcripts[tx_id] = tx_list
            continue
        tx = tx_list[0]
        if not tx.exons:
            continue
        tx.finalize()
        transcripts[tx_id] = tx

    gene_to_transcripts: Dict[str, List[str]] = {}
    for tx_id, tx in transcripts.items():
        gene_to_transcripts.setdefault(tx.gene_id, []).append(tx_id)

    for gene_id in gene_to_transcripts:
        gene_to_transcripts[gene_id].sort()

    return transcripts, gene_to_transcripts, ambiguous_transcripts


def load_queries(bed_path: str, strip_version: bool = False) -> List[Query]:
    """Load transcript-coordinate BED-like input.

    Standard BED is tab-delimited and normally has no header. For convenience,
    this reader also accepts runs of spaces as delimiters and skips a common
    header such as: transcript_id start end name.
    """
    queries: List[Query] = []
    with open_text(bed_path) as handle:
        for line_no, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            # Prefer true BED tab separation, but tolerate whitespace-delimited
            # files produced by tables, copied text, or awk without OFS="\t".
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                fields = stripped.split()

            if len(fields) < 3:
                raise ValueError(
                    f"BED line {line_no} has fewer than 3 columns: {stripped!r}"
                )

            first3 = [field.strip().lower() for field in fields[:3]]
            header_col1 = {"transcript", "transcript_id", "transcriptid", "tx", "tx_id"}
            header_col2 = {"start", "tx_start", "transcript_start"}
            header_col3 = {"end", "stop", "tx_end", "transcript_end"}
            if (
                first3[0] in header_col1
                and first3[1] in header_col2
                and first3[2] in header_col3
            ):
                continue

            tx_id = normalize_id(fields[0].strip(), strip_version)
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                raise ValueError(
                    f"BED line {line_no} has non-integer start/end: {stripped!r}. "
                    "Expected: transcript_id start end [name]"
                )

            if start < 0 or end <= start:
                raise ValueError(
                    f"Invalid BED interval at line {line_no}: start={start}, end={end}"
                )

            name = (
                fields[3].strip()
                if len(fields) >= 4 and fields[3].strip()
                else f"query_{line_no}"
            )
            queries.append(Query(tx_id, start, end, name, line_no))
    return queries


def transcript_interval_to_genome(
    tx: Transcript, start: int, end: int
) -> List[Tuple[int, int, int]]:
    """Return [(genomic_start, genomic_end, source_block_index), ...]."""
    if start < 0 or end > tx.length or start >= end:
        raise ValueError(
            f"Interval [{start}, {end}) is outside transcript {tx.transcript_id} "
            f"length {tx.length}"
        )

    blocks: List[Tuple[int, int, int]] = []
    block_index = 0
    for exon in tx.exons:
        ov_start = max(start, exon.tx_start)
        ov_end = min(end, exon.tx_end)
        if ov_start >= ov_end:
            continue

        rel_start = ov_start - exon.tx_start
        rel_end = ov_end - exon.tx_start
        if tx.strand == "+":
            g_start = exon.start + rel_start
            g_end = exon.start + rel_end
        else:
            g_start = exon.end - rel_end
            g_end = exon.end - rel_start

        block_index += 1
        blocks.append((g_start, g_end, block_index))

    mapped = sum(e - s for s, e, _ in blocks)
    if mapped != end - start:
        raise RuntimeError(
            f"Internal mapping error for {tx.transcript_id}: mapped {mapped} of {end-start} bases"
        )
    return blocks


def genomic_interval_to_transcript_blocks(
    tx: Transcript,
    genomic_blocks: Iterable[Tuple[int, int, int]],
) -> List[Tuple[int, int, int, int, int]]:
    """
    Project genomic blocks to a transcript.

    Return tuples:
        (target_tx_start, target_tx_end, genomic_start, genomic_end, source_block_index)
    """
    out: List[Tuple[int, int, int, int, int]] = []

    for g_start, g_end, source_block_index in genomic_blocks:
        for exon in tx.exons:
            ov_start = max(g_start, exon.start)
            ov_end = min(g_end, exon.end)
            if ov_start >= ov_end:
                continue

            if tx.strand == "+":
                t_start = exon.tx_start + (ov_start - exon.start)
                t_end = exon.tx_start + (ov_end - exon.start)
            else:
                t_start = exon.tx_start + (exon.end - ov_end)
                t_end = exon.tx_start + (exon.end - ov_start)

            out.append((t_start, t_end, ov_start, ov_end, source_block_index))

    out.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    return out


def merge_half_open(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def cds_boundary_to_offset(tx: Transcript, position: int) -> Optional[int]:
    """Return the 0-based CDS offset of a transcript-coordinate boundary.

    A boundary at the end of a CDS interval is valid and returns the offset just
    after the last CDS base in that interval. None means that the position is not
    a boundary within a CDS interval.
    """
    offset = 0
    for start, end in tx.cds_tx_intervals:
        if start <= position <= end:
            return offset + (position - start)
        offset += end - start
    return None


def classify_tx_segment(
    tx: Transcript, start: int, end: int
) -> Tuple[str, Optional[int], Optional[int]]:
    """Classify a target transcript segment and return CDS-relative coordinates."""
    if not tx.cds_tx_intervals:
        return "NONCODING", None, None

    for cds_start, cds_end in tx.cds_tx_intervals:
        if start >= cds_start and end <= cds_end:
            c_start = cds_boundary_to_offset(tx, start)
            c_end = cds_boundary_to_offset(tx, end)
            return "CDS", c_start, c_end

    first_cds_start = tx.cds_tx_intervals[0][0]
    last_cds_end = tx.cds_tx_intervals[-1][1]
    if end <= first_cds_start:
        return "UTR5", None, None
    if start >= last_cds_end:
        return "UTR3", None, None
    return "NON_CDS_INTERNAL", None, None


def split_mapped_block_by_cds(
    tx: Transcript,
    block: Tuple[int, int, int, int, int],
) -> List[Tuple[int, int, int, int, int, str, Optional[int], Optional[int]]]:
    """Split one projected block at CDS boundaries.

    Return tuples:
      target_tx_start, target_tx_end, genomic_start, genomic_end,
      source_block_index, region, target_cds_start, target_cds_end
    """
    t_start, t_end, g_start, g_end, source_block_index = block
    boundaries = {t_start, t_end}
    for cds_start, cds_end in tx.cds_tx_intervals:
        if t_start < cds_start < t_end:
            boundaries.add(cds_start)
        if t_start < cds_end < t_end:
            boundaries.add(cds_end)

    ordered = sorted(boundaries)
    out = []
    for seg_start, seg_end in zip(ordered, ordered[1:]):
        if seg_start >= seg_end:
            continue
        region, cds_start, cds_end = classify_tx_segment(tx, seg_start, seg_end)

        if tx.strand == "+":
            seg_g_start = g_start + (seg_start - t_start)
            seg_g_end = g_start + (seg_end - t_start)
        else:
            seg_g_start = g_end - (seg_end - t_start)
            seg_g_end = g_end - (seg_start - t_start)

        out.append(
            (
                seg_start,
                seg_end,
                seg_g_start,
                seg_g_end,
                source_block_index,
                region,
                cds_start,
                cds_end,
            )
        )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Map transcript-coordinate BED intervals to the genome and to all other "
            "transcripts of the same gene using GTF exon annotations."
        )
    )
    parser.add_argument("--gtf", required=True, help="GTF or GTF.gz file")
    parser.add_argument(
        "--transcript-bed",
        required=True,
        help="Transcript-coordinate BED: transcript_id, start, end, optional name",
    )
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    parser.add_argument(
        "--strip-version",
        action="store_true",
        help="Strip trailing version suffixes such as .1 from gene/transcript IDs",
    )
    parser.add_argument(
        "--include-source",
        action="store_true",
        help=(
            "Always include projection back to the source transcript. Without "
            "this option, the source is included only as a fallback when no "
            "other transcript overlaps the query interval."
        ),
    )
    args = parser.parse_args()

    transcripts, gene_to_transcripts, ambiguous_transcripts = load_gtf(
        args.gtf, args.strip_version
    )
    queries = load_queries(args.transcript_bed, args.strip_version)

    prefix = Path(args.out_prefix)
    genome_path = str(prefix) + ".genome.bed"
    source_path = str(prefix) + ".source_coordinates.tsv"
    blocks_path = str(prefix) + ".other_transcript_blocks.tsv"
    summary_path = str(prefix) + ".other_transcript_summary.tsv"
    error_path = str(prefix) + ".errors.tsv"
    conflict_path = str(prefix) + ".annotation_conflicts.tsv"

    n_ok = 0
    n_error = 0

    with open(conflict_path, "w", newline="") as conflict_f:
        conflict_w = csv.writer(conflict_f, delimiter="\t", lineterminator="\n")
        conflict_w.writerow(
            [
                "transcript_id",
                "gene_id",
                "chrom",
                "strand",
                "exon_count",
                "total_exon_bases",
                "genomic_start",
                "genomic_end",
            ]
        )
        for tx_id in sorted(ambiguous_transcripts):
            for tx in sorted(
                ambiguous_transcripts[tx_id],
                key=lambda x: (x.gene_id, x.chrom, x.strand),
            ):
                exon_keys = {(e.chrom, e.start, e.end, e.strand) for e in tx.exons}
                starts = [e[1] for e in exon_keys]
                ends = [e[2] for e in exon_keys]
                conflict_w.writerow(
                    [
                        tx_id,
                        tx.gene_id,
                        tx.chrom,
                        tx.strand,
                        len(exon_keys),
                        sum(e - s for _, s, e, _ in exon_keys),
                        min(starts) if starts else "",
                        max(ends) if ends else "",
                    ]
                )

    with open(genome_path, "w", newline="") as genome_f, open(
        source_path, "w", newline=""
    ) as source_f, open(blocks_path, "w", newline="") as blocks_f, open(
        summary_path, "w", newline=""
    ) as summary_f, open(
        error_path, "w", newline=""
    ) as error_f:
        genome_w = csv.writer(genome_f, delimiter="\t", lineterminator="\n")
        source_w = csv.writer(source_f, delimiter="\t", lineterminator="\n")
        blocks_w = csv.writer(blocks_f, delimiter="\t", lineterminator="\n")
        summary_w = csv.writer(summary_f, delimiter="\t", lineterminator="\n")
        error_w = csv.writer(error_f, delimiter="\t", lineterminator="\n")

        source_w.writerow(
            [
                "query_name",
                "gene_id",
                "source_transcript",
                "query_tx_start",
                "query_tx_end",
                "source_tx_start",
                "source_tx_end",
                "source_region",
                "source_cds_start",
                "source_cds_end",
                "source_transcript_length",
                "source_cds_length",
                "chrom",
                "genomic_start",
                "genomic_end",
                "strand",
                "source_genomic_block",
            ]
        )
        blocks_w.writerow(
            [
                "query_name",
                "gene_id",
                "source_transcript",
                "query_tx_start",
                "query_tx_end",
                "source_tx_start",
                "source_tx_end",
                "source_region",
                "source_cds_start",
                "source_cds_end",
                "source_transcript_length",
                "source_cds_length",
                "chrom",
                "genomic_start",
                "genomic_end",
                "strand",
                "source_genomic_block",
                "target_transcript",
                "target_tx_start",
                "target_tx_end",
                "target_region",
                "target_cds_start",
                "target_cds_end",
                "target_transcript_length",
                "target_cds_length",
            ]
        )
        summary_w.writerow(
            [
                "query_name",
                "gene_id",
                "source_transcript",
                "source_tx_start",
                "source_tx_end",
                "query_length",
                "source_transcript_length",
                "source_cds_length",
                "source_regions",
                "source_cds_blocks",
                "target_transcript",
                "target_transcript_length",
                "target_cds_length",
                "mapped_bases",
                "unmapped_bases",
                "mapping_status",
                "target_regions",
                "target_blocks",
                "target_cds_blocks",
            ]
        )
        error_w.writerow(
            ["line_no", "query_name", "transcript_id", "start", "end", "error"]
        )

        for query in queries:
            try:
                source = transcripts.get(query.transcript_id)
                if source is None:
                    if query.transcript_id in ambiguous_transcripts:
                        loci = ", ".join(
                            f"{tx.gene_id}:{tx.chrom}:{tx.strand}"
                            for tx in ambiguous_transcripts[query.transcript_id]
                        )
                        raise ValueError(
                            f"Ambiguous transcript_id in GTF: {query.transcript_id}; "
                            f"found at multiple loci/strands ({loci}). See {conflict_path}"
                        )
                    raise KeyError(
                        f"Transcript not found in GTF: {query.transcript_id}"
                    )

                genomic_blocks = transcript_interval_to_genome(
                    source, query.start, query.end
                )

                # BED6 + metadata. One row per exon block.
                for g_start, g_end, block_index in genomic_blocks:
                    genome_w.writerow(
                        [
                            source.chrom,
                            g_start,
                            g_end,
                            query.name,
                            0,
                            source.strand,
                            source.gene_id,
                            source.transcript_id,
                            query.start,
                            query.end,
                            block_index,
                        ]
                    )

                # Split the source interval at source-CDS boundaries and report
                # exact source transcript, CDS, and genomic coordinates.
                source_mapped_blocks = genomic_interval_to_transcript_blocks(
                    source, genomic_blocks
                )
                source_split_blocks = []
                for source_block in source_mapped_blocks:
                    source_split_blocks.extend(
                        split_mapped_block_by_cds(source, source_block)
                    )

                source_regions = ";".join(sorted({b[5] for b in source_split_blocks}))
                source_cds_block_text = ";".join(
                    f"{b[6]}-{b[7]}"
                    for b in source_split_blocks
                    if b[5] == "CDS" and b[6] is not None and b[7] is not None
                )

                for (
                    s_tx_start,
                    s_tx_end,
                    s_g_start,
                    s_g_end,
                    source_block_index,
                    source_region,
                    source_cds_start,
                    source_cds_end,
                ) in source_split_blocks:
                    source_w.writerow(
                        [
                            query.name,
                            source.gene_id,
                            source.transcript_id,
                            query.start,
                            query.end,
                            s_tx_start,
                            s_tx_end,
                            source_region,
                            "" if source_cds_start is None else source_cds_start,
                            "" if source_cds_end is None else source_cds_end,
                            source.length,
                            source.cds_length,
                            source.chrom,
                            s_g_start,
                            s_g_end,
                            source.strand,
                            source_block_index,
                        ]
                    )

                all_target_ids = gene_to_transcripts.get(source.gene_id, [])
                other_target_ids = [
                    target_id
                    for target_id in all_target_ids
                    if target_id != source.transcript_id
                ]

                # Default behavior: project to other transcripts. If none of the
                # other transcripts contains any queried genomic base, append the
                # source transcript as a FULL self-projection fallback. This keeps
                # source-only sites in the complete blocks output while avoiding
                # redundant self rows when a real alternative transcript matches.
                if args.include_source:
                    target_ids = list(all_target_ids)
                else:
                    has_other_overlap = False
                    for other_target_id in other_target_ids:
                        other_target = transcripts[other_target_id]
                        if (
                            other_target.chrom != source.chrom
                            or other_target.strand != source.strand
                        ):
                            continue
                        if genomic_interval_to_transcript_blocks(
                            other_target, genomic_blocks
                        ):
                            has_other_overlap = True
                            break

                    target_ids = list(other_target_ids)
                    if not has_other_overlap:
                        target_ids.append(source.transcript_id)

                for target_id in target_ids:
                    target = transcripts[target_id]

                    # Same gene should normally be on the same chromosome and strand.
                    if target.chrom != source.chrom or target.strand != source.strand:
                        summary_w.writerow(
                            [
                                query.name,
                                source.gene_id,
                                source.transcript_id,
                                query.start,
                                query.end,
                                query.length,
                                source.length,
                                source.cds_length,
                                source_regions,
                                source_cds_block_text,
                                target.transcript_id,
                                target.length,
                                target.cds_length,
                                0,
                                query.length,
                                "INCOMPATIBLE_LOCUS",
                                "",
                                "",
                                "",
                            ]
                        )
                        continue

                    mapped_blocks = genomic_interval_to_transcript_blocks(
                        target, genomic_blocks
                    )

                    # Calculate unique mapped genomic bases to avoid double counting in malformed annotations.
                    merged_genomic = merge_half_open(
                        (b[2], b[3]) for b in mapped_blocks
                    )
                    mapped_bases = sum(e - s for s, e in merged_genomic)
                    unmapped_bases = query.length - mapped_bases
                    status = (
                        "FULL"
                        if mapped_bases == query.length
                        else ("PARTIAL" if mapped_bases > 0 else "NONE")
                    )

                    split_blocks = []
                    for mapped_block in mapped_blocks:
                        split_blocks.extend(
                            split_mapped_block_by_cds(target, mapped_block)
                        )

                    target_block_text = ";".join(
                        f"{b[0]}-{b[1]}" for b in mapped_blocks
                    )
                    target_regions = ";".join(sorted({b[5] for b in split_blocks}))
                    target_cds_block_text = ";".join(
                        f"{b[6]}-{b[7]}"
                        for b in split_blocks
                        if b[5] == "CDS" and b[6] is not None and b[7] is not None
                    )
                    summary_w.writerow(
                        [
                            query.name,
                            source.gene_id,
                            source.transcript_id,
                            query.start,
                            query.end,
                            query.length,
                            source.length,
                            source.cds_length,
                            source_regions,
                            source_cds_block_text,
                            target.transcript_id,
                            target.length,
                            target.cds_length,
                            mapped_bases,
                            unmapped_bases,
                            status,
                            target_regions,
                            target_block_text,
                            target_cds_block_text,
                        ]
                    )

                    # Build a complete source-to-target row for each exact shared
                    # genomic segment. First use target-CDS boundaries, then split
                    # again at source-CDS boundaries so both CDS coordinate systems
                    # are valid for the exact same bases in every output row.
                    for (
                        _,
                        _,
                        g_start,
                        g_end,
                        source_block_index,
                        _,
                        _,
                        _,
                    ) in split_blocks:
                        exact_source_blocks = genomic_interval_to_transcript_blocks(
                            source, [(g_start, g_end, source_block_index)]
                        )
                        exact_source_segments = []
                        for exact_source_block in exact_source_blocks:
                            exact_source_segments.extend(
                                split_mapped_block_by_cds(source, exact_source_block)
                            )

                        for (
                            s_tx_start,
                            s_tx_end,
                            s_g_start,
                            s_g_end,
                            exact_source_block_index,
                            source_region,
                            source_cds_start,
                            source_cds_end,
                        ) in exact_source_segments:
                            exact_target_blocks = genomic_interval_to_transcript_blocks(
                                target,
                                [(s_g_start, s_g_end, exact_source_block_index)],
                            )
                            exact_target_segments = []
                            for exact_target_block in exact_target_blocks:
                                exact_target_segments.extend(
                                    split_mapped_block_by_cds(
                                        target, exact_target_block
                                    )
                                )

                            for (
                                t_start,
                                t_end,
                                tg_start,
                                tg_end,
                                final_source_block_index,
                                target_region,
                                target_cds_start,
                                target_cds_end,
                            ) in exact_target_segments:
                                blocks_w.writerow(
                                    [
                                        query.name,
                                        source.gene_id,
                                        source.transcript_id,
                                        query.start,
                                        query.end,
                                        s_tx_start,
                                        s_tx_end,
                                        source_region,
                                        (
                                            ""
                                            if source_cds_start is None
                                            else source_cds_start
                                        ),
                                        (
                                            ""
                                            if source_cds_end is None
                                            else source_cds_end
                                        ),
                                        source.length,
                                        source.cds_length,
                                        source.chrom,
                                        tg_start,
                                        tg_end,
                                        source.strand,
                                        final_source_block_index,
                                        target.transcript_id,
                                        t_start,
                                        t_end,
                                        target_region,
                                        (
                                            ""
                                            if target_cds_start is None
                                            else target_cds_start
                                        ),
                                        (
                                            ""
                                            if target_cds_end is None
                                            else target_cds_end
                                        ),
                                        target.length,
                                        target.cds_length,
                                    ]
                                )

                n_ok += 1
            except Exception as exc:  # continue processing other queries
                n_error += 1
                error_w.writerow(
                    [
                        query.line_no,
                        query.name,
                        query.transcript_id,
                        query.start,
                        query.end,
                        str(exc),
                    ]
                )

    print(
        f"Loaded {len(transcripts)} unambiguous transcripts from {args.gtf}",
        file=sys.stderr,
    )
    print(
        f"Excluded {len(ambiguous_transcripts)} ambiguous transcript IDs; "
        f"details: {conflict_path}",
        file=sys.stderr,
    )
    print(f"Processed {n_ok} queries; {n_error} queries had errors", file=sys.stderr)
    print(f"Genome BED: {genome_path}", file=sys.stderr)
    print(f"Source coordinates: {source_path}", file=sys.stderr)
    print(
        "Self-fallback: enabled (source is emitted only when no other "
        "transcript overlaps), unless --include-source was used",
        file=sys.stderr,
    )
    print(f"Complete source-to-target mapping: {blocks_path}", file=sys.stderr)
    print(f"Projection summary: {summary_path}", file=sys.stderr)
    print(f"Errors: {error_path}", file=sys.stderr)
    return 0 if n_error == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
