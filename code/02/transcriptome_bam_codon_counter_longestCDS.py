#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
import re
import multiprocessing as mp
from pathlib import Path
from collections import Counter, defaultdict

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
SIMPLE_CIGAR_RE = re.compile(r"^(?:(\d+)S)?(\d+)([M=X])(?:(\d+)S)?$")

G = {}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "RNA-seq transcriptome BAM CDS-codon counter. "
            "Can restrict to one representative transcript per gene using txinfo, "
            "default selecting the longest CDS transcript per gene."
        )
    )

    parser.add_argument("--bam", required=True, help="Input STAR transcriptome BAM")

    parser.add_argument(
        "--no-auto-bam-resolve",
        action="store_true",
        help=(
            "Do not automatically try alternative BAM names. "
            "By default, if the provided BAM does not exist, the script tries "
            "Aligned.toTranscriptome.out_sorted.bam and other common sorted-BAM names."
        ),
    )

    parser.add_argument(
        "--cds-pos", required=True, help="CDS position file: transcript start end CDS"
    )

    parser.add_argument(
        "--txinfo",
        default=None,
        help=(
            "Optional txinfo file. If provided, select one transcript per gene. "
            "Default columns: gene_id=1, transcript_id=2, transcript_len=5, cds_len=6. "
            "The selected transcript is the longest CDS transcript per gene by default."
        ),
    )

    parser.add_argument(
        "--select-mode",
        choices=["longest_cds", "longest_transcript"],
        default="longest_cds",
        help="Representative transcript selection mode. Default: longest_cds",
    )

    parser.add_argument(
        "--txinfo-gene-col",
        type=int,
        default=1,
        help="1-based gene_id column in txinfo. Default: 1",
    )
    parser.add_argument(
        "--txinfo-transcript-col",
        type=int,
        default=2,
        help="1-based transcript_id column in txinfo. Default: 2",
    )
    parser.add_argument(
        "--txinfo-transcript-len-col",
        type=int,
        default=5,
        help="1-based transcript length column in txinfo. Default: 5",
    )
    parser.add_argument(
        "--txinfo-cds-len-col",
        type=int,
        default=6,
        help="1-based CDS length column in txinfo. Default: 6",
    )

    parser.add_argument(
        "--use-txinfo-bed",
        action="store_true",
        help=(
            "Use selected transcript BED with samtools view -L to reduce scanning. "
            "Requires transcript lengths in txinfo and a BAM that can be queried efficiently. "
            "If this causes problems, omit this option."
        ),
    )

    parser.add_argument(
        "--cds-start-base",
        choices=["0", "1", "auto"],
        default="1",
        help="CDS coordinate base. Use 1 for GTF-derived CDS_position.tsv. Default: 1",
    )

    parser.add_argument("--out-prefix", required=True, help="Output prefix")

    parser.add_argument(
        "--threads", type=int, default=4, help="samtools view threads. Default: 4"
    )

    parser.add_argument(
        "--workers", type=int, default=1, help="Python worker processes. Default: 1"
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=20000,
        help="SAM records per multiprocessing chunk. Default: 20000",
    )

    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ. bcftools mpileup default is 0. Default: 0",
    )

    parser.add_argument(
        "--min-baseq",
        type=int,
        default=13,
        help="Minimum base quality for each codon base. bcftools mpileup default is 13. Default: 13",
    )

    parser.add_argument(
        "--skip-secondary",
        action="store_true",
        help="Skip secondary alignments, FLAG 0x100. Default: keep secondary.",
    )

    parser.add_argument(
        "--skip-supplementary",
        action="store_true",
        help="Skip supplementary alignments, FLAG 0x800. Default: keep supplementary.",
    )

    parser.add_argument(
        "--write-old4",
        action="store_true",
        help="Write old 4-column per-read codon output",
    )

    parser.add_argument(
        "--write-records",
        action="store_true",
        help="Detailed per-read output is not supported in this multiprocessing version",
    )

    parser.add_argument(
        "--old4-only",
        action="store_true",
        help="Only write old4; skip codon_counts and codon_depth aggregation",
    )

    return parser.parse_args()


def safe_int(x, default=None):
    try:
        return int(x)
    except Exception:
        return default


def load_selected_transcripts_from_txinfo(
    path,
    mode="longest_cds",
    gene_col=1,
    transcript_col=2,
    transcript_len_col=5,
    cds_len_col=6,
):
    """
    Select one transcript per gene.

    Default txinfo columns:
      1 gene_id
      2 transcript_id
      5 transcript_len
      6 CDS_len

    For longest_cds:
      primary score = CDS_len
      tie-breaker 1 = transcript_len
      tie-breaker 2 = transcript_id lexicographic

    For longest_transcript:
      primary score = transcript_len
      tie-breaker 1 = CDS_len
      tie-breaker 2 = transcript_id lexicographic
    """

    gene_idx = gene_col - 1
    tx_idx = transcript_col - 1
    txlen_idx = transcript_len_col - 1
    cdslen_idx = cds_len_col - 1

    max_idx = max(gene_idx, tx_idx, txlen_idx, cdslen_idx)

    gene_best = {}
    tx_len_by_tx = {}
    cds_len_by_tx = {}

    n_total = 0
    n_valid = 0
    n_skip_short = 0
    n_skip_bad_len = 0

    with open(path, "r") as f:
        for line_no, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            n_total += 1

            if len(fields) <= max_idx:
                n_skip_short += 1
                continue

            gene_id = fields[gene_idx]
            tx_id = fields[tx_idx]

            tx_len = safe_int(fields[txlen_idx])
            cds_len = safe_int(fields[cdslen_idx])

            if tx_len is None or cds_len is None:
                # allows header line
                n_skip_bad_len += 1
                continue

            if cds_len <= 0:
                n_skip_bad_len += 1
                continue

            tx_len_by_tx[tx_id] = tx_len
            cds_len_by_tx[tx_id] = cds_len

            if mode == "longest_cds":
                key = (cds_len, tx_len, tx_id)
            else:
                key = (tx_len, cds_len, tx_id)

            n_valid += 1

            if gene_id not in gene_best:
                gene_best[gene_id] = (key, tx_id, tx_len, cds_len)
            else:
                old_key = gene_best[gene_id][0]
                # largest numeric scores; for transcript_id, choose lexicographically smaller in tie
                old_score1, old_score2, old_tx = old_key
                new_score1, new_score2, new_tx = key

                replace = False
                if new_score1 > old_score1:
                    replace = True
                elif new_score1 == old_score1 and new_score2 > old_score2:
                    replace = True
                elif (
                    new_score1 == old_score1
                    and new_score2 == old_score2
                    and new_tx < old_tx
                ):
                    replace = True

                if replace:
                    gene_best[gene_id] = (key, tx_id, tx_len, cds_len)

    selected = {}
    selected_tx = set()

    for gene_id, (_, tx_id, tx_len, cds_len) in gene_best.items():
        selected[gene_id] = {
            "transcript_id": tx_id,
            "transcript_len": tx_len,
            "cds_len": cds_len,
        }
        selected_tx.add(tx_id)

    sys.stderr.write(f"[INFO] txinfo records total: {n_total}\n")
    sys.stderr.write(f"[INFO] txinfo valid coding records: {n_valid}\n")
    sys.stderr.write(f"[INFO] txinfo skipped short lines: {n_skip_short}\n")
    sys.stderr.write(f"[INFO] txinfo skipped bad/noncoding lengths: {n_skip_bad_len}\n")
    sys.stderr.write(f"[INFO] selected genes: {len(selected)}\n")
    sys.stderr.write(f"[INFO] selected transcripts: {len(selected_tx)}\n")
    sys.stderr.write(f"[INFO] representative selection mode: {mode}\n")

    if not selected_tx:
        raise RuntimeError("No transcripts selected from txinfo.")

    return selected, selected_tx, tx_len_by_tx, cds_len_by_tx


def write_selected_transcript_files(
    out_prefix, selected, selected_tx, tx_len_by_tx, cds_len_by_tx
):
    ids_out = out_prefix + ".selected_transcripts.txt"
    table_out = out_prefix + ".selected_transcripts.tsv"
    bed_out = out_prefix + ".selected_transcripts.bed"

    with open(ids_out, "w") as f:
        for tx in sorted(selected_tx):
            f.write(tx + "\n")

    with open(table_out, "w") as f:
        f.write("gene_id\tselected_transcript\ttranscript_len\tcds_len\n")
        for gene_id in sorted(selected.keys()):
            rec = selected[gene_id]
            f.write(
                f"{gene_id}\t{rec['transcript_id']}\t"
                f"{rec['transcript_len']}\t{rec['cds_len']}\n"
            )

    with open(bed_out, "w") as f:
        for tx in sorted(selected_tx):
            tx_len = tx_len_by_tx.get(tx)
            if tx_len is None or tx_len <= 0:
                continue
            # BED: 0-based half-open, transcript contig whole length
            f.write(f"{tx}\t0\t{tx_len}\n")

    return ids_out, table_out, bed_out


def load_cds_regions(
    path, cds_start_base="1", keep_transcripts=None, selected_cds_pos_out=None
):
    raw_records = []
    n_total = 0
    n_skip = 0

    selected_rows = []

    with open(path, "r") as f:
        for line_no, line in enumerate(f, 1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 3:
                n_skip += 1
                continue

            if line_no == 1:
                try:
                    int(fields[1])
                    int(fields[2])
                except ValueError:
                    continue

            transcript = fields[0]

            if keep_transcripts is not None and transcript not in keep_transcripts:
                continue

            try:
                start_raw = int(fields[1])
                end_raw = int(fields[2])
            except ValueError:
                n_skip += 1
                continue

            region = fields[3] if len(fields) >= 4 else "CDS"
            n_total += 1

            selected_rows.append((transcript, start_raw, end_raw, region))

            if region != "CDS":
                continue

            raw_records.append((transcript, start_raw, end_raw, line_no))

    if selected_cds_pos_out is not None:
        with open(selected_cds_pos_out, "w") as f:
            for row in selected_rows:
                f.write("\t".join(map(str, row)) + "\n")

    if not raw_records:
        raise RuntimeError(f"No CDS records loaded from: {path}")

    if cds_start_base == "auto":
        base = "0" if any(start == 0 for _, start, _, _ in raw_records) else "1"
    else:
        base = cds_start_base

    cds = {}
    n_invalid = 0

    for transcript, start_raw, end_raw, line_no in raw_records:
        if base == "0":
            start = start_raw + 1
            end = end_raw
        else:
            start = start_raw
            end = end_raw

        if start < 1 or end < start:
            n_invalid += 1
            continue

        cds[transcript] = (start, end)

    sys.stderr.write(
        f"[INFO] CDS position records after transcript filter: {n_total}\n"
    )
    sys.stderr.write(f"[INFO] CDS transcripts loaded: {len(cds)}\n")
    sys.stderr.write(f"[INFO] CDS invalid skipped: {n_invalid}\n")
    sys.stderr.write(f"[INFO] CDS coordinate base used: {base}\n")

    if not cds:
        raise RuntimeError("No valid CDS regions after coordinate conversion.")

    return cds, base, n_invalid


def phred33(c):
    return ord(c) - 33


def cigar_query_len(cigar):
    total = 0
    parsed_len = 0

    for m in CIGAR_RE.finditer(cigar):
        n = int(m.group(1))
        op = m.group(2)
        parsed_len += len(m.group(0))

        if op in ("M", "I", "S", "=", "X"):
            total += n

    if parsed_len != len(cigar):
        return None

    return total


def parse_simple_cigar(cigar):
    m = SIMPLE_CIGAR_RE.match(cigar)

    if not m:
        return None

    left_s = int(m.group(1)) if m.group(1) else 0
    match_len = int(m.group(2))

    return left_s, match_len


def codon_start_ge(cds_pos):
    return ((cds_pos - 1 + 2) // 3) * 3 + 1


def pass_baseq3(qual, q_start, min_baseq):
    if min_baseq <= 0:
        return True

    if qual == "*":
        return False

    if q_start < 0 or q_start + 2 >= len(qual):
        return False

    return (
        phred33(qual[q_start]) >= min_baseq
        and phred33(qual[q_start + 1]) >= min_baseq
        and phred33(qual[q_start + 2]) >= min_baseq
    )


def iter_codons_simple(
    transcript,
    aln_start_tx,
    left_soft_clip,
    match_len,
    seq,
    qual,
    cds_start_tx,
    cds_end_tx,
    min_baseq,
    stats,
):
    tx_match_start = aln_start_tx
    tx_match_end = aln_start_tx + match_len - 1

    ov_start_tx = max(tx_match_start, cds_start_tx)
    ov_end_tx = min(tx_match_end, cds_end_tx)

    if ov_start_tx > ov_end_tx:
        return

    first_cds = ov_start_tx - cds_start_tx + 1
    last_cds = ov_end_tx - cds_start_tx + 1

    first_codon_start_cds = codon_start_ge(first_cds)
    last_allowed_start_cds = last_cds - 2

    if first_codon_start_cds > last_allowed_start_cds:
        return

    for codon_start_cds in range(first_codon_start_cds, last_allowed_start_cds + 1, 3):
        codon_end_cds = codon_start_cds + 2

        codon_start_tx = cds_start_tx + codon_start_cds - 1
        codon_end_tx = codon_start_tx + 2

        q_start = left_soft_clip + (codon_start_tx - tx_match_start)

        if q_start < 0 or q_start + 2 >= len(seq):
            stats["skip_simple_qpos_out_of_range"] += 1
            continue

        stats["candidate_codons_total"] += 1

        if not pass_baseq3(qual, q_start, min_baseq):
            stats["candidate_codons_failed_baseq"] += 1
            continue

        codon = seq[q_start : q_start + 3]

        yield (
            transcript,
            codon_start_cds,
            codon_end_cds,
            codon_start_tx,
            codon_end_tx,
            codon,
        )


def build_txpos_to_base_complex(
    pos_1based, cigar, seq, qual, min_baseq, cds_start_tx, cds_end_tx
):
    if cigar == "*" or seq == "*":
        return None

    ops = CIGAR_RE.findall(cigar)

    if not ops:
        return None

    reconstructed = "".join(n + op for n, op in ops)

    if reconstructed != cigar:
        return None

    expected_len = cigar_query_len(cigar)

    if expected_len is None:
        return None

    if expected_len != len(seq):
        return None

    if min_baseq > 0:
        if qual == "*" or len(qual) != len(seq):
            return None

    tx_pos = pos_1based
    q_pos = 0
    tx_to_base = {}

    aligned_match_bases = 0
    low_baseq_bases = 0

    for n_str, op in ops:
        n = int(n_str)

        if op in ("M", "=", "X"):
            for i in range(n):
                this_tx = tx_pos + i
                q_index = q_pos + i

                if this_tx < cds_start_tx or this_tx > cds_end_tx:
                    continue

                aligned_match_bases += 1

                if min_baseq > 0:
                    bq = phred33(qual[q_index])
                    if bq < min_baseq:
                        low_baseq_bases += 1
                        continue

                tx_to_base[this_tx] = seq[q_index]

            tx_pos += n
            q_pos += n

        elif op == "I":
            q_pos += n

        elif op == "S":
            q_pos += n

        elif op in ("D", "N"):
            tx_pos += n

        elif op in ("H", "P"):
            pass

        else:
            return None

    return tx_to_base, aligned_match_bases, low_baseq_bases


def iter_codons_complex(transcript, tx_to_base, cds_start_tx, cds_end_tx, stats):
    if not tx_to_base:
        return

    min_tx = min(tx_to_base.keys())
    max_tx = max(tx_to_base.keys())

    first_cds = max(1, min_tx - cds_start_tx + 1)
    last_cds = min(cds_end_tx - cds_start_tx + 1, max_tx - cds_start_tx + 1)

    first_codon_start_cds = codon_start_ge(first_cds)
    last_allowed_start_cds = last_cds - 2

    if first_codon_start_cds > last_allowed_start_cds:
        return

    for codon_start_cds in range(first_codon_start_cds, last_allowed_start_cds + 1, 3):
        codon_end_cds = codon_start_cds + 2

        codon_start_tx = cds_start_tx + codon_start_cds - 1
        codon_end_tx = codon_start_tx + 2

        p1 = codon_start_tx
        p2 = codon_start_tx + 1
        p3 = codon_start_tx + 2

        stats["candidate_codons_total"] += 1

        if p1 in tx_to_base and p2 in tx_to_base and p3 in tx_to_base:
            codon = tx_to_base[p1] + tx_to_base[p2] + tx_to_base[p3]

            yield (
                transcript,
                codon_start_cds,
                codon_end_cds,
                codon_start_tx,
                codon_end_tx,
                codon,
            )
        else:
            stats["candidate_codons_not_fully_covered"] += 1


def init_worker(cds_regions, min_baseq, old4_only):
    G["cds_regions"] = cds_regions
    G["min_baseq"] = min_baseq
    G["old4_only"] = old4_only


def process_chunk(lines):
    cds_regions = G["cds_regions"]
    min_baseq = G["min_baseq"]
    old4_only = G["old4_only"]

    stats = Counter()
    codon_counts = Counter()
    old4_lines = []

    for line in lines:
        stats["sam_records_after_samtools_filter"] += 1

        line = line.rstrip("\n")
        if not line:
            continue

        fields = line.split("\t")

        if len(fields) < 11:
            stats["skip_malformed_sam"] += 1
            continue

        try:
            transcript = fields[2]
            pos = int(fields[3])
        except ValueError:
            stats["skip_malformed_sam"] += 1
            continue

        cigar = fields[5]
        seq = fields[9]
        qual = fields[10]

        if transcript == "*" or transcript not in cds_regions:
            stats["skip_transcript_without_cds"] += 1
            continue

        if seq == "*" or cigar == "*":
            stats["skip_no_seq_or_cigar"] += 1
            continue

        cds_start_tx, cds_end_tx = cds_regions[transcript]

        emitted = 0
        simple = parse_simple_cigar(cigar)

        if simple is not None:
            stats["fast_simple_cigar_records"] += 1

            left_s, match_len = simple

            if left_s + match_len > len(seq):
                stats["skip_simple_cigar_seq_length_mismatch"] += 1
                continue

            iterator = iter_codons_simple(
                transcript=transcript,
                aln_start_tx=pos,
                left_soft_clip=left_s,
                match_len=match_len,
                seq=seq,
                qual=qual,
                cds_start_tx=cds_start_tx,
                cds_end_tx=cds_end_tx,
                min_baseq=min_baseq,
                stats=stats,
            )

        else:
            stats["slow_complex_cigar_records"] += 1

            mapped = build_txpos_to_base_complex(
                pos_1based=pos,
                cigar=cigar,
                seq=seq,
                qual=qual,
                min_baseq=min_baseq,
                cds_start_tx=cds_start_tx,
                cds_end_tx=cds_end_tx,
            )

            if mapped is None:
                stats["skip_invalid_cigar_seq_or_qual"] += 1
                continue

            tx_to_base, aligned_match_bases, low_baseq_bases = mapped

            stats["complex_aligned_match_bases_total"] += aligned_match_bases
            stats["complex_aligned_match_bases_filtered_low_baseq"] += low_baseq_bases

            iterator = iter_codons_complex(
                transcript=transcript,
                tx_to_base=tx_to_base,
                cds_start_tx=cds_start_tx,
                cds_end_tx=cds_end_tx,
                stats=stats,
            )

        for (
            tx,
            codon_start_cds,
            codon_end_cds,
            codon_start_tx,
            codon_end_tx,
            codon,
        ) in iterator:

            emitted += 1

            if not old4_only:
                key = (
                    tx,
                    codon_start_cds,
                    codon_end_cds,
                    codon_start_tx,
                    codon_end_tx,
                    codon,
                )
                codon_counts[key] += 1

            old4_lines.append(f"{tx}\t{codon_start_cds}\t{codon_end_cds}\t{codon}\n")

        if emitted > 0:
            stats["records_with_at_least_one_full_codon"] += 1
            stats["codon_records_emitted"] += emitted
        else:
            stats["records_with_zero_full_codon"] += 1

    return "".join(old4_lines), codon_counts, stats


def line_chunks(stream, chunk_size):
    chunk = []

    for line in stream:
        chunk.append(line)

        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []

    if chunk:
        yield chunk


def resolve_bam_path(bam_path, auto_resolve=True):
    p = Path(bam_path)

    if p.is_file() and p.stat().st_size > 0:
        return str(p), []

    tried = []

    if not auto_resolve:
        return str(p), tried

    s = str(p)
    candidates = []

    if s.endswith("Aligned.toTranscriptome.out.bam"):
        candidates.append(
            s[: -len("Aligned.toTranscriptome.out.bam")]
            + "Aligned.toTranscriptome.out_sorted.bam"
        )

    if s.endswith(".out.bam"):
        candidates.append(s[: -len(".out.bam")] + ".out_sorted.bam")

    if s.endswith(".bam"):
        candidates.append(s[: -len(".bam")] + "_sorted.bam")
        candidates.append(s[: -len(".bam")] + ".sorted.bam")
        candidates.append(s[: -len(".bam")] + ".sort.bam")

    if s.endswith("Aligned.toTranscriptome.out_sorted.bam"):
        candidates.append(
            s.replace(
                "Aligned.toTranscriptome.out_sorted.bam",
                "Aligned.toTranscriptome.out.bam",
            )
        )

    seen = set()
    candidates_unique = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            candidates_unique.append(c)

    for c in candidates_unique:
        tried.append(c)
        cp = Path(c)
        if cp.is_file() and cp.stat().st_size > 0:
            return str(cp), tried

    return str(p), tried


def open_samtools_view(
    bam,
    threads,
    min_mapq,
    skip_secondary,
    skip_supplementary,
    selected_bed=None,
):
    flag_filter = 0x4

    if skip_secondary:
        flag_filter |= 0x100

    if skip_supplementary:
        flag_filter |= 0x800

    cmd = [
        "samtools",
        "view",
        "-@",
        str(threads),
        "-q",
        str(min_mapq),
        "-F",
        str(flag_filter),
    ]

    if selected_bed is not None:
        cmd.extend(["-L", selected_bed])

    cmd.append(bam)

    sys.stderr.write("[INFO] running: " + " ".join(cmd) + "\n")

    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )


def write_count_outputs(out_counts, out_depth, codon_counts):
    with open(out_counts, "w") as f:
        f.write(
            "transcript\tcodon_start_cds\tcodon_end_cds\t"
            "codon_start_tx\tcodon_end_tx\tobserved_codon\tcount\n"
        )

        for key in sorted(codon_counts.keys(), key=lambda x: (x[0], x[1], x[5])):
            transcript, cds_s, cds_e, tx_s, tx_e, codon = key
            count = codon_counts[key]

            f.write(
                f"{transcript}\t{cds_s}\t{cds_e}\t"
                f"{tx_s}\t{tx_e}\t{codon}\t{count}\n"
            )

    depth = defaultdict(int)
    obs = defaultdict(Counter)

    for key, count in codon_counts.items():
        transcript, cds_s, cds_e, tx_s, tx_e, codon = key
        pos_key = (transcript, cds_s, cds_e, tx_s, tx_e)
        depth[pos_key] += count
        obs[pos_key][codon] += count

    with open(out_depth, "w") as f:
        f.write(
            "transcript\tcodon_start_cds\tcodon_end_cds\t"
            "codon_start_tx\tcodon_end_tx\tdepth\tobserved_codons\n"
        )

        for pos_key in sorted(depth.keys(), key=lambda x: (x[0], x[1])):
            transcript, cds_s, cds_e, tx_s, tx_e = pos_key

            obs_str = ";".join(
                f"{codon}:{cnt}" for codon, cnt in obs[pos_key].most_common()
            )

            f.write(
                f"{transcript}\t{cds_s}\t{cds_e}\t"
                f"{tx_s}\t{tx_e}\t{depth[pos_key]}\t{obs_str}\n"
            )


def main():
    args = parse_args()

    if args.min_mapq < 0:
        raise ValueError("--min-mapq must be >= 0")

    if args.min_baseq < 0:
        raise ValueError("--min-baseq must be >= 0")

    if args.workers < 1:
        raise ValueError("--workers must be >= 1")

    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be >= 1")

    if args.write_records:
        raise ValueError(
            "This multiprocessing version currently does not support --write-records. "
            "Use --write-old4 instead."
        )

    selected = None
    selected_tx = None
    tx_len_by_tx = {}
    cds_len_by_tx = {}
    ids_out = ""
    table_out = ""
    bed_out = ""
    selected_cds_pos_out = ""

    if args.txinfo:
        selected, selected_tx, tx_len_by_tx, cds_len_by_tx = (
            load_selected_transcripts_from_txinfo(
                args.txinfo,
                mode=args.select_mode,
                gene_col=args.txinfo_gene_col,
                transcript_col=args.txinfo_transcript_col,
                transcript_len_col=args.txinfo_transcript_len_col,
                cds_len_col=args.txinfo_cds_len_col,
            )
        )

        ids_out, table_out, bed_out = write_selected_transcript_files(
            args.out_prefix,
            selected,
            selected_tx,
            tx_len_by_tx,
            cds_len_by_tx,
        )

        selected_cds_pos_out = args.out_prefix + ".selected_CDS_position.tsv"
    else:
        selected_tx = None

    cds_regions, cds_base_used, n_cds_invalid = load_cds_regions(
        args.cds_pos,
        args.cds_start_base,
        keep_transcripts=selected_tx,
        selected_cds_pos_out=selected_cds_pos_out if args.txinfo else None,
    )

    resolved_bam, bam_candidates_tried = resolve_bam_path(
        args.bam,
        auto_resolve=(not args.no_auto_bam_resolve),
    )

    if resolved_bam != args.bam:
        sys.stderr.write(f"[INFO] input BAM not found or empty: {args.bam}\n")
        sys.stderr.write(f"[INFO] auto-resolved BAM: {resolved_bam}\n")
    elif not Path(resolved_bam).is_file():
        if bam_candidates_tried:
            sys.stderr.write(f"[ERROR] input BAM not found: {args.bam}\n")
            sys.stderr.write("[ERROR] also tried:\n")
            for x in bam_candidates_tried:
                sys.stderr.write(f"  {x}\n")
        raise FileNotFoundError(f"BAM not found: {args.bam}")

    selected_bed_for_samtools = None
    if args.txinfo and args.use_txinfo_bed:
        if Path(bed_out).is_file() and Path(bed_out).stat().st_size > 0:
            selected_bed_for_samtools = bed_out
        else:
            sys.stderr.write(
                "[WARN] selected BED is empty; fallback to streaming BAM and filtering in Python.\n"
            )

    out_counts = args.out_prefix + ".codon_counts.tsv"
    out_depth = args.out_prefix + ".codon_depth.tsv"
    out_log = args.out_prefix + ".run.log"
    out_old4 = args.out_prefix + ".old4.tsv"

    stats_total = Counter()
    codon_counts_total = Counter()

    old4_fh = None

    if args.write_old4:
        old4_fh = open(out_old4, "w", buffering=1024 * 1024)

    proc = open_samtools_view(
        bam=resolved_bam,
        threads=args.threads,
        min_mapq=args.min_mapq,
        skip_secondary=args.skip_secondary,
        skip_supplementary=args.skip_supplementary,
        selected_bed=selected_bed_for_samtools,
    )

    if args.workers == 1:
        init_worker(cds_regions, args.min_baseq, args.old4_only)

        for chunk in line_chunks(proc.stdout, args.chunk_size):
            old4_str, counts, stats = process_chunk(chunk)

            if old4_fh is not None:
                old4_fh.write(old4_str)

            if not args.old4_only:
                codon_counts_total.update(counts)

            stats_total.update(stats)

    else:
        with mp.Pool(
            processes=args.workers,
            initializer=init_worker,
            initargs=(cds_regions, args.min_baseq, args.old4_only),
        ) as pool:

            for old4_str, counts, stats in pool.imap_unordered(
                process_chunk,
                line_chunks(proc.stdout, args.chunk_size),
                chunksize=1,
            ):
                if old4_fh is not None:
                    old4_fh.write(old4_str)

                if not args.old4_only:
                    codon_counts_total.update(counts)

                stats_total.update(stats)

    stderr_data = proc.stderr.read()
    ret = proc.wait()

    if old4_fh is not None:
        old4_fh.close()

    if ret != 0:
        sys.stderr.write(stderr_data)
        raise RuntimeError(f"samtools view failed with exit code {ret}")

    if not args.old4_only:
        write_count_outputs(out_counts, out_depth, codon_counts_total)

    with open(out_log, "w") as f:
        f.write("[INPUT]\n")
        f.write(f"bam_argument\t{args.bam}\n")
        f.write(f"bam_resolved\t{resolved_bam}\n")
        if bam_candidates_tried:
            f.write("bam_candidates_tried\t" + ";".join(bam_candidates_tried) + "\n")
        f.write(f"cds_pos\t{args.cds_pos}\n")
        f.write(f"cds_start_base_argument\t{args.cds_start_base}\n")
        f.write(f"cds_start_base_used\t{cds_base_used}\n")
        f.write(f"loaded_cds_transcripts\t{len(cds_regions)}\n")
        f.write(f"invalid_cds_records_skipped\t{n_cds_invalid}\n")
        f.write(f"out_prefix\t{args.out_prefix}\n")
        f.write(f"samtools_threads\t{args.threads}\n")
        f.write(f"workers\t{args.workers}\n")
        f.write(f"chunk_size\t{args.chunk_size}\n")
        f.write(f"min_mapq\t{args.min_mapq}\n")
        f.write(f"min_baseq\t{args.min_baseq}\n")
        f.write(f"skip_secondary\t{args.skip_secondary}\n")
        f.write(f"skip_supplementary\t{args.skip_supplementary}\n")
        f.write(f"old4_only\t{args.old4_only}\n")

        f.write("\n[TRANSCRIPT_SELECTION]\n")
        f.write(f"txinfo\t{args.txinfo if args.txinfo else 'NA'}\n")
        f.write(f"select_mode\t{args.select_mode if args.txinfo else 'NA'}\n")
        f.write(f"txinfo_gene_col\t{args.txinfo_gene_col}\n")
        f.write(f"txinfo_transcript_col\t{args.txinfo_transcript_col}\n")
        f.write(f"txinfo_transcript_len_col\t{args.txinfo_transcript_len_col}\n")
        f.write(f"txinfo_cds_len_col\t{args.txinfo_cds_len_col}\n")
        f.write(
            f"selected_transcript_count\t{len(selected_tx) if selected_tx is not None else 'NA'}\n"
        )
        f.write(f"selected_transcript_ids\t{ids_out if ids_out else 'NA'}\n")
        f.write(f"selected_transcript_table\t{table_out if table_out else 'NA'}\n")
        f.write(f"selected_transcript_bed\t{bed_out if bed_out else 'NA'}\n")
        f.write(
            f"selected_CDS_position\t{selected_cds_pos_out if selected_cds_pos_out else 'NA'}\n"
        )
        f.write(f"use_txinfo_bed\t{args.use_txinfo_bed}\n")
        f.write(
            f"selected_bed_used_by_samtools\t{selected_bed_for_samtools if selected_bed_for_samtools else 'NA'}\n"
        )

        f.write("\n[STATS]\n")
        for k in sorted(stats_total.keys()):
            f.write(f"{k}\t{stats_total[k]}\n")

        f.write("\n[OUTPUT]\n")
        if not args.old4_only:
            f.write(f"codon_counts\t{out_counts}\n")
            f.write(f"codon_depth\t{out_depth}\n")
        if args.write_old4:
            f.write(f"old4\t{out_old4}\n")

    if not args.old4_only:
        sys.stderr.write(f"[DONE] codon counts: {out_counts}\n")
        sys.stderr.write(f"[DONE] codon depth : {out_depth}\n")

    sys.stderr.write(f"[DONE] log         : {out_log}\n")

    if args.txinfo:
        sys.stderr.write(f"[DONE] selected transcripts: {table_out}\n")
        sys.stderr.write(f"[DONE] selected CDS position: {selected_cds_pos_out}\n")

    if args.write_old4:
        sys.stderr.write(f"[DONE] old4 format : {out_old4}\n")


if __name__ == "__main__":
    main()
