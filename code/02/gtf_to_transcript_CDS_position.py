#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
from collections import defaultdict, Counter

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate transcript-coordinate CDS positions directly from GTF. "
            "Output is compatible with transcriptome BAM codon counting."
        )
    )

    parser.add_argument(
        "--gtf", required=True, help="Input GTF file, e.g. dmel-all-r6.12.gtf"
    )

    parser.add_argument("--out-prefix", required=True, help="Output prefix")

    parser.add_argument(
        "--mixed-order",
        choices=["gtf", "genomic"],
        default="gtf",
        help=(
            "How to order exons for mixed-strand or multi-contig transcripts. "
            "Default: gtf, using exon line order in GTF. "
            "This is important for trans-spliced transcripts such as mod(mdg4)."
        ),
    )

    parser.add_argument(
        "--with-header",
        action="store_true",
        help="Add header to CDS_position.tsv. Default: no header, matching your old file.",
    )

    return parser.parse_args()


def parse_attrs(attr_text):
    attrs = {}
    for key, value in ATTR_RE.findall(attr_text):
        attrs[key] = value
    return attrs


def safe_int(x):
    try:
        return int(x)
    except Exception:
        return None


def dedup_features(features):
    """
    Remove exact duplicate features while preserving first occurrence.
    Feature key includes seqid, start, end, strand, feature.
    """
    seen = set()
    out = []

    for f in features:
        key = (
            f["seqid"],
            f["start"],
            f["end"],
            f["strand"],
            f["feature"],
        )

        if key in seen:
            continue

        seen.add(key)
        out.append(f)

    return out


def order_exons(exons, mixed_order="gtf"):
    """
    For normal transcripts:
      + strand: genomic ascending
      - strand: genomic descending

    For mixed-strand / multi-contig transcripts:
      default uses GTF line order.
      This is safer for FlyBase trans-spliced annotations.
    """
    exons = dedup_features(exons)

    seqids = {e["seqid"] for e in exons}
    strands = {e["strand"] for e in exons if e["strand"] in ("+", "-")}

    is_mixed = False
    reason = "normal"

    if len(seqids) != 1:
        is_mixed = True
        reason = "multi_seqid"
    elif len(strands) != 1:
        is_mixed = True
        reason = "mixed_or_unknown_strand"

    if is_mixed:
        if mixed_order == "gtf":
            ordered = sorted(exons, key=lambda x: x["line_no"])
        else:
            ordered = sorted(
                exons,
                key=lambda x: (
                    x["seqid"],
                    x["start"],
                    x["end"],
                    x["line_no"],
                ),
            )

        return ordered, True, reason

    strand = list(strands)[0] if strands else "+"

    if strand == "+":
        ordered = sorted(
            exons, key=lambda x: (x["seqid"], x["start"], x["end"], x["line_no"])
        )
    else:
        ordered = sorted(
            exons, key=lambda x: (x["seqid"], -x["end"], -x["start"], x["line_no"])
        )

    return ordered, False, reason


def assign_transcript_coordinates(exons):
    """
    Add tx_start / tx_end to ordered exons.
    Transcript coordinate is 1-based, 5' -> 3' transcript order.
    """
    tx_pos = 1
    out = []

    for e in exons:
        length = e["end"] - e["start"] + 1
        if length <= 0:
            continue

        e2 = dict(e)
        e2["tx_start"] = tx_pos
        e2["tx_end"] = tx_pos + length - 1
        e2["length"] = length

        out.append(e2)
        tx_pos += length

    tx_len = tx_pos - 1

    return out, tx_len


def map_interval_to_exon_tx(cds, exon):
    """
    Map overlap between one genomic CDS interval and one exon to transcript coordinates.

    For exon on +:
      genomic start -> transcript start

    For exon on -:
      genomic end -> transcript start

    For mixed-strand trans-spliced transcripts, each exon carries its own strand.
    """
    if cds["seqid"] != exon["seqid"]:
        return None

    if cds["strand"] in ("+", "-") and exon["strand"] in ("+", "-"):
        if cds["strand"] != exon["strand"]:
            return None

    ov_start = max(cds["start"], exon["start"])
    ov_end = min(cds["end"], exon["end"])

    if ov_start > ov_end:
        return None

    if exon["strand"] == "-":
        tx_start = exon["tx_start"] + (exon["end"] - ov_end)
        tx_end = exon["tx_start"] + (exon["end"] - ov_start)
    else:
        tx_start = exon["tx_start"] + (ov_start - exon["start"])
        tx_end = exon["tx_start"] + (ov_end - exon["start"])

    if tx_start > tx_end:
        tx_start, tx_end = tx_end, tx_start

    return {
        "tx_start": tx_start,
        "tx_end": tx_end,
        "genomic_start": ov_start,
        "genomic_end": ov_end,
        "seqid": cds["seqid"],
        "strand": cds["strand"],
        "cds_line_no": cds["line_no"],
        "exon_line_no": exon["line_no"],
    }


def merge_tx_blocks(blocks):
    """
    Merge overlapping or adjacent transcript-coordinate blocks.
    """
    if not blocks:
        return []

    blocks = sorted(blocks, key=lambda x: (x["tx_start"], x["tx_end"]))

    merged = []
    cur = dict(blocks[0])

    for b in blocks[1:]:
        if b["tx_start"] <= cur["tx_end"] + 1:
            cur["tx_end"] = max(cur["tx_end"], b["tx_end"])
        else:
            merged.append(cur)
            cur = dict(b)

    merged.append(cur)
    return merged


def parse_gtf(gtf_path):
    tx_data = defaultdict(
        lambda: {
            "exon": [],
            "CDS": [],
            "transcript_feature": [],
            "gene_id": "",
            "gene_symbol": "",
            "transcript_symbol": "",
        }
    )

    stats = Counter()

    with open(gtf_path, "r") as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                stats["skip_malformed_gtf"] += 1
                continue

            seqid, source, feature, start_s, end_s, score, strand, frame, attrs_s = (
                fields
            )

            if feature not in ("exon", "CDS", "transcript", "mRNA"):
                continue

            start = safe_int(start_s)
            end = safe_int(end_s)

            if start is None or end is None or end < start:
                stats["skip_invalid_coordinate"] += 1
                continue

            attrs = parse_attrs(attrs_s)
            tid = attrs.get("transcript_id", "")

            if not tid:
                stats["skip_no_transcript_id"] += 1
                continue

            rec = {
                "seqid": seqid,
                "source": source,
                "feature": feature,
                "start": start,
                "end": end,
                "score": score,
                "strand": strand,
                "frame": frame,
                "line_no": line_no,
                "attrs": attrs,
            }

            if feature == "exon":
                tx_data[tid]["exon"].append(rec)
                stats["exon_records"] += 1
            elif feature == "CDS":
                tx_data[tid]["CDS"].append(rec)
                stats["cds_records"] += 1
            else:
                tx_data[tid]["transcript_feature"].append(rec)
                stats["transcript_or_mrna_records"] += 1

            if not tx_data[tid]["gene_id"]:
                tx_data[tid]["gene_id"] = attrs.get("gene_id", "")

            if not tx_data[tid]["gene_symbol"]:
                tx_data[tid]["gene_symbol"] = attrs.get("gene_symbol", "")

            if not tx_data[tid]["transcript_symbol"]:
                tx_data[tid]["transcript_symbol"] = attrs.get("transcript_symbol", "")

    return tx_data, stats


def main():
    args = parse_args()

    out_pos = args.out_prefix + ".CDS_position.tsv"
    out_blocks = args.out_prefix + ".CDS_blocks.tsv"
    out_problem = args.out_prefix + ".problem_transcripts.tsv"
    out_mixed = args.out_prefix + ".mixed_strand_transcripts.tsv"
    out_log = args.out_prefix + ".log"

    tx_data, stats = parse_gtf(args.gtf)

    n_total_tx = 0
    n_ok = 0
    n_no_exon = 0
    n_no_cds = 0
    n_no_cds_mapped = 0
    n_mixed = 0
    n_cds_gap = 0
    n_cds_not_multiple3 = 0

    with open(out_pos, "w") as pos_fh, open(out_blocks, "w") as block_fh, open(
        out_problem, "w"
    ) as prob_fh, open(out_mixed, "w") as mixed_fh:

        if args.with_header:
            pos_fh.write("transcript\tstart\tend\tregion\n")

        block_fh.write(
            "transcript\tblock_index\tcds_tx_start\tcds_tx_end\t"
            "seqid\tgenomic_start\tgenomic_end\tstrand\t"
            "cds_line_no\texon_line_no\tgene_id\tgene_symbol\ttranscript_symbol\n"
        )

        prob_fh.write(
            "transcript\treason\tgene_id\tgene_symbol\ttranscript_symbol\tmessage\n"
        )

        mixed_fh.write(
            "transcript\treason\tgene_id\tgene_symbol\ttranscript_symbol\t"
            "exon_count\tcds_count\n"
        )

        for tid in sorted(tx_data.keys()):
            n_total_tx += 1

            data = tx_data[tid]
            exons = data["exon"]
            cds_list = data["CDS"]

            gene_id = data["gene_id"]
            gene_symbol = data["gene_symbol"]
            transcript_symbol = data["transcript_symbol"]

            if not exons:
                n_no_exon += 1
                prob_fh.write(
                    f"{tid}\tno_exon\t{gene_id}\t{gene_symbol}\t{transcript_symbol}\t"
                    f"transcript has CDS or mRNA records but no exon records\n"
                )
                continue

            if not cds_list:
                n_no_cds += 1
                continue

            ordered_exons, is_mixed, mixed_reason = order_exons(
                exons, mixed_order=args.mixed_order
            )

            if is_mixed:
                n_mixed += 1
                mixed_fh.write(
                    f"{tid}\t{mixed_reason}\t{gene_id}\t{gene_symbol}\t{transcript_symbol}\t"
                    f"{len(exons)}\t{len(cds_list)}\n"
                )

            tx_exons, tx_len = assign_transcript_coordinates(ordered_exons)

            cds_blocks = []

            for cds in cds_list:
                mapped_any = False

                for exon in tx_exons:
                    mapped = map_interval_to_exon_tx(cds, exon)

                    if mapped is not None:
                        mapped_any = True
                        cds_blocks.append(mapped)

                if not mapped_any:
                    prob_fh.write(
                        f"{tid}\tcds_not_mapped_to_exon\t{gene_id}\t{gene_symbol}\t{transcript_symbol}\t"
                        f"CDS line {cds['line_no']} did not overlap any exon with compatible strand\n"
                    )

            if not cds_blocks:
                n_no_cds_mapped += 1
                prob_fh.write(
                    f"{tid}\tno_cds_mapped\t{gene_id}\t{gene_symbol}\t{transcript_symbol}\t"
                    f"no CDS block could be mapped to transcript coordinates\n"
                )
                continue

            merged_blocks = merge_tx_blocks(cds_blocks)

            cds_start_tx = min(b["tx_start"] for b in merged_blocks)
            cds_end_tx = max(b["tx_end"] for b in merged_blocks)

            cds_len = cds_end_tx - cds_start_tx + 1

            # Check if merged CDS blocks are continuous in transcript coordinate
            has_gap = False
            if len(merged_blocks) > 1:
                for i in range(len(merged_blocks) - 1):
                    if (
                        merged_blocks[i + 1]["tx_start"]
                        > merged_blocks[i]["tx_end"] + 1
                    ):
                        has_gap = True
                        break

            if has_gap:
                n_cds_gap += 1
                prob_fh.write(
                    f"{tid}\tcds_internal_gap\t{gene_id}\t{gene_symbol}\t{transcript_symbol}\t"
                    f"CDS blocks are not continuous in transcript coordinates\n"
                )

            if cds_len % 3 != 0:
                n_cds_not_multiple3 += 1
                prob_fh.write(
                    f"{tid}\tcds_length_not_multiple_of_3\t{gene_id}\t{gene_symbol}\t{transcript_symbol}\t"
                    f"CDS transcript-coordinate length is {cds_len}\n"
                )

            # Output block-level coordinates
            for bi, b in enumerate(
                sorted(cds_blocks, key=lambda x: (x["tx_start"], x["tx_end"])), 1
            ):
                block_fh.write(
                    f"{tid}\t{bi}\t{b['tx_start']}\t{b['tx_end']}\t"
                    f"{b['seqid']}\t{b['genomic_start']}\t{b['genomic_end']}\t{b['strand']}\t"
                    f"{b['cds_line_no']}\t{b['exon_line_no']}\t"
                    f"{gene_id}\t{gene_symbol}\t{transcript_symbol}\n"
                )

            # Match your old output style: transcript, start, end, region
            # All transcript coordinates are 1-based inclusive.
            pos_fh.write(f"{tid}\t{cds_start_tx}\t{cds_end_tx}\tCDS\n")

            if cds_start_tx > 1:
                pos_fh.write(f"{tid}\t1\t{cds_start_tx - 1}\tUTR5\n")

            if cds_end_tx < tx_len:
                pos_fh.write(f"{tid}\t{cds_end_tx + 1}\t{tx_len}\tUTR3\n")

            pos_fh.write(f"{tid}\t1\t{tx_len}\tTranscript\n")

            n_ok += 1

    with open(out_log, "w") as log_fh:
        log_fh.write("[INPUT]\n")
        log_fh.write(f"gtf\t{args.gtf}\n")
        log_fh.write(f"out_prefix\t{args.out_prefix}\n")
        log_fh.write(f"mixed_order\t{args.mixed_order}\n")
        log_fh.write("\n[COUNTS]\n")
        log_fh.write(f"total_transcripts_seen\t{n_total_tx}\n")
        log_fh.write(f"ok_transcripts_with_cds_position\t{n_ok}\n")
        log_fh.write(f"no_exon\t{n_no_exon}\n")
        log_fh.write(f"no_cds\t{n_no_cds}\n")
        log_fh.write(f"no_cds_mapped\t{n_no_cds_mapped}\n")
        log_fh.write(f"mixed_or_complex_transcripts\t{n_mixed}\n")
        log_fh.write(f"cds_internal_gap\t{n_cds_gap}\n")
        log_fh.write(f"cds_length_not_multiple_of_3\t{n_cds_not_multiple3}\n")

        log_fh.write("\n[GTF_STATS]\n")
        for k in sorted(stats.keys()):
            log_fh.write(f"{k}\t{stats[k]}\n")

        log_fh.write("\n[OUTPUT]\n")
        log_fh.write(f"CDS_position\t{out_pos}\n")
        log_fh.write(f"CDS_blocks\t{out_blocks}\n")
        log_fh.write(f"problem_transcripts\t{out_problem}\n")
        log_fh.write(f"mixed_strand_transcripts\t{out_mixed}\n")

    sys.stderr.write(f"[DONE] CDS position: {out_pos}\n")
    sys.stderr.write(f"[DONE] CDS blocks:   {out_blocks}\n")
    sys.stderr.write(f"[DONE] Problems:     {out_problem}\n")
    sys.stderr.write(f"[DONE] Mixed-strand: {out_mixed}\n")
    sys.stderr.write(f"[DONE] Log:          {out_log}\n")
    sys.stderr.write(f"[STATS] ok transcripts: {n_ok}\n")
    sys.stderr.write(f"[STATS] mixed/complex transcripts: {n_mixed}\n")
    sys.stderr.write(f"[STATS] CDS length not multiple of 3: {n_cds_not_multiple3}\n")


if __name__ == "__main__":
    main()
