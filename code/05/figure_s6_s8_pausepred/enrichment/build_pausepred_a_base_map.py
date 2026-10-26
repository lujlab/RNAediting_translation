#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Build a transcript-coordinate map of reference A bases in CDS regions.

PausePred uses full-transcript, 1-based coordinates.  The CDS interval table
provides the CDS coordinates and the reference transcriptome FASTA provides
the sequence.  The output is a nucleotide-position map: every reference A
within a complete CDS is one eligible unedited-background position.  No
codon-level aggregation is performed.
"""

import argparse
import csv
import gzip
import re
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a full-transcript-coordinate map of reference A bases in CDS."
    )
    parser.add_argument("--cds-intervals", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def read_cds_intervals(path):
    raw = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        comment="#",
        dtype=str,
        engine="python",
    )
    if raw.shape[1] < 3:
        raise RuntimeError(
            "CDS interval file must contain transcript_id, start, and end: " + str(path)
        )

    raw = raw.iloc[:, :4].copy()
    while raw.shape[1] < 4:
        raw[raw.shape[1]] = "CDS"
    raw.columns = [
        "transcript_id",
        "cds_start_1based",
        "cds_end_1based",
        "feature",
    ]
    raw["transcript_id"] = raw["transcript_id"].astype(str).str.strip()
    raw["feature"] = raw["feature"].astype(str).str.strip().str.upper()
    raw["cds_start_1based"] = pd.to_numeric(raw["cds_start_1based"], errors="coerce")
    raw["cds_end_1based"] = pd.to_numeric(raw["cds_end_1based"], errors="coerce")
    cds = raw[
        raw["transcript_id"].str.startswith("FBtr", na=False)
        & raw["feature"].eq("CDS")
        & raw["cds_start_1based"].notna()
        & raw["cds_end_1based"].notna()
    ].copy()
    cds["cds_start_1based"] = cds["cds_start_1based"].astype(int)
    cds["cds_end_1based"] = cds["cds_end_1based"].astype(int)
    cds = cds[
        (cds["cds_start_1based"] >= 1)
        & (cds["cds_end_1based"] >= cds["cds_start_1based"])
    ].copy()

    conflicts = cds.groupby("transcript_id").agg(
        n_start=("cds_start_1based", "nunique"),
        n_end=("cds_end_1based", "nunique"),
    )
    conflicts = conflicts[(conflicts["n_start"] > 1) | (conflicts["n_end"] > 1)]
    if not conflicts.empty:
        raise RuntimeError(
            "Conflicting CDS intervals found; first examples:\n"
            + conflicts.head(10).to_string()
        )
    return cds.drop_duplicates("transcript_id", keep="first").set_index("transcript_id")


def fasta_records(path):
    name = None
    chunks = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper().replace("U", "T")
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper().replace("U", "T")


def main():
    args = parse_args()
    cds = read_cds_intervals(args.cds_intervals)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_transcripts_seen = 0
    n_transcripts_written = 0
    n_a_bases_written = 0
    seen_transcripts = set()

    with gzip.open(out_path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "transcript_id",
                "cds_start_1based",
                "cds_end_1based",
                "cds_pos_1based",
                "transcript_pos_1based",
                "reference_base",
            ]
        )

        for transcript_id, sequence in fasta_records(args.reference_fasta):
            if transcript_id not in cds.index:
                continue
            n_transcripts_seen += 1
            seen_transcripts.add(transcript_id)
            row = cds.loc[transcript_id]
            start = int(row["cds_start_1based"])
            end = int(row["cds_end_1based"])
            if end > len(sequence):
                raise RuntimeError(
                    f"CDS end {end} exceeds FASTA length {len(sequence)} for {transcript_id}"
                )

            cds_sequence = sequence[start - 1 : end]
            wrote_this_transcript = False
            for offset, base in enumerate(cds_sequence):
                if base != "A":
                    continue
                cds_pos = offset + 1
                transcript_pos = start + offset
                writer.writerow(
                    [
                        transcript_id,
                        start,
                        end,
                        cds_pos,
                        transcript_pos,
                        base,
                    ]
                )
                n_a_bases_written += 1
                wrote_this_transcript = True
            if wrote_this_transcript:
                n_transcripts_written += 1

    missing = sorted(set(cds.index) - seen_transcripts)
    print(f"[INFO] CDS transcripts: {len(cds):,}")
    print(f"[INFO] FASTA transcripts matched: {n_transcripts_seen:,}")
    print(f"[INFO] Transcripts with CDS A bases: {n_transcripts_written:,}")
    print(f"[INFO] CDS A bases written: {n_a_bases_written:,}")
    if missing:
        print(
            f"[WARN] CDS transcripts absent from FASTA: {len(missing):,}; "
            f"first examples: {', '.join(missing[:5])}"
        )
    print(f"[INFO] Output: {out_path}")


if __name__ == "__main__":
    main()
