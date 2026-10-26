#!/usr/bin/env python3
"""Select one representative transcript per gene using the gmRC counter rule.

Tie-breaking is identical to transcriptome_bam_codon_counter_longestCDS_NHweighted.py:
largest CDS length, then largest transcript length, then lexicographically
smallest transcript ID.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--txinfo", type=Path, required=True)
    parser.add_argument("--ids-out", type=Path, required=True)
    parser.add_argument("--table-out", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    best: dict[str, tuple[int, int, str]] = {}
    with args.txinfo.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            try:
                gene, transcript = fields[0], fields[1]
                transcript_length, cds_length = int(fields[4]), int(fields[5])
            except ValueError:
                continue
            if cds_length <= 0:
                continue
            candidate = (cds_length, transcript_length, transcript)
            previous = best.get(gene)
            if previous is None:
                best[gene] = candidate
            elif (candidate[0], candidate[1]) > (previous[0], previous[1]) or (
                candidate[0] == previous[0]
                and candidate[1] == previous[1]
                and candidate[2] < previous[2]
            ):
                best[gene] = candidate

    args.ids_out.parent.mkdir(parents=True, exist_ok=True)
    args.table_out.parent.mkdir(parents=True, exist_ok=True)
    with args.ids_out.open("w") as handle:
        for transcript in sorted(value[2] for value in best.values()):
            handle.write(f"{transcript}\n")
    with args.table_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["gene_id", "selected_transcript", "transcript_length", "cds_length"]
        )
        for gene in sorted(best):
            cds_length, transcript_length, transcript = best[gene]
            writer.writerow([gene, transcript, transcript_length, cds_length])


if __name__ == "__main__":
    main()
