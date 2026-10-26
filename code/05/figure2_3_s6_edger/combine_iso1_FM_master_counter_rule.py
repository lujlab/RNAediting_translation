#!/usr/bin/env python3
"""Pool ISO1 female/male flank tables after rebuilding their codon annotation.

The generic master-table builder adds ``site_id`` before replicate columns,
whereas raw per-library flank tables do not.  This variant adds that stable
identifier while combining the two ISO1 libraries into the canonical rep1/2
layout required by the rigorous-v2 edgeR runner.
"""

from __future__ import annotations

import argparse
import csv
from itertools import zip_longest
from pathlib import Path

RAW_META = [
    "chrom",
    "genomic_codon_first_pos_1based",
    "transcript",
    "cds_block_index",
    "strand",
    "cds_codon_first_pos_1based",
    "codon_index_1based",
]
META = RAW_META + ["site_id"]
COUNTS = [
    "rna_count",
    "rpf_count",
    "flank_rna_count",
    "flank_rpf_count",
    "flank_n",
    "library_RNA_sum",
    "library_RPF_sum",
]


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--female", required=True, type=Path)
    parser.add_argument("--male", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    return parser.parse_args()


def run(ns: argparse.Namespace) -> None:
    ns.output.parent.mkdir(parents=True, exist_ok=True)
    ns.manifest.parent.mkdir(parents=True, exist_ok=True)
    expected = RAW_META + COUNTS
    with ns.female.open(newline="") as fh, ns.male.open(
        newline=""
    ) as mh, ns.output.open("w", newline="") as oh:
        fr, mr = csv.reader(fh, delimiter="\t"), csv.reader(mh, delimiter="\t")
        if next(fr) != expected or next(mr) != expected:
            raise ValueError("Unexpected raw ISO1 flank-table columns")
        out = csv.writer(oh, delimiter="\t", lineterminator="\n")
        out.writerow(
            META + [f"rep1__{x}" for x in COUNTS] + [f"rep2__{x}" for x in COUNTS]
        )
        n = 0
        for female, male in zip_longest(fr, mr):
            if (
                female is None
                or male is None
                or female[: len(RAW_META)] != male[: len(RAW_META)]
            ):
                raise ValueError(f"Female/male metadata mismatch at row {n + 2}")
            meta = female[: len(RAW_META)]
            site_id = f"{meta[2]}_{meta[6]}"
            out.writerow(
                meta + [site_id] + female[len(RAW_META) :] + male[len(RAW_META) :]
            )
            n += 1
    fields = [
        "group_id",
        "species",
        "genotype",
        "sex",
        "flank_size",
        "n_replicates",
        "replicates",
        "sample_ids",
        "head_replicates",
        "rows",
        "output_path",
    ]
    with ns.manifest.open("w", newline="") as h:
        out = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerow(
            {
                "group_id": "dmel_ISO1_FM",
                "species": "dmel",
                "genotype": "ISO1",
                "sex": "F+M",
                "flank_size": 30,
                "n_replicates": 2,
                "replicates": "rep1,rep2",
                "sample_ids": "Dmel_ISO1_female,Dmel_ISO1_male",
                "head_replicates": "rep1,rep2",
                "rows": n,
                "output_path": str(ns.output.resolve()),
            }
        )
    print(f"[write] {ns.output} ({n:,} rows)")


if __name__ == "__main__":
    run(args())
