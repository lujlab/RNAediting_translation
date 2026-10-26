#!/usr/bin/env python3
"""Validate the original FASTA editor/codon extractor against real input sequences."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
import re

def fasta_records(path):
    head = None
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(seq).upper()
                head = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
    if head is not None:
        yield head, "".join(seq).upper()

def txid(head):
    m = re.search(r"FBtr\d+", head)
    return m.group() if m else None


def main():
    ap = argparse.ArgumentParser()
    for name in ("reference", "edited", "sites", "codons"):
        ap.add_argument("--" + name, required=True)
    a = ap.parse_args()
    counts = Counter()
    reference = list(fasta_records(a.reference))
    edited = list(fasta_records(a.edited))
    edits = defaultdict(set)
    with open(a.sites) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            counts["input_rows"] += 1
            if row["target_region"] != "CDS":
                counts["non_CDS_rows"] += 1
                continue
            edits[row["new_transcript"]].add(int(row["new_cds_pos_1based"]))
    details = []
    seen = set()
    codon_records = 0
    if len(reference) != len(edited):
        raise ValueError("Reference/edited FASTA record counts differ")
    with open(a.codons) as codons:
        for (rh, seq), (eh, outseq) in zip(reference, edited):
            tx = txid(rh)
            counts["fasta_records"] += 1
            if rh != eh:
                raise ValueError("FASTA header/order changed: " + rh)
            if tx in seen:
                raise ValueError("Duplicate transcript in FASTA: " + str(tx))
            seen.add(tx)
            expected = list(seq)
            for p in sorted(edits.get(tx, ())):
                if not 1 <= p <= len(seq):
                    status = "outside_fasta"
                elif seq[p - 1] != "A":
                    status = "reference_not_A"
                else:
                    expected[p - 1] = "G"
                    status = "applied_A_to_G"
                counts[status] += 1
                details.append(
                    (tx, p, seq[p - 1] if 1 <= p <= len(seq) else "", status)
                )
            if outseq != "".join(expected):
                raise ValueError("Unexpected sequence change: " + str(tx))
            if len(outseq) % 3:
                counts["partial_terminal_codon_transcripts"] += 1
            for start in range(0, len(outseq), 3):
                actual = codons.readline().rstrip("\n").split("\t")
                expected_codon = [f"{tx}_{start//3+1}", outseq[start : start + 3]]
                if actual != expected_codon:
                    raise ValueError(
                        f"Codon mismatch: expected {expected_codon}, observed {actual}"
                    )
                codon_records += 1
        if codons.readline():
            raise ValueError("Extra records in codon table")
    for tx in sorted(set(edits) - seen):
        for p in sorted(edits[tx]):
            details.append((tx, p, "", "missing_transcript"))
            counts["missing_transcript"] += 1
    counts["codon_records_verified"] = codon_records
    warnings = sum(
        counts[k] for k in ("outside_fasta", "reference_not_A", "missing_transcript")
    )
    result = dict(
        status="WARN" if warnings else "PASS",
        counts=dict(counts),
        reference=a.reference,
        edited=a.edited,
        sites=a.sites,
    )
    # Stop the wrapper if any intended coding edit was silently skipped.
    return int(bool(warnings))


if __name__ == "__main__":
    raise SystemExit(main())
