#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
import argparse
import re
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract unedited CDS windows using coordinates from edited window FASTA."
    )
    parser.add_argument("-w", "--window_fasta", required=True)
    parser.add_argument("-r", "--reference_fasta", required=True)
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args()


def parse_record_id(record_id):
    protein = record_id.split("|")[0]

    # Parse window coordinates from the identifier.
    m = re.search(r"window_(\d+)_(\d+)", record_id)
    if m:
        start_coord = int(m.group(1))
        end_coord = int(m.group(2))
        return protein, start_coord, end_coord

    # Use the minimum and maximum CDS positions when window coordinates are absent.
    nums = list(map(int, re.findall(r"cdspos_(\d+)", record_id)))
    if len(nums) >= 2:
        start_coord = min(nums)
        end_coord = max(nums)
        return protein, start_coord, end_coord
    elif len(nums) == 1:
        # A single coordinate does not define a window.
        raise ValueError(
            f"Only one cdspos found in ID: {record_id}, cannot determine window."
        )

    raise ValueError(f"Cannot parse window coordinates from ID: {record_id}")


def extract_protein_and_coordinates(window_fasta):
    protein_coordinates = []
    for record in SeqIO.parse(window_fasta, "fasta"):
        try:
            protein, start_coord, end_coord = parse_record_id(record.id)
        except ValueError as e:
            print(f"Warning: {e}", file=sys.stderr)
            continue
        protein_coordinates.append((protein, start_coord, end_coord))
    return protein_coordinates


def load_reference_fasta(reference_fasta):
    """
    Extract the transcript ID from the FlyBase CDS FASTA parent attribute.
    """
    fasta_dict = {}
    for record in SeqIO.parse(reference_fasta, "fasta"):
        desc = record.description
        fbtr = None

        for field in desc.split(";"):
            field = field.strip()
            if field.startswith("parent="):
                val = field.replace("parent=", "").strip()
                for v in val.split(","):
                    v = v.strip()
                    if v.startswith("FBtr"):
                        fbtr = v
                        break
                if fbtr:
                    break

        if fbtr is None:
            # Otherwise search the full description for a transcript ID.
            m = re.search(r"(FBtr\d+)", desc)
            if m:
                fbtr = m.group(1)

        if fbtr is None:
            print(f"Warning: skip record without FBtr: {desc[:80]}", file=sys.stderr)
            continue

        if fbtr in fasta_dict:
            raise ValueError(f"Duplicate FBtr in reference FASTA: {fbtr}")

        fasta_dict[fbtr] = record

    return fasta_dict


def extract_sequences(reference_fasta, protein_coordinates):
    ref_dict = load_reference_fasta(reference_fasta)

    extracted = []
    seen = set()
    missing = set()
    out_of_range = []

    for protein, start_coord, end_coord in protein_coordinates:
        if protein not in ref_dict:
            missing.add(protein)
            continue

        seq = ref_dict[protein].seq
        L = len(seq)

        if start_coord < 1 or end_coord > L:
            out_of_range.append((protein, start_coord, end_coord, L))
            continue

        window_seq = seq[start_coord - 1 : end_coord]

        name = f"{protein}|unedited_window_{start_coord}_{end_coord}"
        if name not in seen:
            extracted.append((name, window_seq))
            seen.add(name)

    return extracted, missing, out_of_range


def main():
    args = parse_args()

    coords = extract_protein_and_coordinates(args.window_fasta)
    extracted, missing, out_of_range = extract_sequences(args.reference_fasta, coords)

    with open(args.output, "w") as f:
        for name, seq in extracted:
            f.write(f">{name}\n{seq}\n")

    print("Done.")
    print(f"Input window records: {len(coords)}")
    print(f"Output unedited records: {len(extracted)}")
    print(f"Missing FBtr: {len(missing)}")
    print(f"Out-of-range windows: {len(out_of_range)}")

    if missing:
        print("First 20 missing FBtr:")
        for x in sorted(missing)[:20]:
            print(x)

    if out_of_range:
        print("First 20 out-of-range windows:")
        for p, s, e, L in out_of_range[:20]:
            print(f"{p} window {s}-{e} exceeds CDS length {L}")


if __name__ == "__main__":
    main()
