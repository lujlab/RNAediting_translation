#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq


def normalize_fbtr(tx_id):
    return re.sub(r"\.\d+$", "", str(tx_id).strip())


def extract_fbtr(record):
    text = "{} {}".format(record.id, record.description)
    match = re.search(r"(FBtr\d+(?:\.\d+)?)", text)
    return normalize_fbtr(match.group(1)) if match else None


def is_missing(value):
    text = str(value).strip()
    return text == "" or text.upper() in {"NA", "NAN", "NONE", "NULL", "."}


def read_edit_sites(
    edit_file,
    coord_column="new_cds_pos_1based",
    transcript_column="new_transcript",
    region_column="target_region",
):
    transcript_edits = defaultdict(set)

    n_input = 0
    n_valid = 0
    n_nonpositive = 0
    n_invalid = 0
    n_non_cds = 0

    with open(edit_file, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(
                "Cannot read a header from edit file: {}".format(edit_file)
            )

        reader.fieldnames = [
            x.strip() if x is not None else x for x in reader.fieldnames
        ]

        required_columns = {transcript_column, coord_column}
        if coord_column == "new_cds_pos_1based":
            required_columns.add(region_column)

        missing_columns = required_columns - set(reader.fieldnames)

        if missing_columns:
            raise ValueError(
                "Missing required columns: "
                + ", ".join(sorted(missing_columns))
                + "\nColumns found: "
                + ", ".join(reader.fieldnames)
            )

        for row in reader:
            n_input += 1

            tx_text = row.get(transcript_column, "")
            coordinate_text = str(row.get(coord_column, "")).strip()

            if is_missing(tx_text):
                n_invalid += 1
                continue

            tx_id = normalize_fbtr(tx_text)

            if coord_column == "new_cds_pos_1based":
                target_region = str(row.get(region_column, "")).strip()
                if target_region != "CDS":
                    n_non_cds += 1
                    continue

            if is_missing(coordinate_text) or not re.fullmatch(
                r"-?\d+", coordinate_text
            ):
                n_invalid += 1
                continue

            coordinate = int(coordinate_text)

            if coordinate <= 0:
                n_nonpositive += 1
                continue

            transcript_edits[tx_id].add(coordinate)
            n_valid += 1

    print("[INFO] Input editing rows: {}".format(n_input), file=sys.stderr)
    print("[INFO] Valid positive coordinates: {}".format(n_valid), file=sys.stderr)
    print("[INFO] Non-CDS rows skipped: {}".format(n_non_cds), file=sys.stderr)
    print(
        "[INFO] Non-positive coordinates skipped: {}".format(n_nonpositive),
        file=sys.stderr,
    )
    print(
        "[INFO] Invalid or missing rows skipped: {}".format(n_invalid),
        file=sys.stderr,
    )
    print(
        "[INFO] Transcripts with edits: {}".format(len(transcript_edits)),
        file=sys.stderr,
    )

    return transcript_edits


def modify_fasta(
    fasta_file,
    edit_file,
    output_file,
    coord_column="new_cds_pos_1based",
    force=False,
):
    transcript_edits = read_edit_sites(
        edit_file,
        coord_column=coord_column,
    )

    n_fasta_records = 0
    n_matched_transcripts = 0
    n_edits_applied = 0
    n_out_of_range = 0
    n_non_a = 0
    mismatch_examples = []

    with open(output_file, "w", encoding="utf-8") as output_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            n_fasta_records += 1
            tx_id = extract_fbtr(record)

            if tx_id is None or tx_id not in transcript_edits:
                SeqIO.write(record, output_handle, "fasta")
                continue

            n_matched_transcripts += 1
            sequence = list(str(record.seq).upper())

            for coordinate in sorted(transcript_edits[tx_id]):
                index = coordinate - 1

                if index < 0 or index >= len(sequence):
                    n_out_of_range += 1
                    if len(mismatch_examples) < 20:
                        mismatch_examples.append(
                            "{}: coordinate {} outside sequence length {}".format(
                                tx_id, coordinate, len(sequence)
                            )
                        )
                    continue

                original_base = sequence[index]

                if original_base != "A":
                    n_non_a += 1
                    if len(mismatch_examples) < 20:
                        mismatch_examples.append(
                            "{}: coordinate {}, expected A but found {}".format(
                                tx_id, coordinate, original_base
                            )
                        )
                    if not force:
                        continue

                sequence[index] = "G"
                n_edits_applied += 1

            record.seq = Seq("".join(sequence))
            SeqIO.write(record, output_handle, "fasta")

    print("[INFO] FASTA records: {}".format(n_fasta_records), file=sys.stderr)
    print(
        "[INFO] FASTA transcripts matched to edit table: {}".format(
            n_matched_transcripts
        ),
        file=sys.stderr,
    )
    print("[INFO] Edits applied: {}".format(n_edits_applied), file=sys.stderr)
    print(
        "[INFO] Coordinates outside sequence: {}".format(n_out_of_range),
        file=sys.stderr,
    )
    print(
        "[INFO] Edited positions not containing A: {}".format(n_non_a),
        file=sys.stderr,
    )
    print("[INFO] Output written to: {}".format(output_file), file=sys.stderr)

    if mismatch_examples:
        print("\n[WARN] First mismatch examples:", file=sys.stderr)
        for message in mismatch_examples:
            print("  {}".format(message), file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Apply A-to-G RNA-editing substitutions to transcript or CDS "
            "FASTA sequences using the new mapped editing-site table."
        )
    )

    parser.add_argument("fasta_file", help="Input FASTA file.")
    parser.add_argument(
        "edit_file",
        help=(
            "Headered editing table containing new_transcript, target_region, "
            "new_cds_pos_1based and/or new_tx_pos_1based."
        ),
    )
    parser.add_argument("output_file", help="Output edited FASTA file.")
    parser.add_argument(
        "--coord-column",
        default="new_cds_pos_1based",
        choices=["new_cds_pos_1based", "new_tx_pos_1based"],
        help=(
            "Use new_cds_pos_1based for CDS FASTA or "
            "new_tx_pos_1based for full-transcript FASTA."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Replace the selected nucleotide with G even when the original "
            "base is not A."
        ),
    )

    args = parser.parse_args()

    modify_fasta(
        fasta_file=args.fasta_file,
        edit_file=args.edit_file,
        output_file=args.output_file,
        coord_column=args.coord_column,
        force=args.force,
    )


if __name__ == "__main__":
    main()
