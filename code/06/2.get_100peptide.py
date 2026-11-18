# -*- coding: utf-8 -*-

from Bio import SeqIO
import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract local CDS windows around edited positions from edited CDS FASTA."
    )

    parser.add_argument("-i", "--input", required=True, help="Input edited CDS FASTA.")

    parser.add_argument(
        "-o", "--output", required=True, help="Output window CDS FASTA."
    )

    parser.add_argument(
        "--flank",
        type=int,
        default=150,
        help="Flanking nucleotide length on each side. Default: 150 nt.",
    )

    parser.add_argument(
        "--max_len",
        type=int,
        default=600,
        help="Maximum extracted sequence length. Default: 600 nt.",
    )

    return parser.parse_args()


def parse_edit_positions(record_id):
    """
    Parse editing positions from FASTA header.

    Supported new format:
    FBtr0310279|FBtr0310279_cluster1|cdspos_2053
    FBtr0310279|FBtr0310279_cluster1|cdspos_2053_2056_2062

    Supported old format:
    FBtr0310279_2053
    FBtr0310279_2053_2056_2062
    """

    # New format: extract numbers after "cdspos_"
    if "cdspos_" in record_id:
        m = re.search(r"cdspos_([0-9_]+)", record_id)
        if m is None:
            raise ValueError("Cannot parse cdspos from record ID: {}".format(record_id))

        pos_string = m.group(1)
        positions = [int(x) for x in pos_string.split("_") if x != ""]
        return positions

    # Old format: extract numeric parts after underscores
    parts = record_id.split("_")
    positions = []

    for x in parts[1:]:
        if x.isdigit():
            positions.append(int(x))

    if len(positions) == 0:
        raise ValueError(
            "Cannot parse editing positions from record ID: {}".format(record_id)
        )

    return positions


def get_codon_aligned_window(min_pos, max_pos, seq_len, flank=150):
    """
    Get codon-aligned CDS window.

    CDS coordinates are 1-based.
    start_pos should satisfy start_pos % 3 == 1.
    end_pos should satisfy end_pos % 3 == 0.
    """

    start_pos = min_pos - flank

    while start_pos % 3 != 1:
        start_pos -= 1

    if start_pos < 1:
        start_pos = 1

    end_pos = max_pos + flank

    while end_pos % 3 != 0:
        end_pos += 1

    if end_pos > seq_len:
        end_pos = seq_len

        # Trim the end to the nearest complete codon boundary if possible
        while end_pos > start_pos and end_pos % 3 != 0:
            end_pos -= 1

    return start_pos, end_pos


def process_sequences(input_file, output_file, flank=150, max_len=600):
    seen_sequences = set()

    total_records = 0
    kept_records = 0
    skipped_too_long = 0
    skipped_duplicate = 0
    skipped_parse_error = 0

    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(input_file, "fasta"):
            total_records += 1

            try:
                edit_positions = parse_edit_positions(record.id)
            except ValueError:
                skipped_parse_error += 1
                continue

            min_pos = min(edit_positions)
            max_pos = max(edit_positions)

            seq_len = len(record.seq)

            start_pos, end_pos = get_codon_aligned_window(
                min_pos=min_pos, max_pos=max_pos, seq_len=seq_len, flank=flank
            )

            window_len = end_pos - start_pos + 1

            if window_len > max_len:
                skipped_too_long += 1
                continue

            window_seq = str(record.seq[start_pos - 1 : end_pos])

            if window_seq in seen_sequences:
                skipped_duplicate += 1
                continue

            seen_sequences.add(window_seq)

            new_seq_name = "{}|window_{}_{}".format(record.id, start_pos, end_pos)

            out_f.write(">{}\n{}\n".format(new_seq_name, window_seq))
            kept_records += 1

    print("Done.")
    print("Input FASTA:", input_file)
    print("Output FASTA:", output_file)
    print("Total input records:", total_records)
    print("Kept records:", kept_records)
    print("Skipped too long:", skipped_too_long)
    print("Skipped duplicated sequence:", skipped_duplicate)
    print("Skipped parse error:", skipped_parse_error)


def main():
    args = parse_args()

    process_sequences(
        input_file=args.input,
        output_file=args.output,
        flank=args.flank,
        max_len=args.max_len,
    )


if __name__ == "__main__":
    main()
