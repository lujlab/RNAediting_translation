# -*- coding: utf-8 -*-

from Bio import SeqIO
import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract unedited CDS windows using coordinates from edited window FASTA."
    )

    parser.add_argument(
        "-w",
        "--window_fasta",
        required=True,
        help="Edited window FASTA containing window coordinates in sequence IDs.",
    )

    parser.add_argument(
        "-r", "--reference_fasta", required=True, help="Reference unedited CDS FASTA."
    )

    parser.add_argument(
        "-o", "--output", required=True, help="Output unedited window FASTA."
    )

    return parser.parse_args()


def parse_record_id(record_id):
    """
    Parse transcript ID and window coordinates from FASTA record ID.

    Supported new format:
    FBtr0310279|FBtr0310279_cluster1|cdspos_2053|window_1903_2202

    Supported old format:
    FBtr0310279_2053_1903_2202
    or
    FBtr0310279_2053_2056_1903_2202
    """

    # New format with "|window_start_end"
    if "window_" in record_id:
        first_part = record_id.split("|")[0]
        m = re.search(r"window_([0-9]+)_([0-9]+)", record_id)

        if m is None:
            raise ValueError(
                "Cannot parse window coordinates from ID: {}".format(record_id)
            )

        start_coord = int(m.group(1))
        end_coord = int(m.group(2))

        return first_part, start_coord, end_coord

    # Old underscore format
    parts = record_id.split("_")

    if len(parts) < 3:
        raise ValueError("Cannot parse old-format ID: {}".format(record_id))

    protein = parts[0]
    start_coord = int(parts[-2])
    end_coord = int(parts[-1])

    return protein, start_coord, end_coord


def extract_protein_and_coordinates(window_fasta):
    """
    Extract transcript ID and window coordinates from edited window FASTA.
    """
    protein_coordinates = []

    for record in SeqIO.parse(window_fasta, "fasta"):
        try:
            protein, start_coord, end_coord = parse_record_id(record.id)
        except ValueError as e:
            print("Warning:", e)
            continue

        protein_coordinates.append((protein, start_coord, end_coord))

    return protein_coordinates


def load_reference_fasta(reference_fasta):
    """
    Load reference CDS FASTA into a dictionary.
    """
    return SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))


def find_reference_sequence(ref_dict, protein):
    """
    Find reference CDS sequence by transcript/protein ID.

    First try exact FASTA ID match.
    If not found, search in FASTA description.
    """
    if protein in ref_dict:
        return ref_dict[protein]

    for seq_id, record in ref_dict.items():
        if protein in record.description:
            return record

    return None


def extract_sequences(reference_fasta, protein_coordinates):
    """
    Extract unedited sequences from reference FASTA.
    """
    ref_dict = load_reference_fasta(reference_fasta)

    extracted_sequences = []
    processed_sequence_names = set()

    missing_proteins = set()
    out_of_range = []

    for protein, start_coord, end_coord in protein_coordinates:
        record = find_reference_sequence(ref_dict, protein)

        if record is None:
            missing_proteins.add(protein)
            continue

        seq_len = len(record.seq)

        if start_coord < 1 or end_coord > seq_len:
            out_of_range.append((protein, start_coord, end_coord, seq_len))
            continue

        sequence = record.seq[start_coord - 1 : end_coord]

        sequence_name = "{}|unedited_window_{}_{}".format(
            protein, start_coord, end_coord
        )

        if sequence_name not in processed_sequence_names:
            extracted_sequences.append((sequence_name, sequence))
            processed_sequence_names.add(sequence_name)

    return extracted_sequences, missing_proteins, out_of_range


def main():
    args = parse_args()

    protein_coordinates = extract_protein_and_coordinates(args.window_fasta)

    extracted_sequences, missing_proteins, out_of_range = extract_sequences(
        reference_fasta=args.reference_fasta, protein_coordinates=protein_coordinates
    )

    with open(args.output, "w") as f_out:
        for sequence_name, sequence in extracted_sequences:
            f_out.write(">{}\n{}\n".format(sequence_name, sequence))

    print("Done.")
    print("Input edited window FASTA:", args.window_fasta)
    print("Reference unedited FASTA:", args.reference_fasta)
    print("Output unedited window FASTA:", args.output)
    print("Input window records:", len(protein_coordinates))
    print("Output unedited records:", len(extracted_sequences))
    print("Missing reference IDs:", len(missing_proteins))
    print("Out-of-range windows:", len(out_of_range))

    if missing_proteins:
        print("First 20 missing reference IDs:")
        for x in sorted(missing_proteins)[:20]:
            print(x)

    if out_of_range:
        print("First 20 out-of-range windows:")
        for item in out_of_range[:20]:
            protein, start_coord, end_coord, seq_len = item
            print(
                "{} window {}-{} exceeds reference length {}".format(
                    protein, start_coord, end_coord, seq_len
                )
            )


if __name__ == "__main__":
    main()
