# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
import itertools
import argparse
import csv
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate edited CDS/protein FASTA from RNA editing sites."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input editing site file. Supports header format with FBtr/cds_pos or old 6-column format.",
    )

    parser.add_argument("-f", "--fasta", required=True, help="Input CDS FASTA file.")

    parser.add_argument(
        "-o", "--output_cds", required=True, help="Output edited CDS FASTA file."
    )

    parser.add_argument(
        "-p",
        "--output_protein",
        default=None,
        help="Output edited protein FASTA file. If not set, protein FASTA will not be written.",
    )

    parser.add_argument(
        "--fbtr_col",
        type=int,
        default=4,
        help="0-based FBtr column index for no-header input. Default: 4.",
    )

    parser.add_argument(
        "--cdspos_col",
        type=int,
        default=5,
        help="0-based cds_pos column index for no-header input. Default: 5.",
    )

    parser.add_argument(
        "--fbtr_name",
        default="FBtr",
        help="FBtr column name for header input. Default: FBtr.",
    )

    parser.add_argument(
        "--cdspos_name",
        default="cds_pos",
        help="cds_pos column name for header input. Default: cds_pos.",
    )

    parser.add_argument(
        "--max_distance",
        type=int,
        default=900,
        help="Maximum CDS distance for clustering sites. Default: 900 nt.",
    )

    parser.add_argument("--edit_base", default="G", help="Edited base. Default: G.")

    return parser.parse_args()


def has_header(file_path, fbtr_name="FBtr", cdspos_name="cds_pos"):
    """
    Check whether the input file has a header.

    If the first line contains both FBtr and cds_pos,
    the file is treated as a header-format file.
    """
    with open(file_path, "r") as f:
        first_line = f.readline().rstrip("\n").split("\t")

    return fbtr_name in first_line and cdspos_name in first_line


def read_site_file(
    file_path, fbtr_col=4, cdspos_col=5, fbtr_name="FBtr", cdspos_name="cds_pos"
):
    """
    Read the editing site file.

    Supported formats:

    1. Header format:
       site_id gene FBgn FBtr refAA altAA trans_pos trans_len cds_pos ...

    2. No-header old format:
       2L 10533212 G FBtr0080094 FBtr0080094 905

       By default:
       column 5 is FBtr
       column 6 is cds_pos

       The column indexes are 0-based in this script.
    """
    rows = []

    if has_header(file_path, fbtr_name, cdspos_name):
        with open(file_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")

            for row in reader:
                fbtr = row[fbtr_name]
                cds_pos = int(row[cdspos_name])

                rows.append({"FBtr": fbtr, "cds_pos": cds_pos, "raw": row})

    else:
        with open(file_path, "r") as f:
            for line in f:
                line = line.rstrip("\n")

                if not line:
                    continue

                fields = line.split("\t")

                fbtr = fields[fbtr_col]
                cds_pos = int(fields[cdspos_col])

                rows.append({"FBtr": fbtr, "cds_pos": cds_pos, "raw": fields})

    return rows


def group_sites_by_fbtr_and_distance(rows, max_distance=900):
    """
    Group sites by FBtr first.

    Within each FBtr, sort sites by cds_pos.
    Adjacent sites with distance <= max_distance are assigned
    to the same cluster.
    """
    rows_by_fbtr = defaultdict(list)

    for row in rows:
        rows_by_fbtr[row["FBtr"]].append(row)

    groups = {}

    for fbtr, fbtr_rows in rows_by_fbtr.items():
        fbtr_rows = sorted(fbtr_rows, key=lambda x: x["cds_pos"])

        cluster_id = 1
        current_cluster = []

        for row in fbtr_rows:
            if not current_cluster:
                current_cluster.append(row)
            else:
                previous_pos = current_cluster[-1]["cds_pos"]
                current_pos = row["cds_pos"]

                if abs(current_pos - previous_pos) <= max_distance:
                    current_cluster.append(row)
                else:
                    group_key = "{}_cluster{}".format(fbtr, cluster_id)
                    groups[group_key] = current_cluster

                    cluster_id += 1
                    current_cluster = [row]

        if current_cluster:
            group_key = "{}_cluster{}".format(fbtr, cluster_id)
            groups[group_key] = current_cluster

    return groups


def load_fasta(fasta_file):
    return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))


def get_sequence_by_fbtr(fasta_dict, fbtr):
    """
    Find the CDS sequence for a given FBtr.

    First try to match the FASTA record ID.
    If not found, search the FASTA description.
    """
    if fbtr in fasta_dict:
        return str(fasta_dict[fbtr].seq)

    for seq_id, seq_record in fasta_dict.items():
        if fbtr in seq_record.description:
            return str(seq_record.seq)

    return None


def edit_cds_sequence(cds_sequence, positions, edit_base="G"):
    """
    Edit the CDS sequence.

    positions are 1-based CDS coordinates.
    Each selected position is replaced by edit_base.
    """
    seq_list = list(cds_sequence)

    for pos in positions:
        pos = int(pos)

        if pos < 1 or pos > len(seq_list):
            raise ValueError(
                "cds_pos {} out of range. CDS length = {}".format(pos, len(seq_list))
            )

        seq_list[pos - 1] = edit_base

    return "".join(seq_list)


def generate_combinations(positions):
    """
    Generate all non-empty editing combinations.
    """
    positions = sorted(set([int(x) for x in positions]))

    combinations = []

    for i in range(1, len(positions) + 1):
        for subset in itertools.combinations(positions, i):
            combinations.append(subset)

    return combinations


def translate_cds_to_protein(cds_sequence):
    """
    Translate CDS sequence to protein sequence.

    A terminal stop codon is removed.
    Internal stop codons are replaced by X.
    """
    protein = str(Seq(cds_sequence).translate(to_stop=False))

    if protein.endswith("*"):
        protein = protein[:-1]

    if "*" in protein:
        protein = protein.replace("*", "X")

    return protein


def main():
    args = parse_args()

    rows = read_site_file(
        file_path=args.input,
        fbtr_col=args.fbtr_col,
        cdspos_col=args.cdspos_col,
        fbtr_name=args.fbtr_name,
        cdspos_name=args.cdspos_name,
    )

    groups = group_sites_by_fbtr_and_distance(rows, max_distance=args.max_distance)

    fasta_dict = load_fasta(args.fasta)

    edited_cds_records = []
    edited_protein_records = []

    missing_fbtr = []
    total_groups = 0
    total_edited_records = 0

    for group_key, group_rows in groups.items():
        total_groups += 1

        fbtr = group_rows[0]["FBtr"]
        cds_sequence = get_sequence_by_fbtr(fasta_dict, fbtr)

        if cds_sequence is None:
            missing_fbtr.append(fbtr)
            continue

        positions = [row["cds_pos"] for row in group_rows]
        combinations = generate_combinations(positions)

        for combination in combinations:
            edited_cds = edit_cds_sequence(
                cds_sequence, combination, edit_base=args.edit_base
            )

            combination_name = "_".join([str(x) for x in combination])

            record_name = "{}|{}|cdspos_{}".format(fbtr, group_key, combination_name)

            edited_cds_records.append((record_name, edited_cds))

            if args.output_protein is not None:
                edited_protein = translate_cds_to_protein(edited_cds)
                edited_protein_records.append((record_name, edited_protein))

            total_edited_records += 1

    with open(args.output_cds, "w") as f_out:
        for record_name, sequence in edited_cds_records:
            f_out.write(">{}\n{}\n".format(record_name, sequence))

    if args.output_protein is not None:
        with open(args.output_protein, "w") as f_out:
            for record_name, sequence in edited_protein_records:
                f_out.write(">{}\n{}\n".format(record_name, sequence))

    print("Done.")
    print("Input site file:", args.input)
    print("Input CDS FASTA:", args.fasta)
    print("Output edited CDS FASTA:", args.output_cds)

    if args.output_protein is not None:
        print("Output edited protein FASTA:", args.output_protein)

    print("Total input sites:", len(rows))
    print("Total groups:", total_groups)
    print("Total edited records:", total_edited_records)

    if missing_fbtr:
        missing_fbtr = sorted(set(missing_fbtr))
        print("Warning: {} FBtr not found in CDS FASTA.".format(len(missing_fbtr)))
        print("First 20 missing FBtr:")
        print("\n".join(missing_fbtr[:20]))


if __name__ == "__main__":
    main()
