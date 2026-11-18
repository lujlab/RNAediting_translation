#!/usr/bin/env python3
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

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-o", "--output_cds", required=True)
    parser.add_argument("-p", "--output_protein", default=None)

    parser.add_argument("--fbtr_col", type=int, default=4)
    parser.add_argument("--cdspos_col", type=int, default=5)
    parser.add_argument("--fbtr_name", default="FBtr")
    parser.add_argument("--cdspos_name", default="cds_pos")
    parser.add_argument("--max_distance", type=int, default=900)
    parser.add_argument("--edit_base", default="G")

    return parser.parse_args()


def has_header(file_path, fbtr_name, cdspos_name):
    with open(file_path, "r") as f:
        first_line = f.readline().rstrip("\n").split("\t")
    return fbtr_name in first_line and cdspos_name in first_line


def read_site_file(file_path, fbtr_col, cdspos_col, fbtr_name, cdspos_name):
    rows = []

    if has_header(file_path, fbtr_name, cdspos_name):
        with open(file_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rows.append(
                    {
                        "FBtr": row[fbtr_name],
                        "cds_pos": int(row[cdspos_name]),
                        "raw": row,
                    }
                )
    else:
        with open(file_path, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                rows.append(
                    {
                        "FBtr": fields[fbtr_col],
                        "cds_pos": int(fields[cdspos_col]),
                        "raw": fields,
                    }
                )
    return rows


def group_sites_by_fbtr_and_distance(rows, max_distance):
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
                if abs(row["cds_pos"] - current_cluster[-1]["cds_pos"]) <= max_distance:
                    current_cluster.append(row)
                else:
                    groups[f"{fbtr}_cluster{cluster_id}"] = current_cluster
                    cluster_id += 1
                    current_cluster = [row]

        if current_cluster:
            groups[f"{fbtr}_cluster{cluster_id}"] = current_cluster

    return groups


def load_fasta(fasta_file):
    fasta_dict = {}
    bad_parent = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        desc = record.description

        fbtr = None

        # Extract the transcript ID from the parent attribute.
        for field in desc.split(";"):
            field = field.strip()
            if not field.startswith("parent="):
                continue

            # parent=FBgn0165004,FBtr0193038
            val = field.replace("parent=", "").strip()
            # Parent identifiers are comma separated.
            for v in val.split(","):
                v = v.strip()
                if v.startswith("FBtr"):
                    fbtr = v
                    break
            if fbtr:
                break

        # Otherwise search the complete description for a transcript ID.
        if fbtr is None:
            import re

            m = re.search(r"(FBtr\d+)", desc)
            if m:
                fbtr = m.group(1)
            else:
                bad_parent.append(desc[:120])
                continue

        if fbtr in fasta_dict:
            raise ValueError(f"Duplicate FBtr in FASTA: {fbtr}")

        fasta_dict[fbtr] = record

    if bad_parent:
        print(
            f"Warning: {len(bad_parent)} records cannot find FBtr:",
            file=__import__("sys").stderr,
        )
        for b in bad_parent[:5]:
            print(f"  {b}", file=__import__("sys").stderr)

    return fasta_dict


def edit_cds_sequence(cds_sequence, positions, edit_base):
    seq_list = list(cds_sequence)
    for pos in positions:
        pos = int(pos)
        seq_list[pos - 1] = edit_base
    return "".join(seq_list)


def generate_combinations(positions):
    positions = sorted(set([int(x) for x in positions]))
    combos = []
    for i in range(1, len(positions) + 1):
        combos.extend(itertools.combinations(positions, i))
    return combos


def translate_cds_to_protein(cds_sequence):
    protein = str(Seq(cds_sequence).translate(to_stop=False))
    if protein.endswith("*"):
        protein = protein[:-1]
    protein = protein.replace("*", "X")
    return protein


def main():
    args = parse_args()

    rows = read_site_file(
        args.input, args.fbtr_col, args.cdspos_col, args.fbtr_name, args.cdspos_name
    )

    groups = group_sites_by_fbtr_and_distance(rows, args.max_distance)
    fasta_dict = load_fasta(args.fasta)

    edited_cds_records = []
    edited_protein_records = []

    missing_fbtr = []

    for group_key, group_rows in groups.items():
        fbtr = group_rows[0]["FBtr"]
        if fbtr not in fasta_dict:
            missing_fbtr.append(fbtr)
            continue

        cds_sequence = str(fasta_dict[fbtr].seq)
        positions = [row["cds_pos"] for row in group_rows]
        combinations = generate_combinations(positions)

        for combination in combinations:
            edited_cds = edit_cds_sequence(cds_sequence, combination, args.edit_base)

            # Use the established sequence identifier format.
            combination_name = "_".join([str(x) for x in combination])
            record_name = "{}|{}|cdspos_{}".format(fbtr, group_key, combination_name)

            edited_cds_records.append((record_name, edited_cds))

            if args.output_protein is not None:
                edited_protein = translate_cds_to_protein(edited_cds)
                edited_protein_records.append((record_name, edited_protein))

    # Write CDS sequences.
    with open(args.output_cds, "w") as f:
        for name, seq in edited_cds_records:
            f.write(f">{name}\n{seq}\n")

    # Write protein sequences.
    if args.output_protein is not None:
        with open(args.output_protein, "w") as f:
            for name, seq in edited_protein_records:
                f.write(f">{name}\n{seq}\n")

    print("Done.")
    print(f"Total groups: {len(groups)}")
    print(f"Total edited records: {len(edited_cds_records)}")

    if missing_fbtr:
        print(f"Warning: {len(set(missing_fbtr))} FBtr not found:")
        for f in sorted(set(missing_fbtr))[:10]:
            print(f"  {f}")


if __name__ == "__main__":
    main()
