# -*- coding: utf-8 -*-

import argparse
import re
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rename FASTA headers from pipe-style format to old underscore-style format."
    )

    parser.add_argument("-i", "--input", required=True, help="Input FASTA file.")

    parser.add_argument(
        "-o", "--output", required=True, help="Output FASTA file with renamed headers."
    )

    parser.add_argument(
        "--line_width",
        type=int,
        default=60,
        help="FASTA sequence line width. Default: 60.",
    )

    return parser.parse_args()


def wrap_sequence(seq, width=60):
    seq = str(seq)
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def rename_record_id(record_id):
    """
    Convert new FASTA header format to old-style format.

    Edited example:
    FBtr0310279|FBtr0310279_cluster1|cdspos_2053|window_1903_2205
    -> FBtr0310279_1903_2053_2205

    Edited multiple-site example:
    FBtr0310279|FBtr0310279_cluster1|cdspos_2053_2056_2062|window_1903_2205
    -> FBtr0310279_1903_2053_2056_2062_2205

    Unedited example:
    FBtr0310279|unedited_window_1903_2205
    -> FBtr0310279_1903_2205

    The edited naming convention is:
    transcript_window_start_edit_site(s)_window_end
    """

    # Transcript ID is before the first pipe
    protein = record_id.split("|")[0]

    # Extract window start and end
    m_window = re.search(r"(?:window|unedited_window)_([0-9]+)_([0-9]+)", record_id)

    if m_window is not None:
        start = m_window.group(1)
        end = m_window.group(2)
    else:
        start = None
        end = None

    # Extract edited CDS positions
    m_cdspos = re.search(r"cdspos_([0-9_]+)", record_id)

    # Edited sequence
    if m_cdspos is not None:
        cdspos = m_cdspos.group(1)

        if start is not None and end is not None:
            return "{}_{}_{}_{}".format(protein, start, cdspos, end)
        else:
            return "{}_{}".format(protein, cdspos)

    # Unedited sequence
    # Do not add "unedited" in the new name
    if "unedited" in record_id:
        if start is not None and end is not None:
            return "{}_{}_{}".format(protein, start, end)
        else:
            return protein

    # Already old-style or unknown format
    cleaned_id = record_id.replace("|", "_").replace("+", "_")
    return cleaned_id


def rename_fasta(input_fasta, output_fasta, line_width=60):
    seen_ids = {}

    total_records = 0
    renamed_records = 0
    duplicated_ids = 0

    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            total_records += 1

            old_id = record.id
            new_id = rename_record_id(old_id)

            # Avoid duplicated FASTA IDs
            if new_id in seen_ids:
                seen_ids[new_id] += 1
                new_id_final = "{}_dup{}".format(new_id, seen_ids[new_id])
                duplicated_ids += 1
            else:
                seen_ids[new_id] = 1
                new_id_final = new_id

            out_f.write(
                ">{}\n{}\n".format(
                    new_id_final, wrap_sequence(record.seq, width=line_width)
                )
            )

            if new_id_final != old_id:
                renamed_records += 1

    print("Done.")
    print("Input FASTA:", input_fasta)
    print("Output FASTA:", output_fasta)
    print("Total records:", total_records)
    print("Renamed records:", renamed_records)
    print("Duplicated IDs fixed:", duplicated_ids)


def main():
    args = parse_args()

    rename_fasta(
        input_fasta=args.input, output_fasta=args.output, line_width=args.line_width
    )


if __name__ == "__main__":
    main()
