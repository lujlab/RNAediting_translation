# -*- coding: utf-8 -*-
from Bio import SeqIO
import argparse


def is_multiple_of_three(sequence):
    return len(sequence) % 3 == 0


def pad_sequence(sequence):
    remainder = len(sequence) % 3
    if remainder != 0:
        padding_length = 3 - remainder
        padded_sequence = sequence + "N" * padding_length
        return padded_sequence
    else:
        return sequence


def filter_and_pad_fasta_sequences(fasta_file, output_file):
    with open(output_file, "w") as out_f:
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if is_multiple_of_three(record.seq):
                    SeqIO.write(record, out_f, "fasta")
                else:
                    print(
                        f"Sequence '{record.id}' length is not a multiple of three. Padding with 'N'."
                    )
                    padded_sequence = pad_sequence(record.seq)
                    record.seq = padded_sequence
                    SeqIO.write(record, out_f, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and pad FASTA sequences.")
    parser.add_argument("fasta_file", type=str, help="Input FASTA file path")
    parser.add_argument("output_file", type=str, help="Output FASTA file path")

    args = parser.parse_args()
    filter_and_pad_fasta_sequences(args.fasta_file, args.output_file)
