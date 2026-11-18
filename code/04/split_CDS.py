# -*- coding: utf-8 -*-
import re
import sys

# Read the CDS input filename.
cds_file = sys.argv[1]

# Read CDS sequences.
with open(cds_file, "r") as file:
    sequences = file.read().split(">")

# Process each sequence.
for seq in sequences:
    seq_lines = seq.strip().split("\n")
    if len(seq_lines) > 1:  # Skip empty entries.
        header = seq_lines[0]  # Read the FASTA header.
        sequence = "".join(seq_lines[1:])  # Join sequence lines.
        # Extract the parent identifier.
        parent_match = re.search(r"parent=(\S+)", header)
        if parent_match:
            parent = parent_match.group(1)
            # Split the sequence into comma-separated codons.
            codons = [sequence[i : i + 3] for i in range(0, len(sequence), 3)]
            codon_str = ",".join(codons)
            print(f"{parent}\t{codon_str}")
