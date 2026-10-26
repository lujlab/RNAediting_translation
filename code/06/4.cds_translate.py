from Bio import SeqIO
import sys


def translate_cds_sequences(cds_fasta_file, output_fasta_file):
    """
    Translate CDS sequences into Protein sequences
    """
    cds_sequences = list(SeqIO.parse(cds_fasta_file, "fasta"))

    protein_sequences = []

    for cds_sequence in cds_sequences:
        protein_sequence = cds_sequence.seq.translate()

        protein_record = SeqIO.SeqRecord(
            protein_sequence, id=cds_sequence.id, description=""
        )
        protein_sequences.append(protein_record)

    with open(output_fasta_file, "w+") as output_file:
        SeqIO.write(protein_sequences, output_file, "fasta")


translate_cds_sequences(sys.argv[1], sys.argv[2])
