import sys
import pickle

import click
import pandas as pd
import numpy as np
import pysam
from utils import read_fasta, chunk_iter, rev_comp, CLICK_CS


@click.command(context_settings=CLICK_CS)
@click.argument("path_ref", type=click.STRING)
@click.argument("path_bam", type=click.STRING)
@click.argument("path_model", type=click.STRING)
@click.argument("path_txstart", type=click.STRING)
@click.argument("path_out", type=click.STRING)
@click.option(
    "-i",
    "--ignore_txversion",
    is_flag=True,
    default=False,
    help=r'ignore transcript version in ".\d+" format',
)
@click.option(
    "-l",
    "--rlen_min",
    type=click.INT,
    default=None,
    help="lower bound for mapped read length",
)
@click.option(
    "-u",
    "--rlen_max",
    type=click.INT,
    default=None,
    help="upper bound for mapped read length",
)
@click.option(
    "-c",
    "--chunk_size",
    type=click.INT,
    default=100000,
    help="chunk size for prediction batch",
)
def predict(
    path_ref,
    path_bam,
    path_model,
    path_txstart,
    path_out,
    ignore_txversion=False,
    chunk_size=100000,
    rlen_min=None,
    rlen_max=None,
):
    """
    load pre-trained model and predict P-site offsets with frame
    """

    print("...load ref and model", file=sys.stderr)
    ref = read_fasta(path_ref, ignore_version=ignore_txversion)
    model = pickle.load(open(path_model, "rb"))
    nts = model.nts

    print("...load start_codon info", file=sys.stderr)
    tx_start = pd.read_csv(path_txstart, sep="\t")
    tx_start = tx_start.set_index("tx_name")["start_codon"].to_dict()

    if rlen_min is None:
        rlen_min = model.qwidth_range.min()
    if rlen_max is None:
        rlen_max = model.qwidth_range.max()

    print("...Parse bam by chunk and predict", file=sys.stderr)
    nuc_cols = [f"s{i}" for i in range(2 * nts)]

    with pysam.AlignmentFile(path_bam) as bam:
        output = pysam.AlignmentFile(path_out, "wb", header=bam.header)

        for chunk in chunk_iter(bam, size=chunk_size):
            aligns = list(chunk)
            aligns = [
                i
                for i in aligns
                if (
                    i.reference_name in ref
                    and i.query_alignment_length >= rlen_min
                    and i.query_alignment_length <= rlen_max
                )
            ]
            if len(aligns) == 0:
                continue

            qlen = np.zeros(len(aligns), dtype=int)
            nucs = [None] * len(aligns)

            for i, align in enumerate(aligns):
                seq = ref[align.reference_name]
                qlen[i] = align.query_alignment_length
                if align.is_reverse:
                    if align.reference_end + nts <= len(seq):
                        flank_5p = seq[
                            (align.reference_end - nts) : (align.reference_end + nts)
                        ]
                    else:
                        flank_5p = seq[(align.reference_end - nts) :].ljust(
                            2 * nts, "-"
                        )
                    nucs[i] = rev_comp(flank_5p)
                else:
                    if align.reference_start - nts >= 0:
                        flank_5p = seq[
                            (align.reference_start - nts) : (
                                align.reference_start + nts
                            )
                        ]
                    else:
                        flank_5p = seq[: (align.reference_start + nts)].rjust(
                            2 * nts, "-"
                        )
                    nucs[i] = flank_5p

            nucs = pd.DataFrame([list(i) for i in nucs])
            batch = pd.DataFrame({"qwidth": qlen}, index=range(len(qlen)))
            batch = pd.concat([batch, nucs], axis=1)
            batch.columns = ["qwidth"] + nuc_cols
            for i in nuc_cols:
                batch[i] = pd.Categorical(batch[i], categories=["A", "C", "G", "T"])
            batch = pd.get_dummies(batch, drop_first=True)

            p_sites = model.predict(batch)

            for p, align in zip(p_sites, aligns):
                tx = align.reference_name
                tstart = align.reference_start + 1  # 1-based

                if tx not in tx_start:
                    output.write(align)
                    continue

                start_codon = tx_start[tx]
                frame = (tstart + p - start_codon) % 3

                align.set_tag("PS", p, "i")
                align.set_tag("XF", frame, "i")
                output.write(align)

        output.close()


if __name__ == "__main__":
    predict()
