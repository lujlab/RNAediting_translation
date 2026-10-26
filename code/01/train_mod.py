#!/usr/bin/env python

import sys
from pickle import dump
from random import random
from re import sub

import click
import numpy as np
import pandas as pd
import seaborn as sns
import pysam
import matplotlib.pyplot as plt

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.utils.multiclass import unique_labels

try:
    from psite.utils import CLICK_CS, read_fasta, read_txinfo, strip_version
except ModuleNotFoundError:
    from utils import CLICK_CS, read_fasta, read_txinfo, strip_version


class DistArgMax(BaseEstimator, ClassifierMixin):
    """
    Assign offset based on the most frequent distance
    to start codon in the training data.
    """

    def __init__(self):
        self.offsets = dict()

    def fit(self, X, y):
        X, y = check_X_y(X, y)
        self.classes_ = unique_labels(y)

        X_qwidth = X[:, 0]

        for qw in unique_labels(X_qwidth):
            val, cnt = np.unique(y[X_qwidth == qw], return_counts=True)
            self.offsets[qw] = val[cnt.argmax()]

        return self

    def predict(self, X):
        check_is_fitted(self)
        X = check_array(X)

        return np.array([self.offsets[i] for i in X[:, 0]])


def frame_test(test, psite):
    """
    Calculate proportions of P-sites in frame 0/1/2.

    frame = (tstart + psite - start_codon) % 3
    """
    frame = (
        test["tstart"].to_numpy() + np.asarray(psite) - test["start_codon"].to_numpy()
    ) % 3

    out = pd.Series(frame).value_counts(normalize=True).reindex([0, 1, 2], fill_value=0)

    out.index.name = "frame"
    return out


def get_txrep(txinfo, type_rep="longest", path_exp=None, ignore_version=False):
    """
    Get representative transcript per gene.

    longest:
        original behavior, choose longest transcript by tx_len.

    cds_longest:
        choose transcript with longest CDS by cds_len.
        If tied, choose longer tx_len.
        If still tied, choose smaller tx_name alphabetically.

    principal / salmon / kallisto:
        original psite behavior.
    """

    txinfo = txinfo.copy()

    if "gene_biotype" in txinfo.columns:
        txinfo = txinfo.loc[txinfo["gene_biotype"] == "protein_coding"]

    if "transcript_biotype" in txinfo.columns:
        txinfo = txinfo.loc[txinfo["transcript_biotype"] == "protein_coding"]

    if txinfo.empty:
        print("No protein-coding transcript remains.", file=sys.stderr)
        sys.exit(1)

    if type_rep == "longest":
        txinfo["tx_len"] = pd.to_numeric(txinfo["tx_len"], errors="coerce").fillna(0)

        txrep = txinfo.sort_values(
            ["gene_id", "tx_len", "tx_name"], ascending=[True, False, True]
        )

        txrep = txrep.groupby("gene_id").first().reset_index(drop=True)

    elif type_rep == "cds_longest":
        txinfo["cds_len"] = pd.to_numeric(txinfo["cds_len"], errors="coerce").fillna(0)
        txinfo["tx_len"] = pd.to_numeric(txinfo["tx_len"], errors="coerce").fillna(0)

        txinfo = txinfo.loc[txinfo["cds_len"] > 0].copy()

        if txinfo.empty:
            print("No transcript with cds_len > 0 remains.", file=sys.stderr)
            sys.exit(1)

        txrep = txinfo.sort_values(
            ["gene_id", "cds_len", "tx_len", "tx_name"],
            ascending=[True, False, False, True],
        )

        txrep = txrep.groupby("gene_id").first().reset_index(drop=True)

    elif type_rep == "principal":
        try:
            tx_extra = pd.read_table(path_exp)

            tx_extra.rename(
                columns={
                    "Gene stable ID": "gene_id",
                    "Transcript stable ID": "tx_name",
                    "APPRIS annotation": "appris",
                    "Ensembl Canonical": "canonical",
                    "RefSeq match transcript (MANE Select)": "mane_select",
                },
                inplace=True,
            )

            tx_extra = tx_extra[~tx_extra.appris.isna()]

            tx_extra = tx_extra.assign(
                appris=pd.Categorical(
                    tx_extra.appris,
                    categories=[
                        "principal1",
                        "principal2",
                        "principal3",
                        "principal4",
                        "principal5",
                        "alternative1",
                        "alternative2",
                    ],
                )
            )

            tx_extra = tx_extra.sort_values(by="appris").reset_index(drop=True)

            appris_prin = [
                "principal1",
                "principal2",
                "principal3",
                "principal4",
                "principal5",
            ]

            tx_extra = tx_extra[tx_extra.appris.isin(appris_prin)]
            tx_extra = tx_extra[~tx_extra.duplicated(subset=["gene_id"])]

            txrep = txinfo[txinfo.tx_name.isin(tx_extra.tx_name)]

        except Exception as e:
            print("input extra txinfo file is incorrect", file=sys.stderr)
            print(e, file=sys.stderr)
            sys.exit(1)

    elif type_rep == "salmon":
        try:
            tx_quant = pd.read_csv(path_exp, sep="\t")

            if ignore_version:
                tx_quant["Name"] = tx_quant["Name"].map(strip_version)

            tx_quant = tx_quant.merge(
                txinfo[["gene_id", "tx_name"]],
                how="inner",
                left_on="Name",
                right_on="tx_name",
            )

            tx_quant = tx_quant.sort_values(
                by=["gene_id", "TPM"], ascending=[True, False]
            )

            tx_quant = tx_quant.groupby("gene_id").first().reset_index()

            txrep = txinfo.loc[txinfo["tx_name"].isin(tx_quant["Name"])]

        except Exception as e:
            print("input salmon quant results is incorrect", file=sys.stderr)
            print(e, file=sys.stderr)
            sys.exit(1)

    elif type_rep == "kallisto":
        try:
            tx_quant = pd.read_csv(path_exp, sep="\t")

            if ignore_version:
                tx_quant["target_id"] = tx_quant["target_id"].map(strip_version)

            tx_quant = tx_quant.merge(
                txinfo[["gene_id", "tx_name"]],
                how="inner",
                left_on="target_id",
                right_on="tx_name",
            )

            tx_quant = tx_quant.sort_values(
                by=["gene_id", "tpm"], ascending=[True, False]
            )

            tx_quant = tx_quant.groupby("gene_id").first().reset_index()

            txrep = txinfo.loc[txinfo["tx_name"].isin(tx_quant["target_id"])]

        except Exception as e:
            print("input kallisto quant results is incorrect", file=sys.stderr)
            print(e, file=sys.stderr)
            sys.exit(1)

    else:
        print("Incorrect txinfo_rep option!", file=sys.stderr)
        sys.exit(1)

    return txrep


def extract_features(path_bam, ref, nts=3, frac=1.0):
    """
    Extract features from transcriptome BAM.
    """

    alignments = []

    with pysam.AlignmentFile(path_bam) as bam:
        for align in bam:
            if frac < 1 and random() > frac:
                continue

            if align.is_unmapped:
                continue

            if align.reference_name not in ref:
                continue

            if align.is_reverse:
                continue

            seq = ref[align.reference_name]

            if align.reference_start - nts >= 0:
                flank_5p = seq[
                    (align.reference_start - nts) : (align.reference_start + nts)
                ]
            else:
                flank_5p = seq[: (align.reference_start + nts)].rjust(2 * nts, "-")

            out = [
                align.query_name,
                align.reference_name,
                align.reference_start + 1,
                align.query_alignment_length,
            ] + list(flank_5p)

            alignments.append(out)

    cols = ["rname", "tx_name", "tstart", "qwidth"] + [f"s{i}" for i in range(2 * nts)]

    return pd.DataFrame(alignments, columns=cols)


@click.command(context_settings=CLICK_CS)
@click.argument("path_ref", type=click.STRING)
@click.argument("path_bam", type=click.STRING)
@click.argument("output_prefix", type=click.STRING)
@click.argument("path_txinfo", type=click.STRING)
@click.option(
    "-t",
    "--type_rep",
    default="longest",
    type=click.Choice(["longest", "cds_longest", "principal", "kallisto", "salmon"]),
    help="type of representative transcripts",
)
@click.option(
    "-e",
    "--path_exp",
    type=click.STRING,
    default=None,
    help="transcript expression results or APPRIS info",
)
@click.option(
    "-i",
    "--ignore_txversion",
    is_flag=True,
    default=False,
    help='ignore transcript version in ".\\d+" format',
)
@click.option(
    "-n",
    "--nts",
    type=click.INT,
    default=3,
    help="flanking nucleotides around 5 prime end",
)
@click.option(
    "-f",
    "--frac",
    type=click.FLOAT,
    default=1.0,
    help="fraction of alignments for training",
)
@click.option(
    "--offset_min",
    type=click.INT,
    default=10,
    help="lower bound of distance between RPF 5p and start codon",
)
@click.option(
    "--offset_max",
    type=click.INT,
    default=14,
    help="upper bound of distance between RPF 5p and start codon",
)
@click.option("-d", "--max_depth", type=click.INT, default=3, help="max depth of trees")
@click.option(
    "-m",
    "--min_samples_split",
    type=click.INT,
    default=6,
    help="min samples required to split an internal node",
)
@click.option(
    "-k", "--keep", is_flag=True, default=False, help="keep intermediate training data"
)
def train(
    path_ref,
    path_bam,
    output_prefix,
    path_txinfo,
    type_rep="longest",
    path_exp=None,
    ignore_txversion=False,
    nts=3,
    keep=False,
    frac=1.0,
    offset_min=10,
    offset_max=14,
    max_depth=3,
    min_samples_split=6,
):
    """
    Train P-site offset prediction model.
    """

    log = open(f"{output_prefix}.log", "wt")

    print(f"# TxInfo: {path_txinfo}", file=log)
    print(f"# Transcript quant: {path_exp}", file=log)
    print(f"# Alignment file: {path_bam}", file=log)
    print(f"# Representative transcript type: {type_rep}", file=log)
    print(f"# Ignore transcript version: {ignore_txversion}", file=log)
    print(f"# nts: {nts}", file=log)
    print(f"# train_offset_min: {offset_min}", file=log)
    print(f"# train_offset_max: {offset_max}", file=log)
    print(f"# train_max_depth: {max_depth}", file=log)
    print(f"# train_min_samples_split: {min_samples_split}", file=log)
    print(f"# frac: {frac}", file=log)

    print("# ...load reference fasta", file=log)
    ref = read_fasta(path=path_ref, ignore_version=ignore_txversion)

    print("# ...load txinfo", file=log)
    tx_info = read_txinfo(path_txinfo, "\t")

    if ignore_txversion:
        tx_info["tx_name"] = tx_info["tx_name"].map(strip_version)

    if "chrom" in tx_info.columns:
        tx_info = tx_info[~tx_info.chrom.isin(["MT", "mitochondrion_genome", "chrM"])]

    if (tx_info.gene_id.str.startswith("FBgn")).sum() > 0:
        tx_info = tx_info[tx_info.gene_id.str.startswith("FBgn")]

    print("# ...get representative transcript per gene", file=log)
    txrep = get_txrep(tx_info, type_rep, path_exp, ignore_txversion)

    txrep.to_csv(
        f"{output_prefix}.representative_transcripts.before_filter.tsv",
        sep="\t",
        index=False,
    )

    print("# ...keep transcripts common in fasta and txinfo", file=log)

    txrep = txrep[txrep["tx_name"].isin(ref.keys())].copy()

    txrep["tx_len"] = pd.to_numeric(txrep["tx_len"], errors="coerce").fillna(0)
    txrep["cds_len"] = pd.to_numeric(txrep["cds_len"], errors="coerce").fillna(0)
    txrep["utr5_len"] = pd.to_numeric(txrep["utr5_len"], errors="coerce").fillna(0)

    ref_len = np.array([len(ref[i]) for i in txrep["tx_name"]])
    txrep["ref_fasta_len"] = ref_len

    txrep_mismatch = txrep[txrep["tx_len"] != txrep["ref_fasta_len"]].copy()
    txrep_mismatch.to_csv(f"{output_prefix}.tx_len_mismatch.tsv", sep="\t", index=False)

    txrep = txrep[txrep["tx_len"] == txrep["ref_fasta_len"]].copy()

    if txrep.empty:
        log.close()
        raise ValueError(
            "No representative transcript remains after fasta / tx_len filtering."
        )

    ref = {k: ref[k] for k in txrep["tx_name"]}

    txrep["start_codon"] = txrep["utr5_len"] + 1
    txrep["stop_codon"] = txrep["utr5_len"] + txrep["cds_len"] - 2

    txrep.to_csv(
        f"{output_prefix}.representative_transcripts.tsv", sep="\t", index=False
    )

    print("# Representative transcripts used:", txrep.shape[0], file=log)
    print(sub(r"(^|\n)", r"\1# ", str(txrep.head())), file=log)

    print("# ...parse alignments", file=log)
    alignments = extract_features(path_bam, ref, nts=nts, frac=frac)

    alignments = alignments[alignments.tx_name.isin(txrep.tx_name)].copy()

    if alignments.empty:
        log.close()
        raise ValueError(
            "No alignment remains after representative transcript filtering."
        )

    alignments.drop_duplicates(subset="rname", inplace=True)
    alignments.drop(columns="rname", inplace=True)

    nuc_cols = [f"s{i}" for i in range(2 * nts)]

    for col in nuc_cols:
        alignments[col] = pd.Categorical(
            alignments[col], categories=["A", "C", "G", "T"]
        )

    data = pd.merge(
        alignments, txrep[["tx_name", "start_codon", "stop_codon"]], on="tx_name"
    )

    data = data.assign(
        dist_start=lambda x: x.start_codon - x.tstart,
        dist_stop=lambda x: x.stop_codon - x.tstart - 3,
    )

    del alignments

    if keep:
        data.to_csv(
            f"{output_prefix}.training_data.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    print("# ...draw read length distribution", file=log)

    qwidth_freq = data.qwidth.value_counts().sort_index()
    qwidth_freq = qwidth_freq[(qwidth_freq.index >= 25) & (qwidth_freq.index <= 40)]

    if qwidth_freq.empty:
        log.close()
        raise ValueError("No reads with length between 25 and 40 nt.")

    qwidth_freq = qwidth_freq / qwidth_freq.sum()

    ax = qwidth_freq.plot(kind="bar")
    ax.figure.savefig(f"{output_prefix}.qwidth_distribution.pdf")
    plt.close(ax.figure)

    qwidth_range = qwidth_freq.index[qwidth_freq > 0.05]

    if len(qwidth_range) == 0:
        log.close()
        raise ValueError("No read length with frequency > 0.05.")

    print(f"# Read lengths used: {list(qwidth_range)}", file=log)

    data = data[data.qwidth.isin(qwidth_range)].copy()

    print("# ...draw offset distribution", file=log)

    tmp = data[data.dist_start.between(8, 20)].copy()

    if not tmp.empty:
        g = sns.FacetGrid(
            tmp, col="qwidth", col_wrap=5, sharex=False, sharey=False, height=3
        )

        g.map(sns.histplot, "dist_start", binwidth=0.25)
        g.set(xticks=[8, 10, 12, 14, 16, 18, 20])
        g.savefig(f"{output_prefix}.offset_distribution.pdf")
        plt.close("all")

    del tmp

    print("# ...filter alignments for training", file=log)

    train_start = data[data.dist_start.between(offset_min, offset_max)].copy()
    train_start["dist"] = train_start.dist_start

    train_stop = data[data.dist_stop.between(offset_min, offset_max)].copy()
    train_stop["dist"] = train_stop.dist_stop

    train = pd.concat([train_start, train_stop], axis=0)

    print(f"# Training alignments around start codon: {train_start.shape[0]}", file=log)
    print(f"# Training alignments around stop codon: {train_stop.shape[0]}", file=log)
    print(f"# Training alignments total: {train.shape[0]}", file=log)

    del train_start, train_stop

    if train.empty:
        log.close()
        raise ValueError("No alignment remains for training.")

    print("# ...train DAM and GBT", file=log)

    train_features = ["qwidth"] + nuc_cols

    X = pd.get_dummies(train[train_features], drop_first=True)
    y = train.dist.to_numpy()

    dam = DistArgMax()
    dam.fit(X, y)
    dam.qwidth_range = qwidth_range.to_numpy()
    dam.nts = nts

    dump(dam, open(f"{output_prefix}.dam.pickle", "wb"))

    gbt = GradientBoostingClassifier(
        n_estimators=100,
        learning_rate=0.1,
        min_samples_split=min_samples_split,
        max_depth=max_depth,
    )

    gbt.fit(X, y)
    gbt.qwidth_range = qwidth_range.to_numpy()
    gbt.nts = nts

    dump(gbt, open(f"{output_prefix}.gbt.pickle", "wb"))

    print("# ...frame distribution test", file=log)

    cds_aligns = data[(data.dist_start < 9) & (data.dist_stop > 15)].copy()

    if not cds_aligns.empty:
        cds_X = pd.get_dummies(cds_aligns[train_features], drop_first=True)

        for col in X.columns:
            if col not in cds_X.columns:
                cds_X[col] = 0

        cds_X = cds_X[X.columns]

        dam_cds_yhat = dam.predict(cds_X)
        gbt_cds_yhat = gbt.predict(cds_X)

        frame_df = pd.DataFrame(
            {
                "none": frame_test(cds_aligns, 0),
                "DAM": frame_test(cds_aligns, dam_cds_yhat),
                "GBT": frame_test(cds_aligns, gbt_cds_yhat),
            }
        )

        frame_df = frame_df.reindex([0, 1, 2], fill_value=0)
        frame_df.index.name = "frame"

        print("# Frame distribution of CDS-internal alignments", file=log)
        print(frame_df.to_string(float_format=lambda x: f"{x:.6f}"), file=log)

        frame_df.to_csv(f"{output_prefix}.frame_distribution.tsv", sep="\t", index=True)

    feature_importance = pd.DataFrame(
        {"feat": X.columns, "GBT": gbt.feature_importances_}
    )

    print("# GBT feature importance", file=log)
    print(feature_importance.to_string(index=False), file=log)

    feature_importance.to_csv(
        f"{output_prefix}.feature_importance.tsv", sep="\t", index=False
    )

    txrep[["tx_name", "start_codon"]].to_csv(
        f"{output_prefix}.tx_startcodon.tsv", sep="\t", index=False
    )

    print("# Done!", file=log)

    log.close()


if __name__ == "__main__":
    train()
