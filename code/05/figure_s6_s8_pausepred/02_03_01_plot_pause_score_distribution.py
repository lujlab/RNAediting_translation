#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot Pause_score tail distribution from multiple PausePred CSV files. "
            "Y-axis is fraction of positions with Pause_score >= threshold."
        )
    )

    parser.add_argument(
        "--indir", required=True, help="Directory containing PausePred CSV files."
    )

    parser.add_argument(
        "--pattern",
        default="*.pausepred.csv",
        help="Glob pattern for PausePred CSV files. Default: *.pausepred.csv",
    )

    parser.add_argument(
        "--outprefix",
        required=True,
        help="Output prefix. The script writes <outprefix>.pdf, <outprefix>.png, and <outprefix>.summary.tsv",
    )

    parser.add_argument(
        "--xmax", type=float, default=50, help="Maximum x-axis value. Default: 50"
    )

    parser.add_argument(
        "--xmin", type=float, default=0, help="Minimum x-axis value. Default: 0"
    )

    parser.add_argument(
        "--ymax", type=float, default=0.8, help="Maximum y-axis value. Default: 0.8"
    )

    parser.add_argument(
        "--ymin", type=float, default=0, help="Minimum y-axis value. Default: 0"
    )

    parser.add_argument(
        "--step",
        type=float,
        default=0.5,
        help="Step size for Pause_score thresholds. Default: 0.5",
    )

    parser.add_argument(
        "--min-reads",
        type=float,
        default=None,
        help="Optional filter: keep rows with number_of_reads_mapped >= this value.",
    )

    parser.add_argument(
        "--min-coverage",
        type=float,
        default=None,
        help="Optional filter: keep rows with coverage(%) >= this value.",
    )

    parser.add_argument(
        "--include-empty-position",
        action="store_true",
        help="By default, rows with empty coordinate_position are removed. Use this flag to keep them.",
    )

    parser.add_argument(
        "--dpi", type=int, default=300, help="PNG resolution. Default: 300"
    )

    parser.add_argument(
        "--exclude-samples",
        default="AdarG_F,AdarG_M,AdarS_F,AdarS_M",
        help="Comma-separated exact sample names to omit before reading scores.",
    )
    return parser.parse_args()


def clean_colname(x):
    return str(x).replace("\ufeff", "").strip()


def find_pausepred_header_line(path):
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if line.startswith("gene_name,") or line.startswith("\ufeffgene_name,"):
                return i
    raise RuntimeError(
        f"Cannot find PausePred header line starting with gene_name, in file: {path}"
    )


def read_pausepred_scores(
    path, min_reads=None, min_coverage=None, include_empty_position=False
):
    header_idx = find_pausepred_header_line(path)

    wanted = {
        "gene_name",
        "coordinate_position",
        "number_of_reads_mapped",
        "Pause_score",
        "coverage(%)",
        "coverage",
        "Z-score",
    }

    def usecols(c):
        return clean_colname(c) in wanted

    try:
        df = pd.read_csv(
            path,
            skiprows=header_idx,
            usecols=usecols,
            dtype=str,
            low_memory=False,
            on_bad_lines="skip",
        )
    except TypeError:
        df = pd.read_csv(
            path, skiprows=header_idx, usecols=usecols, dtype=str, low_memory=False
        )

    df.columns = [clean_colname(c) for c in df.columns]

    if "coverage(%)" in df.columns:
        df = df.rename(columns={"coverage(%)": "coverage"})

    required = ["Pause_score"]

    for c in required:
        if c not in df.columns:
            raise RuntimeError(
                f"Missing column {c} in {path}. Columns: {df.columns.tolist()}"
            )

    if not include_empty_position and "coordinate_position" in df.columns:
        df = df[
            df["coordinate_position"].notna()
            & (df["coordinate_position"].astype(str).str.strip() != "")
        ].copy()

    df["Pause_score"] = pd.to_numeric(df["Pause_score"], errors="coerce")

    if "number_of_reads_mapped" in df.columns:
        df["number_of_reads_mapped"] = pd.to_numeric(
            df["number_of_reads_mapped"], errors="coerce"
        )

    if "coverage" in df.columns:
        df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")

    df = df[df["Pause_score"].notna()].copy()
    df = df[df["Pause_score"] >= 0].copy()

    if min_reads is not None and "number_of_reads_mapped" in df.columns:
        df = df[df["number_of_reads_mapped"] >= min_reads].copy()

    if min_coverage is not None and "coverage" in df.columns:
        df = df[df["coverage"] >= min_coverage].copy()

    return df["Pause_score"].to_numpy(dtype=float)


def infer_sample_name(path):
    base = os.path.basename(path)

    for suffix in [".pausepred.csv", ".longestCDS.pausepred.csv", ".csv"]:
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break

    return base


def main():
    args = parse_args()

    indir = Path(args.indir)
    outprefix = Path(args.outprefix)

    outprefix.parent.mkdir(parents=True, exist_ok=True)

    excluded = {s.strip() for s in args.exclude_samples.split(",") if s.strip()}
    files = [
        p
        for p in sorted(glob.glob(str(indir / args.pattern)))
        if infer_sample_name(p) not in excluded
    ]

    if not files:
        raise RuntimeError(f"No files found: {indir}/{args.pattern}")

    print("[INFO] Found files:")
    for f in files:
        print("  " + f)

    thresholds = np.arange(args.xmin, args.xmax + args.step, args.step)

    summary_records = []

    fig, ax = plt.subplots(figsize=(7.2, 5.2))

    for f in files:
        sample = infer_sample_name(f)

        scores = read_pausepred_scores(
            path=f,
            min_reads=args.min_reads,
            min_coverage=args.min_coverage,
            include_empty_position=args.include_empty_position,
        )

        if len(scores) == 0:
            print(f"[WARN] No valid Pause_score values in {f}; skipped.")
            continue

        rel_freq = np.array([(scores >= t).sum() / len(scores) for t in thresholds])

        ax.plot(thresholds, rel_freq, linewidth=1.8, label=sample)

        summary_records.append(
            {
                "sample": sample,
                "file": f,
                "n_positions": len(scores),
                "min_pause_score": np.min(scores),
                "median_pause_score": np.median(scores),
                "mean_pause_score": np.mean(scores),
                "max_pause_score": np.max(scores),
                "fraction_score_ge_5": np.mean(scores >= 5),
                "fraction_score_ge_10": np.mean(scores >= 10),
                "fraction_score_ge_20": np.mean(scores >= 20),
                "fraction_score_ge_50": np.mean(scores >= 50),
            }
        )

    ax.set_xlim(args.xmin, args.xmax)
    ax.set_ylim(args.ymin, args.ymax)

    ax.set_xlabel("Pause Scores")
    ax.set_ylabel("Relative Frequency")
    ax.set_title("Pause_score distribution across samples")

    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()

    pdf_out = str(outprefix) + ".pdf"
    png_out = str(outprefix) + ".png"
    summary_out = str(outprefix) + ".summary.tsv"

    fig.savefig(pdf_out)
    fig.savefig(png_out, dpi=args.dpi)
    plt.close(fig)

    summary_df = pd.DataFrame(summary_records)
    summary_df.to_csv(summary_out, sep="\t", index=False)

    print(f"[INFO] PDF written: {pdf_out}")
    print(f"[INFO] PNG written: {png_out}")
    print(f"[INFO] Summary written: {summary_out}")


if __name__ == "__main__":
    main()
