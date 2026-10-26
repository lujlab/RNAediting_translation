#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Merge per-sample PausePred enrichment tables without BH correction."""

import argparse
import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge PausePred enrichment TSVs and retain raw Fisher P values."
    )
    parser.add_argument("--indir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--offset", type=float, default=None)
    parser.add_argument("--pattern", default="*.pausepred_enrichment.tsv")
    parser.add_argument(
        "--samples", default=None, help="Comma-separated exact sample names to merge"
    )
    return parser.parse_args()


def infer_sample_name(path):
    base = os.path.basename(path)
    for suffix in (
        ".pausepred_enrichment.tsv",
        ".pausepred_editing_enrichment.tsv",
    ):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return re.sub(r"\.tsv$", "", base)


def find_input_files(indir, pattern):
    files = []
    for path in sorted(glob.glob(os.path.join(indir, pattern))):
        base = os.path.basename(path)
        if base.endswith(
            (
                ".pausepred_enrichment.tsv",
                ".pausepred_editing_enrichment.tsv",
            )
        ):
            files.append(path)
    return files


def normalize_table(df, path):
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    aliases = {
        "cutoff": "Pause_score_cutoff",
        "oddsratio": "odds_ratio",
        "pvalue": "p_value",
        "n_universe": "total_sites",
        "n_pause_or_near_pause_positions": "total_sites_in_peak",
        "n_editing_in_universe": "editing_total_used_in_universe",
        "editing_and_pause": "editing_in_peak",
        "editing_and_nonpause": "editing_not_in_peak",
        "nonediting_and_pause": "nonediting_in_peak",
        "nonediting_and_nonpause": "nonediting_not_in_peak",
    }
    rename = {
        source: target
        for source, target in aliases.items()
        if source in df.columns and target not in df.columns
    }
    df = df.rename(columns=rename)

    required = [
        "Pause_score_cutoff",
        "window_nt",
        "total_sites",
        "editing_total_used_in_universe",
        "editing_in_peak",
        "odds_ratio",
        "p_value",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing columns in {path}: {missing}; available={list(df.columns)}"
        )

    numeric = [
        "Pause_score_cutoff",
        "window_nt",
        "total_sites",
        "total_sites_in_peak",
        "total_sites_not_in_peak",
        "editing_total_extracted_from_input",
        "editing_total_extracted_from_snpeff",
        "editing_total_used_in_universe",
        "editing_in_peak",
        "editing_not_in_peak",
        "nonediting_total",
        "nonediting_in_peak",
        "nonediting_not_in_peak",
        "odds_ratio",
        "p_value",
        "p_adj_BH",
    ]
    for column in numeric:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    if "editing_not_in_peak" not in df.columns:
        df["editing_not_in_peak"] = (
            df["editing_total_used_in_universe"] - df["editing_in_peak"]
        )
    if "nonediting_total" not in df.columns:
        df["nonediting_total"] = (
            df["total_sites"] - df["editing_total_used_in_universe"]
        )
    if "nonediting_in_peak" not in df.columns:
        df["nonediting_in_peak"] = df["total_sites_in_peak"] - df["editing_in_peak"]
    if "nonediting_not_in_peak" not in df.columns:
        df["nonediting_not_in_peak"] = df["nonediting_total"] - df["nonediting_in_peak"]
    if "total_sites_in_peak" not in df.columns:
        df["total_sites_in_peak"] = df["editing_in_peak"] + df["nonediting_in_peak"]
    if "total_sites_not_in_peak" not in df.columns:
        df["total_sites_not_in_peak"] = df["total_sites"] - df["total_sites_in_peak"]

    # This analysis deliberately does not use multiple-testing correction.
    df["p_adj_BH"] = np.nan
    df["multiple_testing_correction"] = "none"
    return df


def add_fraction_columns(df):
    df = df.copy()
    df["editing_peak_fraction"] = (
        df["editing_in_peak"] / df["editing_total_used_in_universe"]
    )
    df["nonediting_peak_fraction"] = df["nonediting_in_peak"] / df["nonediting_total"]
    df["editing_peak_percent"] = df["editing_peak_fraction"] * 100
    df["nonediting_peak_percent"] = df["nonediting_peak_fraction"] * 100
    return df


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    files = find_input_files(args.indir, args.pattern)
    if args.samples is not None:
        selected = {s.strip() for s in args.samples.split(",") if s.strip()}
        files = [p for p in files if infer_sample_name(p) in selected]
        missing = selected - {infer_sample_name(p) for p in files}
        if missing:
            raise ValueError(
                "Missing selected sample tables: " + ", ".join(sorted(missing))
            )
    if not files:
        raise RuntimeError(f"No enrichment tables found in {args.indir}")

    all_tables = []
    for path in files:
        table = normalize_table(pd.read_csv(path, sep="\t", low_memory=False), path)
        if args.offset is not None:
            table = table[np.isclose(table["Pause_score_cutoff"], args.offset)].copy()
        if table.empty:
            continue
        table["sample"] = infer_sample_name(path)
        all_tables.append(table)

    if not all_tables:
        raise RuntimeError("No rows remain after cutoff filtering")
    merged = add_fraction_columns(pd.concat(all_tables, ignore_index=True))

    preferred = [
        "sample",
        "Pause_score_cutoff",
        "window_nt",
        "min_reads_for_pause",
        "min_coverage_for_pause",
        "total_sites_original_pausepred_universe",
        "total_sites_original_unedited",
        "a_codon_sites_in_pausepred_universe",
        "unedited_A_codon_sites_removed",
        "total_sites",
        "total_sites_in_peak",
        "total_sites_not_in_peak",
        "editing_total_extracted_from_input",
        "editing_total_extracted_from_snpeff",
        "editing_total_used_in_universe",
        "editing_in_peak",
        "editing_not_in_peak",
        "editing_peak_percent",
        "nonediting_total",
        "nonediting_in_peak",
        "nonediting_not_in_peak",
        "nonediting_peak_percent",
        "odds_ratio",
        "p_value",
        "p_adj_BH",
        "multiple_testing_correction",
        "unedit_background_filter",
    ]
    preferred = [c for c in preferred if c in merged.columns]
    rest = [c for c in merged.columns if c not in preferred]
    merged = merged[preferred + rest]

    output = outdir / "all_samples.pausepred_enrichment.merged.tsv"
    merged.to_csv(output, sep="\t", index=False)
    print(f"[INFO] Input tables: {len(all_tables)}")
    print(f"[INFO] Rows written: {len(merged)}")
    print(f"[INFO] BH correction: none; raw p_value retained")
    print(f"[INFO] Merged output: {output}")


if __name__ == "__main__":
    main()
