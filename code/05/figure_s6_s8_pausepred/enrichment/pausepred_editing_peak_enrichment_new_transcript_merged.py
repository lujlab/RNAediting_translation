#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Test whether editing sites from a mapped editing-site table are enriched in "
            "PausePred-defined RPF peak sites."
        )
    )

    parser.add_argument(
        "--pausepred",
        required=True,
        help="PausePred output CSV. First line may be 'Command used: ...'.",
    )

    parser.add_argument(
        "--edit-sites",
        "--snpeff",
        dest="edit_sites",
        required=True,
        help=(
            "New mapped editing-site TSV. Required columns: "
            "new_transcript and new_tx_pos_1based. "
            "The old --snpeff option name is retained as an alias."
        ),
    )

    parser.add_argument(
        "--target-region",
        default="all",
        choices=["all", "CDS", "UTR5", "UTR3"],
        help=(
            "Optional target_region filter for the new editing-site table. "
            "Default: all rows with valid new transcript coordinates."
        ),
    )

    parser.add_argument(
        "--out", required=True, help="Output Fisher enrichment result TSV."
    )

    parser.add_argument(
        "--min-reads",
        type=float,
        default=2,
        help="Minimum number_of_reads_mapped for defining the background universe. Default: 2.",
    )

    parser.add_argument(
        "--min-coverage",
        type=float,
        default=20,
        help="Minimum coverage(%%) for defining the background universe. Default: 20.",
    )

    parser.add_argument(
        "--pause-cutoffs",
        default="5,10,20",
        help="Comma-separated Pause_score cutoffs. Default: 5,10,20.",
    )

    parser.add_argument(
        "--windows",
        default="0",
        help=(
            "Comma-separated window sizes in nt around PausePred peaks. "
            "0 means exact overlap. Example: 0,9,30. Default: 0."
        ),
    )

    parser.add_argument(
        "--write-intermediate",
        action="store_true",
        help="Write extracted editing sites, universe sites, and editing-in-peak sites.",
    )

    return parser.parse_args()


def clean_col(x):
    return str(x).replace("\ufeff", "").strip()


def find_pausepred_header_line(path):
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if line.startswith("gene_name,") or line.startswith("\ufeffgene_name,"):
                return i, line.rstrip("\n")
    raise RuntimeError(
        "Cannot find PausePred header line starting with 'gene_name,'. "
        "Please check the input PausePred file."
    )


def first_existing(clean_header_cols, candidates):
    col_set = set(clean_header_cols)
    for c in candidates:
        if c in col_set:
            return c
    raise RuntimeError(
        "Cannot find any of these columns: "
        + ",".join(candidates)
        + "\nAvailable columns: "
        + ",".join(clean_header_cols)
    )


def read_pausepred(path, min_reads, min_coverage):
    header_idx, header_line = find_pausepred_header_line(path)
    header_cols = [clean_col(x) for x in header_line.split(",")]

    transcript_col = first_existing(header_cols, ["gene_name", "transcript_id"])
    pos_col = first_existing(
        header_cols, ["coordinate_position", "transcript_pos", "position"]
    )
    reads_col = first_existing(header_cols, ["number_of_reads_mapped", "reads"])
    score_col = first_existing(header_cols, ["Pause_score", "pause_score"])
    cov_col = first_existing(header_cols, ["coverage(%)", "coverage"])

    required_clean_cols = {
        transcript_col,
        pos_col,
        reads_col,
        score_col,
        cov_col,
    }

    def use_column(c):
        return clean_col(c) in required_clean_cols

    try:
        pause = pd.read_csv(
            path,
            skiprows=header_idx,
            usecols=use_column,
            dtype=str,
            low_memory=False,
            on_bad_lines="warn",
        )
    except TypeError:
        # compatibility for older pandas
        pause = pd.read_csv(
            path,
            skiprows=header_idx,
            usecols=use_column,
            dtype=str,
            low_memory=False,
        )

    pause.columns = [clean_col(c) for c in pause.columns]

    rename_map = {
        transcript_col: "transcript_id",
        pos_col: "transcript_pos",
        reads_col: "reads",
        score_col: "Pause_score",
        cov_col: "coverage",
    }

    pause = pause.rename(columns=rename_map)

    required_after = [
        "transcript_id",
        "transcript_pos",
        "reads",
        "Pause_score",
        "coverage",
    ]
    missing = [c for c in required_after if c not in pause.columns]
    if missing:
        raise RuntimeError(
            "Missing required columns after renaming: "
            + ",".join(missing)
            + "\nCurrent columns: "
            + ",".join(pause.columns)
        )

    pause["transcript_id"] = pause["transcript_id"].astype(str).str.strip()
    pause["transcript_pos"] = pd.to_numeric(pause["transcript_pos"], errors="coerce")
    pause["reads"] = pd.to_numeric(pause["reads"], errors="coerce")
    pause["Pause_score"] = pd.to_numeric(pause["Pause_score"], errors="coerce")
    pause["coverage"] = pd.to_numeric(pause["coverage"], errors="coerce")

    valid = (
        pause["transcript_id"].str.startswith("FBtr", na=False)
        & pause["transcript_pos"].notna()
        & pause["reads"].notna()
        & pause["Pause_score"].notna()
        & pause["coverage"].notna()
        & (pause["transcript_pos"] > 0)
        & (pause["reads"] >= min_reads)
        & (pause["coverage"] >= min_coverage)
    )

    pause = pause.loc[
        valid, ["transcript_id", "transcript_pos", "reads", "Pause_score", "coverage"]
    ].copy()

    pause["transcript_pos"] = pause["transcript_pos"].astype(int)

    # Same transcript-position may appear more than once; keep the highest Pause_score.
    pause = pause.sort_values(
        ["transcript_id", "transcript_pos", "Pause_score"],
        ascending=[True, True, False],
    ).drop_duplicates(["transcript_id", "transcript_pos"], keep="first")

    pause["site_key"] = (
        pause["transcript_id"].astype(str) + "_" + pause["transcript_pos"].astype(str)
    )

    return pause


def normalize_transcript_id(x):
    """Normalize transcript IDs for matching to PausePred gene_name."""
    s = str(x).strip()
    if s == "" or s.lower() in {"na", "nan", "none", "null", "."}:
        return None
    # Remove an optional transcript version suffix, e.g. FBtr0000001.1.
    return re.sub(r"\.[0-9]+$", "", s)


def read_new_editing_sites(path, target_region="all"):
    """Read the new mapped editing-site table.

    Matching key:
      new_transcript + "_" + new_tx_pos_1based

    This matches the PausePred key:
      gene_name + "_" + coordinate_position

    Both positions are 1-based transcript nucleotide coordinates including UTR.
    """
    try:
        edit = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            low_memory=False,
        )
    except Exception:
        edit = pd.read_csv(
            path,
            sep=r"\s+",
            dtype=str,
            low_memory=False,
            engine="python",
        )

    edit.columns = [clean_col(c) for c in edit.columns]

    required = [
        "new_transcript",
        "new_tx_pos_1based",
    ]

    missing = [c for c in required if c not in edit.columns]
    if missing:
        raise RuntimeError(
            "New editing-site table is missing required columns: "
            + ",".join(missing)
            + "\nAvailable columns: "
            + ",".join(edit.columns)
        )

    if target_region != "all" and "target_region" not in edit.columns:
        raise RuntimeError(
            "--target-region was set to "
            + target_region
            + ", but the input table has no target_region column."
        )

    edit = edit.copy()
    edit["row_index"] = range(1, len(edit) + 1)

    edit["transcript_id"] = edit["new_transcript"].map(normalize_transcript_id)

    edit["transcript_pos"] = pd.to_numeric(
        edit["new_tx_pos_1based"],
        errors="coerce",
    )

    if "target_region" in edit.columns:
        edit["target_region"] = (
            edit["target_region"].astype(str).str.strip().str.upper()
        )

    if target_region != "all":
        edit = edit.loc[edit["target_region"] == target_region].copy()

    valid = (
        edit["transcript_id"].notna()
        & edit["transcript_id"].str.startswith("FBtr", na=False)
        & edit["transcript_pos"].notna()
        & (edit["transcript_pos"] > 0)
    )

    edit = edit.loc[valid].copy()

    if edit.empty:
        return pd.DataFrame(
            columns=[
                "site_key",
                "transcript_id",
                "transcript_pos",
                "coord_source",
                "row_index",
                "genomic_key",
                "target_region",
                "new_cds_pos_1based",
                "type",
                "sample",
            ]
        )

    edit["transcript_pos"] = edit["transcript_pos"].astype(int)

    edit["site_key"] = (
        edit["transcript_id"].astype(str) + "_" + edit["transcript_pos"].astype(str)
    )

    edit["coord_source"] = "new_transcript+new_tx_pos_1based"

    # Preserve useful source columns when they exist.
    optional_columns = [
        "genomic_key",
        "chr",
        "pos",
        "Ref",
        "Alt",
        "type",
        "gene",
        "sample",
        "sex",
        "target_region",
        "new_cds_pos_1based",
        "new_transcript",
        "new_tx_pos_1based",
    ]

    keep_columns = [
        "site_key",
        "transcript_id",
        "transcript_pos",
        "coord_source",
        "row_index",
    ] + [
        c
        for c in optional_columns
        if c in edit.columns
        and c
        not in {
            "transcript_id",
            "transcript_pos",
            "site_key",
            "coord_source",
            "row_index",
        }
    ]

    edit = edit[keep_columns].copy()

    # Multiple samples or rows may contain the same mapped transcript position.
    # The enrichment analysis counts each transcript nt position once.
    edit = (
        edit.sort_values(["transcript_id", "transcript_pos", "row_index"])
        .drop_duplicates(["site_key"], keep="first")
        .reset_index(drop=True)
    )

    return edit


def parse_num_list(x, cast=float):
    out = []
    for item in str(x).split(","):
        item = item.strip()
        if item == "":
            continue
        out.append(cast(item))
    return out


def make_target_keys_from_peak(pause, peak, universe_set, window_nt):
    if window_nt == 0:
        return set(peak["site_key"]) & universe_set

    target = set()

    # For window overlap, expand around each peak position and keep only positions in universe.
    for row in peak.itertuples(index=False):
        tx = row.transcript_id
        pos = int(row.transcript_pos)
        start = max(1, pos - window_nt)
        end = pos + window_nt

        for p in range(start, end + 1):
            key = f"{tx}_{p}"
            if key in universe_set:
                target.add(key)

    return target


def main():
    args = parse_args()

    pause_cutoffs = parse_num_list(args.pause_cutoffs, cast=float)
    window_list = parse_num_list(args.windows, cast=int)

    print(f"[INFO] Reading PausePred: {args.pausepred}")
    pause = read_pausepred(
        path=args.pausepred,
        min_reads=args.min_reads,
        min_coverage=args.min_coverage,
    )

    print(f"[INFO] PausePred universe after filters: {len(pause)}")
    print(
        f"[INFO] Unique transcripts in PausePred universe: {pause['transcript_id'].nunique()}"
    )

    print(f"[INFO] Reading mapped editing-site table: {args.edit_sites}")
    edit = read_new_editing_sites(
        path=args.edit_sites,
        target_region=args.target_region,
    )

    if edit.empty:
        print("[WARN] No valid mapped editing sites were extracted.")
        editing_set_all = set()
    else:
        editing_set_all = set(edit["site_key"])

    universe_set = set(pause["site_key"])
    editing_in_universe = editing_set_all & universe_set

    print(f"[INFO] Editing sites extracted from mapped input: {len(editing_set_all)}")
    print(f"[INFO] Editing sites within PausePred universe: {len(editing_in_universe)}")

    results = []
    hit_records = []

    for cutoff in pause_cutoffs:
        peak = pause[pause["Pause_score"] >= cutoff].copy()

        for window_nt in window_list:
            target_keys = make_target_keys_from_peak(
                pause=pause,
                peak=peak,
                universe_set=universe_set,
                window_nt=window_nt,
            )

            total_sites = len(universe_set)
            total_sites_in_peak = len(target_keys)
            total_sites_not_in_peak = total_sites - total_sites_in_peak

            editing_total = len(editing_in_universe)
            editing_in_peak = len(editing_in_universe & target_keys)
            editing_not_in_peak = editing_total - editing_in_peak

            nonediting_total = total_sites - editing_total
            nonediting_in_peak = total_sites_in_peak - editing_in_peak
            nonediting_not_in_peak = nonediting_total - nonediting_in_peak

            if (
                min(
                    total_sites,
                    total_sites_in_peak,
                    editing_total,
                    editing_in_peak,
                    nonediting_in_peak,
                    nonediting_not_in_peak,
                )
                < 0
            ):
                raise RuntimeError(
                    "Negative value found in Fisher table. Please check input files."
                )

            table = [
                [editing_in_peak, editing_not_in_peak],
                [nonediting_in_peak, nonediting_not_in_peak],
            ]

            odds_ratio, p_value = fisher_exact(table, alternative="greater")

            results.append(
                {
                    "Pause_score_cutoff": cutoff,
                    "window_nt": window_nt,
                    "min_reads_for_pause": args.min_reads,
                    "min_coverage_for_pause": args.min_coverage,
                    "total_sites": total_sites,
                    "total_sites_in_peak": total_sites_in_peak,
                    "total_sites_not_in_peak": total_sites_not_in_peak,
                    "editing_total_extracted_from_input": len(editing_set_all),
                    "editing_total_extracted_from_snpeff": len(editing_set_all),
                    "editing_total_used_in_universe": editing_total,
                    "editing_in_peak": editing_in_peak,
                    "editing_not_in_peak": editing_not_in_peak,
                    "nonediting_total": nonediting_total,
                    "nonediting_in_peak": nonediting_in_peak,
                    "nonediting_not_in_peak": nonediting_not_in_peak,
                    "odds_ratio": odds_ratio,
                    "p_value": p_value,
                }
            )

            if edit.empty:
                continue

            hit_keys = editing_in_universe & target_keys

            if hit_keys:
                hit_df = edit[edit["site_key"].isin(hit_keys)].copy()
                hit_df["Pause_score_cutoff"] = cutoff
                hit_df["window_nt"] = window_nt
                hit_records.append(hit_df)

    out_df = pd.DataFrame(results)
    out_df["p_adj_BH"] = bh_adjust(out_df["p_value"].tolist())
    out_df.to_csv(args.out, sep="\t", index=False)

    print(f"[INFO] Main output written: {args.out}")

    out_path = Path(args.out)
    prefix = str(out_path)
    if prefix.endswith(".tsv"):
        prefix = prefix[:-4]

    hit_out = prefix + ".editing_sites_in_peak.tsv"

    if hit_records:
        hits = pd.concat(hit_records, ignore_index=True)
        hits.to_csv(hit_out, sep="\t", index=False)
    else:
        pd.DataFrame(
            columns=[
                "site_key",
                "transcript_id",
                "transcript_pos",
                "coord_source",
                "row_index",
                "genomic_key",
                "target_region",
                "new_cds_pos_1based",
                "type",
                "sample",
                "Pause_score_cutoff",
                "window_nt",
            ]
        ).to_csv(hit_out, sep="\t", index=False)

    print(f"[INFO] Editing sites in peak written: {hit_out}")
    # -----------------------------
    # Export unedited controls.
    # -----------------------------
    unedit_in_universe = universe_set - editing_in_universe
    unedit_df = pause[pause["site_key"].isin(unedit_in_universe)].copy()

    unedited_out = prefix + ".unedited_sites.csv"
    unedit_df.to_csv(unedited_out, index=False)
    print(f"[INFO] Unedited sites written: {unedited_out}")

    if args.write_intermediate:
        pause_out = prefix + ".pausepred_universe.tsv"
        edit_out = prefix + ".editing_sites_extracted.tsv"
        edit_used_out = prefix + ".editing_sites_used_in_universe.tsv"

        pause.to_csv(pause_out, sep="\t", index=False)

        if edit.empty:
            pd.DataFrame(
                columns=[
                    "site_key",
                    "transcript_id",
                    "transcript_pos",
                    "coord_source",
                    "row_index",
                    "genomic_key",
                    "target_region",
                    "new_cds_pos_1based",
                    "type",
                    "sample",
                ]
            ).to_csv(edit_out, sep="\t", index=False)
            pd.DataFrame(
                columns=[
                    "site_key",
                    "transcript_id",
                    "transcript_pos",
                    "coord_source",
                    "row_index",
                    "genomic_key",
                    "target_region",
                    "new_cds_pos_1based",
                    "type",
                    "sample",
                ]
            ).to_csv(edit_used_out, sep="\t", index=False)
        else:
            edit.to_csv(edit_out, sep="\t", index=False)
            edit[edit["site_key"].isin(editing_in_universe)].to_csv(
                edit_used_out,
                sep="\t",
                index=False,
            )

        print(f"[INFO] PausePred universe written: {pause_out}")
        print(f"[INFO] Editing sites extracted written: {edit_out}")
        print(f"[INFO] Editing sites used in universe written: {edit_used_out}")


def bh_adjust(pvalues):
    """
    Benjamini-Hochberg adjustment.
    No statsmodels dependency.
    """
    n = len(pvalues)
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])

    adjusted = [None] * n
    prev = 1.0

    for rank_from_end, (idx, p) in enumerate(reversed(indexed), start=1):
        rank = n - rank_from_end + 1
        adj = p * n / rank
        if adj > prev:
            adj = prev
        prev = adj
        adjusted[idx] = min(adj, 1.0)

    return adjusted


if __name__ == "__main__":
    main()
