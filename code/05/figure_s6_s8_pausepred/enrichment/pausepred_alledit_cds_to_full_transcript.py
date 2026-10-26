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
            "Test whether CDS editing sites from an alledit cdspos table are "
            "enriched in PausePred-defined RPF peak sites."
        )
    )

    parser.add_argument(
        "--pausepred",
        required=True,
        help="PausePred output CSV. First line may be 'Command used: ...'.",
    )

    parser.add_argument(
        "--alledit-cds",
        required=True,
        help=(
            "alledit cdspos TSV containing new_transcript and "
            "new_cds_pos_1based columns."
        ),
    )

    parser.add_argument(
        "--cds-intervals",
        required=True,
        help=(
            "Whitespace-delimited, headerless CDS interval file with four columns: "
            "transcript_id, CDS_start_1based, CDS_end_1based, feature. "
            "Coordinates refer to the UTR-containing transcript used by PausePred."
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


def read_cds_intervals(path):
    """
    Read the CDS interval in the UTR-containing transcript coordinate system.

    Expected columns, without a header:
        transcript_id  CDS_start_1based  CDS_end_1based  feature

    Example:
        FBtr0193027  68  1225  CDS
    """
    cds = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        comment="#",
        dtype=str,
        engine="python",
    )

    if cds.shape[1] < 3:
        raise RuntimeError(
            "CDS interval file must contain at least three columns: "
            "transcript_id, CDS_start_1based, CDS_end_1based."
        )

    cds = cds.iloc[:, :4].copy()
    while cds.shape[1] < 4:
        cds[cds.shape[1]] = "CDS"

    cds.columns = [
        "transcript_id",
        "cds_start_1based",
        "cds_end_1based",
        "feature",
    ]

    cds["transcript_id"] = cds["transcript_id"].astype(str).str.strip()
    cds["cds_start_1based"] = pd.to_numeric(cds["cds_start_1based"], errors="coerce")
    cds["cds_end_1based"] = pd.to_numeric(cds["cds_end_1based"], errors="coerce")

    cds = cds.loc[
        cds["transcript_id"].str.startswith("FBtr", na=False)
        & cds["cds_start_1based"].notna()
        & cds["cds_end_1based"].notna()
        & cds["feature"].astype(str).str.upper().eq("CDS")
    ].copy()

    cds["cds_start_1based"] = cds["cds_start_1based"].astype(int)
    cds["cds_end_1based"] = cds["cds_end_1based"].astype(int)

    invalid = cds[
        (cds["cds_start_1based"] < 1)
        | (cds["cds_end_1based"] < cds["cds_start_1based"])
    ]
    if not invalid.empty:
        raise RuntimeError(
            "Invalid 1-based CDS intervals found. Example:\n"
            + invalid.head(10).to_string(index=False)
        )

    # Allow duplicated identical rows, but reject conflicting intervals.
    conflict = (
        cds.groupby("transcript_id")
        .agg(
            n_start=("cds_start_1based", "nunique"),
            n_end=("cds_end_1based", "nunique"),
        )
        .reset_index()
    )
    conflict = conflict[(conflict["n_start"] > 1) | (conflict["n_end"] > 1)]
    if not conflict.empty:
        raise RuntimeError(
            "Conflicting CDS intervals found for transcript IDs. Example:\n"
            + conflict.head(10).to_string(index=False)
        )

    cds = cds.drop_duplicates("transcript_id", keep="first").copy()
    cds["cds_length_from_interval"] = (
        cds["cds_end_1based"] - cds["cds_start_1based"] + 1
    )

    return cds


def read_alledit_cds_extract_editing_sites(path, cds_intervals):
    """
    Extract the species-level CDS editing sites from the alledit/cdspos table
    and convert CDS-relative coordinates to the UTR-containing transcript
    coordinates used by PausePred.

    The supplied cdspos table stores ``new_cds_pos_1based`` relative to the
    CDS.  The PausePred files store a full-transcript coordinate, therefore:

        transcript_pos = cds_start_1based + new_cds_pos_1based - 1

    ``new_tx_pos_1based`` is deliberately not used here: for transcripts with
    a 5' UTR it is still a CDS-relative coordinate in this table.  The CDS
    interval file supplies the required UTR offset and also validates the
    converted position.
    """
    alledit = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
    )
    alledit.columns = [clean_col(c) for c in alledit.columns]

    required = {"new_transcript", "new_cds_pos_1based"}
    missing = sorted(required - set(alledit.columns))
    if missing:
        raise RuntimeError(
            "Missing required columns in alledit cdspos table: "
            + ",".join(missing)
            + "\nAvailable columns: "
            + ",".join(alledit.columns)
        )

    input_rows = len(alledit)
    if "target_region" in alledit.columns:
        region = alledit["target_region"].astype(str).str.upper().str.strip()
        before_region = len(alledit)
        alledit = alledit[region.eq("CDS")].copy()
        print(
            "[INFO] alledit rows retained with target_region=CDS: "
            f"{len(alledit)} / {before_region}"
        )

    alledit["transcript_id"] = alledit["new_transcript"].astype(str).str.strip()
    alledit["alledit_cds_pos"] = pd.to_numeric(
        alledit["new_cds_pos_1based"], errors="coerce"
    )

    cds_lookup = cds_intervals.set_index("transcript_id").to_dict("index")
    records = []
    skipped_non_fbtr = 0
    skipped_no_cds_pos = 0
    skipped_no_interval = 0
    skipped_outside_interval = 0

    for row_idx, row in alledit.iterrows():
        transcript_id = str(row["transcript_id"]).strip()
        if not transcript_id.startswith("FBtr"):
            skipped_non_fbtr += 1
            continue

        cds_pos_value = row["alledit_cds_pos"]
        if pd.isna(cds_pos_value) or int(cds_pos_value) < 1:
            skipped_no_cds_pos += 1
            continue
        cds_pos = int(cds_pos_value)

        interval = cds_lookup.get(transcript_id)
        if interval is None:
            skipped_no_interval += 1
            continue

        cds_start = int(interval["cds_start_1based"])
        cds_end = int(interval["cds_end_1based"])
        cds_length_from_interval = int(interval["cds_length_from_interval"])
        transcript_pos = cds_start + cds_pos - 1

        if transcript_pos < cds_start or transcript_pos > cds_end:
            skipped_outside_interval += 1
            continue

        records.append(
            {
                "site_key": f"{transcript_id}_{transcript_pos}",
                "transcript_id": transcript_id,
                "transcript_pos": transcript_pos,
                "coord_source": "alledit_new_cds_pos_plus_CDS_start",
                "alledit_cds_pos": cds_pos,
                "cds_start_1based": cds_start,
                "cds_end_1based": cds_end,
                "cds_length_from_interval": cds_length_from_interval,
                "alledit_row_index": int(row_idx) + 2,
                "target_region": row.get("target_region", "CDS"),
                "genomic_key": row.get("genomic_key", ""),
                "variant_effect": row.get("variant_effect", ""),
            }
        )

    edit = pd.DataFrame(records)

    print(f"[INFO] alledit input rows: {input_rows}")
    print(
        f"[INFO] alledit rows skipped because transcript was not FBtr: {skipped_non_fbtr}"
    )
    print(
        f"[INFO] alledit rows skipped because CDS position was absent/invalid: {skipped_no_cds_pos}"
    )
    print(
        f"[INFO] alledit rows skipped because CDS interval was absent: {skipped_no_interval}"
    )
    print(
        f"[INFO] alledit rows skipped because converted position was outside CDS: {skipped_outside_interval}"
    )

    if edit.empty:
        return edit

    edit = edit.sort_values(
        ["transcript_id", "transcript_pos", "alledit_row_index"]
    ).drop_duplicates(["site_key"], keep="first")
    print(
        f"[INFO] Unique alledit CDS editing sites after coordinate conversion: {len(edit)}"
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

    print(f"[INFO] Reading CDS intervals: {args.cds_intervals}")
    cds_intervals = read_cds_intervals(args.cds_intervals)
    print(f"[INFO] CDS intervals loaded: {len(cds_intervals)}")

    print(f"[INFO] Reading alledit CDS table: {args.alledit_cds}")
    edit = read_alledit_cds_extract_editing_sites(
        path=args.alledit_cds,
        cds_intervals=cds_intervals,
    )

    if edit.empty:
        print("[WARN] No editing sites extracted from the alledit CDS table.")
        editing_set_all = set()
    else:
        editing_set_all = set(edit["site_key"])

    universe_set = set(pause["site_key"])
    editing_in_universe = editing_set_all & universe_set

    print(
        f"[INFO] Editing sites extracted from alledit CDS after coordinate conversion: {len(editing_set_all)}"
    )
    print(
        f"[INFO] Converted editing sites overlapping PausePred universe: {len(editing_in_universe)}"
    )

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
                    "editing_total_extracted_from_alledit_cds": len(editing_set_all),
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
                "alledit_cds_pos",
                "cds_start_1based",
                "cds_end_1based",
                "cds_length_from_interval",
                "alledit_row_index",
                "target_region",
                "genomic_key",
                "variant_effect",
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
                    "alledit_cds_pos",
                    "cds_start_1based",
                    "cds_end_1based",
                    "cds_length_from_interval",
                    "alledit_row_index",
                    "target_region",
                    "genomic_key",
                    "variant_effect",
                ]
            ).to_csv(edit_out, sep="\t", index=False)
            pd.DataFrame(
                columns=[
                    "site_key",
                    "transcript_id",
                    "transcript_pos",
                    "coord_source",
                    "alledit_cds_pos",
                    "cds_start_1based",
                    "cds_end_1based",
                    "cds_length_from_interval",
                    "alledit_row_index",
                    "target_region",
                    "genomic_key",
                    "variant_effect",
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
