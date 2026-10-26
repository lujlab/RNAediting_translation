#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Add a longest-transcript flag to a mapped transcript table and "
            "create one deduplicated record per query/genomic site."
        )
    )
    parser.add_argument(
        "--mapped-tsv",
        required=True,
        help="Input mapped.other_transcript_blocks.tsv file.",
    )
    parser.add_argument(
        "--longest-list",
        required=True,
        help="Text file containing one longest-CDS transcript ID per line.",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix.",
    )
    parser.add_argument(
        "--strip-version",
        action="store_true",
        help="Ignore transcript version suffixes, e.g. FBtr0000001.1 -> FBtr0000001.",
    )
    parser.add_argument(
        "--fallback",
        choices=["exclude", "first", "source"],
        default="exclude",
        help=(
            "Action for a multi-row site when no target transcript has longest=1: "
            "'exclude' writes it only to unresolved output; "
            "'first' keeps the first row; "
            "'source' prefers target_transcript == source_transcript, otherwise first. "
            "Default: exclude."
        ),
    )
    return parser.parse_args()


def normalize_transcript_id(value, strip_version):
    value = value.strip()
    if strip_version and "." in value:
        return value.split(".", 1)[0]
    return value


def load_longest_ids(path, strip_version):
    longest_ids = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            transcript_id = line.split()[0]
            longest_ids.add(normalize_transcript_id(transcript_id, strip_version))
    return longest_ids


def main():
    args = parse_args()

    longest_ids = load_longest_ids(args.longest_list, args.strip_version)
    if not longest_ids:
        raise ValueError(
            "No transcript IDs were read from the longest-transcript list."
        )

    out_with_flag = args.out_prefix + ".with_longest.tsv"
    out_dedup = args.out_prefix + ".dedup_longest.tsv"
    out_unresolved = args.out_prefix + ".dedup_unresolved.tsv"
    out_summary = args.out_prefix + ".dedup_summary.tsv"

    with open(args.mapped_tsv, "r", encoding="utf-8", newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("The mapped TSV file has no header.")

        required = {"query_name", "source_transcript", "target_transcript"}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                "Missing required column(s): " + ", ".join(sorted(missing))
            )

        fieldnames = list(reader.fieldnames)
        if "longest" not in fieldnames:
            fieldnames.append("longest")

        grouped = OrderedDict()

        with open(out_with_flag, "w", encoding="utf-8", newline="") as flag_handle:
            writer = csv.DictWriter(
                flag_handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()

            for row in reader:
                target_id = normalize_transcript_id(
                    row["target_transcript"], args.strip_version
                )
                row["longest"] = "1" if target_id in longest_ids else "0"
                writer.writerow(row)
                grouped.setdefault(row["query_name"], []).append(row)

    kept_rows = []
    unresolved_rows = []

    n_single = 0
    n_multi = 0
    n_multi_selected_longest = 0
    n_multi_no_longest = 0
    n_multiple_longest = 0

    for query_name, rows in grouped.items():
        if len(rows) == 1:
            kept_rows.append(rows[0])
            n_single += 1
            continue

        n_multi += 1
        longest_rows = [row for row in rows if row["longest"] == "1"]

        if longest_rows:
            # Normally there should be one longest-CDS transcript per gene.
            # If multiple rows remain, keep all and report the situation.
            kept_rows.extend(longest_rows)
            n_multi_selected_longest += 1
            if len(longest_rows) > 1:
                n_multiple_longest += 1
            continue

        n_multi_no_longest += 1
        for row in rows:
            unresolved = dict(row)
            unresolved["unresolved_reason"] = "MULTIROW_NO_LONGEST_TARGET"
            unresolved_rows.append(unresolved)

        if args.fallback == "first":
            kept_rows.append(rows[0])
        elif args.fallback == "source":
            source_matches = [
                row
                for row in rows
                if row["target_transcript"] == row["source_transcript"]
            ]
            kept_rows.append(source_matches[0] if source_matches else rows[0])
        # fallback=exclude: do not add any row to kept_rows

    with open(out_dedup, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(
            out_handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(kept_rows)

    unresolved_fields = fieldnames + ["unresolved_reason"]
    with open(out_unresolved, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(
            out_handle,
            fieldnames=unresolved_fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(unresolved_rows)

    summary_rows = [
        ("longest_transcript_ids", len(longest_ids)),
        ("input_rows", sum(len(rows) for rows in grouped.values())),
        ("unique_query_names", len(grouped)),
        ("single_row_queries_kept", n_single),
        ("multi_row_queries", n_multi),
        ("multi_row_queries_selected_longest", n_multi_selected_longest),
        ("multi_row_queries_without_longest_target", n_multi_no_longest),
        ("queries_with_multiple_longest_rows", n_multiple_longest),
        ("deduplicated_output_rows", len(kept_rows)),
        ("fallback_mode", args.fallback),
    ]
    with open(out_summary, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(summary_rows)

    print("Longest IDs loaded:", len(longest_ids))
    print("Input rows:", sum(len(rows) for rows in grouped.values()))
    print("Unique query names:", len(grouped))
    print("Deduplicated rows:", len(kept_rows))
    print("With longest flag:", out_with_flag)
    print("Deduplicated table:", out_dedup)
    print("Unresolved multi-row sites:", out_unresolved)
    print("Summary:", out_summary)


if __name__ == "__main__":
    raise SystemExit(main())
