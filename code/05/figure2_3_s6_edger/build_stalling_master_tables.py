#!/usr/bin/env python3
"""Stream replicate flank tables into genotype-by-sex master tables.

The input tables must have the 14-column schema emitted by
codon_flank_pipeline.py and identical codon order within a species.  Files are
processed line-by-line, so memory use is essentially independent of file size.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

KEY_COLUMNS = (
    "chrom",
    "genomic_codon_first_pos_1based",
    "transcript",
    "cds_block_index",
    "strand",
    "cds_codon_first_pos_1based",
    "codon_index_1based",
)
VALUE_COLUMNS = (
    "rna_count",
    "rpf_count",
    "flank_rna_count",
    "flank_rpf_count",
    "flank_n",
    "library_RNA_sum",
    "library_RPF_sum",
)
MANIFEST_COLUMNS = (
    "group_id",
    "species",
    "genotype",
    "sex",
    "replicate",
    "sample_id",
    "is_head",
    "results_dir",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", default="stalling_groups.tsv")
    parser.add_argument("--output-dir", default="stalling_master_tables")
    parser.add_argument("--flanks", nargs="+", type=int, default=(20, 30))
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument(
        "--only-groups",
        nargs="*",
        default=None,
        help="Optional group_id subset (default: every group).",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    return parser.parse_args()


def replicate_number(label: str) -> int:
    match = re.fullmatch(r"rep([1-9][0-9]*)", label)
    if match is None:
        raise ValueError(f"invalid replicate label {label!r}; expected rep1, rep2, ...")
    return int(match.group(1))


def load_groups(
    path: Path, selected: set[str] | None
) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MANIFEST_COLUMNS:
            raise ValueError(
                f"{path}: expected columns {MANIFEST_COLUMNS}, got {reader.fieldnames}"
            )
        for row in reader:
            if selected is None or row["group_id"] in selected:
                row["input_root"] = str((path.parent / row["results_dir"]).resolve())
                grouped[row["group_id"]].append(row)

    if selected is not None:
        missing = selected - set(grouped)
        if missing:
            raise ValueError(f"unknown --only-groups: {', '.join(sorted(missing))}")
    if not grouped:
        raise ValueError("no groups selected")

    for group_id, rows in grouped.items():
        for field in ("species", "genotype", "sex"):
            values = {row[field] for row in rows}
            if len(values) != 1:
                raise ValueError(f"{group_id}: inconsistent {field}: {sorted(values)}")
        reps = [replicate_number(row["replicate"]) for row in rows]
        if len(reps) != len(set(reps)):
            raise ValueError(f"{group_id}: duplicate replicate labels")
        rows.sort(key=lambda row: replicate_number(row["replicate"]))
    return dict(grouped)


def read_and_validate_header(handle, input_path: Path) -> None:
    header = handle.readline().rstrip(b"\r\n").split(b"\t")
    expected = [column.encode() for column in KEY_COLUMNS + VALUE_COLUMNS]
    if header != expected:
        decoded = [item.decode(errors="replace") for item in header]
        raise ValueError(f"{input_path}: unexpected header: {decoded}")


def merge_one(task: dict) -> dict[str, str | int]:
    group_id = task["group_id"]
    flank = task["flank"]
    rows = task["rows"]
    output_path = Path(task["output_path"])
    overwrite = task["overwrite"]
    progress_every = task["progress_every"]

    if output_path.exists() and not overwrite:
        raise FileExistsError(f"{output_path} exists; use --overwrite to replace it")

    input_paths = [
        Path(row["input_root"]) / f'{row["sample_id"]}.flank{flank}.tsv' for row in rows
    ]
    missing = [str(path) for path in input_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing input files: " + ", ".join(missing))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_name(f".{output_path.name}.{os.getpid()}.tmp")
    handles = []
    row_count = 0
    start = time.monotonic()
    try:
        handles = [path.open("rb", buffering=1024 * 1024) for path in input_paths]
        for handle, path in zip(handles, input_paths):
            read_and_validate_header(handle, path)

        output_header = [column.encode() for column in KEY_COLUMNS] + [b"site_id"]
        for row in rows:
            prefix = row["replicate"].encode()
            output_header.extend(
                prefix + b"__" + column.encode() for column in VALUE_COLUMNS
            )

        with temp_path.open("wb", buffering=1024 * 1024) as output:
            output.write(b"\t".join(output_header) + b"\n")
            while True:
                lines = [handle.readline() for handle in handles]
                if not lines[0]:
                    if any(lines[1:]):
                        raise ValueError(
                            f"{group_id} flank{flank}: first input ended early"
                        )
                    break
                if any(not line for line in lines[1:]):
                    ended = [
                        str(path) for path, line in zip(input_paths, lines) if not line
                    ]
                    raise ValueError(
                        f"{group_id} flank{flank}: input ended early: {', '.join(ended)}"
                    )

                fields = [line.rstrip(b"\r\n").split(b"\t") for line in lines]
                if any(len(parts) != 14 for parts in fields):
                    bad = [
                        f"{path}:{row_count + 2} ({len(parts)} columns)"
                        for path, parts in zip(input_paths, fields)
                        if len(parts) != 14
                    ]
                    raise ValueError("malformed rows: " + ", ".join(bad))
                key = fields[0][:7]
                for path, parts in zip(input_paths[1:], fields[1:]):
                    if parts[:7] != key:
                        raise ValueError(
                            f"{group_id} flank{flank}: codon mismatch at row "
                            f"{row_count + 2} in {path}"
                        )

                site_id = key[2] + b"_" + key[6]
                values = [value for parts in fields for value in parts[7:14]]
                output.write(b"\t".join(key + [site_id] + values) + b"\n")
                row_count += 1
                if progress_every and row_count % progress_every == 0:
                    elapsed = time.monotonic() - start
                    print(
                        f"[{group_id} flank{flank}] {row_count:,} rows "
                        f"({elapsed:.1f}s)",
                        file=sys.stderr,
                        flush=True,
                    )
        os.replace(temp_path, output_path)
    except BaseException:
        if temp_path.exists():
            temp_path.unlink()
        raise
    finally:
        for handle in handles:
            handle.close()

    elapsed = time.monotonic() - start
    print(
        f"[done] {group_id} flank{flank}: {row_count:,} rows in {elapsed:.1f}s",
        file=sys.stderr,
        flush=True,
    )
    first = rows[0]
    return {
        "group_id": group_id,
        "species": first["species"],
        "genotype": first["genotype"],
        "sex": first["sex"],
        "flank_size": flank,
        "n_replicates": len(rows),
        "replicates": ",".join(row["replicate"] for row in rows),
        "sample_ids": ",".join(row["sample_id"] for row in rows),
        "head_replicates": ",".join(
            row["replicate"] for row in rows if row["is_head"] == "1"
        ),
        "rows": row_count,
        "output_path": str(output_path.resolve()),
    }


def write_summary(path: Path, results: list[dict[str, str | int]]) -> None:
    columns = (
        "group_id",
        "species",
        "genotype",
        "sex",
        "flank_size",
        "n_replicates",
        "replicates",
        "sample_ids",
        "head_replicates",
        "rows",
        "output_path",
    )
    temp_path = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temp_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=columns, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(
            sorted(
                results, key=lambda row: (str(row["group_id"]), int(row["flank_size"]))
            )
        )
    os.replace(temp_path, path)


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise ValueError("--jobs must be >= 1")
    if any(flank < 1 for flank in args.flanks):
        raise ValueError("--flanks values must be positive")

    group_path = Path(args.groups).resolve()
    selected = set(args.only_groups) if args.only_groups else None
    grouped = load_groups(group_path, selected)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for group_id, rows in sorted(grouped.items()):
        for flank in sorted(set(args.flanks)):
            tasks.append(
                {
                    "group_id": group_id,
                    "flank": flank,
                    "rows": rows,
                    "output_path": str(
                        output_dir / f"{group_id}.flank{flank}.master.tsv"
                    ),
                    "overwrite": args.overwrite,
                    "progress_every": args.progress_every,
                }
            )

    results = []
    with ProcessPoolExecutor(max_workers=min(args.jobs, len(tasks))) as executor:
        futures = [executor.submit(merge_one, task) for task in tasks]
        for future in as_completed(futures):
            results.append(future.result())

    summary_path = output_dir / "master_table_manifest.tsv"
    write_summary(summary_path, results)
    print(f"Wrote {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
