#!/usr/bin/env python3
"""Convert aggregated RNA/Ribo tables to headerless transcript, CDS start, codon, count.
No BAM counting, rounding, allele selection or count aggregation is performed.
"""

import argparse
import csv
import json
import os
import re
from decimal import Decimal
from pathlib import Path


def convert(source, target):
    source, target = Path(source), Path(target)
    if source.resolve() == target.resolve():
        raise ValueError("Input and output must differ")
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_name(target.name + ".partial")
    kept = skipped = 0
    try:
        with source.open() as fi, tmp.open("w") as fo:
            reader = csv.reader(fi, delimiter="\t")
            first = next(reader, None)
            if first is None:
                raise ValueError(f"Empty count file: {source}")
            first[0] = first[0].lstrip("\ufeff")
            if "codon_start_cds" in first and "observed_codon" in first:
                cols = [
                    first.index(n)
                    for n in [
                        "transcript",
                        "codon_start_cds",
                        "observed_codon",
                        "count",
                    ]
                ]
                end_col = first.index("codon_end_cds")
                data = reader
            elif first[:1] == ["transcript_id"] and "count" in first:
                cols = [
                    first.index(n) for n in ["transcript_id", "start", "codon", "count"]
                ]
                end_col = first.index("end")
                data = reader
            else:
                if len(first) != 4 or not re.fullmatch(r"[+-]?\d+", first[1]):
                    raise ValueError(f"Unsupported count format: {source}: {first}")
                import itertools

                cols, end_col = [0, 1, 2, 3], None
                data = itertools.chain([first], reader)
            for row in data:
                if not row:
                    continue
                tx, start, codon, count = [row[i] for i in cols]
                end = row[end_col] if end_col is not None else start
                # +N explicitly denotes 3-prime UTR in 01 output, not CDS +N.
                if (
                    start.startswith("+")
                    or end.startswith("+")
                    or int(start) < 1
                    or int(end) < 1
                ):
                    skipped += 1
                    continue
                value = Decimal(count)
                if not value.is_finite() or value < 0 or len(codon) != 3:
                    raise ValueError(f"Invalid count row in {source}: {row}")
                fo.write(f"{tx}\t{int(start)}\t{codon}\t{count}\n")
                kept += 1
        os.replace(tmp, target)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return {
        "source": str(source),
        "output": str(target),
        "kept_rows": kept,
        "skipped_UTR_rows": skipped,
    }


def first_existing(paths):
    hits = [p for p in paths if p.is_file()]
    if not hits:
        raise FileNotFoundError(
            "No count file found; tried: " + ", ".join(map(str, paths))
        )
    # Prefer historical .num when both representations are present.
    # --manifest explicitly selects a different representation or output tree.
    return hits[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--manifest",
        help="TSV: group,rna_prefix,ribo_prefix,site_type,rna_file,ribo_file; absolute file paths",
    )
    ap.add_argument("--pair-dir", type=Path)
    ap.add_argument("--rna-dir", type=Path)
    ap.add_argument("--ribo-dir", type=Path)
    ap.add_argument(
        "--groups",
        nargs="+",
        default=["as", "wbg", "ag", "oldDmel", "sim", "sec", "ere"],
    )
    ap.add_argument(
        "--site-types",
        nargs="+",
        choices=["Asite", "Psite", "Esite"],
        default=["Asite", "Psite", "Esite"],
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--check-only", action="store_true")
    args = ap.parse_args()
    records = []
    if args.manifest:
        with open(args.manifest) as f:
            records = list(csv.DictReader(f, delimiter="\t"))
        required = {
            "group",
            "rna_prefix",
            "ribo_prefix",
            "site_type",
            "rna_file",
            "ribo_file",
        }
        if not records or not required.issubset(records[0]):
            ap.error("Empty manifest or missing columns: " + ",".join(sorted(required)))
        records = [
            r
            for r in records
            if r["group"] in args.groups and r["site_type"] in args.site_types
        ]
    else:
        if not all([args.pair_dir, args.rna_dir, args.ribo_dir]):
            ap.error("--pair-dir, --rna-dir and --ribo-dir required without --manifest")
        for group in args.groups:
            with (args.pair_dir / f"rna_psite_prefix_pair_{group}.list").open() as f:
                for line in f:
                    fields = line.split()
                    if (
                        not fields
                        or fields[0].startswith("#")
                        or fields[0] == "rna_prefix"
                    ):
                        continue
                    rna, ribo = fields[:2]
                    rf = first_existing(
                        [
                            args.rna_dir / f"{rna}_longestCDS.old4.tsv.num",
                            args.rna_dir / f"{rna}.longestCDS.codon_counts.tsv",
                            args.rna_dir / f"{rna}_longestCDS.codon_counts.tsv",
                        ]
                    )
                    for site in args.site_types:
                        bf = first_existing(
                            [
                                args.ribo_dir / f"{ribo}.codon_{site}.num",
                                args.ribo_dir / f"{ribo}.{site[0]}.tsv",
                            ]
                        )
                        records.append(
                            dict(
                                group=group,
                                rna_prefix=rna,
                                ribo_prefix=ribo,
                                site_type=site,
                                rna_file=str(rf.resolve()),
                                ribo_file=str(bf.resolve()),
                            )
                        )
    if not records:
        ap.error("No selected RNA/Ribo pairs")
    destinations = {}
    for r in records:
        if (
            r["group"] not in ["as", "wbg", "ag", "oldDmel", "sim", "sec", "ere"]
            or r["site_type"] not in args.site_types
        ):
            raise ValueError(f"Unknown group/site: {r}")
        for k in ["rna_prefix", "ribo_prefix"]:
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", r[k]):
                raise ValueError("Invalid sample prefix: " + r[k])
        for k in ["rna_file", "ribo_file"]:
            if not Path(r[k]).is_absolute() or not Path(r[k]).is_file():
                raise FileNotFoundError(r[k])
        r["rna_num"] = str((args.out_dir / "rna" / f"{r['rna_prefix']}.num").resolve())
        r["ribo_num"] = str(
            (
                args.out_dir / "ribo" / f"{r['ribo_prefix']}.{r['site_type']}.num"
            ).resolve()
        )
        for source, target in [
            (r["rna_file"], r["rna_num"]),
            (r["ribo_file"], r["ribo_num"]),
        ]:
            if (
                target in destinations
                and Path(destinations[target]).resolve() != Path(source).resolve()
            ):
                raise ValueError("Conflicting sources for " + target)
            destinations[target] = source
    keys = [(r["ribo_prefix"], r["site_type"]) for r in records]
    if len(set(keys)) != len(keys):
        raise ValueError("Repeated Ribo sample/site; specify each matched pair once")
    if args.check_only:
        print(
            f"Input paths checked: {len(records)} matched sample/site tasks; {len(destinations)} unique count files"
        )
        return
    stats = [convert(source, target) for target, source in destinations.items()]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "count_manifest.tsv").open("w") as f:
        cols = [
            "group",
            "rna_prefix",
            "ribo_prefix",
            "site_type",
            "rna_file",
            "ribo_file",
            "rna_num",
            "ribo_num",
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows({k: r[k] for k in cols} for r in records)
    (args.out_dir / "conversion_summary.json").write_text(
        json.dumps(stats, indent=2) + "\n"
    )
    print(f"Prepared {len(records)} tasks: {args.out_dir}/count_manifest.tsv")


if __name__ == "__main__":
    main()
