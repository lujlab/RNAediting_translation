#!/usr/bin/env python3
"""Build the stage-resolved embryo editing-site lists used by Duan et al. 2017.

The developmental measurements come from Graveley et al. 2011 Supplementary
Table 29.  A site is considered detected in a stage when the reported
``Edited %`` value is greater than zero, matching the per-stage counts plotted
by Duan et al. in Supplementary Figure 30.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

STAGES = [f"em{i}-{i + 2}hr" for i in range(0, 24, 2)]
EARLY_STAGES = set(STAGES[:8])
LATE_STAGES = set(STAGES[8:])


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def parse_gtf_attributes(text: str) -> dict[str, str]:
    result = {}
    for field in text.rstrip(";").split("; "):
        if " " not in field:
            continue
        key, value = field.split(" ", 1)
        result[key] = value.strip('"')
    return result


def load_gtf_intervals(path: Path):
    genes = defaultdict(list)
    exons = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "exon"}:
                continue
            attrs = parse_gtf_attributes(fields[8])
            item = (
                int(fields[3]),
                int(fields[4]),
                attrs.get("gene_id", ""),
                attrs.get("gene_symbol", ""),
                fields[6],
            )
            (genes if fields[2] == "gene" else exons)[fields[0]].append(item)
    return genes, exons


def overlaps(intervals, position: int):
    by_fbgn = {}
    for start, end, fbgn, symbol, strand in intervals:
        if start <= position <= end and fbgn:
            by_fbgn[fbgn] = (fbgn, symbol, strand)
    return list(by_fbgn.values())


def load_name_maps(current_symbols: Path, annotation_ids: Path):
    name_to_fbgn = defaultdict(set)
    fbgn_to_symbol = {}
    for row in read_tsv(current_symbols):
        name_to_fbgn[row["gene_symbol"]].add(row["FBgn"])
        fbgn_to_symbol[row["FBgn"]] = row["gene_symbol"]
    for row in read_tsv(annotation_ids):
        name_to_fbgn[row["annotation_id"]].add(row["FBgn"])
    # Historical names absent from the r6.12 annotation-ID table.
    name_to_fbgn["det"].add("FBgn0260003")  # detached, current symbol Dys
    name_to_fbgn["nAcRbeta-21C"].add("FBgn0031261")  # current nAChRbeta3
    return name_to_fbgn, fbgn_to_symbol


def choose_host(old_name, exon_hits, gene_hits, name_to_fbgn):
    exon_by_id = {x[0]: x for x in exon_hits}
    gene_by_id = {x[0]: x for x in gene_hits}
    mapped = name_to_fbgn.get(old_name, set())
    mapped_exon = [exon_by_id[x] for x in mapped if x in exon_by_id]
    if len(mapped_exon) == 1:
        return mapped_exon[0], "historical_name_to_FBgn_and_exon"
    exact_exon = [x for x in exon_hits if x[1] == old_name]
    if len(exact_exon) == 1:
        return exact_exon[0], "exact_symbol_and_exon"
    if len(exon_hits) == 1:
        return exon_hits[0], "unique_exon_overlap"
    mapped_gene = [gene_by_id[x] for x in mapped if x in gene_by_id]
    if len(mapped_gene) == 1:
        return mapped_gene[0], "historical_name_to_FBgn_and_gene"
    exact_gene = [x for x in gene_hits if x[1] == old_name]
    if len(exact_gene) == 1:
        return exact_gene[0], "exact_symbol_and_gene"
    if len(gene_hits) == 1:
        return gene_hits[0], "unique_gene_overlap"
    raise ValueError(
        f"Could not uniquely map host gene {old_name!r}; "
        f"exon_hits={exon_hits}, gene_hits={gene_hits}, mapped={sorted(mapped)}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--coordinate-map", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--current-symbols", type=Path, required=True)
    parser.add_argument("--annotation-ids", type=Path, required=True)
    parser.add_argument("--wide-output", type=Path, required=True)
    parser.add_argument("--long-output", type=Path, required=True)
    parser.add_argument("--count-output", type=Path, required=True)
    args = parser.parse_args()

    coords = {int(r["source_row"]): r for r in read_tsv(args.coordinate_map)}
    genes, exons = load_gtf_intervals(args.gtf)
    name_to_fbgn, fbgn_to_symbol = load_name_maps(
        args.current_symbols, args.annotation_ids
    )

    with args.source.open(newline="") as handle:
        source_rows = list(csv.DictReader(handle))

    wide_rows = []
    long_rows = []
    stage_counts = {stage: 0 for stage in STAGES}
    mapping_method_counts = defaultdict(int)

    for source_row, row in enumerate(source_rows, 1):
        positive_stages = [
            stage for stage in STAGES if float(row[f"{stage} Edited %"] or 0) > 0
        ]
        if not positive_stages:
            continue
        coord = coords[source_row]
        chrom = coord["chr_r6"]
        pos = int(coord["pos_r6_1based"])
        exon_hits = overlaps(exons[chrom], pos)
        gene_hits = overlaps(genes[chrom], pos)
        host, mapping_method = choose_host(
            row["Gene"], exon_hits, gene_hits, name_to_fbgn
        )
        fbgn, current_symbol, strand = host
        current_symbol = fbgn_to_symbol.get(fbgn, current_symbol)
        mapping_method_counts[mapping_method] += 1

        fractions = [float(row[f"{stage} Edited %"] or 0) for stage in STAGES]
        base = {
            "source_row": source_row,
            "gene_symbol_2011": row["Gene"],
            "FBgn_r6_12": fbgn,
            "gene_symbol_r6_12": current_symbol,
            "host_strand_r6_12": strand,
            "host_mapping_method": mapping_method,
            "chr_r5": coord["chr_r5"],
            "pos_r5_1based": coord["pos_r5_1based"],
            "chr_r6": chrom,
            "pos_r6_1based": pos,
            "functional_annotation_2011": row["Substitution"],
            "known_in_2011": row["Known?"],
            "detected_early_embryo_0_16h": int(
                any(x in EARLY_STAGES for x in positive_stages)
            ),
            "detected_late_embryo_16_24h": int(
                any(x in LATE_STAGES for x in positive_stages)
            ),
            "detected_embryo_0_24h": 1,
            "n_positive_embryo_stages": len(positive_stages),
            "max_embryo_editing_fraction": max(fractions),
        }
        for stage in STAGES:
            depth = int(float(row[stage] or 0))
            fraction = float(row[f"{stage} Edited %"] or 0)
            base[f"{stage}_depth"] = depth
            base[f"{stage}_editing_fraction"] = fraction
            if fraction <= 0:
                continue
            stage_counts[stage] += 1
            edited = int(round(depth * fraction))
            long_rows.append(
                {
                    "gene": fbgn,
                    "gene_symbol_2011": row["Gene"],
                    "gene_symbol_r6_12": current_symbol,
                    "dmel_chr": chrom,
                    "dmel_pos": pos,
                    "lib": stage,
                    "depth": depth,
                    "altcount": edited,
                    "refcount": depth - edited,
                    "ratio": fraction,
                    "functional_annotation_2011": row["Substitution"],
                    "source_row": source_row,
                }
            )
        wide_rows.append(base)

    wide_fields = list(wide_rows[0])
    long_fields = list(long_rows[0])
    args.wide_output.parent.mkdir(parents=True, exist_ok=True)
    with args.wide_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=wide_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(wide_rows)
    with args.long_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=long_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(long_rows)

    count_rows = [
        {"scope": stage, "site_count": stage_counts[stage]} for stage in STAGES
    ]
    count_rows.extend(
        [
            {
                "scope": "early_embryo_0_16h_union",
                "site_count": sum(
                    int(r["detected_early_embryo_0_16h"]) for r in wide_rows
                ),
            },
            {
                "scope": "late_embryo_16_24h_union",
                "site_count": sum(
                    int(r["detected_late_embryo_16_24h"]) for r in wide_rows
                ),
            },
            {"scope": "all_embryo_0_24h_union", "site_count": len(wide_rows)},
            {
                "scope": "all_embryo_0_24h_unique_genes",
                "site_count": len({r["FBgn_r6_12"] for r in wide_rows}),
            },
        ]
    )
    with args.count_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["scope", "site_count"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(count_rows)

    print(f"Wrote {len(wide_rows)} unique embryo sites to {args.wide_output}")
    print(f"Wrote {len(long_rows)} positive site-stage rows to {args.long_output}")
    print(f"Unique r6.12 FBgn genes: {len({r['FBgn_r6_12'] for r in wide_rows})}")
    print("Host mapping methods:")
    for method, count in sorted(mapping_method_counts.items()):
        print(f"  {method}: {count}")


if __name__ == "__main__":
    main()
