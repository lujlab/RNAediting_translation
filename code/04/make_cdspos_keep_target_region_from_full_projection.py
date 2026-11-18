#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv
from pathlib import Path

HEADER = [
    "genomic_key",
    "base",
    "variant_effect",
    "source_transcript",
    "source_tx_pos_1based",
    "source_tx_length",
    "source_cds_pos_1based",
    "source_cds_length",
    "source_aa_pos_1based",
    "source_aa_length",
    "site_key",
    "new_transcript",
    "new_tx_pos_1based",
    "target_region",
    "new_cds_pos_1based",
]
CODING = {
    "synonymous_variant",
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "stop_retained_variant",
    "missense_variant&splice_region_variant",
    "splice_region_variant&synonymous_variant",
}


def clean(x):
    return "" if x is None else str(x).strip()


def split_pair(x):
    x = clean(x)
    if "/" not in x:
        return "NA", "NA"
    a, b = x.split("/", 1)
    return (a if a else "NA"), (b if b else "NA")


def has_nums(d):
    return (
        d["source_transcript"].startswith("FBtr")
        and d["source_tx_pos_1based"].isdigit()
        and d["source_cds_pos_1based"].isdigit()
        and d["source_cds_length"].isdigit()
        and d["source_aa_pos_1based"].isdigit()
        and d["source_aa_length"].isdigit()
    )


def parse_record(rec, row):
    f = rec.split("|") if rec else []
    f += [""] * max(0, 16 - len(f))
    txp, txl = split_pair(f[11])
    cdsp, cdsl = split_pair(f[12])
    aap, aal = split_pair(f[13])
    return {
        "variant_effect": clean(f[1] or row.get("type") or row.get("variant_effect")),
        "source_transcript": clean(
            f[6] or row.get("transcript_id") or row.get("transcript")
        )
        or "NA",
        "source_tx_pos_1based": txp,
        "source_tx_length": txl,
        "source_cds_pos_1based": cdsp,
        "source_cds_length": cdsl,
        "source_aa_pos_1based": aap,
        "source_aa_length": aal,
    }


def parse_source(row):
    row_type = clean(row.get("type") or row.get("variant_effect"))
    recs = [r for r in clean(row.get("ann")).split(",") if r]
    parsed = [parse_record(r, row) for r in recs]
    for d in parsed:
        if (
            d["variant_effect"] == row_type
            and d["variant_effect"] in CODING
            and has_nums(d)
        ):
            return d
    for d in parsed:
        if d["variant_effect"] in CODING and has_nums(d):
            return d
    return None


def coord_key(c):
    try:
        a, b = c.rsplit("_", 1)
        return (a, int(b), c)
    except Exception:
        return (c, 0, c)


def dedup_key(r):
    if r["target_region"] == "CDS" and r["new_cds_pos_1based"].isdigit():
        return (r["new_transcript"], "CDS", r["new_cds_pos_1based"])
    return (r["new_transcript"], r["target_region"], r["new_tx_pos_1based"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-full", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()
    rows = []
    skipped = 0
    unmapped = 0
    with open(a.input_full, encoding="utf-8") as f:
        for raw in csv.DictReader(f, delimiter="\t"):
            row = {k: clean(v) for k, v in raw.items()}
            src = parse_source(row)
            if (
                not row.get("new_transcript")
                or row.get("new_transcript") in {"NA", "."}
                or not row.get("new_tx_pos_1based", "").isdigit()
                or int(row["new_tx_pos_1based"]) < 1
            ):
                unmapped += 1
                continue
            if src is None:
                skipped += 1
                continue
            genomic_key = (
                clean(row.get("species_coord"))
                or clean(row.get("handle"))
                or (clean(row.get("chr")) + "_" + clean(row.get("pos")))
            )
            rows.append(
                {
                    "genomic_key": genomic_key,
                    "base": clean(row.get("Alt") or row.get("base")) or "NA",
                    **src,
                    "site_key": clean(row.get("handle") or row.get("site_key"))
                    or genomic_key,
                    "new_transcript": clean(row.get("new_transcript")) or "NA",
                    "new_tx_pos_1based": clean(row.get("new_tx_pos_1based")) or "NA",
                    "target_region": clean(row.get("target_region")) or "NA",
                    "new_cds_pos_1based": clean(row.get("new_cds_pos_1based")) or "NA",
                }
            )
    rows.sort(
        key=lambda r: (
            coord_key(r["genomic_key"]),
            r["site_key"],
            r["new_transcript"],
            r["new_tx_pos_1based"],
        )
    )
    seen = set()
    kept = []
    for r in rows:
        k = dedup_key(r)
        if k in seen:
            continue
        seen.add(k)
        kept.append(r)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(kept)
    counts = {}
    for r in kept:
        counts[r["target_region"]] = counts.get(r["target_region"], 0) + 1
    with open(a.summary, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"input_full_rows\t{len(rows)+skipped+unmapped}\n")
        f.write(f"parsed_source_cds_rows\t{len(rows)}\n")
        f.write(f"skipped_no_source_cds\t{skipped}\n")
        f.write(f"skipped_unmapped_rows\t{unmapped}\n")
        f.write(f"output_rows\t{len(kept)}\n")
        for k in sorted(counts):
            f.write(f"target_region_{k}\t{counts[k]}\n")


if __name__ == "__main__":
    main()
