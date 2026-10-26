#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Match MS/Lp sites to the NEW min3 Lm/Lr summary table.

This version assumes --lmlr is the newly recalculated table:
    all_species_LmLr_min3.sex_site_mean_LmLr.tsv

That table already contains all cohorts/genotypes:
    Dmel_ISO1, Dmel_WBg, Dmel_AG, Dmel_AS, Dsim, Dsec, Dere

Therefore this script:
  - DOES NOT recalculate AG/AS
  - DOES NOT use the old filtered Ipos table
  - DOES NOT require --alledit-dir or --ag-as-out

Workflow
--------
1. Read the MS Lp wide table.
2. Map MS source_site_ids -> site_key in the species cdspos table.
3. Obtain new_transcript + new_cds_pos_1based for every MS editing site.
4. Read the new Lm/Lr table. Its `handle` is already:
       new_transcript_new_cds_pos
5. Match by:
       new_transcript + new_cds_pos
       + genotype
       + sex
6. Append:
       mean_Lm_F mean_Lr_F
       mean_Lm_M mean_Lr_M

If an AA-level MS site corresponds to multiple nucleotide editing positions,
all matching positions are queried. mean_Lm and mean_Lr are averaged
INDEPENDENTLY across the matched nucleotide positions:
    - non-NA mean_Lm values are averaged
    - non-NA mean_Lr values are averaged
Thus Lm-only sites remain Lm-only and are not discarded.

Output also retains QC columns showing exactly which handles were matched.
"""

import argparse
import csv
import glob
import math
import os
import re
import statistics
from collections import defaultdict

# =============================================================================
# Helpers
# =============================================================================


def clean(x):
    if x is None:
        return ""
    return str(x).strip().lstrip("\ufeff").strip('"').strip("'")


def norm_fbtr(x):
    return re.sub(r"\.\d+$", "", clean(x))


def norm_int(x):
    x = clean(x)
    if not x:
        return ""
    try:
        return str(int(float(x)))
    except Exception:
        return x


def parse_float_or_none(x):
    x = clean(x)
    if x == "" or x.upper() in {"NA", "NAN", "NULL"}:
        return None
    try:
        v = float(x)
    except ValueError:
        raise ValueError(f"Non-numeric value where numeric expected: {x!r}")
    return v if math.isfinite(v) else None


def fmt(x):
    return "NA" if x is None else f"{x:.12g}"


def mean_non_na(vals):
    vals = [x for x in vals if x is not None and math.isfinite(x)]
    return statistics.mean(vals) if vals else None


def split_multi(x):
    x = clean(x)
    if not x:
        return []
    vals = [z.strip() for z in re.split(r"[;,|]+", x) if z.strip()]
    out, seen = [], set()
    for z in vals:
        if z not in seen:
            out.append(z)
            seen.add(z)
    return out


def parse_handle(x):
    x = clean(x)
    m = re.match(r"^(FBtr\d+)(?:\.\d+)?_(\d+)$", x)
    if not m:
        return None, None
    return norm_fbtr(m.group(1)), str(int(m.group(2)))


def canonical(x):
    return re.sub(r"[^A-Za-z0-9]+", "", clean(x)).lower()


def normalize_genotype(x):
    z = canonical(x)

    if z in {"iso", "iso1", "dmeliso1"}:
        return "iso"
    if z in {"wbg", "dmelwbg"}:
        return "WBg"
    if z in {"ag", "adarg", "dmelag"}:
        return "AG"
    if z in {"as", "adars", "dmelas"}:
        return "AS"
    if z in {"sim", "dsim", "dsimulans"}:
        return "sim"
    if z in {"sec", "dsec", "dsechellia"}:
        return "sec"
    if z in {"ere", "dere", "derecta"}:
        return "ere"

    return clean(x)


def genotype_to_species(genotype):
    g = normalize_genotype(genotype)

    if g in {"iso", "WBg", "AG", "AS"}:
        return "Dmel"
    if g == "sim":
        return "Dsim"
    if g == "sec":
        return "Dsec"
    if g == "ere":
        return "Dere"

    return None


def normalize_sex(x):
    z = clean(x).lower()
    if z in {"f", "female"}:
        return "F"
    if z in {"m", "male"}:
        return "M"
    return None


# =============================================================================
# Coordinate maps
# =============================================================================


def discover_coordinate_files(cdspos_dir):
    patterns = {
        "Dmel": "ISO1*.cdspos_for_alledit.tsv",
        "Dsim": "sim*.cdspos_for_alledit.tsv",
        "Dsec": "sec*.cdspos_for_alledit.tsv",
        "Dere": "ere*.cdspos_for_alledit.tsv",
    }

    out = {}
    for species, pattern in patterns.items():
        hits = sorted(glob.glob(os.path.join(cdspos_dir, pattern)))

        if len(hits) == 0:
            raise FileNotFoundError(
                f"No coordinate table for {species}: "
                f"{os.path.join(cdspos_dir, pattern)}"
            )

        if len(hits) > 1:
            raise RuntimeError(
                f"Multiple coordinate tables found for {species}:\n" + "\n".join(hits)
            )

        out[species] = hits[0]

    return out


def load_coordinate_index(files_by_species):
    """
    Index:
        (species, site_key) -> coordinate records
    """
    index = defaultdict(list)
    stats = {}

    required = [
        "site_key",
        "new_transcript",
        "new_tx_pos_1based",
        "new_cds_pos_1based",
    ]

    for species, path in files_by_species.items():
        n_rows = 0
        n_indexed = 0
        seen = set()

        with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"No header in {path}")

            reader.fieldnames = [clean(x) for x in reader.fieldnames]

            missing = [x for x in required if x not in reader.fieldnames]
            if missing:
                raise ValueError(f"{path} missing columns: {','.join(missing)}")

            for raw in reader:
                n_rows += 1
                row = {clean(k): clean(v) for k, v in raw.items()}

                site_key = clean(row["site_key"])
                new_tx = norm_fbtr(row["new_transcript"])
                new_tx_pos = norm_int(row["new_tx_pos_1based"])
                new_cds = norm_int(row["new_cds_pos_1based"])

                if not site_key or not new_tx or not new_cds:
                    continue

                rec = {
                    "species": species,
                    "site_key": site_key,
                    "new_transcript": new_tx,
                    "new_tx_pos_1based": new_tx_pos,
                    "new_cds_pos_1based": new_cds,
                }

                key = (
                    species,
                    site_key,
                    new_tx,
                    new_tx_pos,
                    new_cds,
                )

                if key in seen:
                    continue
                seen.add(key)

                index[(species, site_key)].append(rec)
                n_indexed += 1

        stats[species] = {
            "file": path,
            "rows": n_rows,
            "indexed": n_indexed,
        }

    return index, stats


# =============================================================================
# Read Lp wide
# =============================================================================


def read_lp_wide(path):
    rows = []

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"No header in Lp wide table: {path}")

        fields = [clean(x) for x in reader.fieldnames]
        reader.fieldnames = fields

        required = ["site", "source_site_ids"]
        missing = [x for x in required if x not in fields]

        if missing:
            raise ValueError("Lp wide table missing columns: " + ",".join(missing))

        has_genotype = "genotype" in fields
        has_group = "group" in fields

        base = os.path.basename(path).lower()
        filename_genotype = None

        if not has_genotype and not has_group:
            if "dsim" in base:
                filename_genotype = "sim"
            elif "dsec" in base:
                filename_genotype = "sec"
            elif "dere" in base:
                filename_genotype = "ere"
            else:
                raise ValueError(
                    "Lp file has neither genotype nor group and genotype "
                    f"cannot be inferred from filename: {path}"
                )

        for raw in reader:
            row = {clean(k): clean(v) for k, v in raw.items()}

            if has_genotype:
                genotype = row.get("genotype", "")
            elif has_group:
                genotype = row.get("group", "")
            else:
                genotype = filename_genotype

            row["genotype"] = normalize_genotype(genotype)
            rows.append(row)

    if "genotype" not in fields:
        if "group" in fields:
            i = fields.index("group") + 1
            fields = fields[:i] + ["genotype"] + fields[i:]
        else:
            fields.append("genotype")

    return rows, fields


# =============================================================================
# Read NEW min3 Lm/Lr table
# =============================================================================


def read_lmlr_table(path):
    """
    Index:
        (transcript, cdspos, genotype, sex) -> records
    """
    index = defaultdict(list)
    stats = {
        "rows_total": 0,
        "rows_indexed": 0,
        "rows_bad_handle": 0,
        "rows_bad_sex": 0,
        "rows_mean_Lm_nonNA": 0,
        "rows_mean_Lr_nonNA": 0,
        "rows_Lm_only": 0,
    }

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"No header in Lm/Lr table: {path}")

        reader.fieldnames = [clean(x) for x in reader.fieldnames]

        required = ["species", "cohort", "sex", "handle", "mean_Lm", "mean_Lr"]

        missing = [x for x in required if x not in reader.fieldnames]
        if missing:
            raise ValueError("Lm/Lr table missing columns: " + ",".join(missing))

        has_genotype = "genotype" in reader.fieldnames
        has_sex_short = "sex_short" in reader.fieldnames

        for raw in reader:
            stats["rows_total"] += 1
            row = {clean(k): clean(v) for k, v in raw.items()}

            tx, cdspos = parse_handle(row["handle"])
            if not tx:
                stats["rows_bad_handle"] += 1
                continue

            sex = (
                normalize_sex(row["sex_short"])
                if has_sex_short and clean(row.get("sex_short", ""))
                else normalize_sex(row["sex"])
            )

            if sex is None:
                stats["rows_bad_sex"] += 1
                continue

            if has_genotype and clean(row.get("genotype", "")):
                genotype = normalize_genotype(row["genotype"])
            else:
                genotype = normalize_genotype(row["cohort"])

            lm = parse_float_or_none(row["mean_Lm"])
            lr = parse_float_or_none(row["mean_Lr"])

            if lm is not None:
                stats["rows_mean_Lm_nonNA"] += 1
            if lr is not None:
                stats["rows_mean_Lr_nonNA"] += 1
            if lm is not None and lr is None:
                stats["rows_Lm_only"] += 1

            rec = {
                "species": clean(row["species"]),
                "cohort": clean(row["cohort"]),
                "genotype": genotype,
                "sex": sex,
                "handle": clean(row["handle"]),
                "new_transcript": tx,
                "new_cds_pos_1based": cdspos,
                "mean_Lm": lm,
                "mean_Lr": lr,
                "sex_mean_status": clean(row.get("sex_mean_status", "")),
                "n_Lm_replicates": clean(row.get("n_Lm_replicates", "")),
                "n_Lr_replicates": clean(row.get("n_Lr_replicates", "")),
            }

            index[(tx, cdspos, genotype, sex)].append(rec)
            stats["rows_indexed"] += 1

    return index, stats


# =============================================================================
# MS coordinate remapping and final matching
# =============================================================================


def remap_lp_row(row, coord_index):
    genotype = normalize_genotype(row["genotype"])
    species = genotype_to_species(genotype)
    site_keys = split_multi(row.get("source_site_ids", ""))

    mapped = []

    if species:
        for site_key in site_keys:
            mapped.extend(coord_index.get((species, site_key), []))

    # de-duplicate exact coordinate records
    out, seen = [], set()
    for rec in mapped:
        k = (
            rec["site_key"],
            rec["new_transcript"],
            rec["new_tx_pos_1based"],
            rec["new_cds_pos_1based"],
        )
        if k not in seen:
            out.append(rec)
            seen.add(k)

    return species, site_keys, out


def write_output(
    lp_rows,
    lp_fields,
    coord_index,
    lmlr_index,
    out_path,
):
    qc_path = out_path + ".match_qc.tsv"
    unmapped_path = out_path + ".coordinate_unmatched.tsv"

    add_fields = [
        "coordinate_species",
        "mapped_source_site_ids",
        "new_transcript",
        "new_tx_pos_1based",
        "new_cds_pos_1based",
        "new_handle",
        "mean_Lm_F",
        "mean_Lr_F",
        "mean_Lm_M",
        "mean_Lr_M",
        "n_LmLr_matches_F",
        "n_LmLr_matches_M",
        "matched_LmLr_handles_F",
        "matched_LmLr_handles_M",
        "matched_LmLr_cohorts_F",
        "matched_LmLr_cohorts_M",
        "LmLr_status_F",
        "LmLr_status_M",
    ]

    out_fields = list(lp_fields)
    for x in add_fields:
        if x not in out_fields:
            out_fields.append(x)

    qc_fields = [
        "site",
        "genotype",
        "sex",
        "source_site_ids",
        "new_handle",
        "n_matches",
        "matched_handles",
        "matched_cohorts",
        "mean_Lm",
        "mean_Lr",
        "status",
    ]

    unmapped_fields = [
        "site",
        "genotype",
        "FBtr",
        "source_site_ids",
        "source_cds_pos",
        "reason",
    ]

    stats = {
        "lp_rows": len(lp_rows),
        "coordinate_mapped": 0,
        "coordinate_unmapped": 0,
        "rows_any_lmlr": 0,
        "row_sex_matches": 0,
        "row_sex_Lm_only": 0,
    }

    with open(out_path, "w", encoding="utf-8", newline="") as fo, open(
        qc_path, "w", encoding="utf-8", newline=""
    ) as fq, open(unmapped_path, "w", encoding="utf-8", newline="") as fu:

        w = csv.DictWriter(fo, delimiter="\t", fieldnames=out_fields)
        qw = csv.DictWriter(fq, delimiter="\t", fieldnames=qc_fields)
        uw = csv.DictWriter(fu, delimiter="\t", fieldnames=unmapped_fields)

        w.writeheader()
        qw.writeheader()
        uw.writeheader()

        for row in lp_rows:
            genotype = normalize_genotype(row["genotype"])

            species, site_keys, mapped = remap_lp_row(row, coord_index)

            row["coordinate_species"] = species or "NA"

            if mapped:
                stats["coordinate_mapped"] += 1

                row["mapped_source_site_ids"] = ";".join(x["site_key"] for x in mapped)
                row["new_transcript"] = ";".join(x["new_transcript"] for x in mapped)
                row["new_tx_pos_1based"] = ";".join(
                    x["new_tx_pos_1based"] for x in mapped
                )
                row["new_cds_pos_1based"] = ";".join(
                    x["new_cds_pos_1based"] for x in mapped
                )

                handles = []
                seen_h = set()
                for x in mapped:
                    h = f"{x['new_transcript']}_{x['new_cds_pos_1based']}"
                    if h not in seen_h:
                        handles.append(h)
                        seen_h.add(h)

                row["new_handle"] = ";".join(handles)

            else:
                stats["coordinate_unmapped"] += 1
                row["mapped_source_site_ids"] = ""
                row["new_transcript"] = ""
                row["new_tx_pos_1based"] = ""
                row["new_cds_pos_1based"] = ""
                row["new_handle"] = ""

                reason = (
                    "cannot_resolve_species"
                    if species is None
                    else "no_source_site_ids" if not site_keys else "site_key_not_found"
                )

                uw.writerow(
                    {
                        "site": row.get("site", ""),
                        "genotype": genotype,
                        "FBtr": row.get("FBtr", ""),
                        "source_site_ids": row.get("source_site_ids", ""),
                        "source_cds_pos": row.get("source_cds_pos", ""),
                        "reason": reason,
                    }
                )

            coords = []
            seen_c = set()

            for x in mapped:
                k = (x["new_transcript"], x["new_cds_pos_1based"])
                if k not in seen_c:
                    coords.append(k)
                    seen_c.add(k)

            row_any = False

            for sex in ("F", "M"):
                matches = []

                for tx, cdspos in coords:
                    matches.extend(lmlr_index.get((tx, cdspos, genotype, sex), []))

                lm = mean_non_na([x["mean_Lm"] for x in matches])
                lr = mean_non_na([x["mean_Lr"] for x in matches])

                handles = sorted({x["handle"] for x in matches})
                cohorts = sorted({x["cohort"] for x in matches})

                if lm is not None and lr is not None:
                    status = "Lm_and_Lr"
                elif lm is not None:
                    status = "Lm_only"
                elif lr is not None:
                    status = "Lr_only"
                else:
                    status = "neither"

                row[f"mean_Lm_{sex}"] = fmt(lm)
                row[f"mean_Lr_{sex}"] = fmt(lr)
                row[f"n_LmLr_matches_{sex}"] = len(matches)
                row[f"matched_LmLr_handles_{sex}"] = ";".join(handles)
                row[f"matched_LmLr_cohorts_{sex}"] = ";".join(cohorts)
                row[f"LmLr_status_{sex}"] = status

                if matches:
                    row_any = True
                    stats["row_sex_matches"] += 1

                if status == "Lm_only":
                    stats["row_sex_Lm_only"] += 1

                qw.writerow(
                    {
                        "site": row.get("site", ""),
                        "genotype": genotype,
                        "sex": sex,
                        "source_site_ids": row.get("source_site_ids", ""),
                        "new_handle": row.get("new_handle", ""),
                        "n_matches": len(matches),
                        "matched_handles": ";".join(handles),
                        "matched_cohorts": ";".join(cohorts),
                        "mean_Lm": fmt(lm),
                        "mean_Lr": fmt(lr),
                        "status": status,
                    }
                )

            if row_any:
                stats["rows_any_lmlr"] += 1

            w.writerow({k: row.get(k, "") for k in out_fields})

    return stats, qc_path, unmapped_path


# =============================================================================
# Main
# =============================================================================


def main():
    p = argparse.ArgumentParser()

    p.add_argument("--lp-wide", required=True)
    p.add_argument("--lmlr", required=True)
    p.add_argument("--cdspos-dir", required=True)
    p.add_argument("--out", required=True)

    args = p.parse_args()

    for x in [args.lp_wide, args.lmlr]:
        if not os.path.isfile(x):
            raise FileNotFoundError(x)

    if not os.path.isdir(args.cdspos_dir):
        raise NotADirectoryError(args.cdspos_dir)

    lp_rows, lp_fields = read_lp_wide(args.lp_wide)

    coord_files = discover_coordinate_files(args.cdspos_dir)
    coord_index, coord_stats = load_coordinate_index(coord_files)

    lmlr_index, lmlr_stats = read_lmlr_table(args.lmlr)

    match_stats, qc_path, unmapped_path = write_output(
        lp_rows,
        lp_fields,
        coord_index,
        lmlr_index,
        args.out,
    )

    print("=" * 100)
    print("MS Lp + new min3 Lm/Lr matching finished")
    print("=" * 100)
    print("Lm/Lr source:")
    print(f"  {args.lmlr}")
    print("No AG/AS recalculation is performed in this script.")
    print("Lm and Lr are matched independently; Lm-only sites are retained.")
    print("-" * 100)

    for species in ("Dmel", "Dsim", "Dsec", "Dere"):
        s = coord_stats[species]
        print(
            f"{species:4s} coordinate map: "
            f"rows={s['rows']} indexed={s['indexed']} "
            f"| {s['file']}"
        )

    print("-" * 100)
    print(f"Lm/Lr rows total:              {lmlr_stats['rows_total']}")
    print(f"Lm/Lr rows indexed:            {lmlr_stats['rows_indexed']}")
    print(f"Lm non-NA rows:                {lmlr_stats['rows_mean_Lm_nonNA']}")
    print(f"Lr non-NA rows:                {lmlr_stats['rows_mean_Lr_nonNA']}")
    print(f"Lm-only rows:                  {lmlr_stats['rows_Lm_only']}")
    print("-" * 100)
    print(f"Lp rows total:                 {match_stats['lp_rows']}")
    print(f"Lp rows coordinate-mapped:     {match_stats['coordinate_mapped']}")
    print(f"Lp rows coordinate-unmapped:   {match_stats['coordinate_unmapped']}")
    print(f"Lp rows with >=1 Lm/Lr match:  {match_stats['rows_any_lmlr']}")
    print(f"Matched row-sex combinations:  {match_stats['row_sex_matches']}")
    print(f"Matched Lm-only row-sex:        {match_stats['row_sex_Lm_only']}")
    print("-" * 100)
    print(f"Output:   {args.out}")
    print(f"QC:       {qc_path}")
    print(f"Unmapped: {unmapped_path}")


if __name__ == "__main__":
    main()
