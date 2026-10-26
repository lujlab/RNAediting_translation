#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import glob
import math
import os
import re
import statistics
from collections import defaultdict


def clean(x):
    return "" if x is None else str(x).strip().lstrip("\ufeff").strip('"').strip("'")


def num(x, path, line_no, col):
    x = clean(x)
    if x == "" or x.upper() in {"NA", "NAN", "NULL"}:
        return None
    try:
        v = float(x)
    except ValueError:
        raise ValueError(f"{path}:{line_no}: non-numeric {col}={x!r}")
    if not math.isfinite(v):
        return None
    if v < 0:
        raise ValueError(f"{path}:{line_no}: negative {col}={v}")
    return v


def ratio(alt, ref):
    if alt is None or ref is None:
        return None
    den = ref + alt
    return None if den <= 0 else alt / den


def mean_non_na(vals):
    vals = [x for x in vals if x is not None and math.isfinite(x)]
    return statistics.mean(vals) if vals else None


def fmt(x):
    return "NA" if x is None else f"{x:.12g}"


def norm_handle(x):
    m = re.match(r"^(FBtr\d+)(?:\.\d+)?_(\d+)$", clean(x))
    if not m:
        return None
    return f"{m.group(1)}_{int(m.group(2))}"


def classify_file(path):
    base = os.path.basename(path)
    suffix = ".Asite_df_alledit"
    if not base.endswith(suffix):
        return None
    stem = base[: -len(suffix)]

    # AdarG / AdarS / WBg
    m = re.fullmatch(r"(AdarG|AdarS|WBg)_([FM])([123])", stem)
    if m:
        g, sx, rep = m.groups()
        cohort = {"AdarG": "Dmel_AG", "AdarS": "Dmel_AS", "WBg": "Dmel_WBg"}[g]
        genotype = {"AdarG": "AG", "AdarS": "AS", "WBg": "WBg"}[g]
        return dict(
            include=True,
            species="Dmel",
            cohort=cohort,
            genotype=genotype,
            sex="female" if sx == "F" else "male",
            sex_short=sx,
            replicate=f"{g}_{sx}{rep}",
            replicate_order=int(rep),
            reason="numbered replicate",
        )

    # ISO1 one sample per sex
    if stem == "Dmel_ISO1_head_female_ribo":
        return dict(
            include=True,
            species="Dmel",
            cohort="Dmel_ISO1",
            genotype="iso",
            sex="female",
            sex_short="F",
            replicate=stem,
            replicate_order=1,
            reason="ISO1 female",
        )
    if stem == "Dmel_ISO1_head_male_ribo":
        return dict(
            include=True,
            species="Dmel",
            cohort="Dmel_ISO1",
            genotype="iso",
            sex="male",
            sex_short="M",
            replicate=stem,
            replicate_order=1,
            reason="ISO1 male",
        )

    # Species numbered samples
    m = re.fullmatch(r"(Dsim|Dsec|Dere)_([FM])([123])", stem)
    if m:
        species, sx, rep = m.groups()
        genotype = {"Dsim": "sim", "Dsec": "sec", "Dere": "ere"}[species]
        return dict(
            include=True,
            species=species,
            cohort=species,
            genotype=genotype,
            sex="female" if sx == "F" else "male",
            sex_short=sx,
            replicate=f"{species}_{sx}{rep}",
            replicate_order=int(rep),
            reason="numbered replicate",
        )

    # Extra head-ribo samples: use as additional replicate
    m = re.fullmatch(r"(Dsim|Dsec|Dere)_(female|male)_head_ribo", stem)
    if m:
        species, sxword = m.groups()
        sx = "F" if sxword == "female" else "M"
        genotype = {"Dsim": "sim", "Dsec": "sec", "Dere": "ere"}[species]
        return dict(
            include=True,
            species=species,
            cohort=species,
            genotype=genotype,
            sex=sxword,
            sex_short=sx,
            replicate=stem,
            replicate_order=4,
            reason="additional head-ribo replicate",
        )

    # Adar5G1 does not have sex in filename; exclude from current sex analysis.
    if re.fullmatch(r"Dmel_Adar5G1_rep[12]_ribo", stem):
        return dict(
            include=False,
            species="Dmel",
            cohort="Dmel_Adar5G1",
            genotype="Adar5G1",
            sex="unknown",
            sex_short="U",
            replicate=stem,
            replicate_order=0,
            reason="ignored: no F/M sex label",
        )

    return dict(
        include=False,
        species="",
        cohort="",
        genotype="",
        sex="",
        sex_short="",
        replicate=stem,
        replicate_order=0,
        reason="ignored: unrecognized filename",
    )


REQ = [
    "FBt_cdspos",
    "ref_number_rna",
    "alt_number_rna",
    "ref_number_ribo",
    "alt_number_ribo",
]


def read_sample(path, meta, min_coverage):
    """
    Row ratio first; if duplicate FBt_cdspos occurs within one sample,
    average row-level ratios independently for Lm and Lr.
    """
    by_site = defaultdict(list)
    stats = dict(
        rows_total=0,
        rows_bad_handle=0,
        rows_pass_rna_coverage=0,
        rows_pass_ribo_coverage=0,
        rows_pass_either_coverage=0,
        rows_fail_both_coverage=0,
        rows_Lm_valid=0,
        rows_Lr_valid=0,
        rows_Lm_valid_Lr_NA=0,
    )

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        rd.fieldnames = [clean(x) for x in rd.fieldnames]
        miss = [x for x in REQ if x not in rd.fieldnames]
        if miss:
            raise ValueError(f"{path}: missing columns: {','.join(miss)}")

        for ln, raw in enumerate(rd, start=2):
            stats["rows_total"] += 1
            r = {clean(k): clean(v) for k, v in raw.items()}
            h = norm_handle(r["FBt_cdspos"])
            if h is None:
                stats["rows_bad_handle"] += 1
                continue

            rrna = num(r["ref_number_rna"], path, ln, "ref_number_rna")
            arna = num(r["alt_number_rna"], path, ln, "alt_number_rna")
            rrib = num(r["ref_number_ribo"], path, ln, "ref_number_ribo")
            arib = num(r["alt_number_ribo"], path, ln, "alt_number_ribo")

            # ----------------------------------------------------------
            # Coverage filter is applied INDEPENDENTLY to RNA and Ribo.
            #
            # RNA coverage  = ref_number_rna  + alt_number_rna
            # Ribo coverage = ref_number_ribo + alt_number_ribo
            #
            # Lm is valid only when RNA coverage >= min_coverage.
            # Lr is valid only when Ribo coverage >= min_coverage.
            #
            # Therefore a row with good RNA coverage but no Ribo coverage
            # is retained as Lm-only.
            # ----------------------------------------------------------
            rna_cov = rrna + arna if rrna is not None and arna is not None else None
            ribo_cov = rrib + arib if rrib is not None and arib is not None else None

            pass_rna = rna_cov is not None and rna_cov >= min_coverage
            pass_ribo = ribo_cov is not None and ribo_cov >= min_coverage

            if pass_rna:
                stats["rows_pass_rna_coverage"] += 1
            if pass_ribo:
                stats["rows_pass_ribo_coverage"] += 1
            if pass_rna or pass_ribo:
                stats["rows_pass_either_coverage"] += 1
            else:
                stats["rows_fail_both_coverage"] += 1
                # mRNA OR A-site coverage >= threshold is required.
                continue

            lm = ratio(arna, rrna) if pass_rna else None
            lr = ratio(arib, rrib) if pass_ribo else None

            if lm is not None:
                stats["rows_Lm_valid"] += 1
            if lr is not None:
                stats["rows_Lr_valid"] += 1
            if lm is not None and lr is None:
                stats["rows_Lm_valid_Lr_NA"] += 1

            by_site[h].append((lm, lr, rna_cov, ribo_cov, pass_rna, pass_ribo))

    out = []
    for h, vals in by_site.items():
        lms = [x[0] for x in vals if x[0] is not None]
        lrs = [x[1] for x in vals if x[1] is not None]
        lm = mean_non_na(lms)
        lr = mean_non_na(lrs)

        rna_covs = [x[2] for x in vals if x[2] is not None]
        ribo_covs = [x[3] for x in vals if x[3] is not None]

        out.append(
            dict(
                species=meta["species"],
                cohort=meta["cohort"],
                genotype=meta["genotype"],
                sex=meta["sex"],
                sex_short=meta["sex_short"],
                replicate=meta["replicate"],
                replicate_order=meta["replicate_order"],
                handle=h,
                Lm_rep=lm,
                Lr_rep=lr,
                n_rows_for_handle=len(vals),
                n_Lm_rows_nonNA=len(lms),
                n_Lr_rows_nonNA=len(lrs),
                max_rna_coverage=max(rna_covs) if rna_covs else None,
                max_ribo_coverage=max(ribo_covs) if ribo_covs else None,
                n_rows_pass_rna_coverage=sum(bool(x[4]) for x in vals),
                n_rows_pass_ribo_coverage=sum(bool(x[5]) for x in vals),
                replicate_status=(
                    "Lm_and_Lr"
                    if lm is not None and lr is not None
                    else (
                        "Lm_only"
                        if lm is not None
                        else "Lr_only" if lr is not None else "neither"
                    )
                ),
                source_file=path,
            )
        )

    return out, stats


def summarize(rep_rows):
    buckets = defaultdict(list)
    for r in rep_rows:
        key = (
            r["species"],
            r["cohort"],
            r["genotype"],
            r["sex"],
            r["sex_short"],
            r["handle"],
        )
        buckets[key].append(r)

    out = []
    for key, rows in buckets.items():
        species, cohort, genotype, sex, sex_short, handle = key
        rows = sorted(rows, key=lambda x: (x["replicate_order"], x["replicate"]))

        lms = [r["Lm_rep"] for r in rows if r["Lm_rep"] is not None]
        lrs = [r["Lr_rep"] for r in rows if r["Lr_rep"] is not None]
        mlm = mean_non_na(lms)
        mlr = mean_non_na(lrs)

        out.append(
            dict(
                species=species,
                cohort=cohort,
                genotype=genotype,
                sex=sex,
                sex_short=sex_short,
                handle=handle,
                mean_Lm=mlm,
                mean_Lr=mlr,
                n_input_samples=len(rows),
                n_Lm_replicates=len(lms),
                n_Lr_replicates=len(lrs),
                n_Lm_only_replicates=sum(
                    r["Lm_rep"] is not None and r["Lr_rep"] is None for r in rows
                ),
                n_Lr_only_replicates=sum(
                    r["Lm_rep"] is None and r["Lr_rep"] is not None for r in rows
                ),
                replicate_Lm=";".join(
                    f"{r['replicate']}={fmt(r['Lm_rep'])}" for r in rows
                ),
                replicate_Lr=";".join(
                    f"{r['replicate']}={fmt(r['Lr_rep'])}" for r in rows
                ),
                samples_used_Lm=",".join(
                    r["replicate"] for r in rows if r["Lm_rep"] is not None
                ),
                samples_used_Lr=",".join(
                    r["replicate"] for r in rows if r["Lr_rep"] is not None
                ),
                sex_mean_status=(
                    "Lm_and_Lr"
                    if mlm is not None and mlr is not None
                    else (
                        "Lm_only"
                        if mlm is not None
                        else "Lr_only" if mlr is not None else "neither"
                    )
                ),
            )
        )
    return out


def write_rep(path, rows):
    fields = [
        "species",
        "cohort",
        "genotype",
        "sex",
        "sex_short",
        "replicate",
        "replicate_order",
        "handle",
        "Lm_rep",
        "Lr_rep",
        "n_rows_for_handle",
        "n_Lm_rows_nonNA",
        "n_Lr_rows_nonNA",
        "max_rna_coverage",
        "max_ribo_coverage",
        "n_rows_pass_rna_coverage",
        "n_rows_pass_ribo_coverage",
        "replicate_status",
        "source_file",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            x = dict(r)
            x["Lm_rep"] = fmt(x["Lm_rep"])
            x["Lr_rep"] = fmt(x["Lr_rep"])
            x["max_rna_coverage"] = fmt(x.get("max_rna_coverage"))
            x["max_ribo_coverage"] = fmt(x.get("max_ribo_coverage"))
            w.writerow({k: x.get(k, "") for k in fields})


def write_mean(path, rows):
    # First six fields intentionally compatible with old Ipos-style readers.
    fields = [
        "species",
        "cohort",
        "sex",
        "handle",
        "mean_Lm",
        "mean_Lr",
        "genotype",
        "sex_short",
        "n_input_samples",
        "n_Lm_replicates",
        "n_Lr_replicates",
        "n_Lm_only_replicates",
        "n_Lr_only_replicates",
        "replicate_Lm",
        "replicate_Lr",
        "samples_used_Lm",
        "samples_used_Lr",
        "sex_mean_status",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            x = dict(r)
            x["mean_Lm"] = fmt(x["mean_Lm"])
            x["mean_Lr"] = fmt(x["mean_Lr"])
            w.writerow({k: x.get(k, "") for k in fields})


def write_manifest(path, rows):
    fields = [
        "file",
        "include",
        "species",
        "cohort",
        "genotype",
        "sex",
        "sex_short",
        "replicate",
        "reason",
        "rows_total",
        "rows_bad_handle",
        "rows_pass_rna_coverage",
        "rows_pass_ribo_coverage",
        "rows_pass_either_coverage",
        "rows_fail_both_coverage",
        "rows_Lm_valid",
        "rows_Lr_valid",
        "rows_Lm_valid_Lr_NA",
        "replicate_site_rows",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument(
        "--min-coverage",
        type=int,
        default=3,
        help=(
            "Minimum total reads required independently for each metric. "
            "Lm requires ref_number_rna + alt_number_rna >= threshold; "
            "Lr requires ref_number_ribo + alt_number_ribo >= threshold. "
            "Default: 3."
        ),
    )
    args = ap.parse_args()

    if args.min_coverage < 1:
        raise ValueError("--min-coverage must be >= 1")

    if not os.path.isdir(args.input_dir):
        raise NotADirectoryError(args.input_dir)

    files = sorted(glob.glob(os.path.join(args.input_dir, "*Asite_df_alledit")))
    if not files:
        raise FileNotFoundError(f"No *Asite_df_alledit files in {args.input_dir}")

    rep_rows = []
    manifest = []

    for path in files:
        meta = classify_file(path)
        if meta is None:
            continue

        mr = dict(file=path, **meta)

        if not meta["include"]:
            manifest.append(mr)
            continue

        rows, stats = read_sample(path, meta, args.min_coverage)
        rep_rows.extend(rows)
        mr.update(stats)
        mr["replicate_site_rows"] = len(rows)
        manifest.append(mr)

    sex_rows = summarize(rep_rows)

    cohort_order = {
        "Dmel_ISO1": 1,
        "Dmel_WBg": 2,
        "Dmel_AG": 3,
        "Dmel_AS": 4,
        "Dsim": 5,
        "Dsec": 6,
        "Dere": 7,
    }
    sex_order = {"female": 1, "male": 2}

    rep_rows.sort(
        key=lambda r: (
            cohort_order.get(r["cohort"], 999),
            sex_order.get(r["sex"], 999),
            r["handle"],
            r["replicate_order"],
            r["replicate"],
        )
    )
    sex_rows.sort(
        key=lambda r: (
            cohort_order.get(r["cohort"], 999),
            sex_order.get(r["sex"], 999),
            r["handle"],
        )
    )

    mean_out = args.out_prefix + ".sex_site_mean_LmLr.tsv"
    rep_out = args.out_prefix + ".replicate_site_LmLr.tsv"
    manifest_out = args.out_prefix + ".manifest.tsv"

    write_mean(mean_out, sex_rows)
    write_rep(rep_out, rep_rows)
    write_manifest(manifest_out, manifest)

    print("=" * 100)
    print("Unfiltered Lm/Lr recalculation finished")
    print("=" * 100)
    print("Lm = alt_number_rna / (ref_number_rna + alt_number_rna)")
    print("Lr = alt_number_ribo / (ref_number_ribo + alt_number_ribo)")
    print(f"Coverage threshold: RNA OR A-site/Ribo >= {args.min_coverage} reads.")
    print(
        f"  Lm requires RNA coverage >= {args.min_coverage}; "
        f"Lr requires Ribo coverage >= {args.min_coverage}."
    )
    print("  The two filters are independent: Lm-only sites are retained.")
    print("Lm and Lr are calculated and averaged independently.")
    print("Lm=0/Lr=0 are valid when their own coverage threshold is passed.")
    print("-" * 100)
    print(f"Files discovered:      {len(files)}")
    print(f"Files included:        {sum(bool(r.get('include')) for r in manifest)}")
    print(f"Files ignored:         {sum(not bool(r.get('include')) for r in manifest)}")
    print(f"Replicate-site rows:   {len(rep_rows)}")
    print(f"Sex-site mean rows:    {len(sex_rows)}")
    print(
        f"Lm-only replicate rows:{sum(r['replicate_status']=='Lm_only' for r in rep_rows)}"
    )
    print(
        f"Lm-only sex-site rows: {sum(r['sex_mean_status']=='Lm_only' for r in sex_rows)}"
    )
    print("-" * 100)

    for cohort in [
        "Dmel_ISO1",
        "Dmel_WBg",
        "Dmel_AG",
        "Dmel_AS",
        "Dsim",
        "Dsec",
        "Dere",
    ]:
        for sex in ["female", "male"]:
            s = [r for r in sex_rows if r["cohort"] == cohort and r["sex"] == sex]
            if not s:
                continue
            n_lm = sum(r["mean_Lm"] is not None for r in s)
            n_lr = sum(r["mean_Lr"] is not None for r in s)
            n_lmo = sum(r["mean_Lm"] is not None and r["mean_Lr"] is None for r in s)
            print(
                f"{cohort:12s} {sex:6s}: "
                f"sites={len(s):6d} Lm={n_lm:6d} Lr={n_lr:6d} Lm_only={n_lmo:6d}"
            )

    print("-" * 100)
    print(f"Main output:      {mean_out}")
    print(f"Replicate output: {rep_out}")
    print(f"Manifest:         {manifest_out}")


if __name__ == "__main__":
    main()
