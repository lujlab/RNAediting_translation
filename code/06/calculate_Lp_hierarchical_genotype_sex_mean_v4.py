#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Hierarchical Lp calculation for 100AA and trypsin DIA-NN site-level tables.

Intended intensity column can be selected by --intensity-col:
  - sum_intensity (default)
  - median_intensity

Lp definition at the replicate/sample level:
    ALT_total = ALT1 + ALT2 + ALT3 + ...
    Lp = ALT_total / (REF + ALT_total)

Valid values:
    REF > 0, ALT = 0   -> Lp = 0
    REF = 0, ALT > 0   -> Lp = 1
    REF > 0, ALT > 0   -> 0 < Lp < 1
    REF = 0, ALT = 0   -> Lp = NA

Aggregation hierarchy:
    Step 1. Compute Lp independently for every replicate/sample.
    Step 2. Within each METHOD x SITE x GENOTYPE x SEX:
                Lp_method_sex = mean(non-NA replicate Lp)
            Lp=0 and Lp=1 ARE included; only NA is excluded.
    Step 3. Within each SITE x GENOTYPE x SEX:
                Lp_mean = mean(non-NA [Lp_100AA, Lp_trypsin])
            Thus:
              - both methods available -> average the two method means
              - one method NA -> use the available method mean
              - both methods NA -> final Lp_mean = NA

Site universe:
    UNION of sites from 100AA and trypsin.

Sample name parsing:
    Expected basename ending in _F1, _F2, _F3, _M1, ...
    Examples:
      DXN_20260401_ere_F1.raw -> group=ere, sex=F, rep=1
      XQH_20260320_AG_M2.raw  -> group=AG,  sex=M, rep=2
      DXN_20260401_iso_F3.raw -> group=iso, sex=F, rep=3

Outputs:
  <prefix>.sample_Lp.tsv
      one row per method x site x sample/replicate

  <prefix>.method_sex_Lp.tsv
      one row per method x site x genotype x sex
      replicate Lp values averaged within method

  <prefix>.combined_sex_Lp.tsv
      one row per site x genotype x sex
      contains Lp_100AA, Lp_trypsin, final Lp_mean

  <prefix>.combined_sex_Lp_wide.tsv
      one row per site x genotype
      F/M columns side-by-side
"""

import argparse
import csv
import gzip
import math
import os
import re
import statistics
from collections import defaultdict

META_FIELDS = [
    "site",
    "gene",
    "FBgn",
    "FBtr",
    "FBpp",
    "aa_pos",
    "refAA",
    "known_alt_AAs",
    "source_site_ids",
    "source_cds_pos",
]


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def clean(x):
    return str(x).strip().lstrip("\ufeff").strip('"').strip("'")


def parse_float(x):
    x = clean(x)
    if x == "" or x.upper() in {"NA", "NAN", "NULL"}:
        return 0.0
    try:
        v = float(x)
        if not math.isfinite(v):
            return 0.0
        return v
    except Exception:
        return 0.0


def fmt(x):
    if x is None:
        return "NA"
    return f"{x:.10g}"


def first_nonempty(a, b):
    a = clean(a)
    b = clean(b)
    return a if a else b


def infer_genotype_sex_rep(sample):
    """
    Parse basename such as:
      DXN_20260401_ere_F1.raw
      XQH_20260320_AG_M2.raw

    Returns:
      genotype = the last token immediately before _F1/_M1
      sex      = F or M
      rep      = integer replicate number

    Examples:
      DXN_20260401_iso_F1.raw -> genotype=iso
      XQH_20260320_AG_M2.raw  -> genotype=AG
      XQH_20260320_AS_F3.raw  -> genotype=AS
      XQH_20260320_WBg_M1.raw -> genotype=WBg
      DXN_20260401_ere_F1.raw -> genotype=ere
    """
    base = os.path.basename(clean(sample))
    stem = re.sub(r"\.[^.]+$", "", base)

    m = re.match(r"^(.*?)[_-]([FM])(\d+)$", stem, flags=re.IGNORECASE)
    if not m:
        raise ValueError(
            f"Cannot parse sex/replicate from sample name: {sample}\n"
            "Expected basename ending in _F1, _F2, _M1, etc."
        )

    prefix = m.group(1)
    sex = m.group(2).upper()
    rep = int(m.group(3))

    tokens = [x for x in re.split(r"[_-]+", prefix) if x]
    genotype = tokens[-1] if tokens else prefix

    return genotype, sex, rep, base


def calculate_lp(ref, alt):
    den = ref + alt
    if den <= 0:
        return None
    return alt / den


def read_method_table(path, method, intensity_col):
    """
    Read upstream state-level table and collapse ALT1/ALT2/... within
    each site x sample.
    """
    data = {}
    metadata = {}
    sites = set()
    samples = []

    sample_seen = set()

    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = [clean(x) for x in next(reader)]
        except StopIteration:
            raise ValueError(f"Empty input file: {path}")

        col = {x: i for i, x in enumerate(header)}

        required = ["site", "state", "sample", intensity_col]

        missing = [x for x in required if x not in col]
        if missing:
            raise ValueError(
                f"{method}: missing required columns: {','.join(missing)}\n"
                f"Detected columns: {','.join(header)}"
            )

        has_aa = "AA" in col

        for row in reader:
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))

            site = clean(row[col["site"]])
            state = clean(row[col["state"]]).upper()
            sample = clean(row[col["sample"]])

            if not site or not sample:
                continue

            sites.add(site)

            if sample not in sample_seen:
                samples.append(sample)
                sample_seen.add(sample)

            if site not in metadata:
                metadata[site] = {k: "" for k in META_FIELDS}
                metadata[site]["site"] = site

            for m in META_FIELDS:
                if m in col:
                    metadata[site][m] = first_nonempty(
                        metadata[site].get(m, ""), row[col[m]]
                    )

            key = (site, sample)
            if key not in data:
                data[key] = {
                    "ref_intensity": 0.0,
                    "alt_total_intensity": 0.0,
                    "alt_AAs": set(),
                    "alt_states": set(),
                }

            intensity = parse_float(row[col[intensity_col]])
            rec = data[key]

            if state == "REF":
                rec["ref_intensity"] += intensity
            elif state.startswith("ALT"):
                rec["alt_total_intensity"] += intensity

                if intensity > 0:
                    rec["alt_states"].add(state)
                    if has_aa:
                        aa = clean(row[col["AA"]]).upper()
                        if aa:
                            rec["alt_AAs"].add(aa)

    return data, metadata, sites, samples


def merge_meta(site, meta100, metatry):
    a = meta100.get(site, {})
    b = metatry.get(site, {})

    out = {"site": site}
    for k in META_FIELDS:
        if k == "site":
            continue
        out[k] = first_nonempty(a.get(k, ""), b.get(k, ""))
    return out


def build_sample_rows(method, data, union_sites, samples, meta100, metatry):
    """
    Expand to the UNION site universe. For a site/sample absent from a method:
    REF=0 and ALT=0 -> Lp=NA.
    """
    rows = []

    for site in sorted(union_sites):
        meta = merge_meta(site, meta100, metatry)

        for sample in samples:
            genotype, sex, rep, sample_name = infer_genotype_sex_rep(sample)

            rec = data.get(
                (site, sample),
                {
                    "ref_intensity": 0.0,
                    "alt_total_intensity": 0.0,
                    "alt_AAs": set(),
                    "alt_states": set(),
                },
            )

            ref = rec["ref_intensity"]
            alt = rec["alt_total_intensity"]
            lp = calculate_lp(ref, alt)

            if lp is None:
                detection_class = "NA"
            elif lp == 0:
                detection_class = "REF_ONLY"
            elif lp == 1:
                detection_class = "ALT_ONLY"
            else:
                detection_class = "BOTH"

            rows.append(
                {
                    **meta,
                    "method": method,
                    "sample": sample,
                    "sample_name": sample_name,
                    "genotype": genotype,
                    "sex": sex,
                    "replicate": rep,
                    "ref_intensity": ref,
                    "alt_total_intensity": alt,
                    "Lp": lp,
                    "detection_class": detection_class,
                    "ALT_AAs": ";".join(sorted(rec["alt_AAs"])),
                    "ALT_states": ";".join(sorted(rec["alt_states"])),
                }
            )

    return rows


def method_sex_summary(sample_rows):
    """
    Step 2:
      mean replicate-level Lp within METHOD x SITE x GENOTYPE x SEX.
      NA excluded; 0 and 1 included.
    """
    buckets = defaultdict(list)
    meta = {}

    for r in sample_rows:
        key = (r["method"], r["site"], r["genotype"], r["sex"])
        buckets[key].append(r)

        if key not in meta:
            meta[key] = {k: r.get(k, "") for k in META_FIELDS}

    out = []

    for key, rr in buckets.items():
        method, site, genotype, sex = key

        rr = sorted(rr, key=lambda x: x["replicate"])

        valid_lps = [x["Lp"] for x in rr if x["Lp"] is not None]

        lp_method_sex = statistics.mean(valid_lps) if valid_lps else None

        rep_values = []
        for x in rr:
            rep_values.append(f"R{x['replicate']}={fmt(x['Lp'])}")

        out.append(
            {
                **meta[key],
                "method": method,
                "genotype": genotype,
                "sex": sex,
                "Lp_method_sex": lp_method_sex,
                "n_rep_total": len(rr),
                "n_rep_nonNA": len(valid_lps),
                "n_rep_NA": len(rr) - len(valid_lps),
                "n_rep_Lp0": sum(x["Lp"] == 0 for x in rr if x["Lp"] is not None),
                "n_rep_Lp1": sum(x["Lp"] == 1 for x in rr if x["Lp"] is not None),
                "replicate_Lp": ";".join(rep_values),
            }
        )

    return out


def combine_methods(method_rows):
    """
    Step 3:
      for each SITE x GENOTYPE x SEX, average the method-specific sex means.

      final Lp_mean = mean(non-NA [Lp_100AA, Lp_trypsin])
    """
    buckets = defaultdict(dict)
    meta = {}

    for r in method_rows:
        key = (r["site"], r["genotype"], r["sex"])
        buckets[key][r["method"]] = r

        if key not in meta:
            meta[key] = {k: r.get(k, "") for k in META_FIELDS}

    out = []

    for key, methods in buckets.items():
        site, genotype, sex = key

        r100 = methods.get("100AA")
        rt = methods.get("trypsin")

        lp100 = r100["Lp_method_sex"] if r100 else None
        lpt = rt["Lp_method_sex"] if rt else None

        vals = [x for x in [lp100, lpt] if x is not None]
        final_mean = statistics.mean(vals) if vals else None

        out.append(
            {
                **meta[key],
                "genotype": genotype,
                "sex": sex,
                "Lp_100AA": lp100,
                "Lp_trypsin": lpt,
                "Lp_mean": final_mean,
                "n_methods_nonNA": len(vals),
                "n_rep_nonNA_100AA": r100["n_rep_nonNA"] if r100 else 0,
                "n_rep_nonNA_trypsin": rt["n_rep_nonNA"] if rt else 0,
                "replicate_Lp_100AA": r100["replicate_Lp"] if r100 else "",
                "replicate_Lp_trypsin": rt["replicate_Lp"] if rt else "",
            }
        )

    return out


def write_sample(path, rows):
    fields = META_FIELDS + [
        "method",
        "sample",
        "sample_name",
        "genotype",
        "sex",
        "replicate",
        "ref_intensity",
        "alt_total_intensity",
        "Lp",
        "detection_class",
        "ALT_AAs",
        "ALT_states",
    ]

    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()

        for r in rows:
            x = dict(r)
            x["ref_intensity"] = fmt(x["ref_intensity"])
            x["alt_total_intensity"] = fmt(x["alt_total_intensity"])
            x["Lp"] = fmt(x["Lp"])
            w.writerow({k: x.get(k, "") for k in fields})


def write_method_sex(path, rows):
    fields = META_FIELDS + [
        "method",
        "genotype",
        "sex",
        "Lp_method_sex",
        "n_rep_total",
        "n_rep_nonNA",
        "n_rep_NA",
        "n_rep_Lp0",
        "n_rep_Lp1",
        "replicate_Lp",
    ]

    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()

        for r in rows:
            x = dict(r)
            x["Lp_method_sex"] = fmt(x["Lp_method_sex"])
            w.writerow({k: x.get(k, "") for k in fields})


def write_combined(path, rows):
    fields = META_FIELDS + [
        "genotype",
        "sex",
        "Lp_100AA",
        "Lp_trypsin",
        "Lp_mean",
        "n_methods_nonNA",
        "n_rep_nonNA_100AA",
        "n_rep_nonNA_trypsin",
        "replicate_Lp_100AA",
        "replicate_Lp_trypsin",
    ]

    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()

        for r in rows:
            x = dict(r)
            x["Lp_100AA"] = fmt(x["Lp_100AA"])
            x["Lp_trypsin"] = fmt(x["Lp_trypsin"])
            x["Lp_mean"] = fmt(x["Lp_mean"])
            w.writerow({k: x.get(k, "") for k in fields})


def write_combined_wide(path, rows):
    """
    One row per site x genotype, with F/M side-by-side.
    """
    buckets = defaultdict(dict)
    meta = {}

    for r in rows:
        key = (r["site"], r["genotype"])
        buckets[key][r["sex"]] = r

        if key not in meta:
            meta[key] = {k: r.get(k, "") for k in META_FIELDS}

    fields = META_FIELDS + [
        "genotype",
        "Lp_100AA_F",
        "Lp_trypsin_F",
        "Lp_mean_F",
        "n_methods_nonNA_F",
        "n_rep_nonNA_100AA_F",
        "n_rep_nonNA_trypsin_F",
        "Lp_100AA_M",
        "Lp_trypsin_M",
        "Lp_mean_M",
        "n_methods_nonNA_M",
        "n_rep_nonNA_100AA_M",
        "n_rep_nonNA_trypsin_M",
    ]

    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()

        for key in sorted(buckets):
            site, genotype = key
            bysex = buckets[key]

            f = bysex.get("F", {})
            m = bysex.get("M", {})

            row = {
                **meta[key],
                "genotype": genotype,
                "Lp_100AA_F": fmt(f.get("Lp_100AA")),
                "Lp_trypsin_F": fmt(f.get("Lp_trypsin")),
                "Lp_mean_F": fmt(f.get("Lp_mean")),
                "n_methods_nonNA_F": f.get("n_methods_nonNA", 0),
                "n_rep_nonNA_100AA_F": f.get("n_rep_nonNA_100AA", 0),
                "n_rep_nonNA_trypsin_F": f.get("n_rep_nonNA_trypsin", 0),
                "Lp_100AA_M": fmt(m.get("Lp_100AA")),
                "Lp_trypsin_M": fmt(m.get("Lp_trypsin")),
                "Lp_mean_M": fmt(m.get("Lp_mean")),
                "n_methods_nonNA_M": m.get("n_methods_nonNA", 0),
                "n_rep_nonNA_100AA_M": m.get("n_rep_nonNA_100AA", 0),
                "n_rep_nonNA_trypsin_M": m.get("n_rep_nonNA_trypsin", 0),
            }

            w.writerow({k: row.get(k, "") for k in fields})


def main():
    p = argparse.ArgumentParser(
        description=(
            "Hierarchical Lp: replicate -> method-specific sex mean "
            "-> mean across 100AA/trypsin within sex."
        )
    )
    p.add_argument("--aa100", required=True, help="100AA site-level TSV")
    p.add_argument("--trypsin", required=True, help="Trypsin site-level TSV")
    p.add_argument(
        "--intensity-col",
        default="sum_intensity",
        choices=["sum_intensity", "median_intensity"],
        help="Intensity column to aggregate (default: sum_intensity)",
    )
    p.add_argument("--out-prefix", required=True, help="Output prefix")
    args = p.parse_args()

    for path in [args.aa100, args.trypsin]:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    intensity_col = args.intensity_col

    data100, meta100, sites100, samples100 = read_method_table(
        args.aa100, "100AA", intensity_col
    )

    datat, metat, sitest, samplest = read_method_table(
        args.trypsin, "trypsin", intensity_col
    )

    union_sites = sites100 | sitest

    sample100 = build_sample_rows(
        "100AA", data100, union_sites, samples100, meta100, metat
    )
    samplet = build_sample_rows("trypsin", datat, union_sites, samplest, meta100, metat)

    sample_rows = sample100 + samplet

    method_rows = method_sex_summary(sample_rows)
    combined_rows = combine_methods(method_rows)

    sample_out = args.out_prefix + ".sample_Lp.tsv"
    method_out = args.out_prefix + ".method_sex_Lp.tsv"
    combined_out = args.out_prefix + ".combined_sex_Lp.tsv"
    wide_out = args.out_prefix + ".combined_sex_Lp_wide.tsv"

    write_sample(sample_out, sample_rows)
    write_method_sex(method_out, method_rows)
    write_combined(combined_out, combined_rows)
    write_combined_wide(wide_out, combined_rows)

    print("=" * 94)
    print("Hierarchical Lp calculation finished")
    print("=" * 94)
    print("Step 1: replicate Lp = ALT_total / (REF + ALT_total)")
    print("Step 2: within each method, mean replicate Lp by genotype and sex")
    print("Step 3: within each genotype and sex, mean 100AA and trypsin method means")
    print("Rule: Lp=0 and Lp=1 INCLUDED; only NA EXCLUDED.")
    print("Site universe = UNION of 100AA and trypsin sites.")
    print("-" * 94)
    print(f"100AA original sites:   {len(sites100)}")
    print(f"trypsin original sites: {len(sitest)}")
    print(f"union sites:            {len(union_sites)}")
    print(f"100AA samples:          {len(samples100)}")
    print(f"trypsin samples:        {len(samplest)}")
    print("-" * 94)
    print("Outputs:")
    print(sample_out)
    print(method_out)
    print(combined_out)
    print(wide_out)


if __name__ == "__main__":
    main()
