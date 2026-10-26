#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Summarize allele-level intensities for DIA-NN 100AA tables.

For each site × genotype × replicate, count allele intensity and derive
proportion = allele_intensity / total_intensity.

Outputs both:
  1) sex-split tables (F/M)
  2) non-sex-split tables (pooled by genotype)

Allele mapping:
  - default (mode=encoded): AA/IA/AI/II
    AA: reference
    IA: single-edit A at first edited position (or only edited position)
    AI: single-edit B at second edited position
    II: double-edit outcome
  - mode=raw: legacy behaviour (REF/AA code)
"""

import argparse
import csv
import math
import os
import re
from collections import defaultdict
from statistics import mean

AA3_TO_1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "GLY": "G",
    "TRY": "W",
}


def clean(x):
    return "" if x is None else str(x).strip().strip('"').strip("'").strip("\ufeff")


def parse_float(x):
    s = clean(x)
    if s == "" or s.upper() in {"NA", "NAN", "NULL"}:
        return None
    try:
        v = float(s)
        if math.isfinite(v):
            return v
    except Exception:
        return None
    return None


def _dataset_key(path):
    base = os.path.basename(clean(path))
    if base.endswith(".tsv"):
        stem = base[:-4]
    else:
        stem = base
    stem = stem.replace(".unlisted_alt", "")
    m = re.match(r"^(.*)_site_level_intensity_100AA", stem)
    if m:
        return m.group(1)
    return stem


def infer_rep(sample):
    base = os.path.basename(clean(sample))
    stem = re.sub(r"\.[^.]+$", "", base)

    m = re.match(r"^(.*?)[_-]([FM])(\d+)$", stem, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Cannot parse replicate from sample name: {sample}")

    prefix = m.group(1)
    sex = m.group(2).upper()
    rep = int(m.group(3))
    tokens = [x for x in re.split(r"[_-]+", prefix) if x]
    genotype = tokens[-1] if tokens else prefix

    return genotype, sex, rep, base


def _fmt(x):
    return "NA" if x is None else f"{x:.10g}"


def split_multi(v):
    s = clean(v)
    if not s:
        return []
    return [clean(x) for x in re.split(r"\s*;\s*", s) if clean(x)]


def aa3_to_1(x):
    if not x:
        return ""
    s = clean(x).replace(" ", "")
    if len(s) == 1:
        return s.upper()
    return AA3_TO_1.get(s.upper(), s)


def parse_p_from_ann_token(tok):
    s = clean(tok)
    if not s.startswith("p."):
        return None
    # p.Lys137Glu, p.R137W, p. =, p.* etc
    m = re.search(r"p\.([A-Za-z*]{1,3})(\d+)([A-Za-z*]{1,3}|=)", s, flags=re.IGNORECASE)
    if not m:
        return None

    # skip synonymous variants unless codon change is explicit
    ref_aa = aa3_to_1(m.group(1))
    alt_aa = aa3_to_1(m.group(3))
    return ref_aa, alt_aa


def parse_ann_position(tok):
    s = clean(tok)
    if not s.startswith("c."):
        return None
    m = re.search(r"c\.(\d+)", s)
    return m.group(1) if m else None


def load_ann_map(paths):
    """
    Return map[(FBtr, handle, cds_pos)] = alt_aa (1-letter)
    """
    ann_map = {}
    if not paths:
        return ann_map

    for p in paths:
        if not os.path.isfile(p):
            raise FileNotFoundError(p)

        with open(p, "r", encoding="utf-8", errors="replace", newline="") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            if not r.fieldnames:
                raise ValueError(f"Empty annotation file: {p}")

            for row in r:
                tx = clean(row.get("transcript_id", "")) or clean(row.get("FBtr", ""))
                site_id = clean(row.get("handle", ""))
                ann = clean(row.get("ann", ""))
                if not tx or not ann or not site_id:
                    continue

                pos = None
                alt_aa = None
                for part in [x.strip() for x in ann.split("|") if clean(x)]:
                    if part.startswith("c."):
                        pos = parse_ann_position(part)
                    if part.startswith("p.") and alt_aa is None:
                        pa = parse_p_from_ann_token(part)
                        if pa is not None:
                            _, alt = pa
                            if alt not in {"", "="}:
                                alt_aa = alt

                if not pos or not alt_aa:
                    # fallback from hgvs_c + aa field if available
                    hgvs = clean(row.get("hgvs_c", ""))
                    if not pos and hgvs:
                        pos = parse_ann_position(hgvs)
                    if not alt_aa:
                        p_tok = clean(row.get("ann", "")).split("|")
                        for x in p_tok:
                            if x.startswith("p."):
                                pa = parse_p_from_ann_token(x)
                                if pa is not None:
                                    _, alt = pa
                                    if alt not in {"", "="}:
                                        alt_aa = alt

                if tx and site_id and pos and alt_aa:
                    ann_map[(tx, site_id, pos)] = alt_aa.upper()

    return ann_map


def load_unlisted_sites(paths):
    if not paths:
        return {}
    site_filter = {}
    for p in paths:
        if not os.path.isfile(p):
            raise FileNotFoundError(p)

        k = _dataset_key(p)
        if k not in site_filter:
            site_filter[k] = set()

        with open(p, "r", encoding="utf-8", errors="replace", newline="") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            if not r.fieldnames:
                raise ValueError(f"Empty unlisted file: {p}")
            if "site" not in r.fieldnames:
                raise ValueError(f"unlisted file missing required column `site`: {p}")
            for row in r:
                site = clean(row.get("site", ""))
                if site:
                    site_filter[k].add(site)

    return {k: v for k, v in site_filter.items() if v}


def get_required_columns(header, candidates):
    for c in candidates:
        if c in header:
            return c
    raise ValueError(f"None of required columns found: {', '.join(candidates)}")


def classify_allele(row, ann_map, mode):
    """
    map observed AA to AA/IA/AI/II by ann-derived codon result when available.
    """
    if mode == "raw":
        state = clean(row.get("state", "")).upper()
        aa = clean(row.get("AA", "")).upper()
        if state == "REF":
            return "REF"
        if state.startswith("ALT") and aa:
            return aa
        return "REF"

    state = clean(row.get("state", "")).upper()
    aa = clean(row.get("AA", "")).upper()
    ref_aa = clean(row.get("refAA", "")).upper()
    tx = (
        clean(row.get("FBtr", "")).upper()
        or clean(row.get("transcript_id", "")).upper()
    )
    src_sites = split_multi(row.get("source_site_ids", ""))
    src_pos = split_multi(row.get("source_cds_pos", ""))
    known_alt = split_multi(row.get("known_alt_AAs", ""))

    # prefer ann-derived codon-to-AA for each edited position
    ann_alt = []
    n = min(len(src_sites), len(src_pos))
    for i in range(n):
        sid = src_sites[i]
        pos = src_pos[i]
        key = (tx, sid, pos)
        alt = ann_map.get(
            (tx.upper(), sid.upper(), str(int(pos)) if pos.isdigit() else pos), ""
        )
        if alt:
            ann_alt.append(alt)

    while len(ann_alt) < n and len(known_alt) > len(ann_alt):
        ann_alt.append(known_alt[len(ann_alt)].upper())

    if not src_pos:
        # fallback without position info
        if state == "REF":
            return "AA"
        return aa if aa else ("ALT" if mode == "encoded" else "AA")

    if state == "REF":
        return "AA"

    if n <= 1:
        # single edit site: use IA as the edited class
        if state.startswith("ALT"):
            return "IA"
        return "AA"

    # double/complex edit site
    if n >= 2:
        alt1 = ann_alt[0] if len(ann_alt) > 0 else ""
        alt2 = ann_alt[1] if len(ann_alt) > 1 else ""
        if aa == alt1:
            return "IA"
        if aa == alt2:
            return "AI"
        # If it matches neither single-edit aa, treat as double-edit outcome.
        # This also handles AA that only exists in the unlisted_alt table.
        return "II"

    # fallback
    return aa if aa else "II"


def read_intensity_tables(paths, intensity_col, ann_map, allele_mode, site_filter=None):
    # with sex: (site, genotype, sex, rep, sample) -> {allele: intensity}
    # no sex:  (site, genotype, rep)           -> {allele: intensity}
    per_rep_sex = defaultdict(lambda: defaultdict(float))
    per_rep_nosex = defaultdict(lambda: defaultdict(float))

    for p in paths:
        if not os.path.isfile(p):
            raise FileNotFoundError(p)

        ds_key = _dataset_key(p)
        has_filter = bool(site_filter)
        filter_sites = site_filter.get(ds_key, set()) if has_filter else set()

        with open(p, "r", encoding="utf-8", errors="replace", newline="") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            header = [clean(h) for h in r.fieldnames] if r.fieldnames else []
            if not header:
                raise ValueError(f"Empty file: {p}")

            header_set = set(header)
            req = {"site", "state", "sample", "AA"}
            if not req.issubset(header_set):
                missing = sorted(req - header_set)
                raise ValueError(f"{p} missing required columns: {','.join(missing)}")

            intensity_name = get_required_columns(
                header,
                [intensity_col, "sum_intensity", "samplesum_intensity", "intensity"],
            )

            for row in r:
                site = clean(row.get("site", ""))
                sample = clean(row.get("sample", ""))
                if not site or not sample:
                    continue

                if has_filter:
                    if ds_key not in site_filter:
                        continue
                    if site not in filter_sites:
                        continue

                try:
                    genotype, sex, rep, sample_name = infer_rep(sample)
                except Exception:
                    continue

                intensity = parse_float(row.get(intensity_name, ""))
                if intensity is None:
                    intensity = 0.0

                allele = classify_allele(row, ann_map, allele_mode)
                if not allele:
                    continue

                key_sex = (site, genotype, sex, rep, sample_name)
                key_nosex = (site, genotype, rep)

                per_rep_sex[key_sex][allele] += intensity
                per_rep_nosex[key_nosex][allele] += intensity

    return per_rep_sex, per_rep_nosex


def _allele_set_by_group(per_rep, group_positions):
    out = defaultdict(set)
    for key, allele_map in per_rep.items():
        base = tuple(key[i] for i in group_positions)
        out[base].update(allele_map.keys())
    return out


def write_per_rep_with_sex(out_path, per_rep):
    base_alleles = _allele_set_by_group(per_rep, [0, 1, 2])
    header = [
        "site",
        "genotype",
        "sex",
        "rep",
        "sample",
        "allele",
        "intensity",
        "total_intensity",
        "proportion",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=header)
        w.writeheader()

        for key in sorted(per_rep):
            site, genotype, sex, rep, sample_name = key
            allele_map = per_rep[key]
            total_intensity = sum(allele_map.values())
            alleles = sorted(base_alleles[(site, genotype, sex)])

            for allele in alleles:
                intensity = allele_map.get(allele, 0.0)
                proportion = (
                    intensity / total_intensity if total_intensity > 0 else None
                )
                w.writerow(
                    {
                        "site": site,
                        "genotype": genotype,
                        "sex": sex,
                        "rep": rep,
                        "sample": sample_name,
                        "allele": allele,
                        "intensity": _fmt(intensity),
                        "total_intensity": _fmt(total_intensity),
                        "proportion": _fmt(proportion),
                    }
                )


def write_per_rep_no_sex(out_path, per_rep):
    base_alleles = _allele_set_by_group(per_rep, [0, 1])
    header = [
        "site",
        "genotype",
        "rep",
        "allele",
        "intensity",
        "total_intensity",
        "proportion",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=header)
        w.writeheader()

        for key in sorted(per_rep):
            site, genotype, rep = key
            allele_map = per_rep[key]
            total_intensity = sum(allele_map.values())
            alleles = sorted(base_alleles[(site, genotype)])

            for allele in alleles:
                intensity = allele_map.get(allele, 0.0)
                proportion = (
                    intensity / total_intensity if total_intensity > 0 else None
                )
                w.writerow(
                    {
                        "site": site,
                        "genotype": genotype,
                        "rep": rep,
                        "allele": allele,
                        "intensity": _fmt(intensity),
                        "total_intensity": _fmt(total_intensity),
                        "proportion": _fmt(proportion),
                    }
                )


def write_rep_means_with_sex(out_path, per_rep):
    # Collapse by replicate first, then average replicate proportions per site×genotype×sex×allele
    rep_allele = defaultdict(lambda: defaultdict(float))
    for key, allele_map in per_rep.items():
        site, genotype, sex, rep, _ = key
        rep_key = (site, genotype, sex, rep)
        for allele, intensity in allele_map.items():
            rep_allele[rep_key][allele] += intensity

    base_alleles = _allele_set_by_group(rep_allele, [0, 1, 2])
    reps_by_group = defaultdict(dict)
    for rep_key, allele_map in rep_allele.items():
        site, genotype, sex, rep = rep_key
        reps_by_group[(site, genotype, sex)][rep] = allele_map

    agg = {}
    for base, rep_maps in reps_by_group.items():
        site, genotype, sex = base
        alleles = sorted(base_alleles[base])
        for allele in alleles:
            props = []
            ints = []
            for rep, allele_map in sorted(rep_maps.items()):
                total = sum(allele_map.values())
                if total <= 0:
                    continue
                intensity = allele_map.get(allele, 0.0)
                props.append(intensity / total)
                ints.append(intensity)

            if props:
                agg[(site, genotype, sex, allele)] = (props, ints)

    header = [
        "site",
        "genotype",
        "sex",
        "allele",
        "n_reps_with_allele",
        "mean_proportion",
        "mean_intensity",
        "prop_min",
        "prop_max",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=header)
        w.writeheader()

        for key in sorted(agg):
            site, genotype, sex, allele = key
            props, ints = agg[key]
            w.writerow(
                {
                    "site": site,
                    "genotype": genotype,
                    "sex": sex,
                    "allele": allele,
                    "n_reps_with_allele": len(props),
                    "mean_proportion": _fmt(mean(props)),
                    "mean_intensity": _fmt(mean(ints)),
                    "prop_min": _fmt(min(props)),
                    "prop_max": _fmt(max(props)),
                }
            )


def write_rep_means_no_sex(out_path, per_rep):
    # Collapse by replicate first, then average replicate proportions per site×genotype×allele
    rep_allele = defaultdict(lambda: defaultdict(float))
    for key, allele_map in per_rep.items():
        site, genotype, rep = key
        rep_key = (site, genotype, rep)
        for allele, intensity in allele_map.items():
            rep_allele[rep_key][allele] += intensity

    base_alleles = _allele_set_by_group(rep_allele, [0, 1])
    reps_by_group = defaultdict(dict)
    for rep_key, allele_map in rep_allele.items():
        site, genotype, rep = rep_key
        reps_by_group[(site, genotype)][rep] = allele_map

    agg = {}
    for base, rep_maps in reps_by_group.items():
        site, genotype = base
        alleles = sorted(base_alleles[base])
        for allele in alleles:
            props = []
            ints = []
            for rep, allele_map in sorted(rep_maps.items()):
                total = sum(allele_map.values())
                if total <= 0:
                    continue
                intensity = allele_map.get(allele, 0.0)
                props.append(intensity / total)
                ints.append(intensity)

            if props:
                agg[(site, genotype, allele)] = (props, ints)

    header = [
        "site",
        "genotype",
        "allele",
        "n_reps_with_allele",
        "mean_proportion",
        "mean_intensity",
        "prop_min",
        "prop_max",
    ]

    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=header)
        w.writeheader()

        for key in sorted(agg):
            site, genotype, allele = key
            props, ints = agg[key]
            w.writerow(
                {
                    "site": site,
                    "genotype": genotype,
                    "allele": allele,
                    "n_reps_with_allele": len(props),
                    "mean_proportion": _fmt(mean(props)),
                    "mean_intensity": _fmt(mean(ints)),
                    "prop_min": _fmt(min(props)),
                    "prop_max": _fmt(max(props)),
                }
            )


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize per-allele intensity and proportion for 100AA tables."
    )
    p.add_argument(
        "--in",
        dest="inputs",
        nargs="+",
        required=True,
        help="One or more 100AA site-level tables",
    )
    p.add_argument(
        "--intensity-col",
        default="samplesum_intensity",
        help="Preferred intensity column name. If absent, fallback order: samplesum_intensity, sum_intensity, intensity",
    )
    p.add_argument(
        "--ann",
        nargs="+",
        default=[],
        help="Optional annotation file(s) with an `ann` column for codon-to-AA mapping (e.g. *.annotations.long.tsv).",
    )
    p.add_argument(
        "--unlisted",
        nargs="+",
        default=[],
        help="Optional species-specific unlisted-alt files. Only sites appearing in these files are kept (paired by filename key).",
    )
    p.add_argument(
        "--allele-mode",
        choices=["encoded", "raw"],
        default="encoded",
        help="encoded=AA/IA/AI/II, raw=legacy REF/AA style",
    )
    p.add_argument(
        "--out-rep",
        required=True,
        help="Output path: per-replicate allele table with sex split",
    )
    p.add_argument(
        "--out-mean",
        required=True,
        help="Output path: per-site×genotype×sex×allele mean table",
    )
    p.add_argument(
        "--out-rep-nosex",
        default=None,
        help="Output path: per-replicate allele table without sex split",
    )
    p.add_argument(
        "--out-mean-nosex",
        default=None,
        help="Output path: per-site×genotype×allele mean table without sex split",
    )
    return p.parse_args()


def _derive_output_path(path, suffix):
    dirname = os.path.dirname(path)
    basename = os.path.basename(path)
    stem = basename[:-4] if basename.endswith(".tsv") else basename
    return os.path.join(dirname, f"{stem}{suffix}.tsv")


def main():
    args = parse_args()

    ann_map = load_ann_map(args.ann)
    site_filter = load_unlisted_sites(args.unlisted)
    per_rep_sex, per_rep_nosex = read_intensity_tables(
        args.inputs,
        args.intensity_col,
        ann_map,
        args.allele_mode,
        site_filter=site_filter,
    )

    out_rep_nosex = args.out_rep_nosex or _derive_output_path(args.out_rep, "_nosex")
    out_mean_nosex = args.out_mean_nosex or _derive_output_path(args.out_mean, "_nosex")

    write_per_rep_with_sex(args.out_rep, per_rep_sex)
    write_rep_means_with_sex(args.out_mean, per_rep_sex)
    write_per_rep_no_sex(out_rep_nosex, per_rep_nosex)
    write_rep_means_no_sex(out_mean_nosex, per_rep_nosex)

    n_filtered_species = len(site_filter)
    n_filtered_sites = sum(len(v) for v in site_filter.values())
    print(f"annotation records loaded: {len(ann_map)}")
    print(
        f"site filter by unlisted files: {n_filtered_species} species, {n_filtered_sites} sites"
    )
    print("saved:", args.out_rep)
    print("saved:", args.out_mean)
    print("saved:", out_rep_nosex)
    print("saved:", out_mean_nosex)


if __name__ == "__main__":
    main()
