#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import gzip
import os


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def clean(x):
    return str(x).strip().lstrip("\ufeff").strip('"').strip("'")


def as_float(x):
    x = clean(x)
    if x == "" or x.upper() in {"NA", "NAN", "NULL"}:
        return 0.0
    try:
        return float(x)
    except ValueError:
        return 0.0


def first_nonempty(a, b):
    a = str(a or "").strip()
    b = str(b or "").strip()
    return a if a else b


def read_table(path, method):
    """Collapse sample-level rows into site-level REF/ALT detection."""
    info = {}
    alt_states = set()  # (site, observed ALT AA)
    stats = {"rows_total": 0, "rows_positive": 0}

    with open_text(path) as f:
        r = csv.reader(f, delimiter="\t")
        header = [clean(x) for x in next(r)]
        col = {x: i for i, x in enumerate(header)}

        required = ["site", "state", "sum_intensity"]
        missing = [x for x in required if x not in col]
        if missing:
            raise ValueError(
                f"{method}: missing required columns: {missing}\n"
                f"Detected columns: {header}"
            )

        meta_cols = [
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

        for row in r:
            stats["rows_total"] += 1
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))

            site = clean(row[col["site"]])
            if not site:
                continue

            if site not in info:
                info[site] = {
                    "site": site,
                    "gene": "",
                    "FBgn": "",
                    "FBtr": "",
                    "FBpp": "",
                    "aa_pos": "",
                    "refAA": "",
                    "known_alt_AAs": "",
                    "source_site_ids": "",
                    "source_cds_pos": "",
                    "ref_detected": False,
                    "alt_detected": False,
                    "alt_AAs": set(),
                }

            x = info[site]
            for m in meta_cols:
                if m in col:
                    x[m] = first_nonempty(x[m], row[col[m]])

            intensity = as_float(row[col["sum_intensity"]])
            if intensity <= 0:
                continue
            stats["rows_positive"] += 1

            state = clean(row[col["state"]]).upper()
            if state == "REF":
                x["ref_detected"] = True
            elif state.startswith("ALT"):
                x["alt_detected"] = True
                aa = clean(row[col["AA"]]).upper() if "AA" in col else ""
                if aa:
                    x["alt_AAs"].add(aa)
                    alt_states.add((site, aa))

    return info, alt_states, stats


def sets_from_info(info):
    ref = {s for s, x in info.items() if x["ref_detected"]}
    alt = {s for s, x in info.items() if x["alt_detected"]}
    return ref, alt, ref | alt


def pct(n, d):
    return "NA" if d == 0 else f"{100*n/d:.4f}"


def jaccard(a, b):
    u = a | b
    return "NA" if not u else f"{len(a & b)/len(u):.6f}"


def merged_meta(site, a, b):
    x = a.get(site, {})
    y = b.get(site, {})
    out = {"site": site}
    for k in [
        "gene",
        "FBgn",
        "FBtr",
        "FBpp",
        "aa_pos",
        "refAA",
        "known_alt_AAs",
        "source_site_ids",
        "source_cds_pos",
    ]:
        out[k] = first_nonempty(x.get(k, ""), y.get(k, ""))
    return out


def write_summary(path, categories):
    fields = [
        "category",
        "n_100AA",
        "n_trypsin",
        "n_overlap",
        "n_100AA_only",
        "n_trypsin_only",
        "pct_100AA_also_in_trypsin",
        "pct_trypsin_also_in_100AA",
        "jaccard",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for label, a, b in categories:
            ov = a & b
            w.writerow(
                {
                    "category": label,
                    "n_100AA": len(a),
                    "n_trypsin": len(b),
                    "n_overlap": len(ov),
                    "n_100AA_only": len(a - b),
                    "n_trypsin_only": len(b - a),
                    "pct_100AA_also_in_trypsin": pct(len(ov), len(a)),
                    "pct_trypsin_also_in_100AA": pct(len(ov), len(b)),
                    "jaccard": jaccard(a, b),
                }
            )


def write_presence(path, all_sites, a, b):
    fields = [
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
        "REF_100AA",
        "REF_trypsin",
        "REF_overlap",
        "ALT_100AA",
        "ALT_trypsin",
        "ALT_overlap",
        "ANY_100AA",
        "ANY_trypsin",
        "ANY_overlap",
        "ALT_AAs_100AA",
        "ALT_AAs_trypsin",
    ]

    rows = []
    for site in all_sites:
        meta = merged_meta(site, a, b)
        xa = a.get(site, {})
        xb = b.get(site, {})
        ra = bool(xa.get("ref_detected", False))
        rb = bool(xb.get("ref_detected", False))
        aa = bool(xa.get("alt_detected", False))
        ab = bool(xb.get("alt_detected", False))
        rows.append(
            {
                **meta,
                "REF_100AA": int(ra),
                "REF_trypsin": int(rb),
                "REF_overlap": int(ra and rb),
                "ALT_100AA": int(aa),
                "ALT_trypsin": int(ab),
                "ALT_overlap": int(aa and ab),
                "ANY_100AA": int(ra or aa),
                "ANY_trypsin": int(rb or ab),
                "ANY_overlap": int((ra or aa) and (rb or ab)),
                "ALT_AAs_100AA": ";".join(sorted(xa.get("alt_AAs", set()))),
                "ALT_AAs_trypsin": ";".join(sorted(xb.get("alt_AAs", set()))),
            }
        )

    def sk(x):
        try:
            p = int(x["aa_pos"])
        except Exception:
            p = 10**15
        return (x["FBtr"], p, x["site"])

    rows.sort(key=sk)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def write_site_subset(path, sites, a, b, category):
    fields = [
        "site",
        "gene",
        "FBgn",
        "FBtr",
        "FBpp",
        "aa_pos",
        "refAA",
        "known_alt_AAs",
        "category",
        "REF_100AA",
        "REF_trypsin",
        "ALT_100AA",
        "ALT_trypsin",
        "ALT_AAs_100AA",
        "ALT_AAs_trypsin",
    ]
    rows = []
    for site in sites:
        meta = merged_meta(site, a, b)
        xa = a.get(site, {})
        xb = b.get(site, {})
        rows.append(
            {
                **{
                    k: meta[k]
                    for k in [
                        "site",
                        "gene",
                        "FBgn",
                        "FBtr",
                        "FBpp",
                        "aa_pos",
                        "refAA",
                        "known_alt_AAs",
                    ]
                },
                "category": category,
                "REF_100AA": int(bool(xa.get("ref_detected", False))),
                "REF_trypsin": int(bool(xb.get("ref_detected", False))),
                "ALT_100AA": int(bool(xa.get("alt_detected", False))),
                "ALT_trypsin": int(bool(xb.get("alt_detected", False))),
                "ALT_AAs_100AA": ";".join(sorted(xa.get("alt_AAs", set()))),
                "ALT_AAs_trypsin": ";".join(sorted(xb.get("alt_AAs", set()))),
            }
        )
    rows.sort(key=lambda x: (x["FBtr"], x["aa_pos"], x["site"]))
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def write_alt_state(path, states, a, b, category):
    fields = [
        "site",
        "gene",
        "FBgn",
        "FBtr",
        "FBpp",
        "aa_pos",
        "refAA",
        "known_alt_AAs",
        "ALT_AA",
        "category",
    ]
    rows = []
    for site, alt_aa in states:
        m = merged_meta(site, a, b)
        rows.append(
            {
                "site": site,
                "gene": m["gene"],
                "FBgn": m["FBgn"],
                "FBtr": m["FBtr"],
                "FBpp": m["FBpp"],
                "aa_pos": m["aa_pos"],
                "refAA": m["refAA"],
                "known_alt_AAs": m["known_alt_AAs"],
                "ALT_AA": alt_aa,
                "category": category,
            }
        )
    rows.sort(key=lambda x: (x["FBtr"], x["aa_pos"], x["site"], x["ALT_AA"]))
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(
        description="Compare 100AA vs trypsin REF/ALT site overlap"
    )
    p.add_argument("--aa100", required=True)
    p.add_argument("--trypsin", required=True)
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    for path in [args.aa100, args.trypsin]:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    a, alt_state_a, stat_a = read_table(args.aa100, "100AA")
    b, alt_state_b, stat_b = read_table(args.trypsin, "trypsin")

    ref_a, alt_a, any_a = sets_from_info(a)
    ref_b, alt_b, any_b = sets_from_info(b)

    categories = [
        ("REF", ref_a, ref_b),
        ("ALT", alt_a, alt_b),
        ("ANY_REF_or_ALT", any_a, any_b),
    ]

    write_summary(args.out_prefix + ".summary.tsv", categories)
    write_presence(args.out_prefix + ".site_presence.tsv", any_a | any_b, a, b)

    for label, sa, sb in [
        ("REF", ref_a, ref_b),
        ("ALT", alt_a, alt_b),
        ("ANY", any_a, any_b),
    ]:
        write_site_subset(
            args.out_prefix + f".{label}.overlap.tsv", sa & sb, a, b, f"{label}_overlap"
        )
        write_site_subset(
            args.out_prefix + f".{label}.100AA_only.tsv",
            sa - sb,
            a,
            b,
            f"{label}_100AA_only",
        )
        write_site_subset(
            args.out_prefix + f".{label}.trypsin_only.tsv",
            sb - sa,
            a,
            b,
            f"{label}_trypsin_only",
        )

    # Diagnostic: compare exact ALT amino-acid state as (site, AA)
    ov_state = alt_state_a & alt_state_b
    write_alt_state(
        args.out_prefix + ".ALT_state.overlap.tsv", ov_state, a, b, "ALT_state_overlap"
    )
    write_alt_state(
        args.out_prefix + ".ALT_state.100AA_only.tsv",
        alt_state_a - alt_state_b,
        a,
        b,
        "ALT_state_100AA_only",
    )
    write_alt_state(
        args.out_prefix + ".ALT_state.trypsin_only.tsv",
        alt_state_b - alt_state_a,
        a,
        b,
        "ALT_state_trypsin_only",
    )

    with open(
        args.out_prefix + ".ALT_state.summary.tsv", "w", encoding="utf-8", newline=""
    ) as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "category",
                "n_100AA",
                "n_trypsin",
                "n_overlap",
                "n_100AA_only",
                "n_trypsin_only",
                "jaccard",
            ]
        )
        w.writerow(
            [
                "ALT_site_plus_AA",
                len(alt_state_a),
                len(alt_state_b),
                len(ov_state),
                len(alt_state_a - alt_state_b),
                len(alt_state_b - alt_state_a),
                jaccard(alt_state_a, alt_state_b),
            ]
        )

    print("=" * 90)
    print(
        "100AA vs trypsin overlap (site-level; positive if sum_intensity > 0 in >=1 sample)"
    )
    print("=" * 90)
    for label, sa, sb in categories:
        ov = sa & sb
        print(
            f"{label:15s} 100AA={len(sa):5d}  trypsin={len(sb):5d}  "
            f"overlap={len(ov):5d}  100AA_only={len(sa-sb):5d}  "
            f"trypsin_only={len(sb-sa):5d}  "
            f"100AA->trypsin={pct(len(ov),len(sa))}%  "
            f"trypsin->100AA={pct(len(ov),len(sb))}%"
        )
    print("-" * 90)
    print(
        f"ALT exact state (site+AA): 100AA={len(alt_state_a)}  trypsin={len(alt_state_b)}  "
        f"overlap={len(ov_state)}"
    )
    print(
        f'Positive input rows: 100AA={stat_a["rows_positive"]}/{stat_a["rows_total"]}; '
        f'trypsin={stat_b["rows_positive"]}/{stat_b["rows_total"]}'
    )
    print("Summary:", args.out_prefix + ".summary.tsv")
    print("Presence:", args.out_prefix + ".site_presence.tsv")


if __name__ == "__main__":
    main()
