#!/usr/bin/env python3
"""Compare edited nucleotide sites and edited codons in completed PausePred/edgeR peaks."""

import argparse, csv, gzip, json, math, hashlib, datetime
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

LIBRARIES = [
    ("Dmel_Adar5G1_M", "dmel_Adar5G1_M", "ISO1"),
    ("Dmel_ISO1_head_female_ribo", "dmel_ISO1_FM", "ISO1"),
    ("Dmel_ISO1_head_male_ribo", "dmel_ISO1_FM", "ISO1"),
    ("WBg_F", "dmel_WBg_F", "ISO1"),
    ("WBg_M", "dmel_WBg_M", "ISO1"),
    ("Dsim_F", "dsim_w501_F", "sim"),
    ("Dsim_M", "dsim_w501_M", "sim"),
    ("Dsec_F", "dsec_Dsec_F", "sec"),
    ("Dsec_M", "dsec_Dsec_M", "sec"),
    ("Dere_F", "dere_Dere_F", "ere"),
    ("Dere_M", "dere_Dere_M", "ere"),
]
STEM = ".candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"


def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(path, rows, fields=None, delimiter="\t"):
    rows = list(rows)
    fields = fields or list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        w = csv.DictWriter(
            f, fieldnames=fields, delimiter=delimiter, lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def integer(value):
    return int(float(value))


def edited_map(path):
    sites = {}
    for r in read(path):
        if r["target_region"] != "CDS":
            continue
        tx = r["new_transcript"].strip()
        p = integer(r["new_cds_pos_1based"])
        if p < 1:
            raise ValueError("Invalid CDS position")
        key = f"{tx}|CDS:{p}"
        codon = f"{tx}_{(p-1)//3+1}"
        sites[key] = {
            "site_key": key,
            "transcript": tx,
            "target_region": "CDS",
            "cds_pos_1based": p,
            "transcript_pos_1based": r.get("new_tx_pos_1based", ""),
            "codon_key": codon,
            "genomic_key": r.get("genomic_key", ""),
        }
    return sites


def scan_edgeR(payload):
    group, path, wanted = payload
    retained = {}
    n = 0
    with gzip.open(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            n += 1
            assert (
                r["stalling"] == "TRUE"
                and float(r["FDR_greater"]) <= 0.05
                and float(r["logFC"]) > 0
            ), (group, r["site_id"])
            key = r["site_id"]
            if key not in wanted:
                continue
            assert key not in retained
            retained[key] = {
                "edgeR_logFC": r["logFC"],
                "edgeR_FDR_greater": r["FDR_greater"],
            }
    print(f"[edgeR] {group}: {n} peaks, {len(retained)} edited peak codons", flush=True)
    return group, retained, n


def pause_sites(path, cutoff, window):
    sites = {}
    for r in read(path):
        if float(r["Pause_score_cutoff"]) != cutoff or int(r["window_nt"]) != window:
            continue
        tx = r["transcript_id"]
        region = r["target_region"]
        cds = ""
        codon = ""
        if region == "CDS":
            cds = integer(r.get("new_cds_pos_1based") or r.get("alledit_cds_pos"))
            key = f"{tx}|CDS:{cds}"
            codon = f"{tx}_{(cds-1)//3+1}"
        else:
            key = f'{tx}|TX:{r["transcript_pos"]}'
        assert key not in sites, (path, key)
        sites[key] = {
            "site_key": key,
            "transcript": tx,
            "target_region": region,
            "cds_pos_1based": cds,
            "transcript_pos_1based": r["transcript_pos"],
            "codon_key": codon,
            "genomic_key": r.get("genomic_key")
            or (r.get("chr", "") + "_" + r.get("pos", "")).strip("_"),
        }
    return sites


def summary(sample, group, a, b):
    overlap = a & b
    return {
        "library": sample,
        "edgeR_edited_in_peak": len(a),
        "pausepred_edited_in_peak": len(b),
        "overlap_sites": len(overlap),
        "overlap_percent_of_edgeR": round(100 * len(overlap) / len(a), 6) if a else "",
        "overlap_percent_of_pausepred": (
            round(100 * len(overlap) / len(b), 6) if b else ""
        ),
        "edgeR_group": group,
        "edgeR_scope_note": (
            "ISO1 F/M joint model; same edgeR set used for both ISO1 rows"
            if group == "dmel_ISO1_FM"
            else "matched library group"
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pausepred-dir", type=Path, required=True)
    ap.add_argument("--edgeR-dir", type=Path, required=True)
    ap.add_argument("--editing-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--run-tag", default="seq90_up90_down90_counterRule_newalledit")
    ap.add_argument("--edgeR-run-label", default="rigorous_v2_RNA_div2")
    ap.add_argument("--cutoff", type=float, default=5)
    ap.add_argument("--window", type=int, default=2)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    edits = {
        species: edited_map(args.editing_dir / (species + STEM))
        for species in sorted({x[2] for x in LIBRARIES})
    }
    sources = [args.editing_dir / (species + STEM) for species in edits]
    groups = {group: species for _, group, species in LIBRARIES}
    payload = []
    for group, species in groups.items():
        path = (
            args.edgeR_dir
            / f"{group}.flank30.{args.edgeR_run_label}.significant.tsv.gz"
        )
        assert path.is_file()
        sources.append(path)
        payload.append(
            (group, str(path), {r["codon_key"] for r in edits[species].values()})
        )
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(scan_edgeR, payload))
    edge = {group: rows for group, rows, _ in results}
    total = {group: n for group, _, n in results}
    original_summary = {
        r["group_id"]: integer(r["edited_event"])
        for r in read(
            args.edgeR_dir
            / "edited_unedited_summary_codonA/edited_unedited_stalling_valley_by_library.tsv"
        )
        if r["event"] == "stalling"
    }
    run_summary = {
        r["group_id"]: r for r in read(args.edgeR_dir / "edgeR_run_summary.tsv")
    }
    for group in groups:
        assert len(edge[group]) == original_summary[group], (
            group,
            len(edge[group]),
            original_summary[group],
        )
        assert total[group] == integer(run_summary[group]["n_significant"])
    nucleotide = []
    codon = []
    cds_only = []
    audit = []
    folder = args.pausepred_dir / f"enrichment_{args.run_tag}_alledit_cds"
    for sample, group, species in LIBRARIES:
        path = folder / f"{sample}.pausepred_enrichment.editing_sites_in_peak.tsv"
        sources.append(path)
        pp = pause_sites(path, args.cutoff, args.window)
        ed = {k: r for k, r in edits[species].items() if r["codon_key"] in edge[group]}
        pset = set(pp)
        eset = set(ed)
        pp_cds = {k for k, r in pp.items() if r["codon_key"]}
        pp_codons = {pp[k]["codon_key"] for k in pp_cds}
        e_codons = {r["codon_key"] for r in ed.values()}
        # Reconcile current PausePred counts and CDS identities against source editing maps.
        stats = [
            r
            for r in read(folder / f"{sample}.pausepred_enrichment.tsv")
            if float(r["Pause_score_cutoff"]) == args.cutoff
            and int(r["window_nt"]) == args.window
        ]
        assert len(stats) == 1 and integer(stats[0]["editing_in_peak"]) == len(pset)
        assert pp_cds <= set(edits[species]), (sample, pp_cds - set(edits[species]))
        nucleotide.append(summary(sample, group, eset, pset))
        codon.append(summary(sample, group, e_codons, pp_codons))
        cds_only.append(summary(sample, group, eset, pp_cds))
        members = []
        for key in sorted(eset | pset):
            r = dict(pp.get(key) or ed[key])
            r.update(
                {
                    "edgeR_in_peak": key in eset,
                    "pausepred_in_peak": key in pset,
                    "overlap": key in eset and key in pset,
                }
            )
            r.update(
                edge[group].get(
                    r["codon_key"], {"edgeR_logFC": "", "edgeR_FDR_greater": ""}
                )
            )
            members.append(r)
        fields = list(members[0])
        write(
            out / "site_lists" / f"{sample}.all_called_edited_sites.tsv",
            members,
            fields,
        )
        write(
            out / "site_lists" / f"{sample}.overlap_edited_sites.tsv",
            [r for r in members if r["overlap"]],
            fields,
        )
        audit.append(
            {
                "library": sample,
                "pausepred_nonCDS_peak_sites": len(pset - pp_cds),
                "edgeR_nucleotide_sites": len(eset),
                "edgeR_edited_codons_reproduces_existing_summary": len(e_codons)
                == original_summary[group],
                "pausepred_reproduces_existing_summary": True,
            }
        )
    write(out / "edited_peak_overlap.by_library.tsv", nucleotide)
    write(out / "edited_peak_overlap.by_library.csv", nucleotide, delimiter=",")
    write(out / "edited_peak_overlap.CDS_only.by_library.tsv", cds_only)
    write(out / "edited_peak_overlap.codon.by_library.tsv", codon)
    report = {
        "completed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "samples": 11,
        "cutoff": args.cutoff,
        "window_nt": args.window,
        "edgeR_peak_definition": "stalling TRUE: logFC > 0 and one-sided BH FDR_greater <= 0.05; RNA /2 results",
        "pausepred_definition": "existing edited-site-in-peak calls: reads >= 2, coverage >= 20%, center in eligible universe, qualifying peak within window",
        "primary_unit": "edited nucleotide site, transcript + CDS position; non-CDS PausePred calls retained as method-only sites",
        "denominators": "all calls of each method; not restricted to jointly testable sites",
        "codon_table_unit": "unique transcript + CDS codon; PausePred CDS hits collapsed to edited codons",
        "excluded_samples": ["AdarG_F", "AdarG_M", "AdarS_F", "AdarS_M"],
        "audit": audit,
        "inputs": [
            {
                "path": str(p.resolve()),
                "bytes": p.stat().st_size,
                "mtime_ns": p.stat().st_mtime_ns,
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            }
            for p in sources
        ],
    }
    pass
    lines = [
        "# PausePred and edgeR edited-site overlap",
        "",
        "Sites are defined by transcript and CDS nucleotide position. edgeR peak codons are projected onto edited nucleotides; PausePred uses score >=5 within +/-2 nt.",
        "",
        "The 11 reported sample groups exclude AdarG/AdarS and retain Adar5G1. ISO1 female and male rows share the jointly fitted edgeR model. Overlap fractions use each method-specific detected set as the denominator. The CDS-only table excludes UTR hits.",
        "",
        "Outputs: edited_peak_overlap.by_library.tsv, edited_peak_overlap.by_library.csv, edited_peak_overlap.CDS_only.by_library.tsv, edited_peak_overlap.codon.by_library.tsv, site_lists/, validation.json.",
        "",
    ]
    pass


if __name__ == "__main__":
    main()
