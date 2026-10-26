#!/usr/bin/env python3
"""Independent gmRC recalculation; standard library only, original inputs read-only."""

from __future__ import annotations

import argparse
import os
from array import array
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import gzip
import hashlib
import heapq
import json
import math
from pathlib import Path
import pickle
import re
from statistics import median
import sys

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[2]
PROJECT = Path(os.environ.get("RNA_EDITING_DATA_ROOT", REPO.parent))
PEAK = Path(
    os.environ.get(
        "GMRC_ANNOTATION_ROOT",
        "/path/to/Proj1_editing/codex/02_stalling/3flank/peak_find",
    )
)
SPECIES = {
    "dmel": ("Dmel", "dmelr6.12", "dmel_r6.12", "dmel-all-CDS-r6.12"),
    "dsim": ("Dsim", "dsimr2.02", "dsim_all_r2.02", "dsim-all-CDS-r2.02"),
    "dsec": ("Dsec", "dsecr1.3", "dsec_all_r1.3", "dsec-all-CDS-r1.3"),
    "dere": ("Dere", "derer1.05", "dere_all_r1.05", "dere-all-CDS-r1.05"),
}
EXCLUDED = {"Dmel_w1118_rep1", "Dmel_w1118_rep2"}
PARAMETERS = {
    "minimum_coding_codons_before_trimming": 150,
    "exclude_first_coding_codons": 15,
    "exclude_last_coding_codons": 5,
    "minimum_interior_RPF_coverage_rate": 0.5,
    "RPF_coverage_definition": "among interior codons with unscaled RNA >= 5, fraction with A-site RPF count > 0; missing RPF counts are zero",
    "minimum_RNA_count_per_codon": 5,
    "RNA_scale": 1,
    "RNA_median_filter": None,
    "minimum_76_RNA_covered_codons_filter": False,
    "new_longest_transcript_selection": False,
    "transcript_universe": "original unfiltered 1codon annotations and original raw RNA/A-site count files; no new isoform selection or pooling",
    "m6A_exclusion": False,
    "editing_exclusion": False,
    "normalization": "mRC=A-site/RNA; gmRC=mRC/arithmetic mean(mRC) over interior codons with unscaled RNA >= 5, independently per transcript and library",
    "group_mean": "equal-weight mean of finite library gmRC at each transcript/codon, including zeros and omitting missing values",
    "excluded_libraries": sorted(EXCLUDED),
}
LEGACY_HEADER = [
    "FBtr",
    "pos",
    "chr",
    "no",
    "strand",
    "cdsposstart",
    "codonpos",
    "mRNA",
    "mRNA_median",
    "Asite",
    "Asite_median",
    "Esite",
    "Esite_median",
    "Psite",
    "Psite_median",
    "gmRC",
    "FBtr_codonpos",
    "codon",
    "AA",
]
EXTRA_HEADER = [
    "gene_id",
    "mRC",
    "cds_coding_codons",
    "interior_codon_n",
    "RPF_coverage_rate",
    "normalization_codon_n",
    "normalization_mean_mRC",
]
HEADER = LEGACY_HEADER + EXTRA_HEADER
AA_NAMES = dict(
    zip(
        "ACDEFGHIKLMNPQRSTVWY*",
        [
            "Ala",
            "Cys",
            "Asp",
            "Glu",
            "Phe",
            "Gly",
            "His",
            "Ile",
            "Lys",
            "Leu",
            "Met",
            "Asn",
            "Pro",
            "Gln",
            "Arg",
            "Ser",
            "Thr",
            "Val",
            "Trp",
            "Tyr",
            "Stop",
        ],
    )
)
AA_ORDER = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
AA = {c: AA_NAMES[a] for c, a in zip(CODONS, AA_ORDER)}
REF_CACHE = {}


def tsv_write(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def json_write(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")


def fmt(x):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "NA"
    return format(x, ".15g") if isinstance(x, float) else str(x)


def fingerprint(path):
    p = Path(path)
    return {"path": str(p), "bytes": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}


def load_fasta(path):
    result, names, chunks = {}, [], []

    def save():
        seq = "".join(chunks).upper()
        for tx in names:
            if tx in result and result[tx] != seq:
                raise ValueError(f"Conflicting CDS sequences for {tx}")
            result[tx] = seq

    with Path(path).open() as f:
        for line in f:
            if line.startswith(">"):
                save()
                names, chunks = re.findall(r"FBtr\d+", line), []
            else:
                chunks.append(line.strip())
        save()
    return result


def build_reference(species):
    label, tag, txprefix, fasta_prefix = SPECIES[species]
    annotation = PEAK / f"{tag}.gtf.cds.bed.ordered.1codon"
    fasta = PROJECT / "ref_genome" / f"{fasta_prefix}.fasta"
    txinfo = PROJECT / "05.psite" / f"{txprefix}_txinfo.biotype.tsv"
    genes = {r["tx_name"]: r["gene_id"] for r in read_tsv(txinfo)}
    seqs = load_fasta(fasta)
    refs = {}
    with annotation.open() as f:
        for line in f:
            chrom, pos, tx, no, strand, cds_start, codon_i = line.split()
            i = int(codon_i)
            if tx not in refs:
                refs[tx] = {
                    "gene": genes.get(tx, "NA"),
                    "chrom": chrom,
                    "strand": strand,
                    "no": no,
                    "positions": array("q"),
                    "seq": seqs.get(tx, ""),
                }
            r = refs[tx]
            if i != len(r["positions"]) + 1 or int(cds_start) != 3 * i - 2:
                raise ValueError(
                    f"Noncontiguous original codon annotation: {species} {tx} {i}"
                )
            r["positions"].append(int(pos))
    audit = []
    for tx, r in refs.items():
        if not r["seq"]:
            raise ValueError(f"Missing CDS sequence: {species} {tx}")
        original_n = len(r["positions"])
        # The original annotations normally omit stop. Remove it only if a
        # terminal annotated triplet is explicitly a stop in the reference.
        last_triplet = r["seq"][3 * (original_n - 1) : 3 * original_n]
        removed_stop = last_triplet in {"TAA", "TAG", "TGA"}
        if removed_stop:
            r["positions"].pop()
        r["n"] = len(r["positions"])
        seq_n = len(r["seq"]) // 3 - int(r["seq"][-3:] in {"TAA", "TAG", "TGA"})
        audit.append(
            {
                "FBtr": tx,
                "gene_id": r["gene"],
                "annotated_codon_n": original_n,
                "coding_codon_n_used": r["n"],
                "fasta_length_nt": len(r["seq"]),
                "fasta_complete_coding_triplets": seq_n,
                "terminal_annotated_stop_removed": int(removed_stop),
                "annotation_fasta_length_difference": r["n"] - seq_n,
            }
        )
    cache = ROOT / "reference" / f"{species}.pickle"
    with cache.open("wb") as f:
        pickle.dump(refs, f, protocol=5)
    tsv_write(ROOT / "reference" / f"{species}.transcript_lengths.tsv", audit)
    return {
        "species": species,
        "transcripts": len(refs),
        "annotation_fasta_length_discrepancies": sum(
            x["annotation_fasta_length_difference"] != 0 for x in audit
        ),
        "reference_inputs": [fingerprint(p) for p in [annotation, fasta, txinfo]],
    }


def sample_name(s):
    return (
        s.replace("Dsim_w501_", "Dsim_")
        .replace("Dmel_ISO1_female", "Dmel_ISO1_head_female")
        .replace("Dmel_ISO1_male", "Dmel_ISO1_head_male")
    )


def prepare():
    for d in [
        "reference",
        "per_library",
        "gene_qc",
        "group_sex_mean",
        "manifests",
        "logs",
        "glm_inputs",
    ]:
        (ROOT / d).mkdir(parents=True, exist_ok=True)
    if (ROOT / "manifests" / "parameters.json").exists():
        if (
            json.loads((ROOT / "manifests" / "parameters.json").read_text())
            != PARAMETERS
        ):
            raise ValueError(
                "Existing parameters differ; use a fresh output directory."
            )
    json_write(ROOT / "manifests" / "parameters.json", PARAMETERS)
    for sp in SPECIES:
        build_reference(sp)
    design = {
        r["sample_id"]: r for r in read_tsv(REPO / "code/05/config/stalling_groups.tsv")
    }
    samples = []
    for sp in SPECIES:
        name = "samples.dmel.all.tsv" if sp == "dmel" else f"samples.{sp}.tsv"
        for r in read_tsv(REPO / "code/05/config" / name):
            s = r["sample_id"]
            if s in EXCLUDED:
                continue
            info = design[s]
            label = SPECIES[sp][0]
            cohort = f'Dmel_{info["genotype"]}' if sp == "dmel" else label
            sex = {"F": "female", "M": "male"}[info["sex"]]
            r.update(
                sample=sample_name(s),
                species=sp,
                cohort=cohort,
                sex=sex,
                replicate=info["replicate"],
                is_head=info["is_head"],
                unit_id=f"{cohort}__{sex}",
            )
            r["output"] = str(ROOT / "per_library" / f'{r["sample"]}.codonTE_AEP.AA.gz')
            for col in ["rna_count", "rpf_count"]:
                p = Path(r[col])
                if not p.is_file():
                    raise FileNotFoundError(p)
                r[col + "_bytes"] = p.stat().st_size
                r[col + "_mtime_ns"] = p.stat().st_mtime_ns
            samples.append(r)
    if len(samples) != 44:
        raise AssertionError(f"Expected 44 libraries, got {len(samples)}")
    tsv_write(ROOT / "manifests" / "libraries.tsv", samples)
    groups = defaultdict(list)
    for r in samples:
        groups[r["unit_id"]].append(r["sample"])
        if r["cohort"] in {"Dmel_ISO1", "Dmel_WBg"}:
            groups["Dmel__" + r["sex"]].append(r["sample"])
    tsv_write(
        ROOT / "manifests" / "groups.tsv",
        [
            {
                "unit_id": unit,
                "samples": ",".join(ss),
                "n_libraries": len(ss),
                "output": str(ROOT / "group_sex_mean" / f"{unit}.mean_gmRC.tsv.gz"),
            }
            for unit, ss in sorted(groups.items())
        ],
    )
    print(
        f"Prepared {len(samples)} libraries and {len(groups)} group-sex means.",
        flush=True,
    )


def get_reference(species):
    if species not in REF_CACHE:
        with (ROOT / "reference" / f"{species}.pickle").open("rb") as f:
            REF_CACHE[species] = pickle.load(f)
    return REF_CACHE[species]


def load_counts(path, refs, allowed=None):
    counts = {}
    stats = Counter()
    with Path(path).open() as f:
        for line in f:
            a = line.split()
            if len(a) not in (3, 4):
                raise ValueError(f"Unexpected raw-count format: {path}: {line[:120]}")
            tx = a[0]
            if (
                tx not in refs
                or refs[tx]["n"] < 150
                or (allowed is not None and tx not in allowed)
            ):
                continue
            if not a[1].isdigit():
                continue
            p = int(a[1])
            if p < 1 or (p - 1) % 3:
                continue
            i = (p - 1) // 3
            n = refs[tx]["n"]
            if not 15 <= i < n - 5:
                continue
            value = float(a[-1])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid count in {path}: {value}")
            if not value:
                continue
            if tx not in counts:
                counts[tx] = array("d", [0.0]) * (n - 20)
            counts[tx][i - 15] += value
            stats["accepted_variant_rows"] += 1
    return counts, dict(stats)


def compute_transcript(rna, rpf):
    """Evaluate RPF coverage and normalize only on internal codons with RNA >= 5."""
    eligible = [i for i, x in enumerate(rna) if x >= 5]
    if not eligible:
        return 0.0, eligible, 0.0, "no_interior_codon_with_RNA_ge5", {}
    rate = sum(rpf[i] > 0 for i in eligible) / len(eligible)
    ratios = [rpf[i] / rna[i] for i in eligible]
    mean = math.fsum(ratios) / len(ratios)
    reason = "retained"
    if rate < 0.5:
        reason = "RPF_coverage_rate_below_0.5"
    elif mean <= 0:
        reason = "zero_mean_mRC_on_RNA_eligible_codons"
    gmrc = (
        {i: v / mean for i, v in zip(eligible, ratios)} if reason == "retained" else {}
    )
    return rate, eligible, mean, reason, gmrc


def run_library(info):
    sample, sp = info["sample"], info["species"]
    completed = ROOT / "gene_qc" / f"{sample}.summary.json"
    if completed.exists():
        return json.loads(completed.read_text())
    print(f"START {sample}", flush=True)
    refs = get_reference(sp)
    rpfs, rpf_load = load_counts(info["rpf_count"], refs)
    rnas, rna_load = load_counts(info["rna_count"], refs)
    out = Path(info["output"])
    tmp = out.with_suffix(".partial")
    qc_rows, stats = [], Counter()
    with gzip.open(tmp, "wt", compresslevel=1, newline="") as handle:
        handle.write("\t".join(HEADER) + "\n")
        for tx in sorted(refs):
            ref = refs[tx]
            n = ref["n"]
            q = {
                "sample": sample,
                "gene_id": ref["gene"],
                "FBtr": tx,
                "cds_coding_codons": n,
                "interior_codon_n": max(n - 20, 0),
                "RNA_ge5_interior_codons": 0,
                "RPF_positive_among_RNA_ge5_codons": 0,
                "RPF_coverage_rate": 0,
                "normalization_codon_n": 0,
                "normalization_mean_mRC": "NA",
                "reason": "",
            }
            if n < 150:
                q["reason"] = "CDS_fewer_than_150_coding_codons"
            else:
                rpf = rpfs.get(tx, array("d", [0.0]) * (n - 20))
                rna = rnas.get(tx, array("d", [0.0]) * (n - 20))
                rate, eligible, mean, reason, gmrc = compute_transcript(rna, rpf)
                q.update(
                    RNA_ge5_interior_codons=len(eligible),
                    RPF_positive_among_RNA_ge5_codons=sum(rpf[i] > 0 for i in eligible),
                    RPF_coverage_rate=rate,
                    reason=reason,
                    normalization_codon_n=len(eligible),
                    normalization_mean_mRC=mean,
                )
                if reason == "retained":
                    med_rna = median(rna[i] for i in eligible)
                    med_rpf = median(rpf[i] for i in eligible)
                    # Write all interior codons for retained transcripts so downstream
                    # coverage sums retain low-RNA observations; gmRC is NA below RNA 5.
                    for j in range(n - 20):
                        i = j + 16
                        codon = ref["seq"][3 * (i - 1) : 3 * i]
                        if len(codon) != 3:
                            raise ValueError(
                                f"Incomplete internal reference codon: {sp} {tx} {i}"
                            )
                        row = [
                            tx,
                            ref["positions"][i - 1],
                            ref["chrom"],
                            ref["no"],
                            ref["strand"],
                            3 * i - 2,
                            i,
                            rna[j],
                            med_rna,
                            rpf[j],
                            med_rpf,
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            gmrc.get(j),
                            f"{tx}_{i}",
                            codon,
                            AA.get(codon, "Xaa"),
                            ref["gene"],
                            rpf[j] / rna[j] if j in gmrc else None,
                            n,
                            n - 20,
                            rate,
                            len(eligible),
                            mean,
                        ]
                        handle.write("\t".join(fmt(v) for v in row) + "\n")
                    stats["output_rows"] += n - 20
                    stats["finite_gmRC_rows"] += len(gmrc)
                    stats["zero_gmRC_rows"] += sum(v == 0 for v in gmrc.values())
            stats[q["reason"]] += 1
            qc_rows.append(q)
    tmp.replace(out)
    tsv_write(ROOT / "gene_qc" / f"{sample}.transcript_qc.tsv", qc_rows)
    summary = {
        "sample": sample,
        "species": sp,
        "unit_id": info["unit_id"],
        **dict(stats),
        "raw_count_loading": {"RNA": rna_load, "RPF": rpf_load},
        "output": str(out),
    }
    json_write(completed, summary)
    print(
        f'DONE {sample}: {stats["retained"]} transcripts, {stats["finite_gmRC_rows"]} finite gmRC',
        flush=True,
    )
    return summary


def iter_rows(path):
    with gzip.open(path, "rt") as f:
        next(f)
        for line in f:
            a = line.rstrip("\n").split("\t")
            yield (a[0], int(a[6])), a


GROUP_HEADER = [
    "FBtr_codonpos",
    "FBtr",
    "codon",
    "AA",
    "gmRC",
    "n_gmRC_reps_finite",
    "n_libraries_expected",
    "n_libraries_with_transcript",
    "sample_ids_finite",
    "codonpos",
    "gene_id",
    "chr",
    "pos",
    "strand",
    "cds_coding_codons",
    "mRNA_sum",
    "mRNA_mean_over_present_libraries",
    "Asite_sum",
    "Asite_any_positive",
    "RPF_coverage_rate_min_present",
    "gmRC_rep_min",
    "gmRC_rep_max",
]


def run_group(info):
    unit = info["unit_id"]
    completed = ROOT / "group_sex_mean" / f"{unit}.summary.json"
    if completed.exists():
        return json.loads(completed.read_text())
    samples = info["samples"].split(",")
    iters = [
        iter(iter_rows(ROOT / "per_library" / f"{s}.codonTE_AEP.AA.gz"))
        for s in samples
    ]
    heap = []

    def push(k):
        try:
            key, row = next(iters[k])
            heapq.heappush(heap, (key, k, row))
        except StopIteration:
            pass

    for k in range(len(iters)):
        push(k)
    output = Path(info["output"])
    tmp = output.with_suffix(".partial")
    counts = Counter()
    with gzip.open(tmp, "wt", compresslevel=1) as f:
        f.write("\t".join(GROUP_HEADER) + "\n")
        while heap:
            key = heap[0][0]
            entries = []
            while heap and heap[0][0] == key:
                _, k, row = heapq.heappop(heap)
                entries.append((k, row))
                push(k)
            finite = [
                (k, float(r[15]))
                for k, r in entries
                if r[15] != "NA" and math.isfinite(float(r[15]))
            ]
            if not finite:
                counts["codons_without_finite_gmRC"] += 1
                continue
            r = entries[0][1]
            if any(
                (v[17], v[18], v[19], v[21]) != (r[17], r[18], r[19], r[21])
                for _, v in entries
            ):
                raise ValueError(f"Conflicting annotations: {unit} {key}")
            vals = [v for _, v in finite]
            mrna_sum = math.fsum(float(v[7]) for _, v in entries)
            row = [
                r[16],
                r[0],
                r[17],
                r[18],
                math.fsum(vals) / len(vals),
                len(vals),
                len(samples),
                len(entries),
                ",".join(samples[k] for k, _ in finite),
                r[6],
                r[19],
                r[2],
                r[1],
                r[4],
                r[21],
                mrna_sum,
                mrna_sum / len(entries),
                math.fsum(float(v[9]) for _, v in entries),
                int(any(float(v[10]) > 0 for _, v in entries)),
                min(float(v[23]) for _, v in entries),
                min(vals),
                max(vals),
            ]
            f.write("\t".join(fmt(v) for v in row) + "\n")
            counts["output_codons"] += 1
            counts[f"codons_with_{len(vals)}_finite_libraries"] += 1
    tmp.replace(output)
    summary = {
        "unit_id": unit,
        "samples": samples,
        **dict(counts),
        "output": str(output),
    }
    json_write(completed, summary)
    print(f'MEAN {unit}: {counts["output_codons"]} codons', flush=True)
    return summary


def execute(jobs):
    samples = read_tsv(ROOT / "manifests" / "libraries.tsv")
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        results = [
            f.result()
            for f in as_completed([pool.submit(run_library, r) for r in samples])
        ]
    fields = [
        "sample",
        "species",
        "unit_id",
        "retained",
        "output_rows",
        "finite_gmRC_rows",
        "zero_gmRC_rows",
        "CDS_fewer_than_150_coding_codons",
        "RPF_coverage_rate_below_0.5",
        "no_interior_codon_with_RNA_ge5",
        "zero_mean_mRC_on_RNA_eligible_codons",
        "output",
    ]
    tsv_write(
        ROOT / "manifests" / "library_results.tsv",
        [
            {k: r.get(k, 0) for k in fields}
            for r in sorted(results, key=lambda x: x["sample"])
        ],
        fields,
    )
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        means = [
            f.result()
            for f in as_completed(
                [
                    pool.submit(run_group, r)
                    for r in read_tsv(ROOT / "manifests" / "groups.tsv")
                ]
            )
        ]
    tsv_write(
        ROOT / "manifests" / "group_results.tsv",
        [
            {
                "unit_id": r["unit_id"],
                "n_libraries": len(r["samples"]),
                "output_codons": r["output_codons"],
                "output": r["output"],
            }
            for r in sorted(means, key=lambda x: x["unit_id"])
        ],
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=["prepare", "run"])
    p.add_argument("--jobs", type=int, default=4)
    a = p.parse_args()
    if a.stage == "prepare":
        prepare()
    else:
        execute(a.jobs)


if __name__ == "__main__":
    main()
