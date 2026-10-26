#!/usr/bin/env python3
"""Generate the Table S10 candidate evidence; optionally retain historical summaries."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

SPECIES = ("mel", "sim", "sec", "ere")
LABEL = {
    "mel": "D. melanogaster",
    "sim": "D. simulans",
    "sec": "D. sechellia",
    "ere": "D. erecta",
}
# One biological baseline cohort per species; AG/AS are not pooled into mel.
BASELINE = {"mel": "Dmel_WT", "sim": "Dsim", "sec": "Dsec", "ere": "Dere"}
PROJECTION_FILE = {
    "mel": "ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
    "sim": "sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
    "sec": "sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
    "ere": "ere.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
}


def read_tsv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        out = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def protein_effect(annotation):
    """Return a p.Xaa123Yaa consequence from a SnpEff ANN field, if present."""
    hit = re.search(r"(?:^|[|])p\.([A-Za-z*]+\d+[A-Za-z*]+)(?:[|]|$)", annotation or "")
    return hit.group(1) if hit else ""


def theoretical_ms_candidates(pipeline_root, projection_root):
    """Require two distinct single-edit missense products per double-edit codon."""
    catalog = read_tsv(pipeline_root / "01_counts" / "double_edit_codons.tsv")
    expected = {}
    for row in catalog:
        positions = tuple(sorted(int(x) for x in row["edit_cds_positions"].split(";")))
        expected[(row["species"], row["target_transcript"], positions)] = row

    effects = defaultdict(set)
    for species, filename in PROJECTION_FILE.items():
        for row in read_tsv(projection_root / filename):
            if (
                row.get("target_region") != "CDS"
                or row.get("change2", "").upper() != "A>G"
            ):
                continue
            if row.get("type") != "missense_variant":
                continue
            effect = protein_effect(row.get("ann", ""))
            if not effect:
                continue
            tx, pos = row.get("new_transcript", ""), row.get("new_cds_pos_1based", "")
            if tx and pos:
                effects[(species, tx, int(pos))].add(effect)

    selected, audit = {s: set() for s in SPECIES}, []
    for (species, tx, positions), row in expected.items():
        effect_sets = [effects[(species, tx, pos)] for pos in positions]
        # A position may be repeated across RNA libraries, but its missense
        # product must be unambiguous. Both single-edit products must differ.
        valid = all(len(x) == 1 for x in effect_sets)
        effect1 = (
            next(iter(effect_sets[0]))
            if len(effect_sets[0]) == 1
            else ";".join(sorted(effect_sets[0]))
        )
        effect2 = (
            next(iter(effect_sets[1]))
            if len(effect_sets[1]) == 1
            else ";".join(sorted(effect_sets[1]))
        )
        distinct = valid and effect1 != effect2
        if distinct:
            selected[species].add(row["target_handle"])
        audit.append(
            {
                "species": LABEL[species],
                "target_handle": row["target_handle"],
                "gene": row["gene"],
                "edit_cds_positions": row["edit_cds_positions"],
                "single_edit_effect_1": effect1,
                "single_edit_effect_2": effect2,
                "both_missense_and_distinct": "yes" if distinct else "no",
            }
        )
    return selected, audit


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pipeline-root", required=True)
    p.add_argument("--projection-root", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument(
        "--candidates-only",
        action="store_true",
        help="Generate only the canonical Table S10 theoretical candidate evidence",
    )
    args = p.parse_args()
    root, projection_root, outdir = (
        Path(args.pipeline_root),
        Path(args.projection_root),
        Path(args.outdir),
    )
    outdir.mkdir(parents=True, exist_ok=True)

    if args.candidates_only:
        _, audit = theoretical_ms_candidates(root, projection_root)
        write_tsv(
            outdir / "theoretical_MS_candidates.tsv",
            [
                "species",
                "target_handle",
                "gene",
                "edit_cds_positions",
                "single_edit_effect_1",
                "single_edit_effect_2",
                "both_missense_and_distinct",
            ],
            audit,
        )
        return

    qc = read_tsv(root / "01_counts" / "candidate_QC.tsv")
    all_candidates = {r["species"]: int(r["double_edit_codons"]) for r in qc}

    scores = read_tsv(root / "02_normalized" / "double_edit_scores.cohort_mean.tsv")
    status = {
        s: {
            "rna": set(),
            "ribo": set(),
            "eligible": {},
            "aa_neg": set(),
            "ii_pos": set(),
        }
        for s in SPECIES
    }
    for r in scores:
        s = r["species"]
        if s not in status or r["cohort"] != BASELINE[s]:
            continue
        h = r["target_handle"]
        if r["RNA_edited_detected_in_valid_samples"] == "yes":
            status[s]["rna"].add(h)
        if r["RPF_edited_detected_in_valid_samples"] == "yes":
            status[s]["ribo"].add(h)
        if r["heatmap_eligible"] == "yes":
            status[s]["eligible"][h] = r
            if float(r["AA_mean_log2_RPF_RNA"]) < 0:
                status[s]["aa_neg"].add(h)
            if float(r["II_mean_log2_RPF_RNA"]) > 0:
                status[s]["ii_pos"].add(h)

    ms_candidate, candidate_audit = theoretical_ms_candidates(root, projection_root)

    ms_rows = read_tsv(root / "03_ms" / "MS_allele_summary.all.tsv")
    ms_total, ms_edited = defaultdict(set), defaultdict(set)
    for r in ms_rows:
        key = (r["species"], r["target_handle"])
        if int(r["detected_sample_count"] or 0) > 0:
            ms_total[r["species"]].add(r["target_handle"])
            if r["allele"] != "AA":
                ms_edited[r["species"]].add(r["target_handle"])

    table_rows, direction_rows = [], []
    for s in SPECIES:
        total_detected = len(ms_total[s] & ms_candidate[s])
        edited_detected = len(ms_edited[s] & ms_candidate[s])
        table_rows.append(
            {
                "species": LABEL[s],
                "All candidates": all_candidates.get(s, 0),
                "RNA edited": len(status[s]["rna"]),
                "Ribo edited": len(status[s]["ribo"]),
                "MS candidate": len(ms_candidate[s]),
                "MS detected (edited/total)": f"{edited_detected}/{total_detected}",
            }
        )
        n = len(status[s]["eligible"])
        direction_rows.append(
            {
                "species": LABEL[s],
                "coverage-qualified multi-editing codons": n,
                "AA negative log2(Ribo/RNA)": f"{len(status[s]['aa_neg'])}/{n}",
                "II positive log2(Ribo/RNA)": f"{len(status[s]['ii_pos'])}/{n}",
            }
        )

    write_tsv(
        outdir / "Table_S11_multiediting_codon_counts.tsv",
        [
            "species",
            "All candidates",
            "RNA edited",
            "Ribo edited",
            "MS candidate",
            "MS detected (edited/total)",
        ],
        table_rows,
    )
    write_tsv(
        outdir / "Table_S11_AA_II_direction_summary.tsv",
        [
            "species",
            "coverage-qualified multi-editing codons",
            "AA negative log2(Ribo/RNA)",
            "II positive log2(Ribo/RNA)",
        ],
        direction_rows,
    )
    write_tsv(
        outdir / "theoretical_MS_candidates.tsv",
        [
            "species",
            "target_handle",
            "gene",
            "edit_cds_positions",
            "single_edit_effect_1",
            "single_edit_effect_2",
            "both_missense_and_distinct",
        ],
        candidate_audit,
    )
    retained = ", ".join(str(len(status[s]["eligible"])) for s in SPECIES)
    aa = ", ".join(
        f"{len(status[s]['aa_neg'])}/{len(status[s]['eligible'])}" for s in SPECIES
    )
    ii = ", ".join(
        f"{len(status[s]['ii_pos'])}/{len(status[s]['eligible'])}" for s in SPECIES
    )
    text = (
        "After requiring sufficient RNA-seq and Ribo-seq coverage, we identified "
        f"{retained} multi-editing codons in D. melanogaster, D. simulans, "
        "D. sechellia, and D. erecta, respectively (Table S11). Across these "
        "multi-editing codons, doubly edited (II) alleles consistently showed "
        "increased representation in Ribo-seq relative to RNA-seq, whereas "
        "unedited (AA) alleles showed the opposite trend. AA allele "
        f"log2(Ribo/RNA) values were negative in {aa} codons, whereas II allele "
        f"values were positive in {ii} codons across the four species, respectively."
    )
    pass


if __name__ == "__main__":
    main()
