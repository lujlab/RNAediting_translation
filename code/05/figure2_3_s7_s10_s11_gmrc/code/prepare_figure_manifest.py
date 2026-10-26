#!/usr/bin/env python3
"""Generate legacy-compatible figure manifest and exact new CDS lengths."""

import recompute as r
from pathlib import Path

rows = []
for g in r.read_tsv(r.ROOT / "manifests/groups.tsv"):
    unit = g["unit_id"]
    species = unit[:4]
    sex = unit.split("__")[1]
    main = unit in [
        f"{s}__{x}"
        for s in ["Dmel", "Dsim", "Dsec", "Dere"]
        for x in ["female", "male"]
    ]
    if not main and unit != "Dmel_Adar5G1__male":
        continue
    rows.append(
        dict(
            species=species,
            cohort=unit.split("__")[0],
            sex=sex,
            unit_id=unit,
            analysis_role="main_8" if main else "Adar5G1_reference",
            n_input_replicates=g["n_libraries"],
            input_samples=g["samples"],
            mean_gmRC_file=g["output"],
        )
    )
r.tsv_write(r.ROOT / "group_sex_mean/mean_gmRC_table_manifest.tsv", rows)
refs = []
old = r.read_tsv(r.REPO / "code/05/config/metagene_references.tsv")
for x in old:
    species = x["species"]
    dest = r.ROOT / "reference" / f"{species}.metagene_lengths.tsv"
    with dest.open("w") as f:
        for z in r.read_tsv(
            r.ROOT / "reference" / f"{species.lower()}.transcript_lengths.tsv"
        ):
            # Existing metagene parser expects a length including one stop codon.
            # Synthetic coordinates encode the exact new coding N, without using old lengths.
            f.write(f"{z['FBtr']}\t1\t{3*(int(z['coding_codon_n_used'])+1)}\tCDS\n")
    refs.append({**x, "trans_file": str(dest)})
r.tsv_write(r.ROOT / "reference/metagene_references.tsv", refs)
print("Prepared 9 analysis units and exact coding-length references.")
