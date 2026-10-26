#!/usr/bin/env python3
"""DIA-only four-state audit and EndoA figure inputs; original sources unchanged."""

import os
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

OUT = Path(
    os.environ.get(
        "FS14_OUTPUT_DIR",
        str(Path(__file__).resolve().parents[3] / "results/06/FS14_F4B"),
    )
)
OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path(os.environ.get("MS_PROJECT_ROOT", "/path/to/MS_analysis_project"))
PROJECT = Path(os.environ.get("EDITING_PROJECT_ROOT", str(ROOT / "editing_project")))
COUNTS = Path(
    os.environ.get(
        "DOUBLE_EDIT_COUNTS",
        str(
            PROJECT
            / "04_alledit_new/02_double_edit_and_MS/01_counts/double_edit_allele_counts.by_sample.tsv"
        ),
    )
)
MS = Path(os.environ.get("DIA_NN_ANALYSIS_DIR", str(ROOT / "03_MS/DIA_NN_analysis")))
STATES = ["AA", "IA", "AI", "II"]
SITES = [
    ("Dmel", "EndoA", "FBtr0083698:137", "FBtr0334522_409"),
    ("Dmel", "alpha-Spec", "FBtr0332878:2101", "FBtr0332878_6301"),
    ("Dsim", "EndoA", "FBtr0361051:137", "FBtr0361051_409"),
    ("Dsim", "alpha-Spec", "FBtr0362198:521", "FBtr0362198_1561"),
    ("Dsim", "Rlip", "FBtr0219301:233", "FBtr0219301_697"),
    ("Dsec", "EndoA", "FBtr0200880:137", "FBtr0200880_409"),
    ("Dsec", "alpha-Spec", "FBtr0197154:2101", "FBtr0197154_6301"),
    ("Dsec", "Rlip", "FBtr0198043:233", "FBtr0198043_697"),
    ("Dere", "EndoA", "FBtr0143136:137", "FBtr0143136_409"),
]
META = {r[3]: r for r in SITES}
ALIAS = {r[2]: r[3] for r in SITES}


def read(p):
    with open(p, newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")


def write(name, rows):
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


samples = defaultdict(lambda: dict.fromkeys(STATES, 0.0))
selected = []
for r in read(COUNTS):
    h = r["target_handle"]
    if h not in META:
        continue
    if r["species"] == "mel" and r["genotype"] not in ["WBg"]:
        continue
    if r["species"] != "mel" and "_head_" in r["sample"]:
        continue
    selected.append(r)
    for assay in ["RNA", "RPF"]:
        samples[h, r["sex"], assay, r["sample"]][r["allele"]] = float(
            r[assay + "_count"]
        )
write("RNA_Ribo_source_counts.tsv", selected)
sources = [COUNTS]
seen = set()
ms_members = defaultdict(set)
for sp, stem in [
    ("Dmel", "AG_AS_WBG_iso1"),
    ("Dsim", "dsim"),
    ("Dsec", "dsec"),
    ("Dere", "dere"),
]:
    path = MS / (stem + "_site_level_intensity_100AA.v2.tsv")
    sources.append(path)
    for r in read(path):
        if r["site"] not in ALIAS:
            continue
        sample = Path(r["sample"]).name
        m = re.search(r"_(iso|WBg|AG|AS|sim|sec|ere)_([FM])([123])\.raw$", sample)
        assert m, sample
        geno, sex, rep = m.groups()
        if sp == "Dmel" and geno not in ["WBg"]:
            continue
        h = ALIAS[r["site"]]
        gene = META[h][1]
        state = (
            {"N": "AA", "D": "IA", "S": "AI", "G": "II"}
            if gene == "alpha-Spec"
            else {"K": "AA", "E": "IA", "R": "AI", "G": "II"}
        )[r["AA"]]
        k = (h, sex, "MS", sample)
        assert (k, state) not in seen
        seen.add((k, state))
        samples[k][state] += float(r["sum_intensity"])
        ms_members[sp, geno, sex].add(rep)
for geno in ["WBg"]:
    for sex in ["F", "M"]:
        assert ms_members["Dmel", geno, sex] == {"1", "2", "3"}
sample_rows = []
bucket = defaultdict(list)
for (h, sex, assay, sample), counts in samples.items():
    total = sum(counts.values())
    r = dict(
        species=META[h][0],
        gene=META[h][1],
        target_handle=h,
        sex=sex,
        assay=assay,
        sample=sample,
        total=total,
    )
    for a in STATES:
        r[a + "_value"] = counts[a]
        r[a + "_percent"] = 100 * counts[a] / total if total > 0 else "NA"
    if total > 0:
        bucket[h, sex, assay].append([r[a + "_percent"] for a in STATES])
    sample_rows.append(r)
write("per_sample_counts_intensities_percentages.tsv", sample_rows)
numeric = []
formatted = []
plot = {}
for sp, gene, site, h in SITES:
    for sex in ["F", "M"]:
        r = dict(
            species=sp,
            group="Dmel (WT, WBg)" if sp == "Dmel" else sp,
            gene=gene,
            MS_site=site,
            target_handle=h,
            sex=sex,
        )
        f = dict(r)
        for assay in ["RNA", "RPF", "MS"]:
            b = bucket[h, sex, assay]
            means = [sum(v[i] for v in b) / len(b) for i in range(4)] if b else []
            assert not means or math.isclose(sum(means), 100, abs_tol=1e-7)
            r[assay + "_n"] = len(b)
            for i, a in enumerate(STATES):
                r[assay + "_" + a + "_percent"] = means[i] if b else "NA"
            f[assay + "_AA_IA_AI_II_percent"] = (
                " / ".join("{:.2f}".format(v) for v in means) if b else "NA"
            )
        numeric.append(r)
        formatted.append(f)
        if gene == "EndoA":
            plot[sp + "_" + sex] = dict(row=r)
write("RNA_Ribo_MS_numeric.tsv", numeric)
write("RNA_Ribo_MS_table.tsv", formatted)

# Per-state raw mRC = RPF/RNA; zero state RNA -> undefined, not zero.
# A sample without any observed RPF at this target has no observed Ribo state
# distribution and is excluded from the matched descriptive mRC summary.
mrc_samples = []
mrc_buckets = defaultdict(list)
for r in selected:
    h = r["target_handle"]
    if META[h][1] != "EndoA":
        continue
    rn = float(r["RNA_count"])
    rp = float(r["RPF_count"])
    valid = float(r["RNA_total"]) > 0 and float(r["RPF_total"]) > 0 and rn > 0
    value = rp / rn if valid else None
    z = dict(
        species=META[h][0],
        sex=r["sex"],
        sample=r["sample"],
        target_handle=h,
        state=r["allele"],
        RNA_count=rn,
        RPF_count=rp,
        RNA_total=r["RNA_total"],
        RPF_total=r["RPF_total"],
        mRC=value if valid else "NA",
        status=(
            "valid"
            if valid
            else (
                "no_target_RPF_coverage"
                if float(r["RPF_total"]) == 0
                else "zero_state_RNA"
            )
        ),
    )
    mrc_samples.append(z)
    if valid:
        mrc_buckets[META[h][0], r["sex"], r["allele"]].append(value)
write("EndoA_mRC_by_sample.tsv", mrc_samples)
mrc_summary = []
for sp in ["Dmel", "Dsim", "Dsec", "Dere"]:
    for sex in ["F", "M"]:
        vals = []
        for a in STATES:
            b = mrc_buckets[sp, sex, a]
            value = sum(b) / len(b) if b else None
            vals.append(value)
            mrc_summary.append(
                dict(species=sp, sex=sex, state=a, mRC=value if b else "NA", n=len(b))
            )
        plot[sp + "_" + sex]["mRC"] = vals
write("EndoA_mRC_mean.tsv", mrc_summary)
(OUT / "figure_data.json").write_text(json.dumps(plot, indent=2))
pass
for r in formatted:
    print(
        r["gene"],
        r["group"],
        r["sex"],
        r["RNA_AA_IA_AI_II_percent"],
        "|",
        r["RPF_AA_IA_AI_II_percent"],
        "|",
        r["MS_AA_IA_AI_II_percent"],
    )
print("DIA membership:", dict((str(k), sorted(v)) for k, v in ms_members.items()))
print("mRC:", dict((k, v["mRC"]) for k, v in plot.items()))
