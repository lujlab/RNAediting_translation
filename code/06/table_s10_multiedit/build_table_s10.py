#!/usr/bin/env python3
"""Separate any-edited MS support from target-codon II support, using full DIA tables."""

import os
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

OUT = Path(
    os.environ.get(
        "TABLE_S10_OUT_DIR",
        str(Path(__file__).resolve().parents[3] / "results/06_MS/TableS10"),
    )
)
OUT.mkdir(parents=True, exist_ok=True)
INPUT = Path(os.environ["TABLE_S10_INPUTS"])
MS = Path(os.environ["DIA_NN_ANALYSIS_DIR"])
NAMES = {"mel": "Dmel", "sim": "Dsim", "sec": "Dsec", "ere": "Dere"}
LONG = dict(zip(["D. melanogaster", "D. simulans", "D. sechellia", "D. erecta"], NAMES))
COHORT = {"mel": "Dmel_WT", "sim": "Dsim", "sec": "Dsec", "ere": "Dere"}
AA = dict(
    zip(
        "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val".split(),
        "ARNDCQEGHILKMFPSTWYV",
    )
)


def read(p):
    with open(p, newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")


def write(name, rows):
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


catalog = {}
aliases = defaultdict(set)
candidates = defaultdict(set)
ms_candidates = defaultdict(dict)
sources = [
    INPUT / n
    for n in [
        "double_edit_codons.tsv",
        "theoretical_MS_candidates.tsv",
        "double_edit_scores.cohort_mean.tsv",
    ]
]
for r in read(sources[0]):
    sp, h = r["species"], r["target_handle"]
    catalog[sp, h] = r
    candidates[sp].add(h)
    for a in r["ms_site_aliases"].split(";"):
        aliases[sp, a].add(h)
for r in read(sources[1]):
    if r["both_missense_and_distinct"] != "yes":
        continue
    sp = LONG[r["species"]]
    h = r["target_handle"]
    effects = [
        re.fullmatch(
            r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", r["single_edit_effect_" + str(i)]
        )
        for i in [1, 2]
    ]
    assert all(effects), (sp, h)
    ref = AA[effects[0][1]]
    first = AA[effects[0][3]]
    second = AA[effects[1][3]]
    assert len({ref, first, second}) == 3
    ms_candidates[sp][h] = (ref, first, second)
rna = defaultdict(set)
ribo = defaultdict(set)
for r in read(sources[2]):
    sp, h = r["species"], r["target_handle"]
    if r["cohort"] != COHORT[sp]:
        continue
    if r["RNA_edited_detected_in_valid_samples"] == "yes":
        rna[sp].add(h)
    if r["RPF_edited_detected_in_valid_samples"] == "yes":
        ribo[sp].add(h)
support = defaultdict(set)
evidence = []
for sp, stem in [
    ("mel", "AG_AS_WBG_iso1"),
    ("sim", "dsim"),
    ("sec", "dsec"),
    ("ere", "dere"),
]:
    path = MS / (stem + "_site_level_intensity_100AA.v2.tsv")
    sources.append(path)
    for r in read(path):
        sample = Path(r["sample"]).name
        if sp == "mel" and not re.search(r"_(iso|WBg)_[FM][123]\.raw$", sample):
            continue
        if float(r["sum_intensity"]) <= 0:
            continue
        handles = aliases[sp, r["site"]] & set(ms_candidates[sp])
        if not handles:
            continue
        assert len(handles) == 1, (sp, r["site"], handles)
        h = next(iter(handles))
        ref, first, second = ms_candidates[sp][h]
        aa = r["AA"]
        if aa == ref:
            state = "AA"
        elif aa == first:
            state = "IA"
        elif aa == second:
            state = "AI"
        else:
            # The only positive unlisted product in these candidates is EndoA K137G:
            # AAG -> GAG (E), AGG (R), GGG (G), uniquely assigning G to target II.
            assert (ref, first, second, aa) == ("K", "E", "R", "G"), (
                sp,
                h,
                ref,
                first,
                second,
                aa,
            )
            assert h in {
                "FBtr0334522_409",
                "FBtr0361051_409",
                "FBtr0200880_409",
                "FBtr0143136_409",
            }
            state = "II"
        support[sp, h].add(state)
        evidence.append(
            dict(
                species=NAMES[sp],
                target_handle=h,
                MS_site=r["site"],
                gene=r["gene"],
                state_CDS_order=state,
                refAA=ref,
                observedAA=aa,
                sample=sample,
                intensity=r["sum_intensity"],
                source=str(path),
            )
        )
pass
table = []
site_audit = []
for sp, name in NAMES.items():
    any_set = {h for h in ms_candidates[sp] if support[sp, h] - {"AA"}}
    ii_set = {h for h in ms_candidates[sp] if "II" in support[sp, h]}
    assert ii_set <= any_set
    table.append(
        {
            "Species": name,
            "Codons with two candidate editing sites": len(candidates[sp]),
            "RNA-supported edited codons (IA, AI or II)": len(rna[sp]),
            "Ribo-supported edited codons (IA, AI or II)": len(ribo[sp]),
            "MS candidate codons (distinct single-edit products)": len(
                ms_candidates[sp]
            ),
            "MS-supported edited codons (IA, AI or II)": len(any_set),
            "MS-supported doubly edited codons (II)": len(ii_set),
        }
    )
    for h in sorted(any_set):
        site_audit.append(
            dict(
                species=name,
                target_handle=h,
                gene=catalog[sp, h]["gene"],
                genomic_sites=catalog[sp, h]["genomic_sites"],
                MS_aliases=catalog[sp, h]["ms_site_aliases"],
                detected_states=";".join(
                    a for a in ["AA", "IA", "AI", "II"] if a in support[sp, h]
                ),
                MS_any_edited="yes",
                MS_II="yes" if h in ii_set else "no",
            )
        )
write("TableS10.tsv", table)
pass
assert [list(r.values())[-2:] for r in table] == [[2, 1], [3, 1], [3, 1], [2, 0]]
for r in table:
    print(r)
