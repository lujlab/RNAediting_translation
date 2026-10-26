#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rebuild RNA-editing site sets from ORIGINAL files using species-union candidate selection.

Requested logic implemented exactly:
  1. Adar5G1 raw ann.test2
     - derive change2 from column 10 and keep change2 == A>G
     - ignore depth/altcount
     - rep1 union rep2 -> Adar5G1_rawAtoG

  2. WBG raw ann.test2
     - derive change2 from column 10 and keep change2 == A>G
     - retain all A>G annotations until final representative annotation
     - consequence priority first, transcript feature priority second
     - keep sites present in at least 2 of 6 WBG libraries
     - subtract Adar5G1_rawAtoG
     - -> WBG_min2

  3. ISO1
     - use ISO1 at2 results as ISO1_union
     - subtract Adar5G1_rawAtoG

  4. Published
     - five published .site files + Duan2017 High_Confidence
     - -> Published_union

  5. Final Dmel UpSet
     - WBG_min2
     - ISO1_union
     - Published_union

  6. Candidate and group/species lists AFTER Dmel 3-set construction
     - candidate = ISO1_union AND (WBG_min2 OR Published_union)
     - the representative annotation of the master candidate list is selected
       exclusively from ISO1; WBG/AG/AS never overwrite it
     - scan WBG/AG/AS and sim/sec/ere only for candidate handles
     - for every group, retain exact transcript-oriented A>G annotations first,
       then prefer coding/transcript annotations
     - remove non-1v1 liftover mappings before species union/upset
     - WBG/AG/AS min2; sim/sec/ere union among candidate handles
     - output one annotated candidate list per group, UpSet membership tables,
       and heatmap-ready editing-level long/wide tables.

Fast implementation:
  - Python multiprocessing + line-by-line streaming.
  - No full raw file is loaded into memory.
  - Raw test2 handle = col1_col2.
  - Lifted species handle = col1_col3 (BED end).

Usage:
  python3 rebuild_candidate_corrected_WBG_change2_species_union.py OUTDIR THREADS N_UPSET ISO1_MODE SPECIES_MODE

Args:
  OUTDIR       output directory
  THREADS      parallel workers, default 8
  N_UPSET      top intersections for UpSetR, default 40
  ISO1_MODE    head or all, default all
               head = only 2 ISO1 head at2 files
               all  = brain at2 + 2 head at2 files, matching old script
  SPECIES_MODE dp5mq10 or test2, default test2 (no DP/MQ filtering)
               dp5mq10 uses w501/sec/ere ann_AG_dp5mq10 lifted files including old M4/F4
               test2 uses test2_A2G lifted files plus Dsim/Dsec/Dere head mRNA as available
"""

from __future__ import annotations

import csv
import os
import re
import sys
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Set, Tuple, Iterable, Optional

# -----------------------------
# Paths
# -----------------------------

# Shared with 03_Editing_site_confirmation.sh; paths match the accepted reproduction.
INPUT_ROOT = Path(
    os.environ.get(
        "INPUT_ROOT",
        str(
            Path(__file__).resolve().parents[3]
            / "04_reproduce_original_20260913/inputs"
        ),
    )
)

ISO1_BRAIN_AT2_FILES = {
    "F_control_1": f"{INPUT_ROOT}/old_brain/Dmel_female_control_first.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "F_control_2": f"{INPUT_ROOT}/old_brain/Dmel_female_control_second.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "F_14h_1": f"{INPUT_ROOT}/old_brain/Dmel_female_temperature_14h_first.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "F_14h_2": f"{INPUT_ROOT}/old_brain/Dmel_female_temperature_14h_second.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "F_48h_2": f"{INPUT_ROOT}/old_brain/Dmel_female_temperature_48h_second.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "M_control_1": f"{INPUT_ROOT}/old_brain/Dmel_male_control_first.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "M_control_2": f"{INPUT_ROOT}/old_brain/Dmel_male_control_second.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "M_control_3": f"{INPUT_ROOT}/old_brain/Dmel_male_control_third.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "M_14h_1": f"{INPUT_ROOT}/old_brain/Dmel_male_temperature_14h_first.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "M_14h_2": f"{INPUT_ROOT}/old_brain/Dmel_male_temperature_14h_second.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "M_14h_3": f"{INPUT_ROOT}/old_brain/Dmel_male_temperature_14h_third.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "M_48h_2": f"{INPUT_ROOT}/old_brain/Dmel_male_temperature_48h_second.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
}

ISO1_HEAD_AT2_FILES = {
    "ISO1_head_F": f"{INPUT_ROOT}/old_head/Dmel_ISO1_head_female_mRNA.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
    "ISO1_head_M": f"{INPUT_ROOT}/old_head/Dmel_ISO1_head_male_mRNA.uniq1.vcf.Filter.ann_AG_dp5mq10fdr0.05ratio0.01at2",
}

ADAR5G1_RAW_FILES = {
    "Adar5G1_R1": f"{INPUT_ROOT}/old_head/Dmel_Adar5G1_rep1_mRNA.uniq1.vcf.Filter.ann.test2",
    "Adar5G1_R2": f"{INPUT_ROOT}/old_head/Dmel_Adar5G1_rep2_mRNA.uniq1.vcf.Filter.ann.test2",
}

WBG_RAW_FILES = {
    "WBG_F1": f"{INPUT_ROOT}/new_snpeff/WBG-F1-LJD7606.uniq1.vcf.Filter.ann.test2",
    "WBG_F2": f"{INPUT_ROOT}/new_snpeff/WBG-F2-LJD7607.uniq1.vcf.Filter.ann.test2",
    "WBG_F3": f"{INPUT_ROOT}/new_snpeff/WBG-F3-LJD7608.uniq1.vcf.Filter.ann.test2",
    "WBG_M1": f"{INPUT_ROOT}/new_snpeff/WBG-M1-LJD7603.uniq1.vcf.Filter.ann.test2",
    "WBG_M2": f"{INPUT_ROOT}/new_snpeff/WBG-M2-LJD7604.uniq1.vcf.Filter.ann.test2",
    "WBG_M3": f"{INPUT_ROOT}/new_snpeff/WBG-M3-LJD7605.uniq1.vcf.Filter.ann.test2",
}

AG_RAW_FILES = {
    "AG_F1": f"{INPUT_ROOT}/new_snpeff/AG-F1-LJD7618.uniq1.vcf.Filter.ann.test2",
    "AG_F2": f"{INPUT_ROOT}/new_snpeff/AG-F2-LJD7619.uniq1.vcf.Filter.ann.test2",
    "AG_F3": f"{INPUT_ROOT}/new_snpeff/AG-F3-LJD7620.uniq1.vcf.Filter.ann.test2",
    "AG_M1": f"{INPUT_ROOT}/new_snpeff/AG-M1-LJD7615.uniq1.vcf.Filter.ann.test2",
    "AG_M2": f"{INPUT_ROOT}/new_snpeff/AG-M2-LJD7616.uniq1.vcf.Filter.ann.test2",
    "AG_M3": f"{INPUT_ROOT}/new_snpeff/AG-M3-LJD7617.uniq1.vcf.Filter.ann.test2",
}

AS_RAW_FILES = {
    "AS_F1": f"{INPUT_ROOT}/new_snpeff/AS-F1-LJD7612.uniq1.vcf.Filter.ann.test2",
    "AS_F2": f"{INPUT_ROOT}/new_snpeff/AS-F2-LJD7613.uniq1.vcf.Filter.ann.test2",
    "AS_F3": f"{INPUT_ROOT}/new_snpeff/AS-F3-LJD7614.uniq1.vcf.Filter.ann.test2",
    "AS_M1": f"{INPUT_ROOT}/new_snpeff/AS-M1-LJD7609.uniq1.vcf.Filter.ann.test2",
    "AS_M2": f"{INPUT_ROOT}/new_snpeff/AS-M2-LJD7610.uniq1.vcf.Filter.ann.test2",
    "AS_M3": f"{INPUT_ROOT}/new_snpeff/AS-M3-LJD7611.uniq1.vcf.Filter.ann.test2",
}

# A site is considered detected in a group when it occurs in at least this
# many independent libraries. Keeping this in one place makes the rule easy
# to change without touching the analysis code.
MIN_LIBRARY_SUPPORT = {
    "WBG": 2,
    "AG": 2,
    "AS": 2,
    "sim": 1,  # union across libraries
    "sec": 1,  # union across libraries
    "ere": 1,  # union across libraries
}

PUBLISHED_DIR = os.environ.get(
    "PUBLISHED_EDITING_DIR", "/path/to/editing_sites_PUBlist"
)
PUBLISHED_FILES = [
    f"{PUBLISHED_DIR}/2011.Graveley.site",
    f"{PUBLISHED_DIR}/2012.Rodriguez.site",
    f"{PUBLISHED_DIR}/2013.Laurant.site",
    f"{PUBLISHED_DIR}/2016.LiXuan.site",
    f"{PUBLISHED_DIR}/2017.Billy.site",
]
DUAN2017_FILE = os.environ.get("DUAN2017_FILE", "/path/to/2017edit_all5925.txt")

G_ALIGN = f"{INPUT_ROOT}/new_snpeff/g_align"

SPECIES_DP5MQ10 = {
    "sim": {
        "sim_F1": f"{G_ALIGN}/w501-F1-totalRNA-LJD7636-RL1.uniq1.vcf.Filter.ann_AG_dp5mq10.droSim2ToDm6.sites",
        "sim_F2": f"{G_ALIGN}/w501-F2-totalRNA-LJD7637-RL1.uniq1.vcf.Filter.ann_AG_dp5mq10.droSim2ToDm6.sites",
        "sim_F3": f"{G_ALIGN}/w501-F3-totalRNA-LJD7638-RL1.uniq1.vcf.Filter.ann_AG_dp5mq10.droSim2ToDm6.sites",
        "sim_F4": f"{G_ALIGN}/w501-F4-old.uniq1.vcf.Filter.ann_AG_dp5mq10.droSim2ToDm6.sites",
        "sim_M1": f"{G_ALIGN}/w501-M1-totalRNA-LJD7633-RL1.uniq1.vcf.Filter.ann_AG_dp5mq10.droSim2ToDm6.sites",
        "sim_M2": f"{G_ALIGN}/w501-M2-totalRNA-LJD7634-RL1.uniq1.vcf.Filter.ann_AG_dp5mq10.droSim2ToDm6.sites",
        "sim_M3": f"{G_ALIGN}/w501-M3-totalRNA-LJD7635-RL1.uniq1.vcf.Filter.ann_AG_dp5mq10.droSim2ToDm6.sites",
        "sim_M4": f"{G_ALIGN}/w501-M4-old.uniq1.vcf.Filter.ann_AG_dp5mq10.droSim2ToDm6.sites",
    },
    "sec": {
        "sec_F1": f"{G_ALIGN}/sec-F1-1-LJD7630.uniq1.vcf.Filter.ann_AG_dp5mq10.droSec1ToDm6.sites",
        "sec_F2": f"{G_ALIGN}/sec-F2-1-LJD7631.uniq1.vcf.Filter.ann_AG_dp5mq10.droSec1ToDm6.sites",
        "sec_F3": f"{G_ALIGN}/sec-F3-1-LJD7632.uniq1.vcf.Filter.ann_AG_dp5mq10.droSec1ToDm6.sites",
        "sec_M1": f"{G_ALIGN}/sec-M1-1-LJD7627.uniq1.vcf.Filter.ann_AG_dp5mq10.droSec1ToDm6.sites",
        "sec_M2": f"{G_ALIGN}/sec-M2-1-LJD7628.uniq1.vcf.Filter.ann_AG_dp5mq10.droSec1ToDm6.sites",
        "sec_M3": f"{G_ALIGN}/sec-M3-1-LJD7629.uniq1.vcf.Filter.ann_AG_dp5mq10.droSec1ToDm6.sites",
        "sec_M4": f"{G_ALIGN}/sec-M4-old.uniq1.vcf.Filter.ann_AG_dp5mq10.droSec1ToDm6.sites",
    },
    "ere": {
        "ere_F1": f"{G_ALIGN}/ere-F1-1-LJD7624.uniq1.vcf.Filter.ann_AG_dp5mq10.droEre2ToDm6.sites",
        "ere_F2": f"{G_ALIGN}/ere-F2-1-LJD7625.uniq1.vcf.Filter.ann_AG_dp5mq10.droEre2ToDm6.sites",
        "ere_F3": f"{G_ALIGN}/ere-F3-1-LJD7626.uniq1.vcf.Filter.ann_AG_dp5mq10.droEre2ToDm6.sites",
        "ere_M1": f"{G_ALIGN}/ere-M1-1-LJD7621.uniq1.vcf.Filter.ann_AG_dp5mq10.droEre2ToDm6.sites",
        "ere_M2": f"{G_ALIGN}/ere-M2-1-LJD7622.uniq1.vcf.Filter.ann_AG_dp5mq10.droEre2ToDm6.sites",
        "ere_M3": f"{G_ALIGN}/ere-M3-1-LJD7623.uniq1.vcf.Filter.ann_AG_dp5mq10.droEre2ToDm6.sites",
        "ere_M4": f"{G_ALIGN}/ere-M4-old.uniq1.vcf.Filter.ann_AG_dp5mq10.droEre2ToDm6.sites",
    },
}

SPECIES_TEST2 = {
    "sim": {
        "sim_F1": f"{G_ALIGN}/w501-F1-totalRNA-LJD7636-RL1.uniq1.vcf.Filter.ann.test2_A2G.droSim2ToDm6.sites",
        "sim_F2": f"{G_ALIGN}/w501-F2-totalRNA-LJD7637-RL1.uniq1.vcf.Filter.ann.test2_A2G.droSim2ToDm6.sites",
        "sim_F3": f"{G_ALIGN}/w501-F3-totalRNA-LJD7638-RL1.uniq1.vcf.Filter.ann.test2_A2G.droSim2ToDm6.sites",
        "sim_M1": f"{G_ALIGN}/w501-M1-totalRNA-LJD7633-RL1.uniq1.vcf.Filter.ann.test2_A2G.droSim2ToDm6.sites",
        "sim_M2": f"{G_ALIGN}/w501-M2-totalRNA-LJD7634-RL1.uniq1.vcf.Filter.ann.test2_A2G.droSim2ToDm6.sites",
        "sim_M3": f"{G_ALIGN}/w501-M3-totalRNA-LJD7635-RL1.uniq1.vcf.Filter.ann.test2_A2G.droSim2ToDm6.sites",
        "sim_F4": f"{G_ALIGN}/Dsim_female_head_mRNA.uniq1.vcf.Filter.ann.test2_A2G.droSim2ToDm6.sites",
        "sim_M4": f"{G_ALIGN}/Dsim_male_head_mRNA.uniq1.vcf.Filter.ann.test2_A2G.droSim2ToDm6.sites",
    },
    "sec": {
        "sec_F1": f"{G_ALIGN}/sec-F1-1-LJD7630.uniq1.vcf.Filter.ann.test2_A2G.droSec1ToDm6.sites",
        "sec_F2": f"{G_ALIGN}/sec-F2-1-LJD7631.uniq1.vcf.Filter.ann.test2_A2G.droSec1ToDm6.sites",
        "sec_F3": f"{G_ALIGN}/sec-F3-1-LJD7632.uniq1.vcf.Filter.ann.test2_A2G.droSec1ToDm6.sites",
        "sec_M1": f"{G_ALIGN}/sec-M1-1-LJD7627.uniq1.vcf.Filter.ann.test2_A2G.droSec1ToDm6.sites",
        "sec_M2": f"{G_ALIGN}/sec-M2-1-LJD7628.uniq1.vcf.Filter.ann.test2_A2G.droSec1ToDm6.sites",
        "sec_M3": f"{G_ALIGN}/sec-M3-1-LJD7629.uniq1.vcf.Filter.ann.test2_A2G.droSec1ToDm6.sites",
        "sec_M4": f"{G_ALIGN}/Dsec_male_head_mRNA.uniq1.vcf.Filter.ann.test2_A2G.droSec1ToDm6.sites",
    },
    "ere": {
        "ere_F1": f"{G_ALIGN}/ere-F1-1-LJD7624.uniq1.vcf.Filter.ann.test2_A2G.droEre2ToDm6.sites",
        "ere_F2": f"{G_ALIGN}/ere-F2-1-LJD7625.uniq1.vcf.Filter.ann.test2_A2G.droEre2ToDm6.sites",
        "ere_F3": f"{G_ALIGN}/ere-F3-1-LJD7626.uniq1.vcf.Filter.ann.test2_A2G.droEre2ToDm6.sites",
        "ere_M1": f"{G_ALIGN}/ere-M1-1-LJD7621.uniq1.vcf.Filter.ann.test2_A2G.droEre2ToDm6.sites",
        "ere_M2": f"{G_ALIGN}/ere-M2-1-LJD7622.uniq1.vcf.Filter.ann.test2_A2G.droEre2ToDm6.sites",
        "ere_M3": f"{G_ALIGN}/ere-M3-1-LJD7623.uniq1.vcf.Filter.ann.test2_A2G.droEre2ToDm6.sites",
        "ere_M4": f"{G_ALIGN}/Dere_male_head_mRNA.uniq1.vcf.Filter.ann.test2_A2G.droEre2ToDm6.sites",
    },
}

# -----------------------------
# Utilities
# -----------------------------


def ensure_paths(files: Dict[str, str], label: str):
    missing = [f"{k}\t{v}" for k, v in files.items() if not os.path.exists(v)]
    if missing:
        raise FileNotFoundError(f"Missing {label} files:\n" + "\n".join(missing))


def ensure_list(files: List[str], label: str):
    missing = [v for v in files if not os.path.exists(v)]
    if missing:
        raise FileNotFoundError(f"Missing {label} files:\n" + "\n".join(missing))


def tsv_escape(x) -> str:
    if x is None:
        return ""
    return str(x).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def write_tsv(path: Path, header: List[str], rows: Iterable[Dict[str, object]]):
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(tsv_escape(r.get(h, "")) for h in header) + "\n")


def write_list(path: Path, values: Iterable[str], header: Optional[str] = None):
    with open(path, "w", encoding="utf-8") as out:
        if header:
            out.write(header + "\n")
        for v in sorted(set(values)):
            out.write(v + "\n")


def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return 0.0


def extract_change2(change: str) -> str:
    """
    Convert an HGVS-like change to the terminal SNV substitution.

    Examples:
      c.223A>G       -> A>G
      c.*4176T>A     -> T>A
      n.1010845A>G   -> A>G

    This is the Python equivalent of the R-side change2 field, but it is more
    robust than gsub(".*[0-9]", "", change) for unusual HGVS strings.
    """
    matches = re.findall(r"([ACGT])\s*>\s*([ACGT])", str(change).upper())
    if not matches:
        return ""
    ref, alt = matches[-1]
    return f"{ref}>{alt}"


def parse_ann_meta(ann: str) -> Dict[str, object]:
    """Parse transcript metadata from a standard SnpEff ANN field."""
    text = str(ann or "").replace("\\_", "_")
    first = text.split(",", 1)[0]
    fields = first.split("|")
    fields += [""] * max(0, 16 - len(fields))
    feature_type = fields[5]
    transcript_id = fields[6]
    is_transcript = (
        feature_type.lower() == "transcript"
        or transcript_id.startswith("FBtr")
        or "|transcript|" in text
    )
    return {
        "feature_type": feature_type,
        "transcript_id": transcript_id,
        "transcript_biotype": fields[7],
        "hgvs_c": fields[9],
        "is_transcript": is_transcript,
    }


def annotation_priority(t: str, ann: str) -> int:
    s = f"{t};{ann}"
    if re.search(
        r"frameshift_variant|stop_gained|stop_lost|start_lost|"
        r"protein_altering_variant",
        s,
    ):
        return 1
    if "missense_variant" in s:
        return 2
    if "synonymous_variant" in s:
        return 3
    if re.search(
        r"coding_sequence_variant|start_retained_variant|"
        r"stop_retained_variant|incomplete_terminal_codon_variant",
        s,
    ):
        return 4
    if re.search(r"splice_acceptor_variant|splice_donor_variant", s):
        return 5
    if "splice_region_variant" in s:
        return 10
    if re.search(r"3_prime_UTR_variant|5_prime_UTR_variant|prime_UTR_variant", s):
        return 20
    if "non_coding_transcript_exon_variant" in s:
        return 30
    if "intron_variant" in s:
        return 40
    if re.search(r"upstream_gene_variant|downstream_gene_variant", s):
        return 50
    if re.search(r"intragenic_variant|intergenic_region|intergenic_variant", s):
        return 60
    return 90


def make_class2(type_all: str) -> int:
    s = type_all or ""
    if "missense_variant" in s:
        return 1
    if "synonymous_variant" in s:
        return 2
    if re.search(r"3_prime_UTR_variant|5_prime_UTR_variant|prime_UTR_variant", s):
        return 3
    if "non_coding_transcript_exon_variant" in s:
        return 4
    if "intragenic_variant" in s:
        return 5
    return 0


def choose_best(rows: List[Dict[str, object]]) -> Dict[str, object]:
    if not rows:
        return {}
    return sorted(
        rows,
        key=lambda r: (
            annotation_priority(str(r.get("type", "")), str(r.get("ann", ""))),
            # Within the same consequence class, prefer a real transcript feature
            # over gene_variant/null-circ/intragenic annotations.
            0 if bool(r.get("is_transcript", False)) else 1,
            -safe_float(r.get("depth", 0)),
            -safe_float(r.get("altcount", 0)),
            str(r.get("transcript_id", "")),
        ),
    )[0]


def collapse_annotation(rows: Iterable[Dict[str, object]]):
    by_h = defaultdict(list)
    for r in rows:
        h = str(r.get("handle", ""))
        if h:
            by_h[h].append(r)
    out = {}
    for h, rs in by_h.items():
        best = choose_best(rs)
        types = sorted({str(r.get("type", "")) for r in rs if str(r.get("type", ""))})
        genes = sorted({str(r.get("gene", "")) for r in rs if str(r.get("gene", ""))})
        transcripts = sorted(
            {
                str(r.get("transcript_id", ""))
                for r in rs
                if bool(r.get("is_transcript", False))
                and str(r.get("transcript_id", ""))
            }
        )
        sources = sorted(
            {
                str(r.get("annotation_source", ""))
                for r in rs
                if str(r.get("annotation_source", ""))
            }
        )
        type_all = ";".join(types)
        out[h] = {
            "handle": h,
            "chr": best.get("chr", h.rsplit("_", 1)[0]),
            "pos": best.get("pos", h.rsplit("_", 1)[-1]),
            "dmel_chr": best.get("dmel_chr", h.rsplit("_", 1)[0]),
            "dmel_pos": best.get("dmel_pos", h.rsplit("_", 1)[-1]),
            "Ref": best.get("Ref", ""),
            "Alt": best.get("Alt", ""),
            "type": best.get("type", "no_annotation"),
            "type_all": type_all or best.get("type", "no_annotation"),
            # gene/transcript/ANN below belong to the selected representative
            # row. The *_all columns preserve every retained A>G annotation.
            "gene": best.get("gene", ""),
            "gene_all": ";".join(genes),
            "transcript": best.get("transcript_id", ""),
            "transcript_all": ";".join(transcripts),
            "feature_type": best.get("feature_type", ""),
            "transcript_biotype": best.get("transcript_biotype", ""),
            "hgvs_c": best.get("hgvs_c", ""),
            "is_transcript": best.get("is_transcript", False),
            "ann": best.get("ann", ""),
            # change is genomic orientation; change2 is transcript-oriented.
            "change": best.get(
                "change", f"{best.get('Ref', '')}>{best.get('Alt', '')}"
            ),
            "change2": best.get("change2", ""),
            "strand": best.get("strand", ""),
            "source_scaffold": best.get("source_scaffold", ""),
            "source_pos": best.get("source_pos", ""),
            "species_coord": best.get("species_coord", ""),
            "annotation_source": best.get("annotation_source", ""),
            "annotation_sources": ";".join(sources),
            "class2": make_class2(type_all or str(best.get("type", ""))),
        }
    return out


# -----------------------------
# Load ISO1 and Published
# -----------------------------


def load_iso1_at2(files: Dict[str, str]):
    handles = set()
    ann_rows = []
    for lib, f in files.items():
        print(f"[ISO1 at2] {lib} <- {os.path.basename(f)}", flush=True)
        with open(f, "r", encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            if "handle" not in idx:
                raise RuntimeError(f"No handle column in {f}")
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < len(header):
                    cols += [""] * (len(header) - len(cols))
                h = cols[idx["handle"]]
                if not h or h == "NA":
                    continue
                ann_value = cols[idx["ann"]] if "ann" in idx else ""
                ann_meta = parse_ann_meta(ann_value)
                genomic_change = (
                    cols[idx["change"]]
                    if "change" in idx
                    else (
                        f"{cols[idx['Ref']]}>{cols[idx['Alt']]}"
                        if "Ref" in idx and "Alt" in idx
                        else ""
                    )
                )
                change2 = (
                    cols[idx["change2"]]
                    if "change2" in idx
                    else extract_change2(ann_meta["hgvs_c"])
                )
                # Do not rely only on the filename suffix "_AG": enforce the
                # same exact transcript-oriented A>G rule in every input.
                if change2 != "A>G":
                    continue
                handles.add(h)
                strand = (
                    cols[idx["strand"]]
                    if "strand" in idx
                    else ("+" if genomic_change == change2 else "-")
                )
                row = {
                    "handle": h,
                    "chr": cols[idx["chr"]] if "chr" in idx else h.rsplit("_", 1)[0],
                    "pos": cols[idx["pos"]] if "pos" in idx else h.rsplit("_", 1)[-1],
                    "dmel_chr": (
                        cols[idx["chr"]] if "chr" in idx else h.rsplit("_", 1)[0]
                    ),
                    "dmel_pos": (
                        cols[idx["pos"]] if "pos" in idx else h.rsplit("_", 1)[-1]
                    ),
                    "Ref": cols[idx["Ref"]] if "Ref" in idx else "",
                    "Alt": cols[idx["Alt"]] if "Alt" in idx else "",
                    "type": cols[idx["type"]] if "type" in idx else "",
                    "gene": cols[idx["gene"]] if "gene" in idx else "",
                    "ann": ann_value,
                    "change": genomic_change,
                    "change2": change2,
                    "strand": strand,
                    "lib": lib,
                    "refcount": cols[idx["refcount"]] if "refcount" in idx else "",
                    "depth": cols[idx["dp"]] if "dp" in idx else "",
                    "altcount": cols[idx["altcount"]] if "altcount" in idx else "",
                    "ratio": cols[idx["ratio"]] if "ratio" in idx else "",
                    "source_file": f,
                    "source_format": "dmel_at2",
                    "annotation_source": "ISO1",
                }
                row.update(ann_meta)
                ann_rows.append(row)
    return handles, ann_rows


def load_published(files: List[str], duan_file: str):
    handles = set()
    long_rows = []
    # .site files are BED-like from the input annotation pipeline, use start + 1.
    for f in files:
        name = os.path.basename(f)
        with open(f, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                try:
                    pos = int(parts[1]) + 1
                except Exception:
                    continue
                h = f"{parts[0]}_{pos}"
                handles.add(h)
                long_rows.append({"handle": h, "published_set": name})
    # Duan High_Confidence, 1-based Position.
    if os.path.exists(duan_file):
        with open(duan_file, "r", encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            if {"Chromosome", "Position", "High_Confidence"}.issubset(idx):
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < len(header):
                        cols += [""] * (len(header) - len(cols))
                    if cols[idx["High_Confidence"]] in {
                        "TRUE",
                        "True",
                        "true",
                        "T",
                        "1",
                        "YES",
                        "yes",
                        "Y",
                    }:
                        h = f"{cols[idx['Chromosome']]}_{cols[idx['Position']]}"
                        handles.add(h)
                        long_rows.append(
                            {"handle": h, "published_set": "Duan2017_HighConf"}
                        )
    return handles, long_rows


# -----------------------------
# Streaming scanners
# -----------------------------


def scan_dmel_raw_file(args):
    lib, f, target, full = args
    # Keep every annotation whose normalized transcript-oriented substitution
    # is exactly A>G. Do not collapse by handle here: coding/transcript
    # prioritization is performed only after all WBG libraries are combined.
    kept_rows = []
    kept_handles = set()
    with open(f, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11:
                continue
            if len(p[2]) != 1 or len(p[3]) != 1:
                continue
            change2 = extract_change2(p[9])
            # Match the ISO1 R logic exactly: change2 must equal A>G.
            if change2 != "A>G":
                continue
            h = f"{p[0]}_{p[1]}"
            if target is not None and h not in target:
                continue
            refc = safe_float(p[4])
            altc = safe_float(p[5])
            genomic_change = f"{p[2]}>{p[3]}"
            strand = "+" if genomic_change == change2 else "-"
            ann_meta = parse_ann_meta(p[10])
            row = {
                "handle": h,
                "chr": p[0],
                "pos": p[1],
                "Ref": p[2],
                "Alt": p[3],
                "dmel_chr": p[0],
                "dmel_pos": p[1],
                "refcount": p[4],
                "altcount": p[5],
                "depth": refc + altc,
                "MQ": p[6],
                "type": p[7],
                "gene": p[8],
                # Keep genomic and transcript-oriented changes separately.
                "change": genomic_change,
                "change_raw": p[9],
                "change2": change2,
                "strand": strand,
                "ann": p[10],
                "ratio": (altc / (refc + altc)) if (refc + altc) > 0 else "",
                "source_scaffold": "",
                "source_pos": "",
                "species_coord": "",
                "lib": lib,
                "source_file": f,
                "source_format": "dmel_test2",
                "annotation_source": (
                    "WBG"
                    if str(lib).upper().startswith("WBG")
                    else (
                        "AG"
                        if str(lib).upper().startswith("AG_")
                        else "AS" if str(lib).upper().startswith("AS_") else "Adar5G1"
                    )
                ),
            }
            row.update(ann_meta)
            kept_rows.append(row)
            kept_handles.add(h)
    print(
        f"  [done dmel] {lib} handles={len(kept_handles)} "
        f"A>G_annotation_rows={len(kept_rows)}",
        flush=True,
    )
    return lib, kept_rows


def scan_lifted_file(args):
    lib, f, target = args
    kept_rows = []
    kept_handles = set()
    with open(f, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 13:
                continue
            if len(p[5]) != 1 or len(p[6]) != 1:
                continue
            h = f"{p[0]}_{p[2]}"  # first column + third column, BED end
            if target is not None and h not in target:
                continue
            refc = safe_float(p[7])
            altc = safe_float(p[8])
            species_coord = f"{p[3]}_{p[4]}"
            ann_meta = parse_ann_meta(p[12])
            change2 = extract_change2(ann_meta["hgvs_c"])
            # As for D. melanogaster, annotation selection begins only after
            # exact transcript-oriented A>G filtering.
            if change2 != "A>G":
                continue
            genomic_change = f"{p[5]}>{p[6]}"
            strand = "+" if genomic_change == change2 else "-"
            row = {
                "handle": h,
                # chr/pos are the native coordinate of this species.
                "chr": p[3],
                "pos": p[4],
                # dmel_chr/dmel_pos retain the dm6 coordinate used for joins.
                "dmel_chr": p[0],
                "dmel_pos": p[2],
                "Ref": p[5],
                "Alt": p[6],
                "refcount": p[7],
                "altcount": p[8],
                "depth": refc + altc,
                "ratio": (altc / (refc + altc)) if (refc + altc) > 0 else "",
                "MQ": p[9],
                "type": p[10],
                "gene": p[11],
                "change": genomic_change,
                "change2": change2,
                "strand": strand,
                "ann": p[12],
                "source_scaffold": p[3],
                "source_pos": p[4],
                "species_coord": species_coord,
                "lib": lib,
                "source_file": f,
                "source_format": "lifted_sites",
                "annotation_source": str(lib).split("_", 1)[0],
            }
            row.update(ann_meta)
            kept_rows.append(row)
            kept_handles.add(h)
    print(
        f"  [done lifted] {lib} handles={len(kept_handles)} "
        f"A>G_annotation_rows={len(kept_rows)}",
        flush=True,
    )
    return lib, kept_rows


def scan_many(
    files: Dict[str, str], fmt: str, threads: int, target: Optional[Set[str]] = None
):
    if not files:
        return {}
    func = scan_dmel_raw_file if fmt == "dmel" else scan_lifted_file
    tasks = [
        (lib, f, target, target is None) if fmt == "dmel" else (lib, f, target)
        for lib, f in files.items()
    ]
    out = {}
    if threads <= 1:
        for t in tasks:
            lib, rows = func(t)
            out[lib] = rows
    else:
        with ProcessPoolExecutor(max_workers=min(threads, len(tasks))) as ex:
            futs = [ex.submit(func, t) for t in tasks]
            for fut in as_completed(futs):
                lib, rows = fut.result()
                out[lib] = rows
    return out


def rows_to_min_lib(by_lib: Dict[str, List[Dict[str, object]]], min_lib: int):
    support = defaultdict(set)
    for lib, rows in by_lib.items():
        for r in rows:
            support[r["handle"]].add(lib)
    keep = {h for h, libs in support.items() if len(libs) >= min_lib}
    return keep, support


def flatten(by_lib: Dict[str, List[Dict[str, object]]]):
    return [r for rows in by_lib.values() for r in rows]


GROUP_CANDIDATE_HEADER = [
    "group",
    "handle",
    "present",
    "n_libraries",
    "libraries",
    "dmel_chr",
    "dmel_pos",
    "chr",
    "pos",
    "species_coord",
    "Ref",
    "Alt",
    "change",
    "change2",
    "strand",
    "type",
    "type_all",
    "gene",
    "gene_all",
    "transcript",
    "transcript_all",
    "feature_type",
    "transcript_biotype",
    "hgvs_c",
    "is_transcript",
    "class2",
    "ann",
    "annotation_source",
    "annotation_sources",
    "ISO1_type",
    "ISO1_gene",
    "ISO1_transcript",
]

RAW_ANNOTATION_HEADER = [
    "handle",
    "dmel_chr",
    "dmel_pos",
    "chr",
    "pos",
    "source_scaffold",
    "source_pos",
    "species_coord",
    "Ref",
    "Alt",
    "refcount",
    "altcount",
    "depth",
    "ratio",
    "MQ",
    "type",
    "gene",
    "change",
    "change_raw",
    "change2",
    "strand",
    "feature_type",
    "transcript_id",
    "transcript_biotype",
    "hgvs_c",
    "is_transcript",
    "ann",
    "lib",
    "annotation_source",
    "source_file",
]


def rows_by_library(rows: Iterable[Dict[str, object]]):
    out = defaultdict(list)
    for r in rows:
        lib = str(r.get("lib", ""))
        if lib:
            out[lib].append(r)
    return dict(out)


def build_group_candidate_rows(
    group: str,
    candidate_set: Set[str],
    detected: Set[str],
    support: Dict[str, Set[str]],
    annotation_rows: Iterable[Dict[str, object]],
    master_by_h: Dict[str, Dict[str, object]],
):
    """Build one representative, group-specific annotation per detected site."""
    retained = [
        r
        for r in annotation_rows
        if str(r.get("handle", "")) in candidate_set
        and str(r.get("handle", "")) in detected
        and str(r.get("change2", "")) == "A>G"
    ]
    ann_by_h = collapse_annotation(retained)
    out = []
    for h in sorted(detected & candidate_set):
        a = ann_by_h.get(h, {})
        master = master_by_h.get(h, {})
        out.append(
            {
                "group": group,
                "handle": h,
                "present": 1,
                "n_libraries": len(support.get(h, set())),
                "libraries": ";".join(sorted(support.get(h, set()))),
                "dmel_chr": a.get("dmel_chr", h.rsplit("_", 1)[0]),
                "dmel_pos": a.get("dmel_pos", h.rsplit("_", 1)[-1]),
                "chr": a.get("chr", ""),
                "pos": a.get("pos", ""),
                "species_coord": a.get("species_coord", ""),
                "Ref": a.get("Ref", ""),
                "Alt": a.get("Alt", ""),
                "change": a.get("change", ""),
                "change2": a.get("change2", ""),
                "strand": a.get("strand", ""),
                "type": a.get("type", "no_annotation"),
                "type_all": a.get("type_all", ""),
                "gene": a.get("gene", ""),
                "gene_all": a.get("gene_all", ""),
                "transcript": a.get("transcript", ""),
                "transcript_all": a.get("transcript_all", ""),
                "feature_type": a.get("feature_type", ""),
                "transcript_biotype": a.get("transcript_biotype", ""),
                "hgvs_c": a.get("hgvs_c", ""),
                "is_transcript": a.get("is_transcript", ""),
                "class2": a.get("class2", 0),
                "ann": a.get("ann", ""),
                "annotation_source": a.get("annotation_source", group),
                "annotation_sources": a.get("annotation_sources", group),
                # ISO1 columns are included only as stable master-list context.
                # They do not participate in this group's representative choice.
                "ISO1_type": master.get("type", ""),
                "ISO1_gene": master.get("gene", ""),
                "ISO1_transcript": master.get("transcript", ""),
            }
        )
    return out, retained


def build_heatmap_measurements(
    group: str,
    by_lib: Dict[str, List[Dict[str, object]]],
    candidate_set: Set[str],
    detected: Set[str],
):
    """
    Return one editing-level measurement per library and candidate handle.

    Only observed rows are written. A missing row is left as NA in the wide
    matrix; it is not silently converted to zero because these site-call files
    do not provide coverage for uncalled positions.
    """
    out = []
    for lib in sorted(by_lib):
        by_h = defaultdict(list)
        for r in by_lib[lib]:
            h = str(r.get("handle", ""))
            if h in candidate_set and str(r.get("change2", "")) == "A>G":
                by_h[h].append(r)
        for h, rs in by_h.items():
            best = choose_best(rs)
            depth = safe_float(best.get("depth", 0))
            altcount = safe_float(best.get("altcount", 0))
            ratio = best.get("ratio", "")
            ratio_value = (
                safe_float(ratio)
                if str(ratio) != ""
                else (altcount / depth if depth > 0 else "")
            )
            sex_match = re.search(r"(?:^|_)([FM])(?:\d|_|$)", lib, re.I)
            rep_match = re.search(r"([0-9]+)$", lib)
            out.append(
                {
                    "group": group,
                    "library": lib,
                    "sex": sex_match.group(1).upper() if sex_match else "",
                    "replicate": rep_match.group(1) if rep_match else "",
                    "handle": h,
                    "editing_level": ratio_value,
                    "refcount": best.get("refcount", ""),
                    "altcount": best.get("altcount", ""),
                    "depth": best.get("depth", ""),
                    "group_min_support_present": 1 if h in detected else 0,
                }
            )
    return out


def write_heatmap_matrix(
    path: Path,
    measurements: List[Dict[str, object]],
    handles: Iterable[str],
):
    selected = set(handles)
    samples = sorted(
        {(str(r["group"]), str(r["library"])) for r in measurements},
        key=lambda x: (
            (
                ["ISO1", "WBG", "AG", "AS", "sim", "sec", "ere"].index(x[0])
                if x[0] in ["ISO1", "WBG", "AG", "AS", "sim", "sec", "ere"]
                else 99
            ),
            x[1],
        ),
    )
    sample_names = [f"{g}:{lib}" for g, lib in samples]
    value = {
        (str(r["handle"]), f"{r['group']}:{r['library']}"): r["editing_level"]
        for r in measurements
        if str(r["handle"]) in selected
    }
    rows = []
    for h in sorted(selected):
        row = {"handle": h}
        for sample in sample_names:
            row[sample] = value.get((h, sample), "")
        rows.append(row)
    write_tsv(path, ["handle"] + sample_names, rows)


def one2one_filter(
    by_lib: Dict[str, List[Dict[str, object]]], species: str, outdir: Path
):
    pairs = set()
    for rows in by_lib.values():
        for r in rows:
            if r.get("handle") and r.get("species_coord"):
                pairs.add((str(r["handle"]), str(r["species_coord"])))
    by_h = defaultdict(set)
    by_s = defaultdict(set)
    for h, s in pairs:
        by_h[h].add(s)
        by_s[s].add(h)
    bad_h = {h for h, ss in by_h.items() if len(ss) > 1}
    bad_s = {s for s, hs in by_s.items() if len(hs) > 1}
    bad_from_s = {h for s in bad_s for h in by_s[s]}
    bad_handles = bad_h | bad_from_s
    kept_pairs = [(h, s) for h, s in pairs if h not in bad_handles and s not in bad_s]
    dropped_pairs = [(h, s) for h, s in pairs if h in bad_handles or s in bad_s]
    write_tsv(
        outdir / f"{species}.bad_dmel_handle_to_multiple_species_coord.tsv",
        ["handle", "species_coords"],
        [
            {"handle": h, "species_coords": ",".join(sorted(by_h[h]))}
            for h in sorted(bad_h)
        ],
    )
    write_tsv(
        outdir / f"{species}.bad_species_coord_to_multiple_dmel_handle.tsv",
        ["species_coord", "handles"],
        [
            {"species_coord": s, "handles": ",".join(sorted(by_s[s]))}
            for s in sorted(bad_s)
        ],
    )
    write_tsv(
        outdir / f"{species}.dropped_non1v1_pairs.tsv",
        ["handle", "species_coord"],
        [{"handle": h, "species_coord": s} for h, s in sorted(dropped_pairs)],
    )
    write_tsv(
        outdir / f"{species}.kept_1v1_pairs.tsv",
        ["handle", "species_coord"],
        [{"handle": h, "species_coord": s} for h, s in sorted(kept_pairs)],
    )
    kept = {h for h, _ in kept_pairs}
    clean = {}
    for lib, rows in by_lib.items():
        clean[lib] = [r for r in rows if r["handle"] in kept]
    print(
        f"[1v1] {species}: before={len(by_h)} kept={len(kept)} dropped={len(bad_handles)}",
        flush=True,
    )
    return clean, kept


# -----------------------------
# Plotting via R
# -----------------------------


def draw_plots(outdir: Path, n_upset: int):
    r_script = outdir / "plot_upset_species_overlap.R"
    r_script.write_text(
        r"""
suppressPackageStartupMessages({
  library(data.table)
})
args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
n_upset <- as.integer(args[2])
plot_dir <- file.path(outdir, "plots")
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)

between <- function(row, min, max) {
  (row["class2"] <= max) & (row["class2"] >= min)
}

plot_upset <- function(file, sets, prefix, width=12, height=8) {
  if (!requireNamespace("UpSetR", quietly=TRUE)) return(invisible(NULL))
  library(UpSetR)
  dt <- fread(file)
  dt$Name <- dt$handle
  pdf(file.path(plot_dir, paste0(prefix, ".pdf")), width=width, height=height)
  grid::grid.newpage()
  ok <- TRUE
  tryCatch({
    print(upset(
      as.data.frame(dt), sets=sets, nsets=length(sets), nintersects=n_upset,
      order.by="freq", keep.order=TRUE, mb.ratio=c(0.55,0.45),
      text.scale=1.8, number.angles=0, point.size=4.5, line.size=1.5,
      mainbar.y.label="Edit Sites Intersections", sets.x.label="Sites Per Set",
      queries=list(
        list(query=between, params=list(1,2), color="#89CAF0", active=TRUE),
        list(query=between, params=list(1,1), color="#f38181", active=TRUE)
      )
    ))
  }, error=function(e){
    ok <<- FALSE
    message("[UpSet warning] fallback without queries: ", conditionMessage(e))
  })
  if (!ok) {
    grid::grid.newpage()
    print(upset(
      as.data.frame(dt), sets=sets, nsets=length(sets), nintersects=n_upset,
      order.by="freq", keep.order=TRUE, mb.ratio=c(0.55,0.45),
      text.scale=1.8, number.angles=0, point.size=4.5, line.size=1.5,
      mainbar.y.label="Edit Sites Intersections", sets.x.label="Sites Per Set"
    ))
  }
  dev.off()
  png(file.path(plot_dir, paste0(prefix, ".png")), width=3200, height=2200, res=250)
  grid::grid.newpage()
  print(upset(
    as.data.frame(dt), sets=sets, nsets=length(sets), nintersects=n_upset,
    order.by="freq", keep.order=TRUE, mb.ratio=c(0.55,0.45),
    text.scale=1.8, number.angles=0, point.size=4.5, line.size=1.5,
    mainbar.y.label="Edit Sites Intersections", sets.x.label="Sites Per Set"
  ))
  dev.off()
}

plot_upset(file.path(outdir, "dmel_3set_upset.input.tsv"), c("Published_union","ISO1_union","WBG_min2"), "dmel_3set_upset", width=11, height=7.5)
plot_upset(file.path(outdir, "candidate.species_overlap.input.tsv"), c("Dere","Dsec","Dsim","Dmel"), "candidate_species_overlap", width=12, height=8)
""",
        encoding="utf-8",
    )
    rscript = os.environ.get("PLOT_RSCRIPT", "Rscript")
    if not shutil.which(rscript):
        raise FileNotFoundError(f"Plotting Rscript not found: {rscript}")
    subprocess.run([rscript, str(r_script), str(outdir), str(n_upset)], check=True)
    for prefix in ("dmel_3set_upset", "candidate_species_overlap"):
        for ext in ("pdf", "png"):
            output = outdir / "plots" / f"{prefix}.{ext}"
            if not output.is_file() or output.stat().st_size == 0:
                raise RuntimeError(f"Missing plot output: {output}")


# -----------------------------
# Main
# -----------------------------


def main():
    outdir = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else Path(
            "/path/to/RNA_editing_DATAprocessing/all_code/final_data/rebuild_logic_WBG_ISO1_Published_then_species_FAST"
        )
    )
    threads = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    n_upset = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    iso1_mode = sys.argv[4].lower() if len(sys.argv) > 4 else "all"
    species_mode = sys.argv[5].lower() if len(sys.argv) > 5 else "test2"
    if iso1_mode not in {"all", "head"}:
        raise ValueError("ISO1_MODE must be all or head")
    if species_mode not in {"dp5mq10", "test2"}:
        raise ValueError("SPECIES_MODE must be dp5mq10 or test2")

    dirs = {
        "site": outdir / "site_lists",
        "raw": outdir / "raw_rows_by_group",
        "group": outdir / "candidate_lists_by_group",
        "heatmap": outdir / "heatmap_inputs",
        "o2o": outdir / "one2one_reports",
        "plot": outdir / "plots",
    }
    for d in [outdir] + list(dirs.values()):
        d.mkdir(parents=True, exist_ok=True)

    print(f"[config] outdir={outdir}")
    print(f"[config] threads={threads}")
    print(f"[config] iso1_mode={iso1_mode}")
    print(f"[config] species_mode={species_mode}")

    iso1_files = dict(ISO1_HEAD_AT2_FILES)
    if iso1_mode == "all":
        iso1_files = {**ISO1_BRAIN_AT2_FILES, **ISO1_HEAD_AT2_FILES}

    species_files = SPECIES_DP5MQ10 if species_mode == "dp5mq10" else SPECIES_TEST2

    ensure_paths(iso1_files, "ISO1 at2")
    ensure_paths(ADAR5G1_RAW_FILES, "Adar5G1 raw")
    ensure_paths(WBG_RAW_FILES, "WBG raw")
    ensure_paths(AG_RAW_FILES, "AG raw")
    ensure_paths(AS_RAW_FILES, "AS raw")
    ensure_paths(species_files["sim"], "sim species")
    ensure_paths(species_files["sec"], "sec species")
    ensure_paths(species_files["ere"], "ere species")
    ensure_list(PUBLISHED_FILES, "published")
    if not os.path.exists(DUAN2017_FILE):
        raise FileNotFoundError(DUAN2017_FILE)

    # Step 1: Adar5G1 raw A>G, full raw files, union.
    print("[1] Adar5G1 raw A>G full scan", flush=True)
    adar_by_lib = scan_many(ADAR5G1_RAW_FILES, "dmel", threads, target=None)
    adar_handles = {r["handle"] for r in sum((v for v in adar_by_lib.values()), [])}
    write_list(dirs["site"] / "Adar5G1_rawAtoG.handles.txt", adar_handles)

    # Step 2: WBG raw A>G, full raw files, min2, minus Adar.
    print("[2] WBG raw A>G full scan", flush=True)
    wbg_by_lib = scan_many(WBG_RAW_FILES, "dmel", threads, target=None)
    wbg_support = defaultdict(set)
    for lib, rows in wbg_by_lib.items():
        for r in rows:
            wbg_support[r["handle"]].add(lib)
    wbg_min2_before = {
        h for h, libs in wbg_support.items() if len(libs) >= MIN_LIBRARY_SUPPORT["WBG"]
    }
    wbg_min2 = wbg_min2_before - adar_handles
    write_list(
        dirs["site"] / "WBG_min2.before_minus_Adar5G1.handles.txt", wbg_min2_before
    )
    write_list(dirs["site"] / "WBG_min2.handles.txt", wbg_min2)

    # Step 3: ISO1 at2 union, minus Adar.
    print("[3] ISO1 at2 union", flush=True)
    iso1_before, iso1_ann_rows = load_iso1_at2(iso1_files)
    iso1_union = iso1_before - adar_handles
    write_list(
        dirs["site"] / "ISO1_union.before_minus_Adar5G1.handles.txt", iso1_before
    )
    write_list(dirs["site"] / "ISO1_union.handles.txt", iso1_union)

    # Step 4: Published.
    print("[4] Published union", flush=True)
    published, published_long = load_published(PUBLISHED_FILES, DUAN2017_FILE)
    write_list(dirs["site"] / "Published_union.handles.txt", published)
    write_tsv(
        outdir / "Published_union.long.tsv", ["handle", "published_set"], published_long
    )

    # Step 5: Dmel 3-set candidate definition and ISO1-authoritative annotation.
    # Biological membership:
    #   candidate = ISO1_union AND (WBG_min2 OR Published_union)
    # Annotation authority:
    #   candidate representative = best retained ISO1 A>G annotation only.
    print("[5] Candidate and ISO1-authoritative master annotation", flush=True)
    candidate_handles = sorted(iso1_union & (wbg_min2 | published))
    candidate_set = set(candidate_handles)
    write_list(
        dirs["site"] / "candidate.ISO1_AND_WBG_or_Published.handles.txt",
        candidate_handles,
    )
    write_list(dirs["site"] / "candidate.final.handles.txt", candidate_handles)

    iso1_candidate_annotation_rows = [
        r
        for r in iso1_ann_rows
        if r["handle"] in candidate_set and r.get("change2") == "A>G"
    ]
    iso1_master_ann = collapse_annotation(iso1_candidate_annotation_rows)
    iso1_support = defaultdict(set)
    for r in iso1_candidate_annotation_rows:
        iso1_support[r["handle"]].add(r["lib"])

    candidate_rows = []
    for h in candidate_handles:
        a = iso1_master_ann.get(h, {})
        t = a.get("type", "no_annotation")
        ta = a.get("type_all", t)
        candidate_rows.append(
            {
                "handle": h,
                "WBG_min2": 1 if h in wbg_min2 else 0,
                "ISO1_union": 1,
                "Published_union": 1 if h in published else 0,
                "n_ISO1_libraries": len(iso1_support.get(h, set())),
                "ISO1_libraries": ";".join(sorted(iso1_support.get(h, set()))),
                "chr": a.get("chr", h.rsplit("_", 1)[0]),
                "pos": a.get("pos", h.rsplit("_", 1)[-1]),
                "Ref": a.get("Ref", ""),
                "Alt": a.get("Alt", ""),
                "change": a.get("change", ""),
                "change2": a.get("change2", ""),
                "strand": a.get("strand", ""),
                "type": t,
                "type_all": ta,
                "gene": a.get("gene", ""),
                "gene_all": a.get("gene_all", ""),
                "transcript": a.get("transcript", ""),
                "transcript_all": a.get("transcript_all", ""),
                "feature_type": a.get("feature_type", ""),
                "transcript_biotype": a.get("transcript_biotype", ""),
                "hgvs_c": a.get("hgvs_c", ""),
                "is_transcript": a.get("is_transcript", ""),
                "ann": a.get("ann", ""),
                "annotation_source": "ISO1",
                "annotation_sources": a.get("annotation_sources", "ISO1"),
                "class2": make_class2(ta),
                "candidate_rule": "ISO1_AND_WBG_or_Published",
                "annotation_rule": "ISO1_only_AtoG_then_coding_transcript_priority",
            }
        )

    master_header = [
        "handle",
        "WBG_min2",
        "ISO1_union",
        "Published_union",
        "n_ISO1_libraries",
        "ISO1_libraries",
        "chr",
        "pos",
        "Ref",
        "Alt",
        "change",
        "change2",
        "strand",
        "type",
        "type_all",
        "gene",
        "gene_all",
        "transcript",
        "transcript_all",
        "feature_type",
        "transcript_biotype",
        "hgvs_c",
        "is_transcript",
        "ann",
        "annotation_source",
        "annotation_sources",
        "class2",
        "candidate_rule",
        "annotation_rule",
    ]
    candidate_row_by_h = {r["handle"]: r for r in candidate_rows}
    write_tsv(
        outdir / "candidate.master.ISO1_annotation.tsv", master_header, candidate_rows
    )
    write_tsv(outdir / "candidate.final.tsv", master_header, candidate_rows)
    write_tsv(
        outdir / "candidate.handle_type.tsv",
        ["handle", "type"],
        [{"handle": r["handle"], "type": r["type"]} for r in candidate_rows],
    )

    # Keep the original three D. melanogaster membership columns for the
    # three-set UpSet, but all annotation columns are now ISO1-only.
    write_tsv(
        outdir / "dmel_3set_candidate_upset.input.tsv", master_header, candidate_rows
    )
    write_tsv(outdir / "dmel_3set_upset.input.tsv", master_header, candidate_rows)
    print(f"    candidate n = {len(candidate_handles)}", flush=True)

    # Step 6: Build WBG, AG and AS candidate lists from raw multi-annotation
    # files. A>G filtering happens inside scan_dmel_raw_file before collapse.
    print("[6] WBG/AG/AS candidate lists", flush=True)
    ag_by_lib = scan_many(AG_RAW_FILES, "dmel", threads, target=candidate_set)
    as_by_lib = scan_many(AS_RAW_FILES, "dmel", threads, target=candidate_set)
    ag_min2, ag_support = rows_to_min_lib(ag_by_lib, MIN_LIBRARY_SUPPORT["AG"])
    as_min2, as_support = rows_to_min_lib(as_by_lib, MIN_LIBRARY_SUPPORT["AS"])
    wbg_candidate_set = candidate_set & wbg_min2
    ag_candidate_set = candidate_set & ag_min2
    as_candidate_set = candidate_set & as_min2
    write_list(dirs["site"] / "WBG.candidate.min2.handles.txt", wbg_candidate_set)
    write_list(dirs["site"] / "AG.candidate.min2.handles.txt", ag_candidate_set)
    write_list(dirs["site"] / "AS.candidate.min2.handles.txt", as_candidate_set)

    # Step 7: Orthologous species lists. Exact A>G filtering is performed
    # before annotation selection; 1:1 mapping cleanup precedes the union across libraries.
    print("[7] sim/sec/ere candidate lists", flush=True)
    species_clean_union = {}
    species_original_union = {}
    species_support = {}
    species_clean_by_lib = {}
    for sp in ["sim", "sec", "ere"]:
        print(f"  [species] {sp}", flush=True)
        by_lib0 = scan_many(species_files[sp], "lifted", threads, target=candidate_set)
        clean_by_lib, kept_1v1 = one2one_filter(by_lib0, sp, dirs["o2o"])
        min_original, support0 = rows_to_min_lib(by_lib0, MIN_LIBRARY_SUPPORT[sp])
        min_clean, support = rows_to_min_lib(clean_by_lib, MIN_LIBRARY_SUPPORT[sp])
        species_original_union[sp] = min_original & candidate_set
        species_clean_union[sp] = min_clean & candidate_set
        species_support[sp] = support
        species_clean_by_lib[sp] = clean_by_lib
        write_list(
            dirs["site"] / f"{sp}_overlap_candidate.before_1v1.union.handles.txt",
            species_original_union[sp],
        )
        write_list(
            dirs["site"] / f"{sp}_overlap_candidate.1v1.union.handles.txt",
            species_clean_union[sp],
        )

    # Step 8: Write one candidate list per group. Each representative
    # annotation is selected only from that group's exact A>G annotations.
    print("[8] Annotated candidate lists by group", flush=True)
    group_by_lib = {
        "ISO1": rows_by_library(iso1_candidate_annotation_rows),
        "WBG": {
            lib: [r for r in rows if r["handle"] in candidate_set]
            for lib, rows in wbg_by_lib.items()
        },
        "AG": ag_by_lib,
        "AS": as_by_lib,
        "sim": species_clean_by_lib["sim"],
        "sec": species_clean_by_lib["sec"],
        "ere": species_clean_by_lib["ere"],
    }
    group_detected = {
        "ISO1": candidate_set,
        "WBG": wbg_candidate_set,
        "AG": ag_candidate_set,
        "AS": as_candidate_set,
        "sim": species_clean_union["sim"],
        "sec": species_clean_union["sec"],
        "ere": species_clean_union["ere"],
    }
    group_support = {
        "ISO1": iso1_support,
        "WBG": wbg_support,
        "AG": ag_support,
        "AS": as_support,
        "sim": species_support["sim"],
        "sec": species_support["sec"],
        "ere": species_support["ere"],
    }

    all_heatmap_measurements = []
    for group in ["ISO1", "WBG", "AG", "AS", "sim", "sec", "ere"]:
        group_rows, retained_rows = build_group_candidate_rows(
            group=group,
            candidate_set=candidate_set,
            detected=group_detected[group],
            support=group_support[group],
            annotation_rows=flatten(group_by_lib[group]),
            master_by_h=candidate_row_by_h,
        )
        write_tsv(
            dirs["group"] / f"{group}.candidate.detected.annotated.tsv",
            GROUP_CANDIDATE_HEADER,
            group_rows,
        )
        write_tsv(
            dirs["raw"] / f"{group}.candidate.AtoG.annotations.long.tsv",
            RAW_ANNOTATION_HEADER,
            retained_rows,
        )
        all_heatmap_measurements.extend(
            build_heatmap_measurements(
                group, group_by_lib[group], candidate_set, group_detected[group]
            )
        )

    # Step 9: UpSet-ready binary membership tables.
    species_rows = []
    all_group_rows = []
    for h in candidate_handles:
        base = dict(candidate_row_by_h[h])
        base.update(
            {
                "Dmel": 1,
                "Dsim": 1 if h in species_clean_union["sim"] else 0,
                "Dsec": 1 if h in species_clean_union["sec"] else 0,
                "Dere": 1 if h in species_clean_union["ere"] else 0,
            }
        )
        species_rows.append(base)
        all_group_rows.append(
            {
                "handle": h,
                "ISO1": 1,
                "WBG": 1 if h in wbg_candidate_set else 0,
                "AG": 1 if h in ag_candidate_set else 0,
                "AS": 1 if h in as_candidate_set else 0,
                "Dsim": 1 if h in species_clean_union["sim"] else 0,
                "Dsec": 1 if h in species_clean_union["sec"] else 0,
                "Dere": 1 if h in species_clean_union["ere"] else 0,
                "ISO1_type": candidate_row_by_h[h]["type"],
                "ISO1_gene": candidate_row_by_h[h]["gene"],
                "ISO1_transcript": candidate_row_by_h[h]["transcript"],
            }
        )
    species_header = master_header + ["Dmel", "Dsim", "Dsec", "Dere"]
    write_tsv(
        outdir / "candidate.species_overlap.input.tsv", species_header, species_rows
    )
    write_tsv(
        outdir / "dmel_upset_universe.species_overlap.input.tsv",
        species_header,
        species_rows,
    )
    all_group_header = [
        "handle",
        "ISO1",
        "WBG",
        "AG",
        "AS",
        "Dsim",
        "Dsec",
        "Dere",
        "ISO1_type",
        "ISO1_gene",
        "ISO1_transcript",
    ]
    write_tsv(
        outdir / "candidate.upset.all_groups.input.tsv",
        all_group_header,
        all_group_rows,
    )

    exact_counter = Counter(
        (r["Dmel"], r["Dsim"], r["Dsec"], r["Dere"]) for r in species_rows
    )
    exact_rows = []
    for combo, n in sorted(exact_counter.items(), key=lambda x: -x[1]):
        names = [
            name for name, val in zip(["Dmel", "Dsim", "Dsec", "Dere"], combo) if val
        ]
        exact_rows.append(
            {
                "intersection": "&".join(names),
                "Dmel": combo[0],
                "Dsim": combo[1],
                "Dsec": combo[2],
                "Dere": combo[3],
                "N": n,
            }
        )
    write_tsv(
        outdir / "candidate.species_overlap.exact_intersections.tsv",
        ["intersection", "Dmel", "Dsim", "Dsec", "Dere", "N"],
        exact_rows,
    )
    write_tsv(
        outdir / "dmel_upset_universe.species_overlap.exact_intersections.tsv",
        ["intersection", "Dmel", "Dsim", "Dsec", "Dere", "N"],
        exact_rows,
    )

    # Step 10: Heatmap-ready observed editing levels. Missing measurements
    # remain blank/NA rather than being incorrectly converted to zero.
    heatmap_long_header = [
        "group",
        "library",
        "sex",
        "replicate",
        "handle",
        "editing_level",
        "refcount",
        "altcount",
        "depth",
        "group_min_support_present",
    ]
    write_tsv(
        dirs["heatmap"] / "candidate.editing_level.observed.long.tsv",
        heatmap_long_header,
        all_heatmap_measurements,
    )
    write_heatmap_matrix(
        dirs["heatmap"] / "candidate.editing_level.observed.matrix.tsv",
        all_heatmap_measurements,
        candidate_set,
    )
    shared_four_species = (
        candidate_set
        & species_clean_union["sim"]
        & species_clean_union["sec"]
        & species_clean_union["ere"]
    )
    write_list(
        dirs["site"] / "candidate.shared_Dmel_Dsim_Dsec_Dere.handles.txt",
        shared_four_species,
    )
    shared_measurements = [
        r for r in all_heatmap_measurements if r["handle"] in shared_four_species
    ]
    write_tsv(
        dirs["heatmap"] / "shared_four_species.editing_level.observed.long.tsv",
        heatmap_long_header,
        shared_measurements,
    )
    write_heatmap_matrix(
        dirs["heatmap"] / "shared_four_species.editing_level.observed.matrix.tsv",
        shared_measurements,
        shared_four_species,
    )

    counts = [
        {"item": "Adar5G1_rawAtoG", "n": len(adar_handles)},
        {"item": "WBG_min2_before_minus_Adar5G1", "n": len(wbg_min2_before)},
        {"item": "WBG_min2", "n": len(wbg_min2)},
        {"item": "ISO1_before_minus_Adar5G1", "n": len(iso1_before)},
        {"item": "ISO1_union", "n": len(iso1_union)},
        {"item": "Published_union", "n": len(published)},
        {"item": "candidate_ISO1_AND_WBG_or_Published", "n": len(candidate_handles)},
        {"item": "candidate_WBG_min2", "n": len(wbg_candidate_set)},
        {"item": "candidate_AG_min2", "n": len(ag_candidate_set)},
        {"item": "candidate_AS_min2", "n": len(as_candidate_set)},
        {
            "item": "candidate_overlap_Dsim_union_1v1",
            "n": len(species_clean_union["sim"]),
        },
        {
            "item": "candidate_overlap_Dsec_union_1v1",
            "n": len(species_clean_union["sec"]),
        },
        {
            "item": "candidate_overlap_Dere_union_1v1",
            "n": len(species_clean_union["ere"]),
        },
        {"item": "candidate_shared_Dmel_Dsim_Dsec_Dere", "n": len(shared_four_species)},
    ]
    write_tsv(outdir / "final_count.tsv", ["item", "n"], counts)

    draw_plots(outdir, int(n_upset))

    print("\n[DONE]")
    print(f"Output directory: {outdir}")
    print("Main outputs:")
    for p in [
        outdir / "final_count.tsv",
        outdir / "candidate.master.ISO1_annotation.tsv",
        outdir / "dmel_3set_upset.input.tsv",
        outdir / "candidate.upset.all_groups.input.tsv",
        outdir / "candidate.species_overlap.input.tsv",
        dirs["site"] / "Adar5G1_rawAtoG.handles.txt",
        dirs["site"] / "WBG_min2.handles.txt",
        dirs["site"] / "ISO1_union.handles.txt",
        dirs["site"] / "Published_union.handles.txt",
        dirs["site"] / "candidate.final.handles.txt",
        dirs["site"] / "candidate.shared_Dmel_Dsim_Dsec_Dere.handles.txt",
        dirs["group"] / "WBG.candidate.detected.annotated.tsv",
        dirs["group"] / "AG.candidate.detected.annotated.tsv",
        dirs["group"] / "AS.candidate.detected.annotated.tsv",
        dirs["group"] / "sim.candidate.detected.annotated.tsv",
        dirs["group"] / "sec.candidate.detected.annotated.tsv",
        dirs["group"] / "ere.candidate.detected.annotated.tsv",
        dirs["heatmap"] / "candidate.editing_level.observed.long.tsv",
        dirs["heatmap"] / "shared_four_species.editing_level.observed.matrix.tsv",
        dirs["plot"] / "dmel_3set_upset.pdf",
        dirs["plot"] / "candidate_species_overlap.pdf",
    ]:
        print("  ", p)


if __name__ == "__main__":
    main()
