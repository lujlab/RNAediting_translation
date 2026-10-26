#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Strictly rebuild cdspos outputs from original candidate row tables.

Correct order:
  1. Start from original sample_rows_by_group/*.candidate.rows.tsv inputs.
  2. Keep only rows whose original SnpEff ANN record is a legacy coding CDS effect.
  3. Parse the source transcript/cDNA/CDS/AA coordinates from that original ANN record.
  4. Project source transcript position to the selected longest transcript.
  5. Keep only target positions mapped to CDS and deduplicate by
     new_transcript + new_cds_pos_1based.

No rows are taken from previous *.cdspos*.new outputs.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

APPENDED = [
    "new_transcript",
    "new_tx_pos_1based",
    "target_region",
    "new_cds_pos_1based",
]
CDSPOS_HEADER = [
    "genomic_key",
    "base",
    "variant_effect",
    "source_transcript",
    "source_tx_pos_1based",
    "source_tx_length",
    "source_cds_pos_1based",
    "source_cds_length",
    "source_aa_pos_1based",
    "source_aa_length",
    "site_key",
    "new_transcript",
    "new_tx_pos_1based",
    "target_region",
    "new_cds_pos_1based",
]

# This is the same coding-effect gate used by the legacy change_tsv_to_CDSpos.sh.
CODING_EFFECTS = {
    "synonymous_variant",
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "stop_retained_variant",
    "missense_variant&splice_region_variant",
    "splice_region_variant&synonymous_variant",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--group", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--longest-list", required=True)
    p.add_argument("--output-cdspos", required=True)
    p.add_argument("--output-full", required=True)
    p.add_argument("--summary", default=None)
    p.add_argument("--tmp-dir", default="/tmp")
    p.add_argument("--missing-value", default="NA")
    p.add_argument("--keep-temp", action="store_true")
    p.add_argument(
        "--fallback", choices=["exclude", "first", "source"], default="exclude"
    )
    p.add_argument(
        "--map-script",
        default="/path/to/RNA_editing_DATAprocessing/all_code/map_transcript_coordinates_v6.py",
    )
    p.add_argument(
        "--longest-script",
        default="/path/to/RNA_editing_DATAprocessing/all_code/add_longest_and_deduplicate.py",
    )
    return p.parse_args()


def clean(value: Optional[str]) -> str:
    return "" if value is None else str(value).strip()


def parse_pos_len(value: str) -> Tuple[str, str]:
    value = clean(value)
    if not value or "/" not in value:
        return "", ""
    pos, length = value.split("/", 1)
    if pos.isdigit() and length.isdigit():
        return pos, length
    return "", ""


def parse_hgvs_c_pos(value: str) -> str:
    value = clean(value)
    if not value.startswith("c."):
        return ""
    rest = value[2:]
    match = re.match(r"^(\d+)([+-](\d+))?", rest)
    if not match:
        return ""
    base = int(match.group(1))
    if match.group(3) is None:
        return str(base)
    delta = int(match.group(3))
    sign = rest[len(match.group(1))]
    return str(base + delta if sign == "+" else max(base - delta, 1))


def parse_ann_record(record: str, row: Dict[str, str]) -> Dict[str, str]:
    fields = record.split("|") if record else []
    fields += [""] * max(0, 16 - len(fields))

    tx_pos, tx_len = parse_pos_len(fields[11])
    cds_pos, cds_len = parse_pos_len(fields[12])
    aa_pos, aa_len = parse_pos_len(fields[13])

    tx = clean(fields[6])
    if not tx or tx in {".", "null", "null_circ"}:
        tx = clean(
            row.get("transcript_id")
            or row.get("transcript")
            or row.get("source_transcript")
        )

    hgvs_c = clean(fields[9] or row.get("hgvs_c"))
    if not tx_pos:
        tx_pos = clean(
            row.get("source_tx_pos_1based")
        )  # c. is CDS-relative, never a cDNA-position fallback

    if not cds_pos:
        cds_pos = clean(row.get("source_cds_pos_1based"))
    if not cds_len:
        cds_len = clean(row.get("source_cds_length"))

    info = {
        "source_transcript": tx,
        "source_tx_pos_1based": tx_pos,
        "source_tx_length": tx_len,
        "source_cds_pos_1based": cds_pos,
        "source_cds_length": cds_len,
        "source_aa_pos_1based": aa_pos,
        "source_aa_length": aa_len,
        "variant_effect": clean(
            fields[1] or row.get("variant_effect") or row.get("type")
        ),
        "hgvs_c": hgvs_c,
    }
    info["has_source_cds"] = (
        "1"
        if (
            info["source_transcript"].startswith("FBtr")
            and info["source_tx_pos_1based"].isdigit()
            and info["source_cds_pos_1based"].isdigit()
            and info["source_cds_length"].isdigit()
            and info["source_aa_pos_1based"].isdigit()
            and info["source_aa_length"].isdigit()
        )
        else "0"
    )
    return info


def parse_matching_source_cds_ann(row: Dict[str, str]) -> Optional[Dict[str, str]]:
    row_type = clean(row.get("type") or row.get("variant_effect"))
    if row_type not in CODING_EFFECTS:
        return None

    ann = clean(row.get("ann"))
    records = [r for r in ann.split(",") if r] if ann else []
    for record in records:
        info = parse_ann_record(record, row)
        if info["variant_effect"] == row_type and info["has_source_cds"] == "1":
            return info
    return None


def edited_base(row: Dict[str, str]) -> str:
    change2 = clean(row.get("change2"))
    if ">" in change2:
        return change2.split(">", 1)[1]
    change = clean(row.get("change"))
    if ">" in change:
        return change.split(">", 1)[1]
    return clean(row.get("Alt") or row.get("base"))


def run(cmd: List[str]) -> Tuple[int, str]:
    done = subprocess.run(
        cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    if done.returncode:
        raise RuntimeError("Command failed: " + " ".join(cmd) + "\n" + done.stdout)
    return done.returncode, done.stdout


def load_mapped(path: Path) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            return out
        for row in reader:
            q = clean(row.get("query_name"))
            if q and q not in out:
                out[q] = row
    return out


def to_1based_zero_start(value: str, missing: str) -> str:
    value = clean(value)
    if not value:
        return missing
    try:
        return str(int(value) + 1)
    except ValueError:
        return missing


def source_coordinate(row: Dict[str, str]) -> str:
    species_coord = clean(row.get("species_coord"))
    if species_coord:
        return species_coord
    handle = clean(row.get("handle"))
    if handle:
        return handle
    return clean(row.get("chr")) + "_" + clean(row.get("pos"))


def site_key(row: Dict[str, str]) -> str:
    return clean(row.get("handle")) or source_coordinate(row)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    out_full = Path(args.output_full)
    out_cdspos = Path(args.output_cdspos)
    summary_path = (
        Path(args.summary)
        if args.summary
        else out_full.with_suffix(out_full.suffix + ".summary.tsv")
    )
    out_full.parent.mkdir(parents=True, exist_ok=True)
    out_cdspos.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, str]] = []
    qnames: List[str] = []
    meta: Dict[str, Dict[str, str]] = {}
    bed_rows: List[Tuple[str, int, int, str]] = []
    stats = defaultdict(int)

    with input_path.open("r", encoding="utf-8") as f:
        header_line = f.readline()
        if not header_line:
            raise ValueError(f"Empty input: {input_path}")
        header = header_line.rstrip("\n").split("\t")
        reader = csv.DictReader(f, fieldnames=header, delimiter="\t")
        for raw_i, row in enumerate(reader, start=1):
            stats["input_rows"] += 1
            row = {k: clean(v) for k, v in row.items()}
            info = parse_matching_source_cds_ann(row)
            if info is None:
                stats["skip_not_source_cds_coding_ann"] += 1
                continue

            q = f"{site_key(row) or 'row'}#{raw_i:08d}"
            rows.append(row)
            qnames.append(q)
            meta[q] = {
                **info,
                "genomic_key": source_coordinate(row),
                "base": edited_base(row),
                "site_key": site_key(row),
                "parseable": "1",
            }
            start0 = int(info["source_tx_pos_1based"]) - 1
            bed_rows.append((info["source_transcript"], start0, start0 + 1, q))

    tmp = Path(args.tmp_dir) / f"projection_strict_{args.group}"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True, exist_ok=True)

    bed_path = tmp / "queries.bed"
    with bed_path.open("w", encoding="utf-8") as f:
        f.write("transcript\tstart\tend\tname\n")
        for tx, s, e, q in bed_rows:
            f.write(f"{tx}\t{s}\t{e}\t{q}\n")

    prefix = tmp / "mapped"
    code, log = run(
        [
            sys.executable,
            args.map_script,
            "--gtf",
            args.gtf,
            "--transcript-bed",
            str(bed_path),
            "--out-prefix",
            str(prefix),
            "--include-source",
        ]
    )
    (tmp / "map_transcript_coordinates.log").write_text(log, encoding="utf-8")

    blocks = Path(str(prefix) + ".other_transcript_blocks.tsv")
    dedup_prefix = tmp / "mapped_longest"
    mapped: Dict[str, Dict[str, str]] = {}
    if blocks.exists() and blocks.stat().st_size > 0:
        run(
            [
                sys.executable,
                args.longest_script,
                "--mapped-tsv",
                str(blocks),
                "--longest-list",
                args.longest_list,
                "--out-prefix",
                str(dedup_prefix),
                "--fallback",
                args.fallback,
            ]
        )
        mapped = load_mapped(Path(str(dedup_prefix) + ".dedup_longest.tsv"))

    seen_target = set()
    with out_full.open("w", encoding="utf-8") as full_f, out_cdspos.open(
        "w", encoding="utf-8"
    ) as cds_f:
        full_w = csv.writer(full_f, delimiter="\t", lineterminator="\n")
        cds_w = csv.writer(cds_f, delimiter="\t", lineterminator="\n")
        full_w.writerow(header + APPENDED)
        cds_w.writerow(CDSPOS_HEADER)

        for row, q in zip(rows, qnames):
            info = meta[q]
            m = mapped.get(q)
            if m:
                stats["mapped"] += 1
                new_tx = clean(m.get("target_transcript")) or args.missing_value
                new_tx_pos = to_1based_zero_start(
                    clean(m.get("target_tx_start")), args.missing_value
                )
                target_region = clean(m.get("target_region")) or args.missing_value
                new_cds_pos = (
                    to_1based_zero_start(
                        clean(m.get("target_cds_start")), args.missing_value
                    )
                    if target_region == "CDS"
                    else args.missing_value
                )
            else:
                stats["unmapped"] += 1
                new_tx = new_tx_pos = target_region = new_cds_pos = args.missing_value

            full_w.writerow(
                [row.get(c, "") for c in header]
                + [new_tx, new_tx_pos, target_region, new_cds_pos]
            )

            if not (target_region == "CDS" and new_cds_pos.isdigit()):
                stats["skip_non_target_cds"] += 1
                continue
            target_key = (new_tx, new_cds_pos)
            if target_key in seen_target:
                stats["skip_duplicate_target"] += 1
                continue
            seen_target.add(target_key)

            cds_w.writerow(
                [
                    info["genomic_key"],
                    info["base"],
                    info["variant_effect"],
                    info["source_transcript"],
                    info["source_tx_pos_1based"],
                    info["source_tx_length"] or args.missing_value,
                    info["source_cds_pos_1based"] or args.missing_value,
                    info["source_cds_length"] or args.missing_value,
                    info["source_aa_pos_1based"] or args.missing_value,
                    info["source_aa_length"] or args.missing_value,
                    info["site_key"],
                    new_tx,
                    new_tx_pos,
                    target_region,
                    new_cds_pos,
                ]
            )
            stats["cdspos"] += 1

    with summary_path.open("w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"group\t{args.group}\n")
        f.write(f"input\t{input_path}\n")
        f.write(f"input_rows\t{stats['input_rows']}\n")
        f.write(f"source_cds_coding_rows\t{len(rows)}\n")
        f.write(
            f"skip_not_source_cds_coding_ann\t{stats['skip_not_source_cds_coding_ann']}\n"
        )
        f.write(f"bed_queries\t{len(bed_rows)}\n")
        f.write(f"mapped_rows\t{stats['mapped']}\n")
        f.write(f"unmapped_rows\t{stats['unmapped']}\n")
        f.write(f"skip_non_target_cds\t{stats['skip_non_target_cds']}\n")
        f.write(f"skip_duplicate_target\t{stats['skip_duplicate_target']}\n")
        f.write(f"cdspos_rows\t{stats['cdspos']}\n")
        f.write(f"map_exit_code\t{code}\n")
        f.write(f"tmp_dir\t{tmp}\n")
        f.write(f"output_full\t{out_full}\n")
        f.write(f"output_cdspos\t{out_cdspos}\n")

    print(f"[INFO] group={args.group}")
    print(f"[INFO] input_rows={stats['input_rows']}")
    print(f"[INFO] source_cds_coding_rows={len(rows)}")
    print(f"[INFO] mapped_rows={stats['mapped']}")
    print(f"[INFO] skip_non_target_cds={stats['skip_non_target_cds']}")
    print(f"[INFO] skip_duplicate_target={stats['skip_duplicate_target']}")
    print(f"[INFO] cdspos_rows={stats['cdspos']}")
    print(f"[INFO] output_cdspos={out_cdspos}")
    print(f"[INFO] summary={summary_path}")

    if not args.keep_temp:
        shutil.rmtree(tmp)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
