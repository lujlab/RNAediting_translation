#!/usr/bin/env python3
"""Accept validated mergefix QC or the original safe-runner completion logs."""

import argparse
import csv
import json
from pathlib import Path
import re


def check_sample(csv_dir, sample, qc_dir):
    csv_file = csv_dir / (sample + ".pausepred.csv")
    if not csv_file.is_file() or csv_file.stat().st_size == 0:
        raise ValueError("Missing or empty PausePred output: " + str(csv_file))
    qc_file = qc_dir / (sample + ".json")
    if qc_file.is_file():
        qc = json.loads(qc_file.read_text())
        st = csv_file.stat()
        if (
            qc.get("sample") != sample
            or qc.get("status") not in ("completed", "reused_completed")
            or qc.get("exit_code") != 0
            or qc.get("csv_size_bytes") != st.st_size
            or qc.get("csv_mtime_ns") != st.st_mtime_ns
            or not re.fullmatch(r"[0-9a-f]{64}", str(qc.get("csv_sha256", "")))
        ):
            raise ValueError(
                "CSV does not match successful mergefix QC: " + str(qc_file)
            )
        # The batch already hashed and checked the complete CSV. Size/mtime ensure
        # this still refers to that result; this is not a fresh checksum scan.
        return "mergefix_qc"
    log = csv_dir / "logs" / (sample + ".pausepred.log")
    if log.is_file():
        text = log.read_text(errors="replace")
        done = r"^\[INFO\] Done: " + re.escape(sample) + r"$"
        if re.search(done, text, re.M) or "[INFO] Skip existing output:" in text:
            return "safe_runner_log"
    raise ValueError("No successful QC or safe-runner log for " + sample)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv-dir", type=Path, required=True)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--sample")
    group.add_argument("--sheet-dir", type=Path)
    ap.add_argument("--qc-dir", type=Path)
    ap.add_argument(
        "--exclude-samples",
        default="",
        help="Exact comma-separated sample names excluded in sheet mode",
    )
    args = ap.parse_args()
    qc_dir = args.qc_dir or args.csv_dir.parent / "qc"
    excluded = {s.strip() for s in args.exclude_samples.split(",") if s.strip()}
    samples = []
    if args.sample:
        samples = [args.sample]
    else:
        for species in ("dmel", "dere", "dsec", "dsim"):
            sheet = args.sheet_dir / ("pausepred_sample_list." + species + ".tsv")
            with sheet.open() as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    sample = row["sample"].strip()
                    if sample and not sample.startswith("#"):
                        if sample in samples:
                            raise ValueError(
                                "Duplicate sample in species sheets: " + sample
                            )
                        samples.append(sample)
        actual = {
            p.name[: -len(".pausepred.csv")]
            for p in args.csv_dir.glob("*.pausepred.csv")
        }
        samples = [s for s in samples if s not in excluded]
        actual -= excluded
        if set(samples) != actual:
            raise ValueError(
                "CSV/sample-sheet mismatch: missing="
                + str(sorted(set(samples) - actual))
                + ", extra="
                + str(sorted(actual - set(samples)))
            )
    if not samples:
        raise ValueError("No samples selected")
    results = {sample: check_sample(args.csv_dir, sample, qc_dir) for sample in samples}
    print(
        json.dumps(
            {"validated_samples": len(results), "completion_evidence": results},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
