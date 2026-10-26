#!/usr/bin/env python3
"""Scale RNA allele counts before table filters; preserve the same input pairs."""

import argparse, csv, hashlib, json, re
from decimal import Decimal
from pathlib import Path


def canonical(sample):
    if re.fullmatch(r"Dsim_[FM][123]", sample):
        return sample.replace("Dsim_", "Dsim_w501_", 1)
    m = re.fullmatch(r"Dmel_ISO1_head_(female|male)_ribo", sample)
    if m:
        return "Dmel_ISO1_" + m.group(1)
    return re.sub(r"_ribo$", "", sample)


def main():
    ap = argparse.ArgumentParser()
    for name in ("input", "manifest", "out"):
        ap.add_argument("--" + name, required=True)
    a = ap.parse_args()
    source = Path(a.input).resolve()
    out = Path(a.out).resolve()
    if out == source or out in source.parents:
        raise ValueError("Input cannot be inside output")
    with open(a.manifest) as f:
        manifest = list(csv.DictReader(f, delimiter="\t"))
    scales = {r["sample_id"]: Decimal(r["rna_scale"]) for r in manifest}
    if len(scales) != len(manifest):
        raise ValueError("Duplicate scaling IDs")
    files = sorted(
        p
        for p in source.glob("*site_df_alledit")
        if re.search(r"\.[APE]site_df_alledit$", p.name)
    )
    if not files:
        raise ValueError("No final A/P/E alledit files")
    samples = {re.sub(r"\.[APE]site_df_alledit$", "", p.name) for p in files}
    unknown = sorted(s for s in samples if canonical(s) not in scales)
    if unknown:
        raise ValueError("Unresolved samples: " + ", ".join(unknown))
    for s in samples:
        modes = {
            re.search(r"\.([APE])site_df_alledit$", p.name).group(1)
            for p in files
            if p.name.startswith(s + ".")
        }
        if modes != {"A", "P", "E"}:
            raise ValueError("Missing A/P/E file: " + s)
    folder = out / "RNAdiv2/alledit"
    folder.mkdir(parents=True, exist_ok=True)
    if {p.name for p in folder.glob("*site_df_alledit")} - {p.name for p in files}:
        raise ValueError("Stale alledit inputs in " + str(folder))
    audit = []
    for p in files:
        sample = re.sub(r"\.[APE]site_df_alledit$", "", p.name)
        scale = scales[canonical(sample)]
        if scale not in (Decimal("0.5"), Decimal("1")):
            raise ValueError("Unsupported RNA factor")
        with p.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            cols = reader.fieldnames
            rows = list(reader)
        required = {
            "ref_number_rna",
            "alt_number_rna",
            "ref_number_ribo",
            "alt_number_ribo",
            "FBt_cdspos",
        }
        if not required.issubset(cols):
            raise ValueError("Missing allele count columns in " + str(p))
        summary = dict(
            sample=sample,
            manifest_sample=canonical(sample),
            site_mode=re.search(r"\.([APE])site", p.name).group(1),
            rna_scale=str(scale),
            rows=len(rows),
            rows_with_fractional_RNA=0,
            old_RNA_total="0",
            new_RNA_total="0",
            input_file=str(p),
            input_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        )
        old_total = Decimal(0)
        new_total = Decimal(0)
        for row in rows:
            ref = Decimal(row["ref_number_rna"])
            alt = Decimal(row["alt_number_rna"])
            if not ref.is_finite() or not alt.is_finite() or min(ref, alt) < 0:
                raise ValueError("Invalid RNA count")
            old_total += ref + alt
            ref *= scale
            alt *= scale
            new_total += ref + alt
            summary["rows_with_fractional_RNA"] += int(ref % 1 != 0 or alt % 1 != 0)
            row["ref_number_rna"] = format(ref, "f")
            row["alt_number_rna"] = format(alt, "f")
            if scale != 1:
                ar = (Decimal(row["alt_number_ribo"]) + 1) / (alt + 1)
                rr = (Decimal(row["ref_number_ribo"]) + 1) / (ref + 1)
                for key, val in [
                    ("alt_ratio", ar),
                    ("ref_ratio", rr),
                    ("gmRC", ar / rr),
                    ("diff", ar - rr),
                ]:
                    if key in row:
                        row[key] = format(val, "f")
        summary["old_RNA_total"] = str(old_total)
        summary["new_RNA_total"] = str(new_total)
        audit.append(summary)
        dest = out / "RNAdiv2/alledit" / p.name
        if dest.is_symlink():
            raise ValueError("Refusing to overwrite symlink: " + str(dest))
        with dest.open("w") as f:
            writer = csv.DictWriter(
                f, fieldnames=cols, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
    pass
    summary = {
        "input": str(source),
        "files": len(files),
        "samples": len(samples),
        "halved_samples": sum(scales[canonical(s)] == Decimal("0.5") for s in samples),
        "unchanged_samples": sum(scales[canonical(s)] == 1 for s in samples),
        "fractional_RNA_rows": sum(r["rows_with_fractional_RNA"] for r in audit),
        "scope": "Scale final alledit RNA allele counts; retain original discovered codon pairs. Reapply Table1 and S4/S5 thresholds downstream. Fisher uses explicit R round; GLM/coverage retain fractions.",
    }
    pass
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
