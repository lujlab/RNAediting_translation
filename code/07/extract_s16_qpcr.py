#!/usr/bin/env python3
"""Extract and validate the released 5/30/60-min S16 Ct wells (stdlib only)."""
import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
FIELDS = ["Region", "Gene", "Group", "Time", "Technical_well", "Ct",
          "Source_file", "Source_cell", "Sample_label"]


def extract(path, region):
    records = []
    with zipfile.ZipFile(path) as z:
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            strings = ["".join(t.text or "" for t in s.findall(".//m:t", NS))
                       for s in root.findall("m:si", NS)]
        root = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        for row in root.findall(".//m:row", NS):
            cells = {}
            for cell in row.findall("m:c", NS):
                ref = cell.attrib["r"]
                val = cell.find("m:v", NS)
                value = val.text if val is not None else ""
                if cell.get("t") == "s":
                    value = strings[int(value)]
                elif cell.get("t") == "inlineStr":
                    value = "".join(t.text or "" for t in cell.findall(".//m:t", NS))
                cells[re.sub(r"\d", "", ref)] = (value, ref)
            for label_col, value_cols in [("A", "BCDEFG"), ("N", "HIJKLM")]:
                label = cells.get(label_col, ("", ""))[0].strip()
                match = re.fullmatch(r"R(\d+)-(\d+)", label)
                if match is None:
                    continue
                for j, col in enumerate(value_cols):
                    value, ref = cells.get(col, ("", ""))
                    ct = float(value)  # Fail on a missing or nonnumeric expected Ct.
                    if not math.isfinite(ct):
                        raise ValueError(f"Invalid Ct at {path}:{ref}")
                    records.append(dict(zip(FIELDS, [
                        region, "renilla" if j < 3 else "firefly", "R" + match[1],
                        int(match[2]), j % 3 + 1, ct, path.name, ref, label])))
    return records


def write_tsv(path, records, fields):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for region, name in [("3prime", "20260413-DECAY_data_R_3.xls.xlsx"),
                         ("5prime", "20260410-DECAY_data_R_5.xls.xlsx")]:
        path = args.input_dir / name
        records.extend(extract(path, region))
    counts = Counter((r["Region"], r["Gene"], r["Group"], r["Time"]) for r in records)
    expected = {(a, b, c, d) for a in ("3prime", "5prime")
                for b in ("firefly", "renilla") for c in ("R0", "R10", "R30", "R50")
                for d in (5, 30, 60)}
    if set(counts) != expected or set(counts.values()) != {3} or len(records) != 144:
        raise ValueError("S16 source data must contain exactly 144 wells: "
                         "2 regions x 2 reporters x 4 ITP groups x 3 times (5/30/60 min) x 3 wells.")
    write_tsv(args.output_dir / "raw_Ct_5_30_60min.tsv", records, FIELDS)
    print("Validated and extracted 144 Ct wells at 5/30/60 min.")



if __name__ == "__main__":
    main()
