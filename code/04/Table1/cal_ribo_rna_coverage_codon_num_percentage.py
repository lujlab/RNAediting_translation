#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from pathlib import Path

import pandas as pd


def extract_total(value):
    """
    Extract the leading total count from strings such as:

        634 (525/109)
        443 (363/80)

    Returns:
        634
        443
    """

    if pd.isna(value):
        return None

    value = str(value).strip()

    match = re.match(r"^(\d+)", value)

    if match is None:
        raise ValueError(f"Cannot extract total count from: {value}")

    return int(match.group(1))


def format_percentage(site_total, mrna_total):
    """
    Format as:

        495 (51.91%)

    Trailing zeros are removed:
        48.80% -> 48.8%
        50.00% -> 50%
    """

    if mrna_total is None or mrna_total == 0:
        return f"{site_total} (NA)"

    percentage = site_total / mrna_total * 100

    percentage_text = f"{percentage:.2f}".rstrip("0").rstrip(".")

    return f"{site_total} ({percentage_text}%)"


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Calculate the percentage of mRNA editing sites "
            "covered at the Ribo-seq A, P and E sites."
        )
    )

    parser.add_argument(
        "--input",
        default=(
            "/path/to/"
            "RNA_editing_DATAprocessing/all_code/"
            "riboseq_editing_site_counts.formatted.tsv"
        ),
        help="Input formatted TSV file.",
    )

    parser.add_argument(
        "--output",
        default=("riboseq_editing_site_counts.CDS_percentage.tsv"),
        help="Output TSV file.",
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)

    if not input_file.is_file():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    # The input file is tab-separated.
    table = pd.read_csv(
        input_file,
        sep="\t",
        dtype=str,
    )

    required_columns = [
        "Library",
        "mRNA-Seq (I-allele count >= 1, N/S)",
        "A-site (N/S)",
        "P-site (N/S)",
        "E-site (N/S)",
        "A/P-site (N/S)",
    ]

    missing_columns = [
        column for column in required_columns if column not in table.columns
    ]

    if missing_columns:
        raise ValueError(
            "The following columns are missing:\n" + "\n".join(missing_columns)
        )

    output_rows = []

    for _, row in table.iterrows():

        library = row["Library"]

        mrna_total = extract_total(row["mRNA-Seq (I-allele count >= 1, N/S)"])

        a_total = extract_total(row["A-site (N/S)"])

        p_total = extract_total(row["P-site (N/S)"])

        e_total = extract_total(row["E-site (N/S)"])

        ap_total = extract_total(row["A/P-site (N/S)"])

        output_rows.append(
            {
                "Library": library,
                "A-site (CDS %)": format_percentage(
                    a_total,
                    mrna_total,
                ),
                "P-site (CDS %)": format_percentage(
                    p_total,
                    mrna_total,
                ),
                "E-site (CDS %)": format_percentage(
                    e_total,
                    mrna_total,
                ),
                "A/P-site (CDS %)": format_percentage(
                    ap_total,
                    mrna_total,
                ),
            }
        )

    output_table = pd.DataFrame(output_rows)

    output_table.to_csv(
        output_file,
        sep="\t",
        index=False,
    )

    print(f"Finished: {output_file.resolve()}")


if __name__ == "__main__":
    main()
