#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# ============================================================
# Variant-effect classification
# ============================================================


def classify_variant_effect(effect):
    """
    Classify SnpEff effects using substring matching.

    Returns:
        N: nonsynonymous
        S: synonymous
        U: unclassified

    Nonsynonymous has higher priority when a compound annotation
    contains both synonymous and nonsynonymous effects.
    """

    if effect is None or pd.isna(effect):
        return "U"

    effect = str(effect).lower()

    nonsynonymous_patterns = [
        "missense_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
    ]

    synonymous_patterns = [
        "synonymous_variant",
        "stop_retained_variant",
        "start_retained_variant",
    ]

    if any(pattern in effect for pattern in nonsynonymous_patterns):
        return "N"

    if any(pattern in effect for pattern in synonymous_patterns):
        return "S"

    return "U"


CLASS_PRIORITY = {
    "U": 0,
    "S": 1,
    "N": 2,
}


def merge_class(old_class, new_class):
    """
    Merge classifications for duplicated site keys.

    Priority:
        N > S > U
    """

    if old_class is None:
        return new_class

    if CLASS_PRIORITY.get(new_class, 0) > CLASS_PRIORITY.get(old_class, 0):
        return new_class

    return old_class


# ============================================================
# Site-key normalization
# ============================================================


def normalize_site_key(value):
    """
    Normalize transcript CDS-position keys.

    Examples:
        FBtr0070151_1950
        FBtr0070151.1_1950

    Transcript version suffixes are removed.
    """

    if value is None or pd.isna(value):
        return None

    value = str(value).strip()

    if not value or value.lower() == "nan":
        return None

    match = re.match(r"^(.+?)_(\d+)$", value)

    if not match:
        return value

    transcript = match.group(1)
    cds_position = match.group(2)

    transcript = re.sub(r"\.\d+$", "", transcript)

    return f"{transcript}_{cds_position}"


# ============================================================
# Annotation-file parsing
# ============================================================


def load_annotation(annotation_file):
    """
    Parse the new headered annotation table.

    Required columns:
        variant_effect
        new_transcript
        target_region
        new_cds_pos_1based

    Site keys are constructed as:

        new_transcript_new_cds_pos_1based

    Example:
        new_transcript       = FBtr0305555
        new_cds_pos_1based   = 1732
        site key             = FBtr0305555_1732

    This site key is matched directly to the FBt_cdspos column
    in the *_alledit files.

    Only rows with target_region == "CDS" and a valid positive
    new_cds_pos_1based value are used.
    """

    annotation_file = Path(annotation_file)

    if not annotation_file.is_file():
        raise FileNotFoundError(f"Annotation file does not exist: {annotation_file}")

    table = pd.read_csv(
        annotation_file,
        sep="\t",
        dtype=str,
        engine="python",
    )

    required_columns = {
        "variant_effect",
        "new_transcript",
        "target_region",
        "new_cds_pos_1based",
    }

    missing_columns = required_columns - set(table.columns)

    if missing_columns:
        raise ValueError(
            f"Missing columns in annotation file {annotation_file}:\n"
            f"{', '.join(sorted(missing_columns))}\n"
            f"Columns found:\n"
            f"{', '.join(map(str, table.columns))}"
        )

    site_effects = {}

    n_input = len(table)
    n_non_cds = 0
    n_missing_transcript = 0
    n_invalid_cds_pos = 0
    n_valid_rows = 0

    for _, row in table.iterrows():

        target_region = str(row.get("target_region", "")).strip()

        if target_region != "CDS":
            n_non_cds += 1
            continue

        transcript = row.get("new_transcript")

        if transcript is None or pd.isna(transcript):
            n_missing_transcript += 1
            continue

        transcript = str(transcript).strip()

        if (
            not transcript
            or transcript.lower() == "nan"
            or not transcript.startswith("FBtr")
        ):
            n_missing_transcript += 1
            continue

        transcript = re.sub(r"\.\d+$", "", transcript)

        cds_pos_text = str(row.get("new_cds_pos_1based", "")).strip()

        if not re.fullmatch(r"\d+", cds_pos_text):
            n_invalid_cds_pos += 1
            continue

        cds_position = int(cds_pos_text)

        if cds_position <= 0:
            n_invalid_cds_pos += 1
            continue

        effect = row.get("variant_effect")
        site_class = classify_variant_effect(effect)

        site_key = normalize_site_key(f"{transcript}_{cds_position}")

        site_effects[site_key] = merge_class(
            site_effects.get(site_key),
            site_class,
        )

        n_valid_rows += 1

    n_count = sum(value == "N" for value in site_effects.values())
    s_count = sum(value == "S" for value in site_effects.values())
    u_count = sum(value == "U" for value in site_effects.values())

    print(
        f"Loaded annotation: {annotation_file}\n"
        f"  input rows: {n_input}\n"
        f"  non-CDS rows skipped: {n_non_cds}\n"
        f"  missing/invalid new_transcript rows: "
        f"{n_missing_transcript}\n"
        f"  invalid new_cds_pos_1based rows: "
        f"{n_invalid_cds_pos}\n"
        f"  valid CDS annotation rows: {n_valid_rows}\n"
        f"  unique transcript-CDS sites: "
        f"{len(site_effects)}\n"
        f"  N/S/U: {n_count}/{s_count}/{u_count}",
        file=sys.stderr,
    )

    return site_effects


def resolve_annotation_file(annotation_dir, filename=None, pattern=None):
    """
    Locate an annotation file by exact filename or glob pattern.
    """

    annotation_dir = Path(annotation_dir)

    if filename is not None:

        candidate = Path(filename)

        if not candidate.is_absolute():
            candidate = annotation_dir / candidate

        if not candidate.is_file():
            raise FileNotFoundError(f"Annotation file does not exist: {candidate}")

        return candidate

    matches = sorted(annotation_dir.glob(pattern))

    if len(matches) == 0:
        raise FileNotFoundError(
            f"No annotation file matched: {annotation_dir / pattern}"
        )

    if len(matches) > 1:

        matched_files = "\n".join(f"  {path}" for path in matches)

        raise RuntimeError(
            f"Multiple annotation files matched {pattern}:\n"
            f"{matched_files}\n"
            f"Please use --dmel-annotation to specify one file."
        )

    return matches[0]


# ============================================================
# Library-species assignment
# ============================================================


def get_species(library):
    """
    Assign each library to its coding-site annotation.

    Dsim libraries -> Dsim
    Dsec libraries -> Dsec
    Dere libraries -> Dere
    all other libraries -> Dmel
    """

    name = library.lower()

    if name.startswith("dsim"):
        return "Dsim"

    if name.startswith("dsec"):
        return "Dsec"

    if name.startswith("dere"):
        return "Dere"

    return "Dmel"


# ============================================================
# A/P/E table processing
# ============================================================


def read_site_file(site_file, annotation_dict):
    """
    Read one A-, P- or E-site table.

    mRNA-site filtering:
        alt_number_rna >= 1

    Ribo-site filtering:
        alt_number_rna >= 1
        and
        ref_number_ribo + alt_number_ribo >= 1

    Returns:
        rna_sites
        ribo_sites

    Each object is:
        dictionary[FBt_cdspos] = N/S/U
    """

    site_file = Path(site_file)

    table = pd.read_csv(
        site_file,
        sep=r"\s+",
        dtype=str,
        engine="python",
    )

    required_columns = {
        "alt_number_rna",
        "ref_number_ribo",
        "alt_number_ribo",
        "FBt_cdspos",
    }

    missing_columns = required_columns - set(table.columns)

    if missing_columns:
        raise ValueError(
            f"Missing columns in {site_file}:\n" f"{', '.join(sorted(missing_columns))}"
        )

    table["site_key"] = table["FBt_cdspos"].map(normalize_site_key)

    numeric_columns = [
        "alt_number_rna",
        "ref_number_ribo",
        "alt_number_ribo",
    ]

    for column in numeric_columns:

        table[column] = pd.to_numeric(
            table[column],
            errors="coerce",
        ).fillna(0)

    table["ribo_depth"] = table["ref_number_ribo"] + table["alt_number_ribo"]

    rna_sites = {}
    ribo_sites = {}

    for _, row in table.iterrows():

        site_key = row["site_key"]

        if site_key is None:
            continue

        site_class = annotation_dict.get(
            site_key,
            "U",
        )

        # mRNA editing site
        if row["alt_number_rna"] >= 1:

            rna_sites[site_key] = merge_class(
                rna_sites.get(site_key),
                site_class,
            )

            # Ribo-seq-covered editing site
            if row["ribo_depth"] >= 1:

                ribo_sites[site_key] = merge_class(
                    ribo_sites.get(site_key),
                    site_class,
                )

    return rna_sites, ribo_sites


def merge_site_dicts(*site_dicts):
    """
    Obtain the unique union of multiple site dictionaries.
    """

    merged = {}

    for site_dict in site_dicts:

        for site_key, site_class in site_dict.items():

            merged[site_key] = merge_class(
                merged.get(site_key),
                site_class,
            )

    return merged


# ============================================================
# Counting and formatting
# ============================================================


def summarize_sites(site_dict):
    """
    Count total, nonsynonymous, synonymous and unclassified sites.
    """

    n_count = sum(site_class == "N" for site_class in site_dict.values())

    s_count = sum(site_class == "S" for site_class in site_dict.values())

    u_count = sum(site_class == "U" for site_class in site_dict.values())

    return {
        "total": len(site_dict),
        "N": n_count,
        "S": s_count,
        "U": u_count,
    }


def format_summary(summary):
    """
    Format results as:

        total (N/S)

    Example:
        950 (772/178)
    """

    return f"{summary['total']} " f"({summary['N']}/{summary['S']})"


# ============================================================
# Main
# ============================================================


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Count mRNA-Seq and Ribo-seq-covered editing sites "
            "at the A, P and E sites."
        )
    )

    parser.add_argument(
        "--site-dir",
        default=(
            "/path/to/"
            "RNA_editing_DATAprocessing/"
            "05.psite/psite_2times/new_alledit"
        ),
        help=(
            "Directory containing "
            "*.Asite_df_alledit, *.Psite_df_alledit and "
            "*.Esite_df_alledit files."
        ),
    )

    parser.add_argument(
        "--annotation-dir",
        default=(
            "/path/to/"
            "RNA_editing_DATAprocessing/all_code/final_data/"
            "rebuild_logic_WBG_ISO1_Published_then_species_FAST_v2_candidate/"
            "postprocess_3steps"
        ),
        help="Directory containing the new headered annotation files.",
    )

    parser.add_argument(
        "--dmel-annotation",
        default=None,
        help=(
            "Exact Dmel annotation filename or full path. "
            "When omitted, the script searches for "
            "ISO1.candidate.rows.tsv.filter.cdspos963.new."
        ),
    )

    parser.add_argument(
        "--dsim-annotation",
        default="sim.candidate.rows.tsv.cdspos818.new",
        help="Exact dsim annotation filename or full path; overrides the legacy default.",
    )

    parser.add_argument(
        "--dsec-annotation",
        default="sec.candidate.rows.tsv.cdspos785.new",
        help="Exact dsec annotation filename or full path; overrides the legacy default.",
    )

    parser.add_argument(
        "--dere-annotation",
        default="ere.candidate.rows.tsv.cdspos807.new",
        help="Exact dere annotation filename or full path; overrides the legacy default.",
    )

    parser.add_argument(
        "--output-prefix",
        default="riboseq_editing_site_counts",
        help="Prefix of the output files.",
    )

    args = parser.parse_args()

    site_dir = Path(args.site_dir)
    annotation_dir = Path(args.annotation_dir)

    if not site_dir.is_dir():
        raise FileNotFoundError(f"Site directory does not exist: {site_dir}")

    if not annotation_dir.is_dir():
        raise FileNotFoundError(
            f"Annotation directory does not exist: {annotation_dir}"
        )

    # --------------------------------------------------------
    # Annotation files
    # --------------------------------------------------------

    annotation_files = {
        "Dsim": resolve_annotation_file(
            annotation_dir,
            filename=args.dsim_annotation,
        ),
        "Dsec": resolve_annotation_file(
            annotation_dir,
            filename=args.dsec_annotation,
        ),
        "Dere": resolve_annotation_file(
            annotation_dir,
            filename=args.dere_annotation,
        ),
    }

    if args.dmel_annotation is not None:

        annotation_files["Dmel"] = resolve_annotation_file(
            annotation_dir,
            filename=args.dmel_annotation,
        )

    else:

        annotation_files["Dmel"] = resolve_annotation_file(
            annotation_dir,
            pattern="ISO1.candidate.rows.tsv.filter.cdspos963.new",
        )

    print("\nAnnotation files:", file=sys.stderr)

    for species, annotation_file in annotation_files.items():
        print(
            f"  {species}: {annotation_file}",
            file=sys.stderr,
        )

    # Load each annotation only once
    annotation_by_species = {}

    for species, annotation_file in annotation_files.items():

        annotation_by_species[species] = load_annotation(annotation_file)

    # --------------------------------------------------------
    # Find A/P/E files
    # --------------------------------------------------------

    suffixes = {
        "A": ".Asite_df_alledit",
        "P": ".Psite_df_alledit",
        "E": ".Esite_df_alledit",
    }

    files_by_library = defaultdict(dict)

    for site_type, suffix in suffixes.items():

        for site_file in sorted(site_dir.glob(f"*{suffix}")):

            # Do not include files such as:
            # WBg_F2.codon_Asite_df_alledit
            if ".codon_" in site_file.name:
                continue

            library = site_file.name[: -len(suffix)]

            files_by_library[library][site_type] = site_file

    if len(files_by_library) == 0:
        raise RuntimeError(f"No A/P/E-site files were found in: {site_dir}")

    print(
        f"\nDetected {len(files_by_library)} libraries.",
        file=sys.stderr,
    )

    formatted_rows = []
    detailed_rows = []

    # --------------------------------------------------------
    # Analyze each library
    # --------------------------------------------------------

    for library in sorted(files_by_library):

        species = get_species(library)
        annotation_dict = annotation_by_species[species]

        print(
            f"Processing {library} using {species} annotation",
            file=sys.stderr,
        )

        rna_sites_by_type = {
            "A": {},
            "P": {},
            "E": {},
        }

        ribo_sites_by_type = {
            "A": {},
            "P": {},
            "E": {},
        }

        for site_type in ["A", "P", "E"]:

            site_file = files_by_library[library].get(site_type)

            if site_file is None:

                print(
                    f"Warning: {library} is missing " f"{site_type}-site file.",
                    file=sys.stderr,
                )

                continue

            rna_sites, ribo_sites = read_site_file(
                site_file,
                annotation_dict,
            )

            rna_sites_by_type[site_type] = rna_sites
            ribo_sites_by_type[site_type] = ribo_sites

        # mRNA sites:
        # union of RNA-edited sites appearing in A/P/E input tables
        mrna_sites = merge_site_dicts(
            rna_sites_by_type["A"],
            rna_sites_by_type["P"],
            rna_sites_by_type["E"],
        )

        a_sites = ribo_sites_by_type["A"]
        p_sites = ribo_sites_by_type["P"]
        e_sites = ribo_sites_by_type["E"]

        # Unique union of A-site and P-site sites
        ap_sites = merge_site_dicts(
            a_sites,
            p_sites,
        )

        summaries = {
            "mRNA": summarize_sites(mrna_sites),
            "A": summarize_sites(a_sites),
            "P": summarize_sites(p_sites),
            "E": summarize_sites(e_sites),
            "AP": summarize_sites(ap_sites),
        }

        formatted_rows.append(
            {
                "Library": library,
                "mRNA-Seq (I-allele count >= 1, N/S)": format_summary(
                    summaries["mRNA"]
                ),
                "A-site (N/S)": format_summary(summaries["A"]),
                "P-site (N/S)": format_summary(summaries["P"]),
                "E-site (N/S)": format_summary(summaries["E"]),
                "A/P-site (N/S)": format_summary(summaries["AP"]),
            }
        )

        detailed_row = {
            "Library": library,
            "annotation_species": species,
        }

        for category in ["mRNA", "A", "P", "E", "AP"]:

            detailed_row[f"{category}_total"] = summaries[category]["total"]

            detailed_row[f"{category}_N"] = summaries[category]["N"]

            detailed_row[f"{category}_S"] = summaries[category]["S"]

            detailed_row[f"{category}_unclassified"] = summaries[category]["U"]

        detailed_rows.append(detailed_row)

    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    formatted_output = Path(f"{args.output_prefix}.formatted.tsv")

    detailed_output = Path(f"{args.output_prefix}.detailed.tsv")

    pd.DataFrame(formatted_rows).to_csv(
        formatted_output,
        sep="\t",
        index=False,
    )

    pd.DataFrame(detailed_rows).to_csv(
        detailed_output,
        sep="\t",
        index=False,
    )

    print(
        "\nFinished.\n"
        f"Formatted output: {formatted_output.resolve()}\n"
        f"Detailed output:  {detailed_output.resolve()}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
