#!/usr/bin/env python3
"""Reproduce Figure S19 upstream/downstream peptide-abundance analysis.

The aggregation order is fixed and explicit:

    sample log2(downstream/upstream)
      -> arithmetic mean across available replicates within sex
      -> equal-weight arithmetic mean of the female and male means.

If only one sex has a qualifying observation for a gene, that sex-specific
mean is retained. Missing gene-sample observations are not imputed.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, wilcoxon

import fs19_site_selection as site_utils

UNIT_ORDER = ["Dmel-ISO1", "Dmel-WBg", "Dsim", "Dsec", "Dere"]
GROUP_ORDER = ["Random", "Edited"]
COLORS = {"Random": "#4C78A8", "Edited": "#B2182B"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--ms-run-root",
        type=Path,
        default=Path("/data/dongxn/MS/DIA_NN_RUN"),
        help="Root containing the retained v7 deduplicated peptide tables and controls.",
    )
    parser.add_argument(
        "--cdspos-dir",
        type=Path,
        default=Path(
            "/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/"
            "03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_"
            "sourceCDS_keepTargetUTR3_fixed"
        ),
        help="Directory containing updated counter-rule longest-CDS site tables.",
    )
    parser.add_argument(
        "--ref-dir",
        type=Path,
        default=Path("/path/to/RNA_editing_DATAprocessing_new/ref_genome"),
        help="Directory containing the four original all-CDS FASTA files.",
    )
    parser.add_argument("--random-multiplier", type=int, default=10)
    return parser.parse_args()


def build_config(args: argparse.Namespace) -> dict:
    ms = args.ms_run_root
    sites = args.cdspos_dir
    refs = args.ref_dir
    return {
        "Dmel": {
            "cds": refs / "dmel-all-CDS-r6.12.fasta",
            "sites": sites
            / "ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
            "dedup": ms
            / "Dmel/WBg_ag_as_ISO1/editing_updown_v7.Dmel.deduplicated_peptide_intensity.tsv",
            "controls": ms
            / "Dmel/WBg_ag_as_ISO1/editing_updown_v7.Dmel.control_sites.seed_20260726.tsv",
        },
        "Dsim": {
            "cds": refs / "dsim-all-CDS-r2.02.fasta",
            "sites": sites
            / "sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
            "dedup": ms
            / "Dsim/Dsim_diann/editing_vs_random_updown_v7.Dsim.deduplicated_peptide_intensity.tsv",
            "controls": ms
            / "Dsim/Dsim_diann/editing_vs_random_updown_v7.Dsim.control_sites.seed_20260726.tsv",
        },
        "Dsec": {
            "cds": refs / "dsec-all-CDS-r1.3.fasta",
            "sites": sites
            / "sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
            "dedup": ms
            / "Dsec/Dsec_diann/editing_vs_random_updown_v7.Dsec.deduplicated_peptide_intensity.tsv",
            "controls": ms
            / "Dsec/Dsec_diann/editing_vs_random_updown_v7.Dsec.control_sites.seed_20260726.tsv",
        },
        "Dere": {
            "cds": refs / "dere-all-CDS-r1.05.fasta",
            "sites": sites
            / "ere.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv",
            "dedup": ms
            / "Dere/Dere_diann/editing_vs_random_updown_v7.Dere.deduplicated_peptide_intensity.tsv",
            "controls": ms
            / "Dere/Dere_diann/editing_vs_random_updown_v7.Dere.control_sites.seed_20260726.tsv",
        },
    }


def validate_inputs(config: dict) -> None:
    missing = []
    for species, cfg in config.items():
        for label, path in cfg.items():
            if not Path(path).is_file():
                missing.append(f"{species}\t{label}\t{path}")
    if missing:
        raise FileNotFoundError("Missing Figure S19 inputs:\n" + "\n".join(missing))


def prepare_updated_sites(species: str, cfg: dict):
    selected, selection, _ = site_utils.select_longest_cds(Path(cfg["cds"]))
    tx_to_gene = dict(zip(selection["transcript_id"], selection["gene_id"]))
    raw = site_utils.read_table_auto(Path(cfg["sites"])).copy()
    raw["_resolved_gene_id"] = site_utils.gene_ids_for_rows(
        raw, "gene", "new_transcript", tx_to_gene
    )
    raw["_is_cds"] = raw["target_region"].map(site_utils.is_cds_region)
    all_cds_edit_genes = set(
        raw.loc[
            raw["_is_cds"] & (raw["_resolved_gene_id"] != ""),
            "_resolved_gene_id",
        ].astype(str)
    )
    edited, qc, all_valid = site_utils.prepare_editing_sites(
        raw.loc[raw["_is_cds"]].copy(),
        selected,
        "new_transcript",
        "gene",
        "new_cds_pos_1based",
        "ratio",
        "lib",
        "cds_nt_1based",
        True,
        "one-per-gene",
    )
    edited = edited.copy()
    edited["species"] = species
    edited["group"] = "Edited"
    edited["source"] = "updated_counter_rule_sites"
    return selected, edited, qc, all_valid, all_cds_edit_genes


def find_occurrences(sequence: str, peptide: str):
    occurrences = []
    start = 0
    while True:
        pos = sequence.find(peptide, start)
        if pos < 0:
            return occurrences
        occurrences.append((pos, pos + len(peptide)))
        start = pos + 1


def classify(sequence: str, peptide: str, site_aa_1based: int):
    occurrences = find_occurrences(sequence, peptide)
    if not occurrences:
        return "not_found"
    site0 = int(site_aa_1based) - 1
    regions = set()
    for start0, end0 in occurrences:
        if end0 - 1 < site0:
            regions.add("upstream")
        elif start0 > site0:
            regions.add("downstream")
        else:
            regions.add("spanning")
    return next(iter(regions)) if len(regions) == 1 else "ambiguous"


def sample_meta(species: str, sample: str):
    text = str(sample)
    sex_match = re.search(r"[_-]([FM])(?:\d+)?(?:$|[_.-])", text, flags=re.I)
    sex = sex_match.group(1).upper() if sex_match else "Unknown"
    if species == "Dmel":
        if re.search(r"(?:^|_)iso(?:_|$)", text, flags=re.I):
            unit = "Dmel-ISO1"
        elif re.search(r"WBg", text, flags=re.I):
            unit = "Dmel-WBg"
        else:
            unit = "Exclude-Dmel-AG-AS"
    else:
        unit = species
    return unit, sex


def choose_site_definitions(
    species: str,
    cfg: dict,
    edited: pd.DataFrame,
    all_cds_edit_genes: set,
    dedup: pd.DataFrame,
    multiplier: int,
):
    available_tx = set(dedup["transcript_id"].dropna().astype(str))
    edited = edited.copy()
    edited["old_quant_available"] = (
        edited["transcript_id"].astype(str).isin(available_tx)
    )
    missing = edited.loc[~edited["old_quant_available"]].copy()
    analyzable = edited.loc[edited["old_quant_available"]].copy()

    controls = pd.read_csv(cfg["controls"], sep="\t", low_memory=False)
    controls = controls.loc[
        ~controls["gene_id"].astype(str).isin(all_cds_edit_genes)
    ].copy()
    controls = controls.drop_duplicates("gene_id", keep="first")
    max_edited = len(controls) // multiplier
    if len(analyzable) > max_edited:
        level_col = (
            "editing_level_mean" if "editing_level_mean" in analyzable else "gene_id"
        )
        sort_cols = [level_col, "gene_id"] if level_col != "gene_id" else ["gene_id"]
        sort_ascending = [False, True] if level_col != "gene_id" else [True]
        analyzable = (
            analyzable.sort_values(sort_cols, ascending=sort_ascending)
            .head(max_edited)
            .copy()
        )
    controls = controls.head(len(analyzable) * multiplier).copy()
    controls["species"] = species
    controls["group"] = "Random"
    controls["source"] = "old_seed_20260726_control_site"
    controls["old_quant_available"] = True
    keep = [
        "species",
        "site_id",
        "group",
        "gene_id",
        "transcript_id",
        "site_nt_1based",
        "site_aa_1based",
        "source",
        "old_quant_available",
    ]
    return pd.concat([analyzable[keep], controls[keep]], ignore_index=True), missing


def rebuild_one_species(species: str, cfg: dict, multiplier: int):
    selected, edited, _, _, all_cds_edit_genes = prepare_updated_sites(species, cfg)
    columns = ["gene_id", "transcript_id", "Stripped.Sequence", "sample", "intensity"]
    dedup = pd.read_csv(cfg["dedup"], sep="\t", usecols=columns, low_memory=False)
    dedup["intensity"] = pd.to_numeric(dedup["intensity"], errors="coerce")
    dedup = dedup.loc[np.isfinite(dedup["intensity"]) & (dedup["intensity"] > 0)].copy()
    dedup = dedup.groupby(
        ["gene_id", "transcript_id", "Stripped.Sequence", "sample"],
        as_index=False,
        dropna=False,
    )["intensity"].max()

    sites, missing = choose_site_definitions(
        species, cfg, edited, all_cds_edit_genes, dedup, multiplier
    )
    merged = dedup.merge(
        sites,
        on=["gene_id", "transcript_id"],
        how="inner",
        validate="many_to_one",
    )
    proteins = {
        record.transcript_id: record.protein_sequence for record in selected.values()
    }
    unique_keys = merged[
        ["transcript_id", "Stripped.Sequence", "site_aa_1based"]
    ].drop_duplicates()
    region_map = {}
    for row in unique_keys.itertuples(index=False):
        sequence = proteins.get(str(row.transcript_id), "")
        region_map[(str(row.transcript_id), str(row._1), int(row.site_aa_1based))] = (
            classify(sequence, str(row._1), int(row.site_aa_1based))
        )
    merged["region"] = [
        region_map[(str(tx), str(peptide), int(position))]
        for tx, peptide, position in zip(
            merged["transcript_id"],
            merged["Stripped.Sequence"],
            merged["site_aa_1based"],
        )
    ]
    mapping_qc = merged.groupby(["species", "group", "region"], as_index=False).agg(
        n_peptide_sample_rows=("intensity", "size"),
        n_genes=("gene_id", "nunique"),
    )
    mapped = merged.loc[merged["region"].isin(["upstream", "downstream"])].copy()
    group_keys = [
        "species",
        "site_id",
        "group",
        "gene_id",
        "transcript_id",
        "site_nt_1based",
        "site_aa_1based",
        "sample",
        "region",
    ]
    aggregated = mapped.groupby(group_keys, as_index=False).agg(
        top1_intensity=("intensity", "max"),
        n_peptides=("Stripped.Sequence", "nunique"),
    )
    index_keys = group_keys[:-1]
    values = aggregated.pivot_table(
        index=index_keys,
        columns="region",
        values="top1_intensity",
        aggfunc="first",
    ).reset_index()
    counts = (
        aggregated.pivot_table(
            index=index_keys,
            columns="region",
            values="n_peptides",
            aggfunc="first",
        )
        .reset_index()
        .rename(
            columns={
                "upstream": "n_upstream_peptides",
                "downstream": "n_downstream_peptides",
            }
        )
    )
    values = values.merge(counts, how="outer")
    for column in ["upstream", "downstream"]:
        if column not in values:
            values[column] = np.nan
    for column in ["n_upstream_peptides", "n_downstream_peptides"]:
        if column not in values:
            values[column] = 0
        values[column] = values[column].fillna(0).astype(int)
    values["n_total_peptides"] = (
        values["n_upstream_peptides"] + values["n_downstream_peptides"]
    )
    values = values.loc[
        (values["n_upstream_peptides"] >= 1)
        & (values["n_downstream_peptides"] >= 1)
        & (values["n_total_peptides"] >= 2)
        & (values["upstream"] > 0)
        & (values["downstream"] > 0)
    ].copy()
    values["log2_upstream"] = np.log2(values["upstream"])
    values["log2_downstream"] = np.log2(values["downstream"])
    values["log2_downstream_upstream"] = (
        values["log2_downstream"] - values["log2_upstream"]
    )
    metadata = [sample_meta(species, sample) for sample in values["sample"]]
    values["unit"] = [x[0] for x in metadata]
    values["sex"] = [x[1] for x in metadata]
    values = values.loc[
        ~values["unit"].str.startswith("Exclude") & values["sex"].isin(["F", "M"])
    ].copy()
    qc = {
        "species": species,
        "updated_valid_edited_genes": int(edited["gene_id"].nunique()),
        "updated_edited_genes_with_old_quant": int(
            sites.loc[sites.group.eq("Edited"), "gene_id"].nunique()
        ),
        "updated_edited_genes_missing_old_quant": int(len(missing)),
        "random_control_genes_retained": int(
            sites.loc[sites.group.eq("Random"), "gene_id"].nunique()
        ),
        "random_multiplier_realized": float(
            sites.loc[sites.group.eq("Random"), "gene_id"].nunique()
            / max(1, sites.loc[sites.group.eq("Edited"), "gene_id"].nunique())
        ),
        "sample_rows_passing_up_down_filters": int(len(values)),
    }
    return sites, missing, mapping_qc, values, qc


def hierarchical_merge(sample_values: pd.DataFrame):
    keys = [
        "species",
        "unit",
        "site_id",
        "group",
        "gene_id",
        "transcript_id",
        "site_nt_1based",
        "site_aa_1based",
    ]
    sex_values = sample_values.groupby(keys + ["sex"], as_index=False).agg(
        sex_mean_log2_upstream=("log2_upstream", "mean"),
        sex_mean_log2_downstream=("log2_downstream", "mean"),
        sex_mean_log2_ratio=("log2_downstream_upstream", "mean"),
        n_reps=("sample", "nunique"),
    )
    unit_values = sex_values.groupby(keys, as_index=False).agg(
        mean_log2_upstream=("sex_mean_log2_upstream", "mean"),
        mean_log2_downstream=("sex_mean_log2_downstream", "mean"),
        mean_log2_downstream_upstream=("sex_mean_log2_ratio", "mean"),
        n_sexes=("sex", "nunique"),
        n_rep_sex_values=("n_reps", lambda x: ";".join(map(str, x))),
    )
    unit_values["ratio_identity_error"] = (
        unit_values["mean_log2_downstream"]
        - unit_values["mean_log2_upstream"]
        - unit_values["mean_log2_downstream_upstream"]
    )
    if unit_values.duplicated(["unit", "gene_id"]).any():
        raise RuntimeError(
            "More than one final observation remained for a gene within a unit"
        )
    return sex_values, unit_values


def fmt_p(p_value):
    if not np.isfinite(p_value):
        return "NA"
    return f"{p_value:.2e}" if p_value < 1e-4 else f"{p_value:.4f}"


def paired_p(data: pd.DataFrame):
    if len(data) == 0:
        return np.nan
    differences = (data["mean_log2_downstream"] - data["mean_log2_upstream"]).to_numpy(
        float
    )
    if np.allclose(differences, 0):
        return 1.0
    return float(
        wilcoxon(differences, alternative="two-sided", zero_method="wilcox").pvalue
    )


def compute_stats(unit_values: pd.DataFrame):
    rows = []
    for unit in UNIT_ORDER:
        unit_data = unit_values.loc[unit_values.unit.eq(unit)]
        row = {"unit": unit}
        for group in GROUP_ORDER:
            group_data = unit_data.loc[unit_data.group.eq(group)]
            row[f"{group}_n"] = len(group_data)
            row[f"{group}_paired_Wilcoxon_p"] = paired_p(group_data)
            row[f"{group}_median_log2_ratio"] = group_data[
                "mean_log2_downstream_upstream"
            ].median()
        edited = unit_data.loc[
            unit_data.group.eq("Edited"), "mean_log2_downstream_upstream"
        ].dropna()
        random = unit_data.loc[
            unit_data.group.eq("Random"), "mean_log2_downstream_upstream"
        ].dropna()
        row["edited_vs_random_KS_p"] = (
            float(ks_2samp(edited, random).pvalue)
            if len(edited) and len(random)
            else np.nan
        )
        rows.append(row)
    return pd.DataFrame(rows)


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)


def scatter_panel(ax, data: pd.DataFrame, unit: str, stats: pd.Series, compact=True):
    for group in GROUP_ORDER:
        subset = data.loc[data.group.eq(group)]
        ax.scatter(
            subset.mean_log2_upstream,
            subset.mean_log2_downstream,
            s=9 if group == "Random" else 15,
            c=COLORS[group],
            alpha=0.58 if group == "Random" else 0.82,
            edgecolors="none",
            label=f"{group} (n={len(subset)})",
        )
    values = data[["mean_log2_upstream", "mean_log2_downstream"]].to_numpy(float)
    values = values[np.isfinite(values)]
    if len(values):
        lower, upper = values.min(), values.max()
        pad = max(0.2, (upper - lower) * 0.04)
        ax.plot(
            [lower - pad, upper + pad],
            [lower - pad, upper + pad],
            "--",
            lw=0.7,
            c="#777777",
        )
        ax.set_xlim(lower - pad, upper + pad)
        ax.set_ylim(lower - pad, upper + pad)
    ax.set_title(unit, fontsize=9, weight="bold")
    ax.set_xlabel("Mean log2 upstream top1 intensity", fontsize=7.5)
    ax.set_ylabel("Mean log2 downstream top1 intensity", fontsize=7.5)
    text = (
        f"Edited: paired Wilcoxon P={fmt_p(stats['Edited_paired_Wilcoxon_p'])}\n"
        f"Random: paired Wilcoxon P={fmt_p(stats['Random_paired_Wilcoxon_p'])}"
    )
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        va="top",
        fontsize=6.3,
        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=2),
    )
    style(ax)


def hist_panel(ax, data: pd.DataFrame, unit: str, stats: pd.Series):
    all_values = data.mean_log2_downstream_upstream.dropna()
    bins = (
        np.linspace(all_values.min(), all_values.max(), 25)
        if len(all_values) and all_values.max() > all_values.min()
        else 20
    )
    for group in GROUP_ORDER:
        values = data.loc[
            data.group.eq(group), "mean_log2_downstream_upstream"
        ].dropna()
        ax.hist(
            values,
            bins=bins,
            density=True,
            alpha=0.48,
            color=COLORS[group],
            label=f"{group} (n={len(values)})",
        )
    ax.axvline(0, ls="--", lw=0.8, c="#555555")
    ax.set_title(unit, fontsize=9, weight="bold")
    ax.set_xlabel("Mean log2(downstream/upstream)", fontsize=7.5)
    ax.set_ylabel("Density", fontsize=7.5)
    ax.text(
        0.02,
        0.98,
        f"Two-sample KS P={fmt_p(stats['edited_vs_random_KS_p'])}",
        transform=ax.transAxes,
        va="top",
        fontsize=6.5,
        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=2),
    )
    style(ax)


def save_plots(unit_values: pd.DataFrame, stats: pd.DataFrame, out_dir: Path):
    scatter_dir = out_dir / "individual_scatter"
    histogram_dir = out_dir / "individual_ratio_histogram"
    scatter_dir.mkdir(parents=True, exist_ok=True)
    histogram_dir.mkdir(parents=True, exist_ok=True)
    for unit in UNIT_ORDER:
        data = unit_values.loc[unit_values.unit.eq(unit)]
        stat = stats.loc[stats.unit.eq(unit)].iloc[0]
        fig, ax = plt.subplots(figsize=(4.5, 4.0))
        scatter_panel(ax, data, unit, stat, False)
        ax.legend(frameon=False, fontsize=7, loc="lower right")
        fig.tight_layout()
        fig.savefig(scatter_dir / f"{unit}.upstream_vs_downstream.scatter.pdf")
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(4.5, 3.5))
        hist_panel(ax, data, unit, stat)
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(histogram_dir / f"{unit}.log2_downstream_upstream.histogram.pdf")
        plt.close(fig)

    for _, drawer, filename in [
        (
            "scatter",
            scatter_panel,
            "all_5_units.upstream_vs_downstream.A4_landscape.pdf",
        ),
        (
            "hist",
            hist_panel,
            "all_5_units.log2_downstream_upstream_histogram.A4_landscape.pdf",
        ),
    ]:
        fig, axes = plt.subplots(2, 3, figsize=(11.69, 8.27))
        for ax, unit in zip(axes.ravel(), UNIT_ORDER):
            data = unit_values.loc[unit_values.unit.eq(unit)]
            stat = stats.loc[stats.unit.eq(unit)].iloc[0]
            drawer(ax, data, unit, stat)
        axes.ravel()[-1].axis("off")
        handles, _ = axes.ravel()[0].get_legend_handles_labels()
        fig.legend(
            handles,
            GROUP_ORDER,
            frameon=False,
            fontsize=8,
            loc="lower right",
            bbox_to_anchor=(0.93, 0.08),
        )
        fig.tight_layout(rect=(0.03, 0.04, 0.98, 0.98))
        fig.savefig(out_dir / filename)
        plt.close(fig)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    config = build_config(args)
    validate_inputs(config)
    site_tables = []
    missing_tables = []
    mapping_qc_tables = []
    sample_tables = []
    run_qc_rows = []
    for species, cfg in config.items():
        sites, missing, mapping_qc, sample_values, qc = rebuild_one_species(
            species, cfg, args.random_multiplier
        )
        site_tables.append(sites)
        missing_tables.append(missing)
        mapping_qc_tables.append(mapping_qc)
        sample_tables.append(sample_values)
        run_qc_rows.append(qc)

    sample_values = pd.concat(sample_tables, ignore_index=True)
    sex_values, unit_values = hierarchical_merge(sample_values)
    stats = compute_stats(unit_values)
    pd.concat(site_tables, ignore_index=True).to_csv(
        args.out_dir / "site_definitions_updated_edit_old_seed_controls.tsv",
        sep="\t",
        index=False,
    )
    pd.concat(missing_tables, ignore_index=True).to_csv(
        args.out_dir / "updated_edit_sites_missing_old_quantification.tsv",
        sep="\t",
        index=False,
    )
    pd.concat(mapping_qc_tables, ignore_index=True).to_csv(
        args.out_dir / "peptide_reclassification_QC.tsv", sep="\t", index=False
    )
    pd.DataFrame(run_qc_rows).to_csv(args.out_dir / "run_QC.tsv", sep="\t", index=False)
    sample_values.to_csv(
        args.out_dir / "sample_level_top1_and_log2_ratio.tsv", sep="\t", index=False
    )
    sex_values.to_csv(
        args.out_dir / "sex_level_rep_mean_log2_values.tsv", sep="\t", index=False
    )
    unit_values.to_csv(
        args.out_dir / "species_unit_FM_equal_mean_log2_values.tsv",
        sep="\t",
        index=False,
    )
    stats.to_csv(args.out_dir / "species_unit_statistics.tsv", sep="\t", index=False)
    save_plots(unit_values, stats, args.out_dir)


if __name__ == "__main__":
    main()
