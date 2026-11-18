#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Combine PausePred odds-ratio points and peak-fraction bars.

The two panels share one species-sex x-axis.  The bottom x-axis contains only
the sample labels; edited/unedit counts and the raw Fisher exact-test P value
are annotated above the corresponding bar group.
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

GROUP_ORDER = [
    "Dmel_Adar5G1_M",
    "Dmel_ISO1_head_female_ribo",
    "Dmel_ISO1_head_male_ribo",
    "WBg_F",
    "WBg_M",
    "Dsec_F",
    "Dsec_M",
    "Dere_F",
    "Dere_M",
    "Dsim_F",
    "Dsim_M",
]

GROUP_LABELS = {
    "Dmel_ISO1_head_female_ribo": r"$\it{Dmel}$ (iso-1), F",
    "Dmel_ISO1_head_male_ribo": r"$\it{Dmel}$ (iso-1), M",
    "WBg_F": r"$\it{Dmel}$ (WBg), F",
    "WBg_M": r"$\it{Dmel}$ (WBg), M",
    "Dmel_Adar5G1_M": r"$\it{Dmel}$ (Adar5G1), M",
    "Dsec_F": r"$\it{Dsec}$, F",
    "Dsec_M": r"$\it{Dsec}$, M",
    "Dere_F": r"$\it{Dere}$, F",
    "Dere_M": r"$\it{Dere}$, M",
    "Dsim_F": r"$\it{Dsim}$, F",
    "Dsim_M": r"$\it{Dsim}$, M",
}

GROUP_AXIS_LABELS = {
    "Dmel_ISO1_head_female_ribo": "$\\it{Dmel}$\n(iso-1), F",
    "Dmel_ISO1_head_male_ribo": "$\\it{Dmel}$\n(iso-1), M",
    "WBg_F": "$\\it{Dmel}$\n(WBg), F",
    "WBg_M": "$\\it{Dmel}$\n(WBg), M",
    "Dmel_Adar5G1_M": "$\\it{Dmel}$\n(Adar5G1), M",
    "Dsec_F": "$\\it{Dsec}$\nF",
    "Dsec_M": "$\\it{Dsec}$\nM",
    "Dere_F": "$\\it{Dere}$\nF",
    "Dere_M": "$\\it{Dere}$\nM",
    "Dsim_F": "$\\it{Dsim}$\nF",
    "Dsim_M": "$\\it{Dsim}$\nM",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot combined PausePred odds ratio and fraction panels."
    )
    parser.add_argument("--indir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--windows", type=float, default=2)
    parser.add_argument("--offset", type=float, default=5)
    parser.add_argument(
        "--significance-column",
        choices=["p_value", "p_adj_BH"],
        default="p_value",
        help="Column used for dark-color significance; default is raw Fisher P.",
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def read_selected_table(indir, outdir, window_nt, cutoff):
    merged = Path(outdir) / "all_samples.pausepred_enrichment.merged.tsv"
    if not merged.exists():
        merged = Path(indir) / "all_samples.pausepred_enrichment.merged.tsv"
    if not merged.exists():
        raise FileNotFoundError("Missing merged table: " + str(merged))

    df = pd.read_csv(merged, sep="\t", low_memory=False)
    required = [
        "sample",
        "Pause_score_cutoff",
        "window_nt",
        "editing_in_peak",
        "editing_total_used_in_universe",
        "nonediting_in_peak",
        "nonediting_total",
        "odds_ratio",
        "p_value",
    ]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError("Missing columns: " + ", ".join(missing))

    if "p_adj_BH" not in df.columns:
        df["p_adj_BH"] = np.nan

    numeric = [col for col in required if col not in ["sample"]]
    for col in numeric:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df[
        np.isclose(df["Pause_score_cutoff"], cutoff)
        & np.isclose(df["window_nt"], window_nt)
        & np.isfinite(df["odds_ratio"])
    ].copy()
    if df.empty:
        raise ValueError(
            f"No row for Pause_score_cutoff={cutoff:g}, window_nt={window_nt:g}"
        )

    duplicated = df[df.duplicated("sample", keep=False)]["sample"].unique()
    if len(duplicated):
        raise ValueError("Duplicate selected rows: " + ", ".join(duplicated))

    observed = set(df["sample"].astype(str))
    missing_groups = [group for group in GROUP_ORDER if group not in observed]
    if missing_groups:
        raise ValueError("Missing expected groups: " + ", ".join(missing_groups))

    extras = sorted(observed - set(GROUP_ORDER))
    order = GROUP_ORDER + extras
    df["sample"] = df["sample"].astype(str)
    df["sample_order"] = pd.Categorical(df["sample"], categories=order, ordered=True)
    df = df.sort_values("sample_order").reset_index(drop=True)

    df["editing_peak_fraction"] = (
        df["editing_in_peak"] / df["editing_total_used_in_universe"]
    )
    df["nonediting_peak_fraction"] = df["nonediting_in_peak"] / df["nonediting_total"]
    df["editing_peak_percent"] = df["editing_peak_fraction"] * 100
    df["nonediting_peak_percent"] = df["nonediting_peak_fraction"] * 100
    return df, merged


def fmt_int(value):
    if pd.isna(value):
        return "NA"
    return f"{int(round(float(value))):,}"


def fmt_p_math(value):
    if pd.isna(value):
        return "NA"
    value = float(value)
    if value == 0:
        return "0"
    if value >= 1e-3:
        return f"{value:.3g}"
    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10**exponent)
    return f"{mantissa:.2f}\\times10^{{{exponent}}}"


def point_color(row, significance_column, alpha):
    significant = (
        pd.notna(row[significance_column]) and float(row[significance_column]) < alpha
    )
    if float(row["odds_ratio"]) > 1:
        return "#8B0000" if significant else "#F4A3A3"
    return "#08306B" if significant else "#9ECAE1"


def make_combined_plot(df, outdir, window_nt, cutoff, significance_column, alpha, dpi):
    df = df.copy()
    df["plot_color"] = df.apply(
        point_color,
        axis=1,
        significance_column=significance_column,
        alpha=alpha,
    )
    df["significant"] = pd.to_numeric(df[significance_column], errors="coerce") < alpha

    x = np.arange(len(df))
    ratios = df["odds_ratio"].to_numpy(float)
    max_ratio = float(np.nanmax(ratios))
    or_label_pad = max(0.08, max_ratio * 0.06)
    or_label_y = ratios + or_label_pad
    or_ymax = max(2.2, float(np.nanmax(or_label_y)) + max(0.32, max_ratio * 0.14))

    fig, (ax_or, ax_fraction) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(8.27, 11.69),
        gridspec_kw={"height_ratios": [1.0, 1.18]},
    )

    # Top panel: one OR point per species-sex group.
    ax_or.axhline(1, color="#555555", linestyle="--", linewidth=1.0, zorder=1)
    ax_or.scatter(
        x,
        ratios,
        s=82,
        c=df["plot_color"].tolist(),
        marker="o",
        edgecolors="black",
        linewidths=0.45,
        zorder=3,
    )
    for i, row in df.iterrows():
        ax_or.text(
            x[i],
            or_label_y[i],
            f"OR={row['odds_ratio']:.2f}",
            ha="center",
            va="bottom",
            fontsize=7.5,
        )
    ax_or.set_ylabel("Odds ratio\n(edited / unedited)", fontsize=8)
    ax_or.set_ylim(0, or_ymax)
    ax_or.set_title(
        f"PausePred editing enrichment (window_nt={window_nt:g}, "
        f"Pause_score_cutoff={cutoff:g})",
        fontsize=8,
    )
    ax_or.grid(axis="y", color="#DDDDDD", linewidth=0.6, alpha=0.7)
    ax_or.spines["top"].set_visible(False)
    ax_or.spines["right"].set_visible(False)
    ax_or.tick_params(axis="x", labelbottom=False)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#8B0000",
            markeredgecolor="black",
            markersize=6,
            label="OR > 1, P < 0.05",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#F4A3A3",
            markeredgecolor="black",
            markersize=6,
            label="OR > 1, P ≥ 0.05",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#08306B",
            markeredgecolor="black",
            markersize=6,
            label="OR < 1, P < 0.05",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#9ECAE1",
            markeredgecolor="black",
            markersize=6,
            label="OR < 1, P ≥ 0.05",
        ),
    ]
    ax_or.legend(
        handles=legend_handles, frameon=False, fontsize=7.5, loc="upper left", ncol=2
    )

    # Bottom panel: the same rows as grouped edited/unedit peak-fraction bars.
    width = 0.36
    edit_percent = df["editing_peak_percent"].to_numpy(float)
    unedit_percent = df["nonediting_peak_percent"].to_numpy(float)
    ax_fraction.bar(
        x - width / 2,
        edit_percent,
        width,
        color="#4C78A8",
        label="edited",
        zorder=2,
    )
    ax_fraction.bar(
        x + width / 2,
        unedit_percent,
        width,
        color="#BDBDBD",
        label="unedit",
        zorder=2,
    )
    fraction_max = max(float(np.nanmax(edit_percent)), float(np.nanmax(unedit_percent)))
    annotation_y = max(50.0, fraction_max * 1.22)
    ax_fraction.set_ylim(0, max(60.0, annotation_y + 25.0))
    ax_fraction.set_ylabel("Sites in PausePred peak (%)", fontsize=8)
    ax_fraction.set_xlabel("Species × sex group", fontsize=8, labelpad=5)
    ax_fraction.set_xticks(x)
    tick_labels = []
    stat_labels = []
    for _, row in df.iterrows():
        sample = row["sample"]
        species_label = GROUP_AXIS_LABELS.get(sample, GROUP_LABELS.get(sample, sample))
        # Keep the x-axis readable: sample names only, with sex on a second
        # horizontal line. Counts and P are drawn above the corresponding bars.
        species_label = species_label.replace(", F", "\nF").replace(", M", "\nM")
        tick_labels.append(species_label)
        edit_n = (
            f"edited n=\n{fmt_int(row['editing_in_peak'])}/\n"
            f"{fmt_int(row['editing_total_used_in_universe'])}"
        )
        unedit_n = (
            f"unedit n=\n{fmt_int(row['nonediting_in_peak'])}/\n"
            f"{fmt_int(row['nonediting_total'])}"
        )
        p_label = f"$P=$\n${fmt_p_math(row['p_value'])}$"
        stat_labels.append(f"{edit_n}\n{unedit_n}\n{p_label}")
    ax_fraction.set_xticklabels(
        tick_labels,
        fontsize=7.5,
        rotation=0,
        ha="center",
    )
    ax_fraction.tick_params(axis="x", pad=5, length=3)
    for i, label in enumerate(stat_labels):
        ax_fraction.text(
            x[i],
            annotation_y,
            label,
            ha="center",
            va="bottom",
            fontsize=6.5,
            linespacing=1.0,
            clip_on=False,
        )
    ax_fraction.grid(axis="y", color="#DDDDDD", linewidth=0.6, alpha=0.7)
    ax_fraction.spines["top"].set_visible(False)
    ax_fraction.spines["right"].set_visible(False)
    ax_fraction.legend(frameon=False, fontsize=7.5, loc="upper left", ncol=2)

    fig.subplots_adjust(left=0.10, right=0.98, top=0.94, bottom=0.22, hspace=0.08)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pdf_path = (
        outdir
        / f"odds_ratio_fraction_combined.window_{window_nt:g}.cutoff_{cutoff:g}.pdf"
    )
    # Keep the physical page exactly A4 portrait.
    fig.savefig(pdf_path, dpi=dpi)
    plt.close(fig)

    table = df[
        [
            "sample",
            "Pause_score_cutoff",
            "window_nt",
            "odds_ratio",
            "editing_in_peak",
            "editing_total_used_in_universe",
            "nonediting_in_peak",
            "nonediting_total",
            "editing_peak_percent",
            "nonediting_peak_percent",
            "p_value",
            "p_adj_BH",
            "significant",
            "plot_color",
        ]
    ].copy()
    table["group_label"] = table["sample"].map(GROUP_LABELS).fillna(table["sample"])
    table["significance_column"] = significance_column
    table["alpha"] = alpha
    table_path = (
        outdir
        / f"odds_ratio_fraction_combined.window_{window_nt:g}.cutoff_{cutoff:g}.tsv"
    )
    table.to_csv(table_path, sep="\t", index=False)
    return pdf_path, table_path


def main():
    args = parse_args()
    df, source = read_selected_table(args.indir, args.outdir, args.windows, args.offset)
    pdf_path, table_path = make_combined_plot(
        df,
        args.outdir,
        args.windows,
        args.offset,
        args.significance_column,
        args.alpha,
        args.dpi,
    )
    print(f"[INFO] Source table: {source}")
    print(f"[INFO] Selected rows: {len(df)}")
    print(f"[INFO] Significance: {args.significance_column} < {args.alpha:g}")
    print(f"[INFO] Combined PDF written: {pdf_path}")
    print(f"[INFO] Combined table written: {table_path}")


if __name__ == "__main__":
    main()
