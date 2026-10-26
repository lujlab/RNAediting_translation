#!/usr/bin/env python3
"""Plot combined PausePred editing enrichment for the hairpin analysis.

The input is the summary table written by
pausepred_hairpin_edit_unedit_enrichment.py.  One PDF is written for each
analysis scope and transcript-coordinate window.  The layout follows the
legacy odds_ratio_fraction_combined figure: an odds-ratio panel above a paired
peak-fraction panel.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SAMPLE_ORDER = [
    "Dmel_Adar5G1_M",
    "Dmel_ISO1_head_female_ribo",
    "Dmel_ISO1_head_male_ribo",
    "WBg_F",
    "WBg_M",
    "AdarG_F",
    "AdarG_M",
    "AdarS_F",
    "AdarS_M",
    "Dsec_F",
    "Dsec_M",
    "Dere_F",
    "Dere_M",
    "Dsim_F",
    "Dsim_M",
]

GROUP_LABELS = {
    "AdarG_F": "$\\it{Dmel}$\n(AdarG)\nF",
    "AdarG_M": "$\\it{Dmel}$\n(AdarG)\nM",
    "AdarS_F": "$\\it{Dmel}$\n(AdarS)\nF",
    "AdarS_M": "$\\it{Dmel}$\n(AdarS)\nM",
    "Dmel_Adar5G1_M": "$\\it{Dmel}$\n(Adar5G1)\nM",
    "Dmel_ISO1_head_female_ribo": "$\\it{Dmel}$\n(iso-1)\nF",
    "Dmel_ISO1_head_male_ribo": "$\\it{Dmel}$\n(iso-1)\nM",
    "WBg_F": "$\\it{Dmel}$\n(WBg)\nF",
    "WBg_M": "$\\it{Dmel}$\n(WBg)\nM",
    "Dsec_F": "$\\it{Dsec}$\nF",
    "Dsec_M": "$\\it{Dsec}$\nM",
    "Dere_F": "$\\it{Dere}$\nF",
    "Dere_M": "$\\it{Dere}$\nM",
    "Dsim_F": "$\\it{Dsim}$\nF",
    "Dsim_M": "$\\it{Dsim}$\nM",
}

SCOPE_LABELS = {
    "all_edit_cds": "all-edit CDS background",
    "unedit_A_base": "unedited A-base background",
}

REQUIRED_COLUMNS = {
    "sample",
    "analysis_scope",
    "Pause_score_cutoff",
    "window_nt",
    "edited_total_in_hairpin",
    "edited_in_peak_in_hairpin",
    "unedited_total_in_hairpin",
    "unedited_in_peak_in_hairpin",
    "edited_peak_fraction",
    "unedited_peak_fraction",
    "odds_ratio_editing_in_peak",
    "p_value_one_sided_fisher_greater",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, help="hairpin summary TSV")
    parser.add_argument("--out-dir", required=True, help="PDF/TSV output directory")
    parser.add_argument("--cutoff", type=float, default=5.0)
    parser.add_argument("--windows", default="0,2")
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser.parse_args()


def as_float(row: dict[str, str], column: str) -> float:
    value = row[column]
    if value in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


def as_int(row: dict[str, str], column: str) -> int:
    return int(round(as_float(row, column)))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {', '.join(sorted(missing))}"
            )
        return list(reader)


def is_close(a: float, b: float) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) < 1e-9


def p_text(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    if value == 0:
        return "<1e-300"
    if value < 1e-3:
        return f"{value:.2e}"
    return f"{value:.3g}"


def format_cutoff(value: float) -> str:
    return str(int(value)) if value.is_integer() else str(value)


def color_for(odds_ratio: float, p_value: float, alpha: float) -> str:
    if odds_ratio >= 1:
        return "#8B0000" if p_value < alpha else "#F4A3A3"
    return "#2166AC" if p_value < alpha else "#9ECAE1"


def ordered_subset(
    rows: list[dict[str, str]], scope: str, window: int, cutoff: float
) -> list[dict[str, str]]:
    subset = [
        row
        for row in rows
        if row["analysis_scope"] == scope
        and is_close(as_float(row, "Pause_score_cutoff"), cutoff)
        and int(round(as_float(row, "window_nt"))) == window
    ]
    by_sample = {row["sample"]: row for row in subset}
    expected = [sample for sample in SAMPLE_ORDER if sample in by_sample]
    unknown = sorted(set(by_sample) - set(SAMPLE_ORDER))
    ordered = [by_sample[sample] for sample in expected + unknown]
    if not ordered:
        raise ValueError(
            f"No rows found for scope={scope}, window={window}, cutoff={cutoff}"
        )
    if len(ordered) != len(by_sample):
        raise ValueError("Duplicate sample rows found in the selected subset")
    return ordered


def write_plot_table(
    rows: list[dict[str, str]], output_path: Path, alpha: float
) -> None:
    fieldnames = [
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
        "significant",
        "plot_color",
        "group_label",
        "significance_column",
        "alpha",
        "analysis_scope",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            odds_ratio = as_float(row, "odds_ratio_editing_in_peak")
            p_value = as_float(row, "p_value_one_sided_fisher_greater")
            writer.writerow(
                {
                    "sample": row["sample"],
                    "Pause_score_cutoff": row["Pause_score_cutoff"],
                    "window_nt": row["window_nt"],
                    "odds_ratio": odds_ratio,
                    "editing_in_peak": as_int(row, "edited_in_peak_in_hairpin"),
                    "editing_total_used_in_universe": as_int(
                        row, "edited_total_in_hairpin"
                    ),
                    "nonediting_in_peak": as_int(row, "unedited_in_peak_in_hairpin"),
                    "nonediting_total": as_int(row, "unedited_total_in_hairpin"),
                    "editing_peak_percent": 100.0
                    * as_float(row, "edited_peak_fraction"),
                    "nonediting_peak_percent": 100.0
                    * as_float(row, "unedited_peak_fraction"),
                    "p_value": p_value,
                    "significant": p_value < alpha,
                    "plot_color": color_for(odds_ratio, p_value, alpha),
                    "group_label": GROUP_LABELS.get(row["sample"], row["sample"]),
                    "significance_column": "p_value_one_sided_fisher_greater",
                    "alpha": alpha,
                    "analysis_scope": row["analysis_scope"],
                }
            )


def style_axis(axis: plt.Axes) -> None:
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.6)
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.0)
    axis.spines["bottom"].set_linewidth(1.0)


def make_plot(
    rows: list[dict[str, str]],
    output_pdf: Path,
    scope: str,
    window: int,
    cutoff: float,
    alpha: float,
) -> None:
    x_values = list(range(len(rows)))
    odds_ratios = [as_float(row, "odds_ratio_editing_in_peak") for row in rows]
    p_values = [as_float(row, "p_value_one_sided_fisher_greater") for row in rows]
    edited_fraction = [100.0 * as_float(row, "edited_peak_fraction") for row in rows]
    unedited_fraction = [
        100.0 * as_float(row, "unedited_peak_fraction") for row in rows
    ]
    colors = [
        color_for(odds_ratio, p_value, alpha)
        for odds_ratio, p_value in zip(odds_ratios, p_values)
    ]
    labels = [GROUP_LABELS.get(row["sample"], row["sample"]) for row in rows]

    fig = plt.figure(figsize=(8.27, 11.69))
    grid = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=(1.0, 1.25),
        hspace=0.18,
        left=0.10,
        right=0.98,
        top=0.93,
        bottom=0.08,
    )
    odds_axis = fig.add_subplot(grid[0])
    fraction_axis = fig.add_subplot(grid[1])

    style_axis(odds_axis)
    odds_axis.axhline(1.0, color="#555555", linestyle="--", linewidth=0.8)
    odds_axis.scatter(
        x_values,
        odds_ratios,
        s=175,
        c=colors,
        edgecolors="black",
        linewidths=0.8,
        zorder=3,
    )
    finite_or = [value for value in odds_ratios if math.isfinite(value)]
    max_or = max(finite_or) if finite_or else 2.5
    odds_top = max(2.8, math.ceil((max_or * 1.22) * 10.0) / 10.0)
    odds_axis.set_ylim(0, odds_top)
    odds_axis.set_ylabel("Odds ratio\n(edited / unedited)", fontsize=10)
    odds_axis.set_xticks(x_values)
    odds_axis.set_xticklabels([])
    odds_axis.tick_params(axis="x", length=0)
    odds_axis.tick_params(axis="y", labelsize=9)
    odds_axis.set_title(
        "PausePred editing enrichment\n"
        f"({SCOPE_LABELS.get(scope, scope)}, "
        f"window_nt={window}, Pause_score_cutoff={format_cutoff(cutoff)})",
        fontsize=12,
        pad=12,
    )
    for x_value, odds_ratio in zip(x_values, odds_ratios):
        if not math.isfinite(odds_ratio):
            continue
        odds_axis.text(
            x_value,
            odds_ratio + 0.055 * odds_top,
            f"OR={odds_ratio:.2f}",
            ha="center",
            va="bottom",
            fontsize=8.5,
        )
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor="#8B0000",
            markeredgecolor="black",
            markersize=9,
            label="OR > 1, P < 0.05",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor="#F4A3A3",
            markeredgecolor="black",
            markersize=9,
            label="OR > 1, P ≥ 0.05",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor="#2166AC",
            markeredgecolor="black",
            markersize=9,
            label="OR < 1, P < 0.05",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor="#9ECAE1",
            markeredgecolor="black",
            markersize=9,
            label="OR < 1, P ≥ 0.05",
        ),
    ]
    odds_axis.legend(
        handles=legend_handles,
        loc="upper left",
        ncol=2,
        frameon=False,
        fontsize=8,
        handletextpad=0.4,
        columnspacing=1.0,
        borderaxespad=0.1,
    )

    style_axis(fraction_axis)
    bar_width = 0.34
    fraction_axis.bar(
        [x - bar_width / 2 for x in x_values],
        edited_fraction,
        width=bar_width,
        color="#4C78A8",
        label="edited",
    )
    fraction_axis.bar(
        [x + bar_width / 2 for x in x_values],
        unedited_fraction,
        width=bar_width,
        color="#BDBDBD",
        label="unedit",
    )
    max_fraction = max(edited_fraction + unedited_fraction)
    fraction_top = max(70.0, math.ceil((max_fraction + 25.0) / 5.0) * 5.0)
    fraction_axis.set_ylim(0, fraction_top)
    fraction_axis.set_ylabel("Sites in PausePred peak (%)", fontsize=10)
    fraction_axis.set_xlabel("Species × sex group", fontsize=10)
    fraction_axis.set_xticks(x_values)
    fraction_axis.set_xticklabels(labels, fontsize=8.5)
    fraction_axis.tick_params(axis="y", labelsize=9)
    fraction_axis.legend(
        loc="upper left",
        frameon=False,
        fontsize=8,
        ncol=2,
        handlelength=1.8,
        columnspacing=1.0,
    )

    # Keep the count/P labels below the legend, as in the legacy figure.
    annotation_y = fraction_top - 8.0
    for x_value, row, p_value in zip(x_values, rows, p_values):
        edited_in = as_int(row, "edited_in_peak_in_hairpin")
        edited_total = as_int(row, "edited_total_in_hairpin")
        unedited_in = as_int(row, "unedited_in_peak_in_hairpin")
        unedited_total = as_int(row, "unedited_total_in_hairpin")
        annotation = (
            f"edited n=\n{edited_in:,}/\n{edited_total:,}\n"
            f"unedit n=\n{unedited_in:,}/\n{unedited_total:,}\n"
            f"P =\n{p_text(p_value)}"
        )
        fraction_axis.text(
            x_value,
            annotation_y,
            annotation,
            ha="center",
            va="top",
            fontsize=6.5,
            linespacing=1.0,
        )

    fig.savefig(output_pdf, format="pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    summary_path = Path(args.summary).resolve()
    output_dir = Path(args.out_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(summary_path)
    windows = [int(value) for value in args.windows.split(",") if value.strip()]
    scopes = ["all_edit_cds", "unedit_A_base"]
    manifest_rows: list[dict[str, str]] = []

    for scope in scopes:
        for window in windows:
            subset = ordered_subset(rows, scope, window, args.cutoff)
            scope_slug = scope
            cutoff_slug = format_cutoff(args.cutoff)
            stem = (
                f"odds_ratio_fraction_combined."
                f"{scope_slug}.window_{window}.cutoff_{cutoff_slug}"
            )
            output_pdf = output_dir / f"{stem}.pdf"
            output_tsv = output_dir / f"{stem}.tsv"
            write_plot_table(subset, output_tsv, args.alpha)
            make_plot(
                subset,
                output_pdf,
                scope,
                window,
                args.cutoff,
                args.alpha,
            )
            manifest_rows.append(
                {
                    "analysis_scope": scope,
                    "window_nt": str(window),
                    "cutoff": str(args.cutoff),
                    "n_samples": str(len(subset)),
                    "pdf": str(output_pdf),
                    "tsv": str(output_tsv),
                }
            )
            print(f"[write] {output_pdf}")

    manifest_path = output_dir / "combined_plot_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "analysis_scope",
                "window_nt",
                "cutoff",
                "n_samples",
                "pdf",
                "tsv",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"[write] {manifest_path}")


if __name__ == "__main__":
    main()
