#!/usr/bin/env python3
"""Draw one S10 panel per species from the model's statistics and density polygons."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=root / "results_3N_sex_model")
    parser.add_argument("--outdir", type=Path, default=root / "figure_3N_sex_model")
    args = parser.parse_args()
    stats = read(args.data_dir / "FS10_plot_stats.tsv")
    density = read(args.data_dir / "FS10_density_polygons.tsv")
    models = read(args.data_dir / "RSCU_site_level_AA_adjusted_models.tsv")
    args.outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"pdf.fonttype": 42, "font.size": 9})
    for species in ("Dmel", "Dsim", "Dsec", "Dere"):
        rows = {r["preference"]: r for r in stats if r["species"] == species}
        test = next(r for r in models if r["species"] == species
                    and r["model"] == "AA_and_sex_adjusted_primary"
                    and r["SE_type"] == "cluster_robust_codon_CR1")
        fig, ax = plt.subplots(figsize=(3.3, 3.6))
        fig.subplots_adjust(top=0.70, bottom=0.17, left=0.22, right=0.97)
        for position, preference, color in [(1, "Preferred", "#9B9B9B"),
                                             (2, "Unpreferred", "#E88E94")]:
            row = rows[preference]
            points = sorted((r for r in density if r["species"] == species
                             and r["preference"] == preference), key=lambda r: int(r["polygon_order"]))
            ax.add_patch(Polygon([(float(r["x"]), float(r["y"])) for r in points],
                                 closed=True, facecolor=color, edgecolor="none"))
            q1, q3, median = [float(row[k]) for k in ("q1_gmRC", "q3_gmRC", "median_gmRC")]
            mean, sd = float(row["mean_gmRC"]), float(row["sd_gmRC"])
            ax.add_patch(Rectangle((position - 0.1, q1), 0.2, q3 - q1,
                                   facecolor="white", edgecolor="black", linewidth=0.6))
            ax.plot([position - 0.1, position + 0.1], [median, median], color="black", linewidth=0.6)
            ax.errorbar(position, mean, yerr=[[mean - max(0, mean - sd)], [sd]],
                        fmt="none", ecolor="black", capsize=3, linewidth=0.7)
            ax.text(position, 1.12, "n = {:,}".format(int(row["n_positions"])),
                    transform=ax.get_xaxis_transform(), ha="center", fontsize=8)
            ax.text(position, 1.03, "{:.3f} ± {:.3f}".format(mean, sd),
                    transform=ax.get_xaxis_transform(), ha="center", fontsize=8)
        ax.set_title(species + "  (P = {:.3g})".format(float(test["P_one_sided_unpreferred_greater"])), pad=43)
        ax.set_xlim(0.4, 2.6)
        ax.set_ylim(0, max(float(r["y_upper"]) for r in rows.values()))
        ax.set_xticks([1, 2], ["Preferred", "Unpreferred"])
        ax.set_ylabel("gmRC")
        ax.spines[["top", "right"]].set_visible(False)
        output = args.outdir / ("FigureS10_" + species + ".pdf")
        fig.savefig(output)
        plt.close(fig)
        print(output)


if __name__ == "__main__":
    main()
