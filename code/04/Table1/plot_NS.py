#!/usr/bin/env python3
"""Plot N/S composition from the counts and Fisher tests prepared by prepare_NS_plot.R."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--y-axis-max", type=int, choices=(16, 20), default=16)
    args = parser.parse_args()
    with args.data.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = [sp + "_" + sex for sp in ("mel", "sim", "sec", "ere") for sex in ("F", "M")]
    if [row["sample"] for row in rows] != expected:
        raise ValueError("Expected one female and one male row for each of the four species")
    plt.rcParams.update({"pdf.fonttype": 42, "font.size": 8})
    fig, ax = plt.subplots(figsize=(10, 3.5))
    heights = [[] for _ in rows]
    for offset, prefix, label, color in [(-0.2, "IgtA", "I>A", "#E55C8D"),
                                          (0.2, "Remaining", "Remaining", "#AAAAAA")]:
        values = [round(int(r[prefix + "_N"]) / int(r[prefix + "_S"]), 2) for r in rows]
        if max(values) > args.y_axis_max:
            raise ValueError("An N/S ratio exceeds --y-axis-max")
        ax.bar([i + offset for i in range(len(rows))], values, 0.38, color=color, label=label)
        for i, (row, value) in enumerate(zip(rows, values)):
            ax.text(i + offset, value + 0.15, row[prefix + "_N"] + "/" + row[prefix + "_S"],
                    ha="center", va="bottom", fontsize=6)
            heights[i].append(value)
    for i, row in enumerate(rows):
        ax.text(i, max(heights[i]) + 0.85, "P" + row["P_label"], ha="center", fontsize=7)
    ax.set_xticks(range(len(rows)), ["D" + r["sample"].replace("_", " ") for r in rows])
    ax.set_ylabel("N/S ratio")
    ax.set_ylim(0, args.y_axis_max)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(0, 1.17))
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
