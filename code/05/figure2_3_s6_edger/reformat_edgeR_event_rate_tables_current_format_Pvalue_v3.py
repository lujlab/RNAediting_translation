#!/usr/bin/env python3
"""Write clean edited/un-edited event-rate tables for edgeR v2.

The output is deliberately human-readable and follows the requested layout:

Library | Unedited codon: outside, inside, ratio | Edited codon: outside,
inside, ratio | FDR

Stalling and valley are written as separate tables because the event name is
part of the column meaning.  The original mixed plot-data path is replaced by
the same two tables as two clearly marked blocks; a machine-readable long
backup is retained as ``edgeR_event_rate_plot_data.long.tsv``.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

ORDER = [
    "dmel_Adar5G1_F",
    "dmel_Adar5G1_M",
    "dmel_ISO1_FM",
    "dmel_WBg_F",
    "dmel_WBg_M",
    "dmel_AdarG_F",
    "dmel_AdarG_M",
    "dmel_AdarS_F",
    "dmel_AdarS_M",
    "dsim_w501_F",
    "dsim_w501_M",
    "dsec_Dsec_F",
    "dsec_Dsec_M",
    "dere_Dere_F",
    "dere_Dere_M",
]

DISPLAY_NAME = {
    "dmel_Adar5G1_F": "5G1_F",
    "dmel_Adar5G1_M": "5G1_M",
    "dmel_ISO1_FM": "ISO1_F+M",
    "dmel_WBg_F": "WBG_F",
    "dmel_WBg_M": "WBG_M",
    "dmel_AdarG_F": "AG_F",
    "dmel_AdarG_M": "AG_M",
    "dmel_AdarS_F": "AS_F",
    "dmel_AdarS_M": "AS_M",
    "dsim_w501_F": "sim_F",
    "dsim_w501_M": "sim_M",
    "dsec_Dsec_F": "sec_F",
    "dsec_Dsec_M": "sec_M",
    "dere_Dere_F": "ere_F",
    "dere_Dere_M": "ere_M",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-summary", required=True)
    parser.add_argument("--hairpin-summary", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser.parse_args()


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


def as_int(row: dict[str, str], key: str) -> int:
    value = as_float(row, key)
    return int(round(value)) if math.isfinite(value) else 0


def read_rows(path: Path, scope: str) -> list[dict[str, object]]:
    """
    Read current edgeR summary formats.

    Supported:
    1) all summary:
       group_id, event, edited_total, edited_event, edited_rate,
       unedited_total, unedited_event, unedited_rate

    2) hairpin long summary:
       library, event,
       edited_total_region_tested, edited_inside_event,
       edited_outside_event, edited_event_rate,
       unedited_total_region_tested, unedited_inside_event,
       unedited_outside_event, unedited_event_rate
    """
    with path.open("r", encoding="utf-8", newline="") as handle:
        raw_rows = list(csv.DictReader(handle, delimiter="\t"))

    if not raw_rows:
        raise ValueError(f"No rows in {path}")

    rows = []

    # all summary format
    if "group_id" in raw_rows[0]:
        for raw in raw_rows:
            edited_inside = as_int(raw, "edited_event")
            edited_total = as_int(raw, "edited_total")
            unedited_inside = as_int(raw, "unedited_event")
            unedited_total = as_int(raw, "unedited_total")

            rows.append(
                {
                    "scope": scope,
                    "library": raw["group_id"],
                    "event": raw["event"],
                    "edited_total": edited_total,
                    "edited_inside": edited_inside,
                    "edited_outside": edited_total - edited_inside,
                    "edited_rate": as_float(raw, "edited_rate"),
                    "unedited_total": unedited_total,
                    "unedited_inside": unedited_inside,
                    "unedited_outside": unedited_total - unedited_inside,
                    "unedited_rate": as_float(raw, "unedited_rate"),
                    "p_value": as_float(raw, "fisher_p"),
                }
            )
        return rows

    # hairpin long format
    required = {
        "library",
        "event",
        "edited_total_region_tested",
        "edited_inside_event",
        "edited_outside_event",
        "edited_event_rate",
        "unedited_total_region_tested",
        "unedited_inside_event",
        "unedited_outside_event",
        "unedited_event_rate",
        "event_fisher_p_greater",
    }
    missing = required - set(raw_rows[0])
    if missing:
        raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")

    for raw in raw_rows:
        rows.append(
            {
                "scope": scope,
                "library": raw["library"],
                "event": raw["event"],
                "edited_total": as_int(raw, "edited_total_region_tested"),
                "edited_inside": as_int(raw, "edited_inside_event"),
                "edited_outside": as_int(raw, "edited_outside_event"),
                "edited_rate": as_float(raw, "edited_event_rate"),
                "unedited_total": as_int(raw, "unedited_total_region_tested"),
                "unedited_inside": as_int(raw, "unedited_inside_event"),
                "unedited_outside": as_int(raw, "unedited_outside_event"),
                "unedited_rate": as_float(raw, "unedited_event_rate"),
                "p_value": as_float(raw, "event_fisher_p_greater"),
            }
        )

    return rows


def bh_adjust(values: list[float]) -> list[float]:
    indexed = [
        (index, value) for index, value in enumerate(values) if math.isfinite(value)
    ]
    output = [math.nan] * len(values)
    indexed.sort(key=lambda item: item[1])
    n = len(indexed)
    running = 1.0
    for rank in range(n, 0, -1):
        index, value = indexed[rank - 1]
        running = min(running, value * n / rank)
        output[index] = min(running, 1.0)
    return output


def add_fdr(rows: list[dict[str, object]], alpha: float) -> None:
    groups: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in rows:
        groups.setdefault((str(row["scope"]), str(row["event"])), []).append(row)
    for group_rows in groups.values():
        fdr_values = bh_adjust([float(row["p_value"]) for row in group_rows])
        for row, fdr in zip(group_rows, fdr_values):
            row["fdr"] = fdr
            row["significant"] = math.isfinite(fdr) and fdr < alpha


def ordered_subset(
    rows: list[dict[str, object]], scope: str, event: str
) -> list[dict[str, object]]:
    lookup = {
        str(row["library"]): row
        for row in rows
        if row["scope"] == scope and row["event"] == event
    }
    return [lookup[library] for library in ORDER if library in lookup]


def ratio_text(value: float) -> str:
    return "NA" if not math.isfinite(value) else f"{100.0 * value:.2f}%"


def fdr_text(value: float) -> str:
    return "NA" if not math.isfinite(value) else f"{value:.6g}"


def event_header(event: str) -> list[str]:
    word = "peak" if event == "stalling" else "valley"
    return [
        "Library",
        f"Unedited codon Outside {word} number",
        f"Unedited codon Inside {word} number",
        f"Unedited codon Ratio inside {word}",
        f"Edited codon Outside {word} number",
        f"Edited codon Inside {word} number",
        f"Edited codon Ratio inside {word}",
        "P",
    ]


def event_values(row: dict[str, object], event: str) -> list[object]:
    return [
        DISPLAY_NAME[str(row["library"])],
        row["unedited_outside"],
        row["unedited_inside"],
        ratio_text(float(row["unedited_rate"])),
        row["edited_outside"],
        row["edited_inside"],
        ratio_text(float(row["edited_rate"])),
        fdr_text(float(row["p_value"])),
    ]


def write_event_table(
    rows: list[dict[str, object]], event: str, output_path: Path
) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(event_header(event))
        for row in rows:
            writer.writerow(event_values(row, event))


def write_combined_human_table(
    stalling_rows: list[dict[str, object]],
    valley_rows: list[dict[str, object]],
    output_path: Path,
) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["# all_edgeR_tested; P = one-sided Fisher test"])
        writer.writerow(event_header("stalling"))
        for row in stalling_rows:
            writer.writerow(event_values(row, "stalling"))
        writer.writerow([])
        writer.writerow(["# all_edgeR_tested; P = one-sided Fisher test"])
        writer.writerow(event_header("valley"))
        for row in valley_rows:
            writer.writerow(event_values(row, "valley"))


def write_long_machine_table(rows: list[dict[str, object]], output_path: Path) -> None:
    fields = [
        "scope",
        "event",
        "Library",
        "source_library",
        "Unedited_outside",
        "Unedited_inside",
        "Unedited_ratio_percent",
        "Edited_outside",
        "Edited_inside",
        "Edited_ratio_percent",
        "P",
        "Fisher_P",
        "P_significant",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "scope": row["scope"],
                    "event": row["event"],
                    "Library": DISPLAY_NAME[str(row["library"])],
                    "source_library": row["library"],
                    "Unedited_outside": row["unedited_outside"],
                    "Unedited_inside": row["unedited_inside"],
                    "Unedited_ratio_percent": 100 * float(row["unedited_rate"]),
                    "Edited_outside": row["edited_outside"],
                    "Edited_inside": row["edited_inside"],
                    "Edited_ratio_percent": 100 * float(row["edited_rate"]),
                    "P": row["p_value"],
                    "Fisher_P": row["p_value"],
                    "P_significant": float(row["p_value"]) < 0.05,
                }
            )


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    all_rows = read_rows(Path(args.all_summary).resolve(), "all_edgeR_tested")
    hairpin_rows = read_rows(
        Path(args.hairpin_summary).resolve(), "hairpin_edgeR_tested"
    )
    rows = all_rows + hairpin_rows
    add_fdr(rows, args.alpha)

    all_stalling = ordered_subset(all_rows, "all_edgeR_tested", "stalling")
    all_valley = ordered_subset(all_rows, "all_edgeR_tested", "valley")
    hairpin_stalling = ordered_subset(hairpin_rows, "hairpin_edgeR_tested", "stalling")

    stalling_path = out_dir / "edgeR_stalling_peak_event_rate_table.tsv"
    valley_path = out_dir / "edgeR_valley_event_rate_table.tsv"
    hairpin_path = out_dir / "edgeR_hairpin_stalling_peak_event_rate_table.tsv"
    combined_path = out_dir / "edgeR_event_rate_plot_data.tsv"
    long_path = out_dir / "edgeR_event_rate_plot_data.reformatted_long.tsv"

    write_event_table(all_stalling, "stalling", stalling_path)
    write_event_table(all_valley, "valley", valley_path)
    write_event_table(hairpin_stalling, "stalling", hairpin_path)
    write_combined_human_table(all_stalling, all_valley, combined_path)
    write_long_machine_table(rows, long_path)

    for path in [stalling_path, valley_path, hairpin_path, combined_path, long_path]:
        print(f"[write] {path}")


if __name__ == "__main__":
    main()
