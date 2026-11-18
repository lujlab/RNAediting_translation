#!/usr/bin/env python3
"""Plot RNA /2 edgeR edited/unedited-A-codon enrichment in the PausePred reference style."""

import argparse, csv, importlib.util, math, json, hashlib, datetime
from pathlib import Path
from scipy.stats import fisher_exact

CODE = Path(__file__).resolve().parent
SHARED = (
    Path(__file__).resolve().parents[2] / "figure_s6_s8_pausepred/plot_pausepred_reference_style.py"
)
spec = importlib.util.spec_from_file_location("pausepred_reference_style", SHARED)
style = importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)
ORDER = [
    "dmel_Adar5G1_M",
    "dmel_ISO1_FM",
    "dmel_WBg_F",
    "dmel_WBg_M",
    "dsim_w501_F",
    "dsim_w501_M",
    "dsec_Dsec_F",
    "dsec_Dsec_M",
    "dere_Dere_F",
    "dere_Dere_M",
]
LABELS = [
    "Adar5G1 M",
    "iso-1 F+M",
    "WT F",
    "WT M",
    "Dsim F",
    "Dsim M",
    "Dsec F",
    "Dsec M",
    "Dere F",
    "Dere M",
]
ALIASES = {
    "library": "group_id",
    "edited_total_region_tested": "edited_total",
    "edited_inside_event": "edited_event",
    "edited_event_rate": "edited_rate",
    "unedited_total_region_tested": "unedited_total",
    "unedited_inside_event": "unedited_event",
    "unedited_event_rate": "unedited_rate",
    "event_fisher_p_greater": "fisher_p",
}


def load(path):
    rows = style.read(path)
    return [{ALIASES.get(k, k): v for k, v in r.items()} for r in rows]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--edgeR-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--include-adarg-adars", action="store_true")
    args = ap.parse_args()
    root = args.edgeR_dir
    out = args.out_dir or root / "plots_pausepred_style"
    out.mkdir(parents=True, exist_ok=True)
    order = list(ORDER)
    labels = list(LABELS)
    if args.include_adarg_adars:
        order[4:4] = ["dmel_AdarG_F", "dmel_AdarG_M", "dmel_AdarS_F", "dmel_AdarS_M"]
        labels[4:4] = ["AdarG F", "AdarG M", "AdarS F", "AdarS M"]
    suffix = ".all14" if args.include_adarg_adars else ""
    records = []
    all_rows = []
    scopes = [
        (
            "all_codons",
            root
            / "edited_unedited_summary_codonA/edited_unedited_peaks_by_library.tsv",
        ),
        (
            "hairpin_inside",
            root
            / "edited_unedited_hairpin_summary/edited_unedited_hairpin_in_peaks_by_library.long.tsv",
        ),
    ]
    for scope, source in scopes:
        data = load(source)
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        for event in ["stalling"]:
            picked = [r for r in data if r["group_id"] in order and r["event"] == event]
            assert len(picked) == len(order)
            rows = []
            for r in picked:
                e, et, u, ut = [
                    int(r[k])
                    for k in [
                        "edited_event",
                        "edited_total",
                        "unedited_event",
                        "unedited_total",
                    ]
                ]
                assert 0 <= e <= et and 0 <= u <= ut and min(et, ut) > 0
                assert math.isclose(e / et, float(r["edited_rate"]), abs_tol=1e-12)
                assert math.isclose(u / ut, float(r["unedited_rate"]), abs_tol=1e-12)
                calculated = fisher_exact(
                    [[e, et - e], [u, ut - u]], alternative="greater"
                )[1]
                assert math.isclose(
                    calculated, float(r["fisher_p"]), rel_tol=1e-7, abs_tol=1e-12
                ), (scope, event, r["group_id"])
                row = {
                    "sample": r["group_id"],
                    "event": event,
                    "scope": scope,
                    "edited_event": e,
                    "edited_total": et,
                    "unedited_event": u,
                    "unedited_total": ut,
                    "p_value": float(r["fisher_p"]),
                }
                rows.append(row)
                all_rows.append(row)
            pdf = out / f"{scope}.{event}{suffix}.pdf"
            name = "peak"
            style.draw(
                rows,
                pdf,
                source,
                scope,
                order=order,
                labels=labels,
                ylabel=f"Codons in edgeR {name} (%)",
                min_ymax=30,
                edited_key="edited_event",
                unedited_key="unedited_event",
            )
            records.append(
                {
                    "scope": scope,
                    "event": event,
                    "pdf": str(pdf.resolve()),
                    "source": str(source.resolve()),
                    "source_sha256": source_hash,
                    "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
                    "groups": len(order),
                }
            )
        assert hashlib.sha256(source.read_bytes()).hexdigest() == source_hash
    with (out / f"plot_data{suffix}.tsv").open("w") as f:
        w = csv.DictWriter(
            f, fieldnames=list(all_rows[0]), delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(all_rows)
    report = {
        "completed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "groups": order,
        "figure_count": 2,
        "count_rate_and_one_sided_Fisher_checks": "passed",
        "source_statistics_unchanged": True,
        "background": "unedited reference codons containing A; edited codons retained",
        "ISO1_scope": "F/M jointly fitted edgeR model",
        "records": records,
    }
    pass
    print("[DONE] Two figures; counts, rates and Fisher P values checked.")


if __name__ == "__main__":
    main()
