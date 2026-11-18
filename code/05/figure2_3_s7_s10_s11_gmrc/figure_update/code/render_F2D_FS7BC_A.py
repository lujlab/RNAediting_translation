#!/usr/bin/env python3
"""Render F2D and FS7B/C from publication group means after exact genomic m6A exclusion."""

import csv
import gzip
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import render_figures as rf

FIG = Path(os.environ["FIGURE_ROOT"]).resolve()
GROUP_DIR = Path(
    os.environ.get("PUBLICATION_GROUP_MEAN_DIR", FIG.parent / "group_sex_mean")
).resolve()
M6A_REFERENCE_CONFIG = Path(os.environ["M6A_REFERENCE_CONFIG"])


def load_m6a_handles():
    with M6A_REFERENCE_CONFIG.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {}
    for row in rows:
        path = Path(row["m6a_file"])
        if not path.is_file():
            raise RuntimeError("Missing m6A reference: %s" % path)
        result[row["species"]] = set(path.read_text().split())
    return result


M6A_BY_SPECIES = load_m6a_handles()

F2D = [
    ("Dmel_Adar5G1__male", "FBtr0309056", 661, (560, 760), "Dmel (Adar5G1) M (jef)"),
    ("Dmel__male", "FBtr0309056", 661, (560, 760), "Dmel (ISO1 + WT) M (jef)"),
    ("Dsim__male", "FBtr0355324", 661, (560, 760), "Dsim M (jef)"),
]

FS7B = [
    ("Dmel__female", "FBtr0309056", 661, (560, 760), "Dmel (ISO1 + WT) F (jef)"),
    ("Dsec__male", "FBtr0203176", 664, (560, 760), "Dsec M (jef)"),
    ("Dere__male", "FBtr0415162", 659, (560, 760), "Dere M (jef)"),
    ("Dsim__female", "FBtr0355324", 661, (560, 760), "Dsim F (jef)"),
    ("Dsec__female", "FBtr0203176", 664, (560, 760), "Dsec F (jef)"),
    ("Dere__female", "FBtr0415162", 659, (560, 760), "Dere F (jef)"),
]

FS7C = [
    ("Dmel_Adar5G1__male", "FBtr0330286", 264, (164, 364), "Dmel (Adar5G1) M (slo)"),
    ("Dmel__male", "FBtr0330286", 264, (164, 364), "Dmel (ISO1 + WT) M (slo)"),
    ("Dsec__male", "FBtr0209603", 281, (181, 381), "Dsec M (slo)"),
    ("Dmel__female", "FBtr0330286", 264, (164, 364), "Dmel (ISO1 + WT) F (slo)"),
    ("Dsec__female", "FBtr0209603", 281, (181, 381), "Dsec F (slo)"),
]


def requests_by_group():
    result = {}
    for group, tx, _site, xlim, _title in F2D + FS7B + FS7C:
        result.setdefault(group, {})[tx] = xlim
    return result


def read_group(item):
    group, requests = item
    species = group.split("__", 1)[0].split("_", 1)[0]
    m6a_handles = M6A_BY_SPECIES[species]
    path = GROUP_DIR / (group + ".mean_gmRC.tsv.gz")
    profiles = {(group, tx): {} for tx in requests}
    qc = {
        (group, tx): {
            "group_id": group,
            "species": species,
            "transcript": tx,
            "n_rows_before_m6a": 0,
            "n_m6a_removed": 0,
            "n_rows_after_m6a": 0,
        }
        for tx in requests
    }
    conditions = []
    for tx, (lo, hi) in sorted(requests.items()):
        conditions.append('(tx=="%s" && $10>=%d && $10<=%d)' % (tx, lo, hi))
    program = (
        'NR>1 {tx=$2; gsub(/"/,"",tx)} ('
        + " || ".join(conditions)
        + ") {print tx,$10,$5,$6,$7,$16,$18,$19,$12,$13}"
    )
    unzip = subprocess.Popen(["gzip", "-cd", str(path)], stdout=subprocess.PIPE)
    awk = subprocess.Popen(
        ["awk", "-F", "\t", "-v", "OFS=\t", program],
        stdin=unzip.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    unzip.stdout.close()
    payload, error = awk.communicate()
    unzip_status = unzip.wait()
    if awk.returncode or unzip_status:
        raise RuntimeError(
            "profile extraction failed for %s: %s"
            % (group, error.decode("utf-8", "replace"))
        )
    for line in payload.decode("utf-8").splitlines():
        tx, pos, gmrc, nrep, expected, mrna, asite, any_positive, chrom, genomic_pos = (
            line.split("\t")
        )
        q = qc[(group, tx)]
        q["n_rows_before_m6a"] += 1
        genomic_pos = genomic_pos.strip('"')
        chrom = chrom.strip('"')
        handle = "%s_%s" % (chrom, genomic_pos)
        if handle in m6a_handles:
            q["n_m6a_removed"] += 1
            continue
        q["n_rows_after_m6a"] += 1
        profiles[(group, tx)][int(pos)] = {
            "FBtr": tx,
            "codonpos": pos,
            "gmRC": gmrc,
            "n_gmRC_reps_finite": nrep,
            "n_libraries_expected": expected,
            "mRNA_sum": mrna,
            "Asite_sum": asite,
            "Asite_any_positive": any_positive,
            "genomic_handle": handle,
        }
    return profiles, list(qc.values())


def read_profiles():
    profiles = {}
    qc_rows = []
    requests = requests_by_group()
    with ThreadPoolExecutor(max_workers=len(requests)) as pool:
        futures = [pool.submit(read_group, item) for item in requests.items()]
        for future in as_completed(futures):
            result_profiles, result_qc = future.result()
            profiles.update(result_profiles)
            qc_rows.extend(result_qc)
    for group, tx, site, _xlim, _title in F2D + FS7B + FS7C:
        if site not in profiles[(group, tx)]:
            raise RuntimeError("Missing focal codon: %s %s %s" % (group, tx, site))
    return profiles, qc_rows


def xticks(xlim):
    if xlim == (560, 760):
        return [600, 650, 700, 750]
    return [200, 250, 300, 350]


def draw_legend(c, width, top):
    rf.rect(c, width - 120, top - 2, 7, 7, rf.EDIT)
    rf.text(c, width - 109, top - 1, "Edited", 6.2)
    rf.rect(c, width - 69, top - 2, 7, 7, rf.CONTROL)
    rf.text(c, width - 58, top - 1, "Unedited", 6.2)


def draw_one(c, spec, profile, x, y, w, h, show_ylabel=False):
    group, tx, site, xlim, title = spec
    X, Y = rf.axis(
        c,
        x,
        y,
        w,
        h,
        xlim,
        (0, 10),
        xticks(xlim),
        [0, 2, 4, 6, 8, 10],
        ylabel="gmRC" if show_ylabel else None,
        size=5.5,
    )
    step = w / float(xlim[1] - xlim[0])
    bar_width = max(0.45, min(0.80, step * 0.86))
    selected = []
    for pos in range(xlim[0], xlim[1] + 1):
        row = profile.get(pos)
        if row is None:
            continue
        value = float(row["gmRC"])
        if value > 0:
            clipped = min(value, 10.0)
            rf.rect(
                c,
                X(pos) - bar_width / 2.0,
                Y(0),
                bar_width,
                Y(clipped) - Y(0),
                rf.EDIT if pos == site else rf.CONTROL,
            )
        selected.append(
            {
                "group_id": group,
                "transcript": tx,
                "title": title,
                "codonpos": pos,
                "gmRC": row["gmRC"],
                "is_focal_edited_codon": int(pos == site),
                "n_gmRC_reps_finite": row["n_gmRC_reps_finite"],
                "n_libraries_expected": row["n_libraries_expected"],
                "Asite_any_positive": row["Asite_any_positive"],
                "mRNA_sum": row["mRNA_sum"],
                "Asite_sum": row["Asite_sum"],
                "genomic_handle": row["genomic_handle"],
                "m6A_excluded": 0,
            }
        )
    focal = profile[site]
    focal_value = float(focal["gmRC"])
    anchor_y = Y(min(focal_value, 10.0))
    label_x = min(x + w - 4, X(site) + 33)
    label_y = min(y + h - 7, max(y + h * 0.62, anchor_y + 10))
    rf.line(
        c, X(site) + bar_width / 2, anchor_y, label_x - 3, label_y - 2, rf.EDIT, 0.55
    )
    rf.text(c, label_x, label_y, "%.2f" % focal_value, 6.2, "right", color=rf.EDIT)
    rf.text(c, x + w / 2.0, y + h + 9, title, 6.7, "center")
    return selected


def render_grid(name, panel_letter, specs, profiles, ncol=3):
    width = 510
    nrow = int(math.ceil(len(specs) / float(ncol)))
    height = 148 if nrow == 1 else 260
    c = rf.new(name, width, height)
    rf.text(c, 6, height - 13, panel_letter, 11, font="SansBold")
    draw_legend(c, width, height - 12)
    pw = 143
    ph = 71
    x0 = 42
    top_y = height - 104
    rows = []
    for index, spec in enumerate(specs):
        col = index % ncol
        row = index // ncol
        x = x0 + col * 157
        y = top_y - row * 112
        profile = profiles[(spec[0], spec[1])]
        rows.extend(draw_one(c, spec, profile, x, y, pw, ph, show_ylabel=(col == 0)))
    rf.text(c, width / 2.0, 12, "Transcript codon position", 7.5, "center")
    c.save()
    rf.write(name + "_profile_data.tsv", rows)
    return rows




def main():
    profiles, qc_rows = read_profiles()
    rf.write(
        "F2D_FS7BC_m6A_filter_QC.tsv",
        sorted(qc_rows, key=lambda x: (x["group_id"], x["transcript"])),
    )
    render_grid("F2D", "D", F2D, profiles)
    render_grid("FS7B", "B", FS7B, profiles)
    render_grid("FS7C", "C", FS7C, profiles)
    print("Rendered F2D, FS7B and FS7C from", GROUP_DIR)


if __name__ == "__main__":
    main()
