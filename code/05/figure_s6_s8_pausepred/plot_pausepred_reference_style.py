#!/usr/bin/env python3
"""Redraw two PausePred A-base bar charts in the supplied Figure 10 style."""

import argparse, csv, json, math, hashlib
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.patches import Rectangle
from fontTools.ttLib import TTFont

ORDER = [
    "Dmel_Adar5G1_M",
    "Dmel_ISO1_head_female_ribo",
    "Dmel_ISO1_head_male_ribo",
    "WBg_F",
    "WBg_M",
    "Dsim_F",
    "Dsim_M",
    "Dsec_F",
    "Dsec_M",
    "Dere_F",
    "Dere_M",
]
LABELS = [
    "Adar5G1 M",
    "iso-1 F",
    "iso-1 M",
    "WT F",
    "WT M",
    "Dsim F",
    "Dsim M",
    "Dsec F",
    "Dsec M",
    "Dere F",
    "Dere M",
]
EDITED = (0.8, 0.078, 0.498)
UNEDIT = (0.29, 0.631, 0.710)
INK = (0.102, 0.102, 0.102)
PAGE_W = 595.276
PAGE_H = 841.890
FONTS = Path(__file__).resolve().parent / "reference_style_fonts"
REG = FontProperties(fname=str(FONTS / "Arial-Regular.ttf"))
BOLD = FontProperties(fname=str(FONTS / "Arial-Bold.ttf"))
ITALIC = FontProperties(fname=str(FONTS / "Arial-Italic.ttf"))
FONT = TTFont(str(FONTS / "Arial-Regular.ttf"))
CMAP = FONT.getBestCmap()


def regular_font(text):
    # The reference PDF embeds only the glyphs used in that figure. Use a metric-
    # compatible font for a tick label needing a new digit, such as 60.
    missing = [
        c
        for c in text
        if not c.isspace()
        and (
            ord(c) not in CMAP
            or not getattr(FONT["glyf"][CMAP[ord(c)]], "numberOfContours", 0)
        )
    ]
    return FontProperties(fname=findfont("Liberation Sans")) if missing else REG


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def marks(p):
    if not math.isfinite(p):
        return "NA"
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def draw(
    rows,
    pdf,
    source,
    scope,
    order=None,
    labels=None,
    ylabel="Sites in PausePred peak (%)",
    min_ymax=50,
    edited_key="edited_in_peak",
    unedited_key="unedited_in_peak",
):
    order = ORDER if order is None else order
    labels = LABELS if labels is None else labels
    n = len(order)
    assert len(labels) == n and n > 0
    label_size = 9 if n <= 11 else 7.5
    by_sample = {r["sample"]: r for r in rows}
    assert len(rows) == n and set(by_sample) == set(order)
    rows = [by_sample[s] for s in order]
    edited = [100 * r[edited_key] / r["edited_total"] for r in rows]
    unedit = [100 * r[unedited_key] / r["unedited_total"] for r in rows]
    top = max(min_ymax, 10 * math.ceil((max(edited + unedit) + 3) / 10))
    plt.rcParams.update(
        {"pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": False}
    )
    fig = plt.figure(figsize=(PAGE_W / 72, PAGE_H / 72))
    ax = fig.add_axes(
        [54.4877 / PAGE_W, 656.8821 / PAGE_H, 524 / PAGE_W, 151.5658 / PAGE_H]
    )
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, top)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.65)
    ax.spines["left"].set_linewidth(0.593)
    for spine in ("bottom", "left"):
        ax.spines[spine].set_color(INK)
    ticks = list(range(0, top + 1, 10))
    ax.set_yticks(ticks)
    ax.tick_params(axis="y", labelsize=9, length=2.4, width=0.45, pad=2.6, colors=INK)
    for label in ax.get_yticklabels():
        label.set_fontproperties(regular_font(label.get_text()))
        label.set_fontsize(9)
    # Reference has internal horizontal grid lines, with no line at the top.
    for value in ticks[1:-1]:
        ax.axhline(value, color=(0.898,) * 3, linewidth=0.45, zorder=0)
    ax.set_xticks([])
    width = 13.2 / (524 / n)
    offset = 7.3 / (524 / n)
    ax.bar(
        [i - offset for i in range(n)],
        edited,
        width=width,
        color=EDITED,
        edgecolor="none",
        zorder=2,
    )
    ax.bar(
        [i + offset for i in range(n)],
        unedit,
        width=width,
        color=UNEDIT,
        edgecolor="none",
        zorder=2,
    )
    plot_rows = []
    for i, r in enumerate(rows):
        mark = marks(r["p_value"])
        ax.annotate(
            mark,
            (i, max(edited[i], unedit[i])),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="baseline",
            fontproperties=BOLD,
            fontsize=9,
            color=EDITED,
            annotation_clip=False,
        )
        center = 54.4877 + (i + 0.5) * (524 / n)
        if i == 0:
            shift = center - (54.4877 + 0.5 * (524 / 11))
            scale = label_size / 9
            anchor = 54.4877 + 0.5 * (524 / 11)
            fig.text(
                (anchor + shift + (61.8488 - anchor) * scale) / PAGE_W,
                638.5918 / PAGE_H,
                "Adar",
                fontproperties=ITALIC,
                fontsize=label_size,
                color=INK,
                va="baseline",
            )
            fig.text(
                (anchor + shift + (80.8574 - anchor) * scale) / PAGE_W,
                (638.5918 + (641.5889 - 638.5918) * scale) / PAGE_H,
                "5G1",
                fontproperties=ITALIC,
                fontsize=5.247 * scale,
                color=INK,
                va="baseline",
            )
            fig.text(
                (anchor + shift + (92.2329 - anchor) * scale) / PAGE_W,
                638.5918 / PAGE_H,
                "M",
                fontproperties=REG,
                fontsize=label_size,
                color=INK,
                va="baseline",
            )
        else:
            fig.text(
                center / PAGE_W,
                638.5918 / PAGE_H,
                labels[i],
                fontproperties=regular_font(labels[i]),
                fontsize=label_size,
                color=INK,
                ha="center",
                va="baseline",
            )
        plot_rows.append(
            dict(
                r,
                display_label=labels[i],
                edited_percent=edited[i],
                unedited_percent=unedit[i],
                significance=mark,
            )
        )
    fig.text(
        24.9229 / PAGE_W,
        667.7075 / PAGE_H,
        ylabel,
        fontproperties=regular_font(ylabel),
        fontsize=10,
        color=INK,
        rotation=90,
        ha="left",
        va="baseline",
    )
    for x, color, text in [(261.986, EDITED, "Edited"), (318.986, UNEDIT, "Unedit")]:
        fig.patches.append(
            Rectangle(
                (x / PAGE_W, 817.778 / PAGE_H),
                7 / PAGE_W,
                6 / PAGE_H,
                transform=fig.transFigure,
                facecolor=color,
                edgecolor="none",
            )
        )
        fig.text(
            (x + 11) / PAGE_W,
            817.1782 / PAGE_H,
            text,
            fontproperties=REG,
            fontsize=9,
            color=INK,
            va="baseline",
        )
    pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf, format="pdf")
    fig.savefig(pdf.with_suffix(".png"), dpi=180)
    # A cropped PNG preview is convenient without changing the A4 PDF canvas.
    from matplotlib.transforms import Bbox

    fig.savefig(
        pdf.with_name(pdf.stem + ".preview.png"),
        dpi=180,
        bbox_inches=Bbox.from_extents(0, 625 / 72, PAGE_W / 72, PAGE_H / 72),
    )
    plt.close(fig)
    with pdf.with_suffix(".reference_style.tsv").open("w") as f:
        w = csv.DictWriter(
            f, fieldnames=list(plot_rows[0]), delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(plot_rows)
    manifest = {
        "scope": scope,
        "source": str(source.resolve()),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "output": str(pdf.resolve()),
        "samples": order,
        "ymax": top,
        "ylabel": ylabel,
        "p_value": "original raw one-sided Fisher greater P",
        "significance_thresholds": {
            "***": "p < 0.001",
            "**": "p < 0.01",
            "*": "p < 0.05",
            "ns": "p >= 0.05",
        },
        "edited_rgb": EDITED,
        "unedited_rgb": UNEDIT,
        "canvas_points": [PAGE_W, PAGE_H],
        "axes_points": [54.4877, 656.8821, 524, 151.5658],
        "fonts": "Arial subsets extracted from reference PDF; unsupported tick digits use Liberation Sans",
    }
    pdf.with_suffix(".reference_style.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print("[PLOT]", pdf, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pausepred-dir", required=True, type=Path)
    ap.add_argument("--run-tag", default="seq90_up90_down90_counterRule_newalledit")
    ap.add_argument("--cutoff", default=5, type=float)
    ap.add_argument("--window", default=2, type=int)
    args = ap.parse_args()
    out = args.pausepred_dir
    tag = args.run_tag
    src = (
        out
        / f"enrichment_{tag}_unedit_A_base/plots/all_samples.pausepred_enrichment.merged.tsv"
    )
    selected = [
        r
        for r in read(src)
        if float(r["Pause_score_cutoff"]) == args.cutoff
        and int(r["window_nt"]) == args.window
    ]
    rows = [
        {
            "sample": r["sample"],
            "edited_in_peak": int(r["editing_in_peak"]),
            "edited_total": int(r["editing_total_used_in_universe"]),
            "unedited_in_peak": int(r["nonediting_in_peak"]),
            "unedited_total": int(r["nonediting_total"]),
            "p_value": float(r["p_value"]),
        }
        for r in selected
    ]
    draw(
        rows,
        src.parent
        / f"odds_ratio_fraction_combined.window_{args.window}.cutoff_{args.cutoff:g}.pdf",
        src,
        "unedit_A_base",
    )
    src = (
        out
        / f"hairpin_edit_unedit_enrichment_{tag}/hairpin_pausepred_edit_unedit_enrichment_summary.tsv"
    )
    selected = [
        r
        for r in read(src)
        if r["analysis_scope"] == "unedit_A_base"
        and float(r["Pause_score_cutoff"]) == args.cutoff
        and int(r["window_nt"]) == args.window
    ]
    rows = [
        {
            "sample": r["sample"],
            "edited_in_peak": int(r["edited_in_peak_in_hairpin"]),
            "edited_total": int(r["edited_total_in_hairpin"]),
            "unedited_in_peak": int(r["unedited_in_peak_in_hairpin"]),
            "unedited_total": int(r["unedited_total_in_hairpin"]),
            "p_value": float(r["p_value_one_sided_fisher_greater"]),
        }
        for r in selected
    ]
    draw(
        rows,
        src.parent
        / "plots"
        / f"odds_ratio_fraction_combined.unedit_A_base.window_{args.window}.cutoff_{args.cutoff:g}.pdf",
        src,
        "hairpin_unedit_A_base",
    )


if __name__ == "__main__":
    main()
