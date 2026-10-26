"""Four-panel human/mouse figure matching the supplied reference PDF; fixed per-gene longest CDS transcripts; RNA alt>=2 cohort."""

from pathlib import Path
import csv, math, itertools, warnings, json
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
from scipy.stats import wilcoxon, gaussian_kde

P = Path(__file__).resolve().parent
rows = [
    r
    for r in csv.DictReader(open(str(P / "scatter_points.tsv")), delimiter="\t")
    if r["species"] in ("human", "mouse")
]
# All allowed edited codon states must be synonymous to call a multi-position codon S.
code = dict(
    zip(
        ("".join(c) for c in itertools.product("TCAG", repeat=3)),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
)
for r in rows:
    ref = r["ref_codon"]
    off = [j for j in range(3) if int(r["edit_mask"]) & (1 << j)]
    effects = set()
    for n in range(1, 1 << len(off)):
        a = list(ref)
        for j, o in enumerate(off):
            if n & (1 << j):
                a[o] = "G"
        effects.add(code.get("".join(a), "?"))
    assert "?" not in effects and ref in code
    r["SN_class"] = "S" if effects == {code[ref]} else "N"
    r["SN_rule"] = (
        "all allowed edited codons retain reference amino acid: S; otherwise N; stop changes included"
    )
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "mathtext.default": "it",
    }
)
spec = [
    ("human", "Humans", ["brain1", "brain2", "liver", "testis", "fibroblasts"]),
    ("mouse", "Mice", ["brain", "liver", "testis"]),
]
labels = {
    "brain1": "Brain 1",
    "brain2": "Brain 2",
    "brain": "Brain",
    "liver": "Liver",
    "testis": "Testis",
    "fibroblasts": "Fibroblasts",
}
colors = {
    "brain1": "#1f77b4",
    "brain2": "#9467bd",
    "brain": "#1f77b4",
    "liver": "#ff7f0e",
    "testis": "#2ca02c",
    "fibroblasts": "#d62728",
}
sncolors = {"S": "#56BFE3", "N": "#F17A82"}
stats = {}
statrows = []
sncounts = []
for sp, title, datasets in spec:
    for ds in ["ALL"] + datasets:
        rr = [
            r
            for r in rows
            if r["species"] == sp and (ds == "ALL" or r["dataset"] == ds)
        ]
        x = np.array([float(r["Lm"]) for r in rr])
        y = np.array([float(r["Lr"]) for r in rr])
        d = np.round(y - x, 12)
        nz = d[d != 0]
        if len(d) == 0:
            w, p = float("nan"), float("nan")
            method = "not tested: empty"
        elif len(nz) == 0:
            w, p = 0.0, 1.0
            method = "all paired differences zero"
        else:
            exact = len(d) <= 50 and len(nz) == len(d) and len(set(abs(nz))) == len(nz)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                z = wilcoxon(
                    d,
                    zero_method="wilcox",
                    correction=not exact,
                    alternative="two-sided",
                    mode="exact" if exact else "approx",
                )
            w, p = float(z.statistic), float(z.pvalue)
            method = (
                "exact" if exact else "normal approximation with continuity correction"
            )
        stats[(sp, ds)] = (len(rr), p)
        statrows.append(
            [
                sp,
                ds,
                len(rr),
                len(nz),
                w,
                p,
                method,
                int(sum(y == 0)),
                float(np.median(x)) if len(x) else "",
                float(np.median(y)) if len(y) else "",
            ]
        )
    for ds in datasets:
        for cl in ["S", "N"]:
            sncounts.append(
                [
                    sp,
                    ds,
                    cl,
                    sum(
                        r["species"] == sp
                        and r["dataset"] == ds
                        and r["SN_class"] == cl
                        for r in rows
                    ),
                ]
            )


def pstr(p, mathmode=False):
    if not math.isfinite(p):
        return "NA"
    if p >= 0.001:
        return "{:.3g}".format(p) if p < 1 else "1.00"
    e = int(math.floor(math.log10(p))) if p else -300
    m = p / 10**e
    return (
        "{:.2f} \\times 10^{{{}}}".format(m, e)
        if mathmode
        else "{:.2f}e{}".format(m, e)
    )


def scatter(ax, sp, datasets):
    for ds in datasets:
        rr = [r for r in rows if r["species"] == sp and r["dataset"] == ds]
        ax.scatter(
            [float(r["Lm"]) for r in rr],
            [float(r["Lr"]) for r in rr],
            s=12,
            c=colors[ds],
            alpha=0.72,
            edgecolors="none",
            clip_on=False,
            zorder=3,
        )
    ax.plot([0, 1], [0, 1], color="black", lw=0.8, dashes=(7, 5), zorder=1)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.tick_params(length=2, pad=1.5, labelsize=6.5)
    ax.set_xlabel(r"$L_m$", fontsize=10, labelpad=1)
    ax.set_ylabel(r"$L_r$", fontsize=10, rotation=0, labelpad=10)
    ax.set_aspect("equal", adjustable="box")


def table(ax, sp, datasets):
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.add_patch(
        Rectangle(
            (0.01, 0.02),
            0.98,
            0.95,
            fill=False,
            linestyle=(0, (4, 4)),
            lw=0.45,
            edgecolor="#777777",
        )
    )
    ax.text(0.55, 0.88, r"$n$", fontsize=7, fontweight="bold", ha="center")
    ax.text(0.73, 0.88, r"$P$", fontsize=7, fontweight="bold", ha="center")
    ys = np.linspace(0.72, 0.12, len(datasets))
    for ds, y in zip(datasets, ys):
        n, p = stats[(sp, ds)]
        ax.plot(0.055, y, "o", color=colors[ds], ms=3.2)
        ax.text(0.105, y, labels[ds], fontsize=6.4, va="center")
        ax.text(0.55, y, str(n), fontsize=6.4, ha="center", va="center")
        ax.text(0.65, y, r"$" + pstr(p, True) + "$", fontsize=6.3, va="center")


def violin(ax, sp, datasets):
    for i, ds in enumerate(datasets):
        for cl, delta in [("S", -0.14), ("N", 0.14)]:
            vals = np.array(
                [
                    float(r["Lm"])
                    for r in rows
                    if r["species"] == sp and r["dataset"] == ds and r["SN_class"] == cl
                ]
            )
            y = i + delta
            if len(vals) > 1 and np.ptp(vals) > 1e-12:
                grid = np.linspace(float(vals.min()), float(vals.max()), 180)
                try:
                    den = gaussian_kde(vals)(grid)
                    half = den / den.max() * 0.115
                    ax.fill_between(
                        grid,
                        y - half,
                        y + half,
                        facecolor=sncolors[cl],
                        alpha=0.88,
                        linewidth=0,
                    )
                except np.linalg.LinAlgError:
                    pass
            if len(vals):
                q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
                iq = q3 - q1
                low = min(vals[vals >= q1 - 1.5 * iq])
                high = max(vals[vals <= q3 + 1.5 * iq])
                ax.plot([low, high], [y, y], color="#444444", lw=0.55)
                ax.add_patch(
                    Rectangle(
                        (q1, y - 0.037),
                        max(q3 - q1, 0.001),
                        0.074,
                        facecolor="white",
                        edgecolor="#444444",
                        lw=0.55,
                        zorder=3,
                    )
                )
                ax.plot(
                    [med, med],
                    [y - 0.037, y + 0.037],
                    color="#444444",
                    lw=0.55,
                    zorder=4,
                )
                if len(vals) == 1:
                    ax.plot(vals, [y], "|", color=sncolors[cl], ms=5)
            ax.text(
                1.01,
                y,
                r"$n = {}$".format(len(vals)),
                fontsize=6.7,
                va="center",
                ha="left",
                clip_on=False,
            )
    ax.set_ylim(len(datasets) - 0.45, -0.45)
    ax.set_xlim(-0.01, 1.12)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks(range(len(datasets)))
    ax.set_yticklabels(
        ["Fibro-\nblasts" if d == "fibroblasts" else labels[d] for d in datasets],
        fontsize=7,
    )
    ax.tick_params(axis="x", length=2, pad=1.5, labelsize=6.5)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.spines["bottom"].set_bounds(0, 1)
    ax.set_xlabel(r"$L_m$", fontsize=10, labelpad=1)
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=sncolors[k], edgecolor="none", label=k)
        for k in ["S", "N"]
    ]
    leg = ax.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(1.03, 1.22),
        ncol=2,
        fontsize=7,
        handlelength=0.9,
        handletextpad=0.25,
        columnspacing=0.4,
        borderpad=0.2,
        frameon=True,
        fancybox=False,
    )
    leg.get_frame().set_linewidth(0.4)
    leg.get_frame().set_linestyle((0, (4, 4)))


fig = plt.figure(figsize=(467.717 / 72, 410 / 72))
for i, (sp, title, datasets) in enumerate(spec):
    y = [0.58, 0.12][i]
    ax = fig.add_axes([0.119, y, 0.249, 0.284])
    scatter(ax, sp, datasets)
    tax = fig.add_axes([0.393, y + 0.02, 0.262, 0.112 if len(datasets) == 5 else 0.087])
    table(tax, sp, datasets)
    vx = fig.add_axes([0.713, y, 0.244, 0.284])
    violin(vx, sp, datasets)
    fig.text(0.247, y + 0.345, title, fontsize=8.8, ha="center")
    n, p = stats[(sp, "ALL")]
    fig.text(
        0.247,
        y + 0.315,
        r"$n = " + str(n) + r",\ P = " + pstr(p, True) + "$",
        fontsize=8.5,
        ha="center",
    )
    fig.text(0.84, y + 0.345, title, fontsize=8.8, ha="center")
    fig.text(0.045, y + 0.35, "AC"[i], fontsize=11, fontweight="bold")
    fig.text(0.66, y + 0.35, "BD"[i], fontsize=11, fontweight="bold")
fig.text(0.99, 0.009, "Longest CDS per gene", fontsize=6, ha="right", color="#555555")
for ext in ["pdf", "png", "svg"]:
    fig.savefig(str(P / ("FS28_alt2_scatter." + ext)), dpi=300)
plt.close(fig)
with open(str(P / "scatter_statistics.tsv"), "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(
        [
            "species",
            "dataset",
            "n_codons",
            "n_nonzero_differences",
            "Wilcoxon_W",
            "Wilcoxon_p",
            "method",
            "n_Asite_alt_zero",
            "median_Lm",
            "median_Lr",
        ]
    )
    w.writerows(statrows)
with open(str(P / "SN_violin_counts.tsv"), "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["species", "dataset", "SN_class", "n"])
    w.writerows(sncounts)
with open(str(P / "figure_points_with_SN.tsv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
assert sum(r[3] for r in sncounts) == len(rows)
print("Four-panel human/mouse figure exported. Total points:", len(rows))
for sp, title, ds in spec:
    print(sp, stats[(sp, "ALL")])
