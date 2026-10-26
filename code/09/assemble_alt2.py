from pathlib import Path
import os, csv, json, math, itertools, zipfile, xml.etree.ElementTree as ET
from collections import defaultdict

P = Path(__file__).resolve().parent
baseline_file = Path(
    os.environ.get(
        "MAMMALS_COVERAGE_BASELINE", str(P / "provenance/original_coverage.tsv")
    )
)
BASELINE = (
    {
        (r["species"], r["dataset"], r["library"]): r["recount_CDS_covered"]
        for r in csv.DictReader(open(baseline_file), delimiter="\t")
    }
    if baseline_file.is_file()
    else {}
)


def read(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def write(p, head, rows):
    with open(P / p, "w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(head)
        w.writerows(rows)


def gkey(r):
    return tuple(r[x] for x in ["chromosome", "strand", "base1", "base2", "base3"])


def txkey(r):
    return (r["transcript"], r["codon_start_tx"])


def rpass(r):
    return r["RNA_total"] >= 5 - 1e-9 and r["RNA_alt"] >= 2 - 1e-9


def bpass(r):
    return r["Asite_total"] >= 5 - 1e-9


# Standard genetic code in TCAG order.
code = dict(
    zip(
        ("".join(c) for c in itertools.product("TCAG", repeat=3)),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
)
headers = [
    "Dataset",
    "RNA-covered CDS codons (>0)",
    "RNA candidate editing codons (DP>=5, alt>=2)",
    "A-site-covered CDS codons (>0)",
    "A-site candidate editing codons (DP>=5)",
    "Joint candidate codons / scatter points",
    "RNA accession",
    "Ribo accession",
]
pointhead = [
    "species",
    "dataset",
    "chromosome",
    "strand",
    "base1",
    "base2",
    "base3",
    "transcript",
    "codon_start_tx",
    "ref_codon",
    "edit_mask",
    "editing_positions_1based",
    "RNA_total",
    "RNA_alt",
    "Asite_total",
    "Asite_alt",
    "Lm",
    "Lr",
    "consequence",
    "eligible_transcripts",
]
gene_by_tx = {}
for sp in ["human", "mouse"]:
    for r in read(P / (sp + ".selected_transcripts.tsv")):
        gene_by_tx[r["transcript"]] = r["gene_id"]
tablerows = []
points = []
allcandidate = []
audits = []
jointrows = []
for sp in ["human", "mouse"]:
    for ds in (
        ["brain1", "brain2", "fibroblasts", "liver", "testis"]
        if sp == "human"
        else ["brain", "liver", "testis"]
    ):
        ar = read(P / f"{sp}.{ds}.RNA.transcript_candidate_counts.tsv")
        br = {
            txkey(r): r
            for r in read(P / f"{sp}.{ds}.RPF.transcript_candidate_counts.tsv")
        }
        assert len(ar) == len(br)
        groups = defaultdict(list)
        for r in ar:
            b = br[txkey(r)]
            assert gkey(r) == gkey(b) and r["edit_mask"] == b["edit_mask"]
            v = dict(r)
            v.update(
                RNA_total=float(r["total"]),
                RNA_alt=float(r["alt"]),
                Asite_total=float(b["total"]),
                Asite_alt=float(b["alt"]),
            )
            assert (
                -1e-9 <= v["RNA_alt"] <= v["RNA_total"] + 1e-8
                and -1e-9 <= v["Asite_alt"] <= v["Asite_total"] + 1e-8
            )
            groups[gkey(v)].append(v)
        na = nb = nj = ncross = 0
        local = []
        for k, vv in sorted(groups.items()):
            pa = any(rpass(v) for v in vv)
            pb = any(bpass(v) for v in vv)
            na += pa
            nb += pb
            ncross += pa and pb
            eligible = [v for v in vv if rpass(v) and bpass(v)]
            assert len({gene_by_tx[v["transcript"]] for v in vv}) == len(
                vv
            ), "Only one fixed transcript per gene is allowed"
            rank = lambda v: (
                -min(v["RNA_total"], v["Asite_total"]),
                -(v["RNA_total"] + v["Asite_total"]),
                v["transcript"],
                int(v["codon_start_tx"]),
            )
            v = sorted(eligible or vv, key=rank)[0]
            allcandidate.append(
                [
                    sp,
                    ds,
                    *k,
                    int(pa),
                    int(pb),
                    int(bool(eligible)),
                    v["transcript"],
                    v["RNA_total"],
                    v["RNA_alt"],
                    v["Asite_total"],
                    v["Asite_alt"],
                ]
            )
            if not eligible:
                continue
            nj += 1
            mask = int(v["edit_mask"])
            ref = v["ref_codon"]
            offsets = [j for j in range(3) if mask & (1 << j)]
            classes = set()
            for n in range(1, 1 << len(offsets)):
                alt = list(ref)
                for j, o in enumerate(offsets):
                    if n & (1 << j):
                        alt[o] = "G"
                aa, ab = code.get(ref, "?"), code.get("".join(alt), "?")
                classes.add(
                    "Unknown"
                    if "?" in (aa, ab)
                    else "Stop" if "*" in (aa, ab) else "S" if aa == ab else "N"
                )
            consequence = next(iter(classes)) if len(classes) == 1 else "Mixed"
            r = [
                sp,
                ds,
                *k,
                v["transcript"],
                v["codon_start_tx"],
                ref,
                mask,
                ",".join(k[2 + j] for j in offsets),
                v["RNA_total"],
                v["RNA_alt"],
                v["Asite_total"],
                v["Asite_alt"],
                v["RNA_alt"] / v["RNA_total"],
                v["Asite_alt"] / v["Asite_total"],
                consequence,
                len(eligible),
            ]
            points.append(r)
            local.append(r)
        s = []
        for lib in ["RNA", "RPF"]:
            q = dict(
                csv.reader(open(P / f"{sp}.{ds}.{lib}.summary.tsv"), delimiter="\t")
            )
            if (sp, ds, lib) in BASELINE:
                assert int(q["covered_CDS_codons"]) <= int(BASELINE[(sp, ds, lib)])
            assert int(q["candidate_CDS_codons_reference"]) == len(groups)
            s.append(int(q["covered_CDS_codons"]))
            audits.append(
                [
                    sp,
                    ds,
                    lib,
                    int(q["input_rows"]),
                    s[-1],
                    BASELINE.get((sp, ds, lib), "not_provided"),
                ]
            )
        label = (
            {"human": "Human", "mouse": "Mouse"}[sp]
            + " "
            + {"brain1": "brain 1", "brain2": "brain 2"}.get(ds, ds)
        )
        accRNA = accRPF = "E-MTAB-7247"
        if sp == "human" and ds == "brain1":
            accRNA = "SRR1562544/SRR1562545/SRR1562546"
            accRPF = "SRR1562539/SRR1562540/SRR1562541"
        if sp == "human" and ds == "fibroblasts":
            accRNA = "SRR2064024"
            accRPF = "SRR2064017"
        tablerows.append([label, s[0], na, s[1], nb, nj, accRNA, accRPF])
        jointrows.append([label, s[0], nj, s[1], nj, accRNA, accRPF])
        assert nj <= na <= s[0] and nj <= nb <= s[1]
        write(f"{sp}.{ds}.scatter_points.tsv", pointhead, local)
        print(
            label,
            "RNA",
            na,
            "Asite",
            nb,
            "joint same tx",
            nj,
            "cross tx intersection",
            ncross,
            flush=True,
        )
        assert nj == len(local)
write("TableS12_alt2.tsv", headers, tablerows)
publication_headers = [headers[i] for i in [0, 1, 2, 3, 5, 6, 7]]
publication_headers[4] = "Candidate edited codons covered by both assays"
publication_rows = [[r[i] for i in [0, 1, 2, 3, 5, 6, 7]] for r in tablerows]
write("TableS12_publication.tsv", publication_headers, publication_rows)
write(
    "TableS12_joint_alt2.tsv",
    [
        headers[0],
        headers[1],
        "Joint candidate editing codons, RNA",
        headers[3],
        "Joint candidate editing codons, A-site",
        *headers[-2:],
    ],
    jointrows,
)
write("scatter_points.tsv", pointhead, points)
write(
    "all_candidate_codon_eligibility.tsv",
    [
        "species",
        "dataset",
        "chromosome",
        "strand",
        "base1",
        "base2",
        "base3",
        "RNA_eligible",
        "Asite_eligible",
        "joint_same_transcript_eligible",
        "representative_transcript",
        "RNA_total",
        "RNA_alt",
        "Asite_total",
        "Asite_alt",
    ],
    allcandidate,
)
pass
from xlsx_writer import write_workbook

sheets = [
    ("TableS12", publication_headers, publication_rows),
    ("independent_assay_counts", headers, tablerows),
    (
        "joint_for_scatter",
        [
            headers[0],
            headers[1],
            "Joint candidate RNA",
            headers[3],
            "Joint candidate A-site",
            *headers[-2:],
        ],
        jointrows,
    ),
    ("scatter_points", pointhead, points),
]
write_workbook(P / "TableS12_alt2.xlsx", sheets)
print("TOTAL SCATTER POINTS", len(points))
with open(P / "new_counts_summary.json", "w") as f:
    json.dump(
        {
            "scatter_points": len(points),
            "points_by_species": {
                sp: sum(r[0] == sp for r in points) for sp in ["human", "mouse"]
            },
        },
        f,
        indent=2,
    )
