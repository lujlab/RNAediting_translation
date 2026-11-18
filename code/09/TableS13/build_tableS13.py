from pathlib import Path
import csv, json, zipfile
from collections import Counter
import xlsxwriter

P = Path(__file__).resolve().parent
BASE = P.parent
fig = [
    r
    for r in csv.DictReader(
        open(str(BASE / "figure_points_with_SN.tsv")), delimiter="\t"
    )
    if r["species"] in ("human", "mouse")
]
meta = {}
for sp in ["human", "mouse"]:
    for r in csv.DictReader(
        open(str(BASE / (sp + ".selected_transcripts.tsv"))), delimiter="\t"
    ):
        meta[(sp, r["transcript"])] = r
spec = {"human": ("Homo sapiens", "GRCh38"), "mouse": ("Mus musculus", "GRCm38/mm10")}
rows = []
for sp, ds in sorted({(r["species"], r["dataset"]) for r in fig}):
    read = lambda lib: {
        (r["transcript"], r["codon_start_tx"]): r
        for r in csv.DictReader(
            open(
                str(
                    BASE
                    / ("{}.{}.{}.transcript_candidate_counts.tsv".format(sp, ds, lib))
                )
            ),
            delimiter="\t",
        )
    }
    aa, bb = read("RNA"), read("RPF")
    for r in [r for r in fig if r["species"] == sp and r["dataset"] == ds]:
        k = (r["transcript"], r["codon_start_tx"])
        a, b = aa[k], bb[k]
        m = meta[(sp, r["transcript"])]
        am, im, ar, ir = (
            float(a["allA"]),
            float(a["alt"]),
            float(b["allA"]),
            float(b["alt"]),
        )
        dm, dr = float(a["total"]), float(b["total"])
        assert (
            abs(im - float(r["RNA_alt"])) < 1e-7
            and abs(ir - float(r["Asite_alt"])) < 1e-7
        )
        om, orr = dm - am - im, dr - ar - ir
        om = 0.0 if abs(om) < 1e-8 else om
        orr = 0.0 if abs(orr) < 1e-8 else orr
        assert min(om, orr) >= 0
        coord = ",".join(r["base" + str(j)] for j in [1, 2, 3])
        rows.append(
            dict(
                Species=spec[sp][0],
                Dataset=ds,
                Genome_assembly=spec[sp][1],
                Chromosome="chr" + r["chromosome"],
                Strand=r["strand"],
                Editing_positions_1based=r["editing_positions_1based"],
                Codon_genomic_positions_translation_order=coord,
                Gene_ID=m["gene_id"],
                Gene_symbol=m["gene_name"],
                Transcript_ID=r["transcript"],
                Reference_codon=r["ref_codon"],
                SN_class=r["SN_class"],
                Am=am,
                Im=im,
                Ar=ar,
                Ir=ir,
                RNA_other_count=om,
                Asite_other_count=orr,
                RNA_total_count=dm,
                Asite_total_count=dr,
                Lm_A_I_only=im / (am + im),
                Lr_A_I_only=ir / (ar + ir),
                Lm_plotted=float(r["Lm"]),
                Lr_plotted=float(r["Lr"]),
                FS28_point_ID=":".join([sp, ds, r["chromosome"], r["strand"], coord]),
                CDS_length=int(m["CDS_length"]),
                Transcript_length=int(m["transcript_length"]),
            )
        )
assert len(rows) == len(fig) and len({r["FS28_point_ID"] for r in rows}) == len(rows)
with open(str(P / "FS28_codon_counts.tsv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
columns = [
    ("Species", "Species"),
    ("Dataset", "Dataset"),
    ("Genome_assembly", "Genome_assembly"),
    ("Chromosome", "Chromosome"),
    ("Editing_positions_1based", "Editing_positions_1based"),
    ("Strand", "Strand"),
    ("Gene_ID", "Gene_ID"),
    ("Gene_symbol", "Gene_symbol"),
    ("Transcript_ID", "Transcript_ID"),
    (
        "Codon_genomic_positions_translation_order",
        "Codon_genomic_positions_translation_order",
    ),
    ("Reference_codon", "Reference_codon"),
    ("SN_class", "SN_class"),
    ("Am", "Am"),
    ("Im", "Im"),
    ("Ar", "Ar"),
    ("Ir", "Ir"),
    ("Lm", "Lm_plotted"),
    ("Lr", "Lr_plotted"),
    ("RNA_total_count", "RNA_total_count"),
    ("Asite_total_count", "Asite_total_count"),
    ("RNA_other_count", "RNA_other_count"),
    ("Asite_other_count", "Asite_other_count"),
    ("Lm_A_I_only", "Lm_A_I_only"),
    ("Lr_A_I_only", "Lr_A_I_only"),
    ("FS28_point_ID", "FS28_point_ID"),
    ("CDS_length", "CDS_length"),
    ("Transcript_length", "Transcript_length"),
]
numeric = {
    "CDS_length",
    "Transcript_length",
    "Am",
    "Im",
    "Ar",
    "Ir",
    "Lm",
    "Lr",
    "RNA_total_count",
    "Asite_total_count",
    "RNA_other_count",
    "Asite_other_count",
    "Lm_A_I_only",
    "Lr_A_I_only",
}
wb = xlsxwriter.Workbook(str(P / "TableS13.tmp.xlsx"))
wb.set_properties(
    {"title": "Table S13. Mammalian editing codons in Figure S28: longest-CDS analysis"}
)
hfmt = wb.add_format(
    {
        "bold": True,
        "bg_color": "#DCEAF3",
        "text_wrap": True,
        "valign": "top",
        "border": 1,
    }
)
fnum = wb.add_format({"num_format": "0.##########"})
wrap = wb.add_format({"text_wrap": True, "valign": "top"})
ws = wb.add_worksheet("TableS13")
ws.freeze_panes(1, 2)
ws.set_row(0, 46)
ws.write_row(0, 0, [x[0] for x in columns], hfmt)
ws.set_column(0, 11, 23)
ws.set_column(12, 23, 18, fnum)
ws.set_column(24, 24, 55)
with open(str(P / "TableS13.tsv"), "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow([x[0] for x in columns])
    for i, r in enumerate(rows, 1):
        v = [float(r[src]) if name in numeric else r[src] for name, src in columns]
        ws.write_row(i, 0, v)
        w.writerow(v)
        assert (
            abs(float(r["Lm_plotted"]) - float(r["Im"]) / float(r["RNA_total_count"]))
            < 1e-12
        )
        assert (
            abs(float(r["Lr_plotted"]) - float(r["Ir"]) / float(r["Asite_total_count"]))
            < 1e-12
        )
ws.autofilter(0, 0, len(rows), len(columns) - 1)
wb.close()
with zipfile.ZipFile(str(P / "TableS13.tmp.xlsx")) as z:
    assert z.testzip() is None
(P / "TableS13.tmp.xlsx").replace(P / "TableS13.xlsx")
qc = {
    "rows": len(rows),
    "unique_plot_points": len({r["FS28_point_ID"] for r in rows}),
    "species_counts": dict(Counter(r["Species"] for r in rows)),
    "scope": "longest-CDS FS28 plotted codons only",
    "Lm_Lr_match_figure": True,
}
pass
print(json.dumps(qc, indent=2))
