import os
from pathlib import Path
from collections import defaultdict, Counter
import re, bisect, csv, json

OUT = Path(__file__).resolve().parent
BASE = Path(os.environ.get("MAMMALS_DATA_ROOT", OUT.parent.parent / "11.mammals"))
OLD = Path(
    os.environ.get("MAMMALS_REFERENCE_ROOT", "/path/to/Proj1_RNAediting/06_mammals")
)
CONFIG = {
    "human": (
        "Homo_sapiens.GRCh38.87",
        "GRCh38.87_CDS_position",
        "Homo_sapiens.GRCh38.cdna.simple.fasta",
        "TABLE1_hg38_v3.exonic.txt2",
    ),
    "mouse": (
        "Mus_musculus.GRCm38.87",
        "GRCm38.87_CDS_position",
        "Mus_musculus.GRCm38.cdna.simple.fasta",
        str(BASE / "TABLE1_mm10_v3.coding"),
    ),
}
summary = {}
for species, (gtf, cdsfile, fasta, editfile) in CONFIG.items():
    cds = {}
    for line in open(OLD / cdsfile):
        f = line.split()
        cds[f[0]] = (int(f[1]), int(f[2]))
    edits = defaultdict(set)
    for line in open(OLD / editfile):
        f = line.split("\t")
        chrom = f[0].removeprefix("chr")
        pos = int(f[1])
        strand = f[3]
        edits[(chrom, strand)].add(pos)
    edits = {k: sorted(v) for k, v in edits.items()}
    exons = defaultdict(list)
    info = {}
    for line in open(OLD / (gtf + ".gtf")):
        if line.startswith("#"):
            continue
        f = line.split("\t", 8)
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        tx = m[1]
        if tx not in cds:
            continue
        exons[tx].append((int(f[3]), int(f[4])))
        info[tx] = (f[0].removeprefix("chr"), f[6])
    cand = defaultdict(dict)
    maps = {}
    for tx, blocks in exons.items():
        chrom, strand = info[tx]
        blocks = sorted(set(blocks), reverse=strand == "-")
        maps[tx] = blocks
        positions = edits.get((chrom, strand), [])
        offset = 0
        cs, ce = cds[tx]
        for lo, hi in blocks:
            for pos in positions[
                bisect.bisect_left(positions, lo) : bisect.bisect_right(positions, hi)
            ]:
                tpos = offset + (pos - lo + 1 if strand == "+" else hi - pos + 1)
                if cs <= tpos <= ce:
                    start = cs + (tpos - cs) // 3 * 3
                    if start + 2 <= ce:
                        cand[tx].setdefault(start, {})[tpos - start] = pos
            offset += hi - lo + 1
    qc = Counter()
    candidates = []

    def consume(tx, seq):
        if tx not in cand:
            return
        for start, positions in cand[tx].items():
            ref = seq[start - 1 : start + 2].upper()
            mask = 0
            valid = []
            if len(ref) != 3:
                qc["reference_short_codon"] += 1
                continue
            for j, pos in positions.items():
                if ref[j] == "A":
                    mask |= 1 << j
                    valid.append(str(pos))
                else:
                    qc["annotated_position_not_reference_A"] += 1
            if mask:
                candidates.append((tx, start, mask, ref, ",".join(valid)))

    tx = None
    parts = []
    for line in open(BASE / fasta):
        if line.startswith(">"):
            if tx is not None:
                consume(tx, "".join(parts))
            tx = line[1:].split()[0].split(".")[0]
            parts = []
        elif tx in cand:
            parts.append(line.strip())
    if tx is not None:
        consume(tx, "".join(parts))
    with open(OUT / f"{species}.transcripts.tsv", "w") as f:
        for tx in sorted(maps):
            chrom, strand = info[tx]
            cs, ce = cds[tx]
            f.write(
                "\t".join(
                    [
                        tx,
                        chrom,
                        strand,
                        str(cs),
                        str(ce),
                        ",".join(f"{a}-{b}" for a, b in maps[tx]),
                    ]
                )
                + "\n"
            )
    with open(OUT / f"{species}.candidates.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "transcript",
                "codon_start_tx",
                "edit_mask",
                "ref_codon",
                "editing_genome_positions",
            ]
        )
        w.writerows(candidates)
    summary[species] = {
        "input_genomic_edit_sites": sum(map(len, edits.values())),
        "coding_transcripts_with_exon_mapping": len(maps),
        "candidate_transcript_codons": len(candidates),
        "qc": dict(qc),
    }
    print(species, summary[species], flush=True)
manifest = []
for species in CONFIG:
    version = "v6"
    meta = BASE / "final" / f"{species}_downstream_{version}" / "pair_metadata.used.tsv"
    for r in csv.DictReader(open(meta), delimiter="\t"):
        tissue = r["tissue_group"]
        if species != "human" and tissue not in ("brain", "liver", "testis"):
            continue
        label = (
            {
                "Normal_brain": "brain1",
                "brain": "brain2",
                "Uninfected_Fibroblasts": "fibroblasts",
            }.get(tissue, tissue)
            if species == "human"
            else tissue
        )
        qc = dict(
            csv.reader(
                open(
                    BASE / "final" / (r["output_prefix"] + ".RNA_RPF_A2G_codon.qc.tsv")
                ),
                delimiter="\t",
            )
        )
        for kind, key in [("RNA", "rna_file"), ("RPF", "rpf_file")]:
            p = Path(qc[key])
            assert p.is_file(), p
            manifest.append([species, label, kind, r["biological_replicate"], str(p)])
with open(OUT / "inputs.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["species", "dataset", "library", "replicate", "path"])
    w.writerows(manifest)
