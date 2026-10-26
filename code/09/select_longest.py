import os
import argparse
import hashlib
from pathlib import Path
import csv, json, shutil, itertools
from collections import defaultdict, Counter

P = Path(__file__).resolve().parent
B = Path(os.environ.get("MAMMALS_DATA_ROOT", P.parent.parent / "11.mammals"))
A = P / "reference_all_transcripts"
A.mkdir(exist_ok=True)
config = {"human": "Homo_sapiens.GRCh38.87", "mouse": "Mus_musculus.GRCm38.87"}
parser = argparse.ArgumentParser(description="Select each gene's longest CDS before coverage filtering.")
parser.add_argument('--refresh-reference', action='store_true',
                    help='Snapshot newly generated full inputs; use only immediately after prepare.py.')
args = parser.parse_args()
cache_names = [sp + suffix for sp in config for suffix in ('.transcripts.tsv', '.candidates.tsv')]
state_path = A / 'selection_input_state.json'


def file_hashes(directory):
    return {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
            for name in cache_names}


# A plain rerun must use the same full reference and previously selected outputs.
# Newly prepared inputs require explicit refresh, preventing silent reuse of an
# older candidate cache. Check every species before overwriting any outputs.
if args.refresh_reference or not any((A / name).exists() for name in cache_names):
    file_hashes(P)  # Validate that all new input files exist before copying.
    for name in cache_names:
        shutil.copyfile(P / name, A / name)
else:
    if not state_path.exists():
        raise RuntimeError('Untracked reference cache. Run prepare.py, then select_longest.py '
                           '--refresh-reference, or use a fresh working directory.')
    state = json.loads(state_path.read_text())
    if file_hashes(A) != state['reference_files'] or file_hashes(P) != state['selected_files']:
        raise RuntimeError('Reference/candidate inputs changed since selection. Run prepare.py, '
                           'then select_longest.py --refresh-reference to use the new full inputs.')
qc = {}
for sp, stem in config.items():
    info = list(csv.DictReader(open(B / (stem + "_txinfo.tsv")), delimiter="\t"))
    bygene = defaultdict(list)
    allinfo = {}
    for r in info:
        try:
            c, t = int(r["cds_len"]), int(r["tx_len"])
        except ValueError:
            continue
        if c <= 0 or t <= 0:
            continue
        r["cds_len"] = c
        r["tx_len"] = t
        bygene[r["gene_id"]].append(r)
        allinfo[r["tx_name"]] = r
    rank = lambda r: (-r["cds_len"], -r["tx_len"], r["tx_name"])
    chosen = {g: sorted(rr, key=rank)[0] for g, rr in bygene.items()}
    wanted = {r["tx_name"] for r in chosen.values()}
    meta_path = A / (sp + ".transcripts.tsv")
    meta = {f[0]: f for f in csv.reader(open(meta_path), delimiter="\t")}
    available = wanted & meta.keys()
    missing = wanted - meta.keys()
    with open(P / (sp + ".transcripts.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerows(meta[tx] for tx in sorted(available))

    def gpos(tx, p):
        f = meta[tx]
        off = 0
        for bl in f[5].split(","):
            lo, hi = map(int, bl.split("-"))
            if off < p <= off + hi - lo + 1:
                return lo + p - off - 1 if f[2] == "+" else hi - (p - off - 1)
            off += hi - lo + 1
        raise ValueError((tx, p))

    with open(A / (sp + ".candidates.tsv")) as candidate_file:
        candidate_reader = csv.DictReader(candidate_file, delimiter="\t")
        candidate_fields = candidate_reader.fieldnames
        if not candidate_fields or not {'transcript', 'codon_start_tx', 'editing_genome_positions'}.issubset(candidate_fields):
            raise ValueError(f'{sp}.candidates.tsv: missing required column headers')
        candidates = list(candidate_reader)
    groups = defaultdict(list)
    for r in candidates:
        tx = r["transcript"]
        if tx not in available:
            continue
        p = int(r["codon_start_tx"])
        k = (meta[tx][1], meta[tx][2], *(gpos(tx, p + j) for j in range(3)))
        groups[k].append(r)
    kept = []
    collisions = []
    for k, rr in sorted(groups.items()):
        rr = sorted(
            rr,
            key=lambda r: (*rank(allinfo[r["transcript"]]), int(r["codon_start_tx"])),
        )
        kept.extend(rr)
        if len(rr) > 1:
            collisions.append(
                [*k, rr[0]["transcript"], ";".join(r["transcript"] for r in rr[1:])]
            )
    with open(P / (sp + ".candidates.tsv"), "w") as f:
        w = csv.DictWriter(
            f, fieldnames=candidate_fields, delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(kept)
    with open(P / (sp + ".selected_transcripts.tsv"), "w") as f:
        head = [
            "species",
            "gene_id",
            "gene_name",
            "transcript",
            "CDS_length",
            "transcript_length",
            "available_in_count_reference",
            "selection_rule",
        ]
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(head)
        for g, r in sorted(chosen.items()):
            w.writerow(
                [
                    sp,
                    g,
                    r["gene_name"],
                    r["tx_name"],
                    r["cds_len"],
                    r["tx_len"],
                    r["tx_name"] in available,
                    "max CDS length; max transcript length; smallest transcript ID",
                ]
            )
    with open(P / (sp + ".overlapping_gene_codon_dedup.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "chromosome",
                "strand",
                "base1",
                "base2",
                "base3",
                "annotation_first_transcript",
                "other_gene_longest_transcripts",
            ]
        )
        w.writerows(collisions)
    mapped = set()
    for r in kept:
        for pos in r["editing_genome_positions"].split(","):
            mapped.add((meta[r["transcript"]][1], meta[r["transcript"]][2], pos))
    qc[sp] = {
        "coding_genes_in_txinfo": len(chosen),
        "selected_transcripts_available": len(available),
        "selected_transcripts_unavailable": len(missing),
        "unavailable_transcripts": sorted(missing),
        "previous_transcript_candidate_codons": len(candidates),
        "fixed_candidate_genomic_codons": len(groups),
        "fixed_candidate_transcript_codons": len(kept),
        "fixed_candidate_reference_A_sites": len(mapped),
        "overlapping_gene_codon_collisions": len(collisions),
        "input_txinfo": str(B / (stem + "_txinfo.tsv")),
    }
    print(
        sp,
        {k: v for k, v in qc[sp].items() if k != "unavailable_transcripts"},
        flush=True,
    )
state_path.write_text(json.dumps({'reference_files': file_hashes(A),
                                  'selected_files': file_hashes(P)}, indent=2) + '\n')
