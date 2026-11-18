#!/usr/bin/env python3
"""Index complete non-stop CDS codons and resolve an existing unscaled master manifest."""

import argparse, csv, json, hashlib
from pathlib import Path

SPECS = {
    "dmel": ("dmelr6.12", "m6A_list.dmel_longest_cds.tsv"),
    "dere": ("derer1.05", "m6A_list.droEre2ToDm6.longest_cds.tsv"),
    "dsec": ("dsecr1.3", "m6A_list.droSec1ToDm6.longest_cds.tsv"),
    "dsim": ("dsimr2.02", "m6A_list.droSim1.longest_cds.tsv"),
}


def write(path, rows, fields=None):
    with path.open("w") as f:
        w = csv.DictWriter(
            f, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotation-dir", type=Path, required=True)
    ap.add_argument("--master-manifest", type=Path, required=True)
    ap.add_argument(
        "--master-dir",
        type=Path,
        help="Relocate historical manifest files to this directory by basename",
    )
    ap.add_argument("--output-dir", type=Path, required=True)
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    refs = []
    qc = []
    for sp, (tag, m6name) in SPECS.items():
        annotation = a.annotation_dir / f"{tag}.gtf.cds.bed.ordered.1codon.longest.tsv"
        m6 = a.annotation_dir / m6name
        if not m6.exists():
            m6 = a.annotation_dir / f"{sp}.m6a.tsv"
        assert m6.is_file(), m6
        lengths = {}
        sha = hashlib.sha256()
        with annotation.open("rb") as f:
            for line in f:
                sha.update(line)
                v = line.rstrip().split(b"\t")
                t = v[2].decode()
                idx = int(v[6])
                start = int(v[5])
                assert idx == lengths.get(t, 0) + 1 and start == 3 * idx - 2, (
                    annotation,
                    t,
                    idx,
                    start,
                )
                lengths[t] = idx
        rows = [dict(transcript=t, codonlength=L) for t, L in sorted(lengths.items())]
        out = a.output_dir / f"{sp}.cds_lengths.tsv"
        write(out, rows)
        refs.append(
            dict(species=sp, cds_lengths=str(out.resolve()), m6a_map=str(m6.resolve()))
        )
        qc.append(
            dict(
                species=sp,
                annotation=str(annotation.resolve()),
                annotation_sha256=sha.hexdigest(),
                transcripts=len(lengths),
                codons=sum(lengths.values()),
                length_definition="complete annotated CDS codons, excluding separately annotated stop_codon",
            )
        )
        print(sp, len(lengths), sum(lengths.values()), flush=True)
    write(a.output_dir / "focal_reference_manifest.tsv", refs)
    with a.master_manifest.open() as f:
        masters = list(csv.DictReader(f, delimiter="\t"))
    selected = []
    for r in masters:
        if int(r["flank_size"]) != 30 or int(r["n_replicates"]) < 2:
            continue
        p = Path(r["output_path"])
        if a.master_dir:
            p = a.master_dir / p.name
        assert p.is_file(), p
        r["output_path"] = str(p.resolve())
        selected.append(r)
    assert len(selected) == len({r["group_id"] for r in selected})
    write(a.output_dir / "master_manifest.tsv", selected)
    pass
    print("masters", len(selected), flush=True)


if __name__ == "__main__":
    main()
