#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import argparse
from collections import defaultdict

AA3_TO_AA1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Sec": "U",
    "Pyl": "O",
    "Ter": "*",
    "Stop": "*",
}


def aa_to_one(x: str) -> str:
    x = x.strip()
    if len(x) == 1:
        return x.upper()
    if x in AA3_TO_AA1:
        return AA3_TO_AA1[x]
    raise ValueError(f"Unknown amino acid code: {x}")


def read_tsv(path):
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield row


def load_fbtr_to_fbpp(map_tsv):
    """
    Input columns:
    flybase_name, FBgn, FBtr, FBpp, length
    """
    fbtr_to_fbpp = defaultdict(list)

    for row in read_tsv(map_tsv):
        fbtr = row["FBtr"].strip()
        fbpp = row["FBpp"].strip()
        if fbtr and fbpp:
            fbtr_to_fbpp[fbtr].append(fbpp)

    cleaned = {}
    for fbtr, vals in fbtr_to_fbpp.items():
        seen = set()
        uniq = []
        for v in vals:
            if v not in seen:
                uniq.append(v)
                seen.add(v)
        cleaned[fbtr] = uniq

    return cleaned


def load_peptides(peptide_tsv):
    """
    dmel_r6.12_tryptic_peptides.tsv
    columns:
    source_id, peptide_seq, aa_start_1based, aa_end_1based, missed_cleavages, ...
    """
    peptides_by_fbpp = defaultdict(list)

    for row in read_tsv(peptide_tsv):
        fbpp = row["source_id"].strip()
        pep = row["peptide_seq"].strip()
        aa_start = int(row["aa_start_1based"])
        aa_end = int(row["aa_end_1based"])
        mc = row.get("missed_cleavages", "")
        enzyme = row.get("enzyme", "")

        peptides_by_fbpp[fbpp].append(
            {
                "peptide_seq": pep,
                "aa_start": aa_start,
                "aa_end": aa_end,
                "missed_cleavages": mc,
                "enzyme": enzyme,
            }
        )

    return peptides_by_fbpp


def make_headers(fbpp, aa_pos, pep_start, pep_end):
    """
    WT:   >FBppxxxxxxx_start_end
    EDIT: >FBppxxxxxxx_start_editpos_end
    """
    wt_header = f">{fbpp}_{pep_start}_{pep_end}"
    edit_header = f">{fbpp}_{pep_start}_{aa_pos}_{pep_end}"
    return wt_header, edit_header


def main():
    parser = argparse.ArgumentParser(
        description="Generate WT and edited peptide FASTA using complete FlyBase FBtr->FBpp mapping."
    )
    parser.add_argument("--map", required=True, help="flybase_fbtr_fbpp_complete.tsv")
    parser.add_argument(
        "--edits", required=True, help="mel_both_3103_editingsites_874mis.pos"
    )
    parser.add_argument(
        "--peptides", required=True, help="dmel_r6.12_tryptic_peptides.tsv"
    )
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    parser.add_argument(
        "--strict-ref-check",
        action="store_true",
        help="Skip peptide if residue at aa_pos does not match refAA",
    )
    args = parser.parse_args()

    fbtr_to_fbpp = load_fbtr_to_fbpp(args.map)
    peptides_by_fbpp = load_peptides(args.peptides)

    wt_fasta = args.out_prefix + ".wt.fasta"
    edit_fasta = args.out_prefix + ".edit.fasta"
    summary_tsv = args.out_prefix + ".summary.tsv"
    skipped_tsv = args.out_prefix + ".skipped.tsv"

    n_wt = 0
    n_edit = 0
    n_skip = 0

    with open(wt_fasta, "w", encoding="utf-8") as wt_out, open(
        edit_fasta, "w", encoding="utf-8"
    ) as edit_out, open(
        summary_tsv, "w", encoding="utf-8", newline=""
    ) as sum_out, open(
        skipped_tsv, "w", encoding="utf-8", newline=""
    ) as skip_out:

        sum_writer = csv.writer(sum_out, delimiter="\t")
        skip_writer = csv.writer(skip_out, delimiter="\t")

        sum_writer.writerow(
            [
                "site_id",
                "gene",
                "FBgn",
                "FBtr",
                "FBpp",
                "aa_pos",
                "refAA",
                "altAA",
                "pep_start",
                "pep_end",
                "local_pos_in_peptide",
                "wt_header",
                "edit_header",
                "wt_peptide",
                "edit_peptide",
                "missed_cleavages",
                "enzyme",
            ]
        )

        skip_writer.writerow(
            [
                "site_id",
                "gene",
                "FBgn",
                "FBtr",
                "FBpp",
                "aa_pos",
                "refAA",
                "altAA",
                "reason",
                "details",
            ]
        )

        for row in read_tsv(args.edits):
            site_id = row["site_id"].strip()
            gene = row["gene"].strip()
            fbgn = row["FBgn"].strip()
            fbtr = row["FBtr"].strip()
            aa_pos = int(row["aa_pos"])
            ref1 = aa_to_one(row["refAA"])
            alt1 = aa_to_one(row["altAA"])

            fbpps = fbtr_to_fbpp.get(fbtr, [])
            if not fbpps:
                n_skip += 1
                skip_writer.writerow(
                    [
                        site_id,
                        gene,
                        fbgn,
                        fbtr,
                        "",
                        aa_pos,
                        ref1,
                        alt1,
                        "FBtr_not_found_in_map",
                        "",
                    ]
                )
                continue

            any_hit = False

            for fbpp in fbpps:
                peps = peptides_by_fbpp.get(fbpp, [])
                if not peps:
                    n_skip += 1
                    skip_writer.writerow(
                        [
                            site_id,
                            gene,
                            fbgn,
                            fbtr,
                            fbpp,
                            aa_pos,
                            ref1,
                            alt1,
                            "FBpp_not_found_in_peptides",
                            "",
                        ]
                    )
                    continue

                for pep_info in peps:
                    pep = pep_info["peptide_seq"]
                    pep_start = pep_info["aa_start"]
                    pep_end = pep_info["aa_end"]

                    if not (pep_start <= aa_pos <= pep_end):
                        continue

                    any_hit = True
                    local_pos = aa_pos - pep_start + 1
                    idx = local_pos - 1

                    observed_ref = pep[idx]

                    if observed_ref != ref1:
                        reason = "refAA_mismatch"
                        detail = f"peptide_has={observed_ref}"
                        if args.strict_ref_check:
                            n_skip += 1
                            skip_writer.writerow(
                                [
                                    site_id,
                                    gene,
                                    fbgn,
                                    fbtr,
                                    fbpp,
                                    aa_pos,
                                    ref1,
                                    alt1,
                                    reason,
                                    detail,
                                ]
                            )
                            continue

                    edit_pep = pep[:idx] + alt1 + pep[idx + 1 :]

                    wt_header, edit_header = make_headers(
                        fbpp=fbpp, aa_pos=aa_pos, pep_start=pep_start, pep_end=pep_end
                    )

                    wt_out.write(wt_header + "\n")
                    wt_out.write(pep + "\n")

                    edit_out.write(edit_header + "\n")
                    edit_out.write(edit_pep + "\n")

                    n_wt += 1
                    n_edit += 1

                    sum_writer.writerow(
                        [
                            site_id,
                            gene,
                            fbgn,
                            fbtr,
                            fbpp,
                            aa_pos,
                            ref1,
                            alt1,
                            pep_start,
                            pep_end,
                            local_pos,
                            wt_header[1:],
                            edit_header[1:],
                            pep,
                            edit_pep,
                            pep_info["missed_cleavages"],
                            pep_info["enzyme"],
                        ]
                    )

            if not any_hit:
                n_skip += 1
                skip_writer.writerow(
                    [
                        site_id,
                        gene,
                        fbgn,
                        fbtr,
                        ",".join(fbpps),
                        aa_pos,
                        ref1,
                        alt1,
                        "no_covering_peptide",
                        "",
                    ]
                )

    print(f"WT FASTA: {wt_fasta} ({n_wt} records)")
    print(f"EDIT FASTA: {edit_fasta} ({n_edit} records)")
    print(f"SUMMARY: {summary_tsv}")
    print(f"SKIPPED: {skipped_tsv} ({n_skip} rows)")


if __name__ == "__main__":
    main()
