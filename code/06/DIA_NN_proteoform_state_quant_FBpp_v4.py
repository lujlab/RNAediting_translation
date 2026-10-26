#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DIA_NN_proteoform_state_quant_FBpp_v4.py

Purpose
-------
Quantify REF and multiple ALT amino-acid states from a DIA-NN
report.pr_matrix.tsv when the search FASTA / DIA-NN IDs are FBpp protein IDs
whose trailing numeric fields encode a protein-AA window plus zero or more
edited AA positions.

The RNA-editing site table is keyed by FBtr, therefore an FBtr<->FBpp mapping
table is used before matching against DIA-NN / FASTA IDs.

Core rules
----------
1. FBpp coordinate IDs: same logic as the previous FBtr script, except ALL
   numeric coordinates are already protein AA coordinates and are NOT divided
   by 3.

   For all trailing numeric fields:
     - remove ONE smallest value as window start;
     - remove ONE largest value as window end;
     - every remaining value is an edited AA position.

   Examples:
     FBppXXXX_83_135         -> start=83,  end=135, edits=[]
     FBppXXXX_83_100_135     -> start=83,  end=135, edits=[100]
     FBppXXXX_83_100_120_135 -> start=83,  end=135, edits=[100,120]
     FBppXXXX_83_83_135      -> start=83,  end=135, edits=[83]

2. DIA-NN Protein.Ids handling:
   - use Protein.Ids preferentially;
   - within one matrix row, de-duplicate by FBpp prefix;
   - for each FBpp keep ONLY the first listed coordinate-bearing FASTA ID;
   - do not fall through to a later same-FBpp ID if the first one fails;
   - Stripped.Sequence must actually be a substring of that exact FASTA entry.

3. Site matching:
   - merge site-table rows sharing FBtr + aa_pos + refAA (same AA/codon site);
   - convert FBtr -> FBpp using --id-map;
   - map Stripped.Sequence back to the retained FASTA sequence;
   - compute the peptide's true protein-AA interval;
   - only peptides that physically cover the target aa_pos are informative;
   - peptides not covering the target aa_pos are ignored completely.

4. REF / multiple-ALT state classification:
   - read the actual residue from Stripped.Sequence at target aa_pos;
   - observed residue == refAA -> REF;
   - observed residue != refAA -> ALT1, ALT2, ...;
   - actual observed residue is retained in AA;
   - known altAA values from the site table are numbered first in site-file
     order; newly observed AA states are marked UNLISTED_ALT;
   - the edited positions parsed from the FBpp ID are retained in the detail
     output as an independent audit of which AA positions that FASTA entry
     encodes as edited.

5. Intensity aggregation:
   - for site + AA-state + sample, precursor rows belonging to the same
     Stripped.Sequence are first SUMMED to one peptide-level intensity;
   - sum_intensity    = sum of positive peptide-level intensities;
   - median_intensity = median of positive peptide-level intensities;
   - n_peptides       = number of unique Stripped.Sequence peptides with
     intensity > 0 in that sample.
"""

import argparse
import csv
import gzip
import os
import re
from collections import defaultdict, OrderedDict
from statistics import median

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

SITE_COLUMNS = [
    "site_id",
    "gene",
    "FBgn",
    "FBtr",
    "refAA",
    "altAA",
    "trans_pos",
    "trans_len",
    "cds_pos",
    "cds_len",
    "aa_pos",
    "aa_len",
]


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def open_output(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return open(path, "w", encoding="utf-8", newline="")


def aa_to_one(x):
    x = str(x).strip()
    if len(x) == 1:
        return x.upper()
    if x in AA3_TO_AA1:
        return AA3_TO_AA1[x]
    raise ValueError(f"Unknown amino acid code: {x}")


def parse_float(x):
    if x is None:
        return 0.0
    x = str(x).strip()
    if x == "" or x.upper() in {"NA", "NAN", "NULL"}:
        return 0.0
    try:
        return float(x)
    except ValueError:
        return 0.0


def format_float(x):
    return f"{float(x):.10g}"


def normalize_seq(x):
    return str(x).strip().strip('"').strip("'").replace(" ", "").upper()


def clean_colname(x):
    return x.strip().lstrip("\ufeff").strip('"').strip("'").rstrip("\r")


def split_auto(line):
    line = line.rstrip("\n").rstrip("\r")
    if "\t" in line:
        return line.split("\t")
    return line.split()


def normalize_id(x):
    x = str(x).strip().strip('"').strip("'")
    if x.endswith("_DROME"):
        x = x[:-6]
    x = re.sub(r"^(FBtr\d+)\.\d+$", r"\1", x)
    x = re.sub(r"^(FBpp\d+)\.\d+$", r"\1", x)
    return x


def derive_unlisted_alt_path(out_path):
    if out_path.endswith(".tsv.gz"):
        return out_path[:-7] + ".unlisted_alt.tsv.gz"
    if out_path.endswith(".tsv"):
        return out_path[:-4] + ".unlisted_alt.tsv"
    return out_path + ".unlisted_alt.tsv"


# -----------------------------------------------------------------------------
# FBpp coordinate-ID parsing / FASTA
# -----------------------------------------------------------------------------


def parse_fbpp_coord_id(x):
    """
    Parse an FBpp coordinate-bearing ID using the same min/max/remnant rule as
    the FBtr version, but with direct protein-AA coordinates.

    One occurrence of the minimum is the window start, one occurrence of the
    maximum is the window end, and all remaining numeric fields are edited AA
    positions. Duplicate boundary values are therefore preserved as edits.
    """
    original = str(x).strip().strip('"').strip("'")
    clean = original[:-6] if original.endswith("_DROME") else original
    toks = clean.split("_")

    numeric = []
    i = len(toks) - 1
    while i >= 0 and toks[i].isdigit():
        numeric.append(int(toks[i]))
        i -= 1
    numeric = list(reversed(numeric))

    if len(numeric) < 2:
        raise ValueError(f"Cannot parse FBpp coordinate ID: {original}")

    prefix = normalize_id("_".join(toks[: len(toks) - len(numeric)]))
    if not re.fullmatch(r"FBpp\d+", prefix):
        raise ValueError(f"Not an FBpp coordinate ID: {original}")

    sorted_numeric = sorted(numeric)
    start = sorted_numeric[0]
    end = sorted_numeric[-1]
    # Slice removes exactly one minimum and one maximum occurrence.
    edited_positions = sorted_numeric[1:-1]

    if start > end:
        raise ValueError(f"Invalid FBpp coordinate window in {original}: {start}>{end}")

    return {
        "prefix": prefix,
        "start": start,
        "end": end,
        "edited_positions": edited_positions,
        "raw_numeric": numeric,
        "raw_numeric_sorted": sorted_numeric,
        "raw_start": start,
        "raw_end": end,
        "raw_edited_positions": list(edited_positions),
        "coordinate_system": "FBpp_protein_AA",
        "id_layout": "min_start_max_end_remaining_edits_AA",
    }


def read_fasta(path):
    header = None
    seq_chunks = []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, normalize_seq("".join(seq_chunks))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, normalize_seq("".join(seq_chunks))


def extract_fasta_ids(header):
    """Return likely ID strings from plain or UniProt-like FASTA headers."""
    ids = []
    m = re.match(r"^(sp|tr)\|([^|]+)\|([^ ]+)\s*(.*)$", header)
    if m:
        ids.extend([m.group(2).strip(), m.group(3).strip()])
        rest = m.group(4)
        gm = re.search(r"\bGN=([^ ]+)", rest)
        if gm:
            ids.append(gm.group(1).strip())
    else:
        ids.append(header.split()[0].strip())

    out = []
    seen = set()
    for x in ids:
        if x and x not in seen:
            out.append(x)
            seen.add(x)
    return out


def load_fasta_coord_map(fasta_files):
    if isinstance(fasta_files, (str, bytes)):
        fasta_files = [fasta_files]

    expanded = []
    for item in fasta_files:
        for p in str(item).split(","):
            p = p.strip()
            if p:
                expanded.append(p)
    if not expanded:
        raise ValueError("No FASTA file supplied")

    id_to_records = defaultdict(list)
    seen_by_id = defaultdict(set)
    total_entries = 0
    total_indexed = 0
    length_mismatch = 0

    for fasta_file in expanded:
        if not os.path.isfile(fasta_file):
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

        file_entries = 0
        file_indexed = 0
        for header, seq in read_fasta(fasta_file):
            file_entries += 1
            total_entries += 1

            for one_id in extract_fasta_ids(header):
                try:
                    parsed = parse_fbpp_coord_id(one_id)
                except ValueError:
                    continue

                expected_len = parsed["end"] - parsed["start"] + 1
                if len(seq) != expected_len:
                    length_mismatch += 1

                key = (
                    one_id,
                    seq,
                    parsed["prefix"],
                    parsed["start"],
                    parsed["end"],
                    tuple(parsed["edited_positions"]),
                )
                if key in seen_by_id[one_id]:
                    continue
                seen_by_id[one_id].add(key)

                rec = {
                    "fasta_id": one_id,
                    "fbpp": parsed["prefix"],
                    "fasta_seq": seq,
                    "fasta_start": parsed["start"],
                    "fasta_end": parsed["end"],
                    "fasta_edited_positions": list(parsed["edited_positions"]),
                    "raw_fasta_start": parsed["raw_start"],
                    "raw_fasta_end": parsed["raw_end"],
                    "raw_fasta_edited_positions": list(parsed["raw_edited_positions"]),
                    "raw_numeric": list(parsed["raw_numeric"]),
                    "coordinate_system": parsed["coordinate_system"],
                    "id_layout": parsed["id_layout"],
                    "fasta_header": header,
                    "source_fasta": fasta_file,
                }
                id_to_records[one_id].append(rec)
                file_indexed += 1
                total_indexed += 1

        print(
            f"Loaded FASTA: {fasta_file} | entries={file_entries} | "
            f"indexed_records={file_indexed}"
        )

    print(
        f"Merged FASTA files: {len(expanded)} | entries={total_entries} | "
        f"unique ID keys={len(id_to_records)} | indexed_records={total_indexed}"
    )
    print(
        f"FASTA window/sequence length mismatches: {length_mismatch} "
        f"(expected length = end-start+1)"
    )
    return id_to_records


def find_all_offsets(haystack, needle):
    offsets = []
    start = 0
    while True:
        idx = haystack.find(needle, start)
        if idx == -1:
            break
        offsets.append(idx)
        start = idx + 1
    return offsets


def make_substring_candidates(rec, peptide_seq):
    """Map Stripped.Sequence to one exact FBpp FASTA entry by substring."""
    peptide_seq = normalize_seq(peptide_seq)
    candidates = []
    for offset in find_all_offsets(rec["fasta_seq"], peptide_seq):
        real_start = rec["fasta_start"] + offset
        real_end = real_start + len(peptide_seq) - 1
        effective_edits = [
            p
            for p in rec.get("fasta_edited_positions", [])
            if real_start <= p <= real_end
        ]
        candidates.append(
            {
                "candidate_id": rec["fasta_id"],
                "fbpp": rec["fbpp"],
                "peptide_seq": peptide_seq,
                "real_start": real_start,
                "real_end": real_end,
                "substring_offset_0based": offset,
                "fasta_start": rec["fasta_start"],
                "fasta_end": rec["fasta_end"],
                "fasta_edited_positions": list(rec.get("fasta_edited_positions", [])),
                "effective_edited_positions": effective_edits,
                "raw_fasta_start": rec.get("raw_fasta_start", ""),
                "raw_fasta_end": rec.get("raw_fasta_end", ""),
                "raw_fasta_edited_positions": list(
                    rec.get("raw_fasta_edited_positions", [])
                ),
                "raw_numeric": list(rec.get("raw_numeric", [])),
                "coordinate_system": rec["coordinate_system"],
                "id_layout": rec["id_layout"],
                "fasta_header": rec["fasta_header"],
                "source_fasta": rec["source_fasta"],
            }
        )
    return candidates


# -----------------------------------------------------------------------------
# Site table + FBtr<->FBpp map
# -----------------------------------------------------------------------------


def load_sites(site_file):
    rows = []
    with open_text(site_file) as f:
        header_line = f.readline()
        header = [clean_colname(x) for x in split_auto(header_line)]
        missing = [c for c in SITE_COLUMNS if c not in header]
        if missing:
            raise ValueError(
                "Missing required columns in site file: " + ",".join(missing)
            )

        for line_no, line in enumerate(f, start=2):
            if not line.strip():
                continue
            vals = split_auto(line)
            if len(vals) < len(header):
                vals += [""] * (len(header) - len(vals))
            row_all = dict(zip(header, vals))
            row = {
                k: row_all.get(k, "").strip().strip('"').rstrip("\r")
                for k in SITE_COLUMNS
            }
            row["FBtr"] = normalize_id(row["FBtr"])
            row["aa_pos"] = int(row["aa_pos"])
            row["_ref1"] = aa_to_one(row["refAA"])
            row["_alt1"] = aa_to_one(row["altAA"])
            row["_line_no"] = line_no
            rows.append(row)
    return rows


def group_sites_by_amino_acid(sites):
    """Merge nucleotide editing rows affecting the same FBtr AA position."""
    groups = OrderedDict()
    for row in sites:
        key = (row["FBtr"], int(row["aa_pos"]), row["_ref1"])
        if key not in groups:
            groups[key] = {
                "site": f"{row['FBtr']}:{row['aa_pos']}",
                "gene": row.get("gene", ""),
                "FBgn": row.get("FBgn", ""),
                "FBtr": row["FBtr"],
                "aa_pos": int(row["aa_pos"]),
                "refAA": row.get("refAA", ""),
                "refAA1": row["_ref1"],
                "aa_len": row.get("aa_len", ""),
                "source_site_ids": [],
                "source_cds_pos": [],
                "source_trans_pos": [],
                "known_alt_AAs": [],
            }
        g = groups[key]
        for field, outkey in [
            ("site_id", "source_site_ids"),
            ("cds_pos", "source_cds_pos"),
            ("trans_pos", "source_trans_pos"),
        ]:
            v = row.get(field, "")
            if v and v not in g[outkey]:
                g[outkey].append(v)
        if row["_alt1"] not in g["known_alt_AAs"]:
            g["known_alt_AAs"].append(row["_alt1"])
    return list(groups.values())


def load_fbtr_fbpp_map(path):
    """
    Read a FlyBase mapping table containing columns FBtr and FBpp.

    Example columns:
      flybase_name  FBgn  FBtr  FBpp  length
    """
    fbtr_to_fbpp = OrderedDict()
    rows_total = 0
    rows_with_pair = 0

    with open_text(path) as f:
        header_line = f.readline()
        header = [clean_colname(x) for x in split_auto(header_line)]
        if "FBtr" not in header or "FBpp" not in header:
            raise ValueError(
                "ID map must contain columns named FBtr and FBpp. Detected: "
                + ",".join(header)
            )
        i_fbtr = header.index("FBtr")
        i_fbpp = header.index("FBpp")

        for line in f:
            if not line.strip():
                continue
            rows_total += 1
            vals = split_auto(line)
            if len(vals) <= max(i_fbtr, i_fbpp):
                continue
            fbtr = normalize_id(vals[i_fbtr])
            fbpp = normalize_id(vals[i_fbpp])
            if not re.fullmatch(r"FBtr\d+", fbtr):
                continue
            if not re.fullmatch(r"FBpp\d+", fbpp):
                continue
            rows_with_pair += 1
            if fbtr not in fbtr_to_fbpp:
                fbtr_to_fbpp[fbtr] = []
            if fbpp not in fbtr_to_fbpp[fbtr]:
                fbtr_to_fbpp[fbtr].append(fbpp)

    n_multi = sum(1 for v in fbtr_to_fbpp.values() if len(v) > 1)
    print(
        f"Loaded FBtr-FBpp map: {path} | rows={rows_total} | "
        f"valid_pairs={rows_with_pair} | FBtr_keys={len(fbtr_to_fbpp)} | "
        f"FBtr_with_multiple_FBpp={n_multi}"
    )
    return fbtr_to_fbpp


# -----------------------------------------------------------------------------
# DIA-NN matrix
# -----------------------------------------------------------------------------


def parse_col_list(x):
    if x is None:
        return []
    x = str(x).strip().strip('"')
    if not x:
        return []
    return [
        item.strip().strip('"').strip("'")
        for item in x.split(";")
        if item.strip().strip('"').strip("'")
    ]


def first_id_per_fbpp(raw_ids):
    """
    Retain only the first coordinate-bearing ID for each FBpp prefix.

    Example:
      FBpp0306255_588_643;FBpp0306255_593_643;FBpp0306255_621_643
    -> FBpp0306255_588_643
    """
    kept = OrderedDict()
    for rid in raw_ids:
        variants = [rid]
        if rid.endswith("_DROME"):
            variants.append(rid[:-6])

        parsed = None
        for v in variants:
            try:
                parsed = parse_fbpp_coord_id(v)
                break
            except ValueError:
                pass
        if parsed is None:
            continue

        fbpp = parsed["prefix"]
        if fbpp not in kept:
            kept[fbpp] = rid
    return list(kept.values())


def load_matrix(
    matrix_file,
    fasta_coord_map,
    candidate_id_cols_1based,
    peptide_seq_col_1based,
    sample_start_col_1based,
):
    rows_by_fbpp = defaultdict(list)
    stats = {
        "matrix_rows_total": 0,
        "matrix_rows_with_candidate_ids": 0,
        "matrix_rows_without_candidate_ids": 0,
        "matrix_rows_with_fasta_substring_candidates": 0,
        "matrix_rows_first_ids_failed_fasta_substring": 0,
    }

    with open_text(matrix_file) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        header_clean = [clean_colname(x) for x in header]

        candidate_idxs = [x - 1 for x in candidate_id_cols_1based]
        seq_idx = peptide_seq_col_1based - 1
        sample_start_idx = sample_start_col_1based - 1

        if seq_idx >= len(header):
            raise ValueError("peptide_seq_col is out of range")
        if sample_start_idx >= len(header):
            raise ValueError("sample_start_col is out of range")

        protein_ids_idx = (
            header_clean.index("Protein.Ids") if "Protein.Ids" in header_clean else None
        )
        modified_seq_idx = (
            header_clean.index("Modified.Sequence")
            if "Modified.Sequence" in header_clean
            else None
        )
        precursor_charge_idx = (
            header_clean.index("Precursor.Charge")
            if "Precursor.Charge" in header_clean
            else None
        )
        precursor_id_idx = (
            header_clean.index("Precursor.Id")
            if "Precursor.Id" in header_clean
            else None
        )

        sample_names = header[sample_start_idx:]

        for line_no, row in enumerate(reader, start=2):
            stats["matrix_rows_total"] += 1
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))

            peptide_seq = normalize_seq(row[seq_idx])
            if not peptide_seq:
                stats["matrix_rows_without_candidate_ids"] += 1
                continue

            raw_ids = []
            if protein_ids_idx is not None:
                raw_ids = parse_col_list(row[protein_ids_idx])
            if not raw_ids:
                for idx in candidate_idxs:
                    if 0 <= idx < len(row):
                        raw_ids.extend(parse_col_list(row[idx]))

            retained_ids = first_id_per_fbpp(raw_ids)
            if not retained_ids:
                stats["matrix_rows_without_candidate_ids"] += 1
                continue
            stats["matrix_rows_with_candidate_ids"] += 1

            candidates = []
            for rid in retained_ids:
                lookup_ids = [rid]
                if rid.endswith("_DROME"):
                    lookup_ids.append(rid[:-6])

                # Crucial rule: use only this first retained ID for the FBpp.
                # Never fall through to a later same-FBpp ID from Protein.Ids.
                found_this_id = False
                for lid in lookup_ids:
                    recs = fasta_coord_map.get(lid, [])
                    if not recs:
                        continue
                    for rec in recs:
                        mapped = make_substring_candidates(rec, peptide_seq)
                        if mapped:
                            candidates.extend(mapped)
                            found_this_id = True
                    if found_this_id:
                        break

            unique_candidates = []
            seen = set()
            for cand in candidates:
                key = (
                    cand["fbpp"],
                    cand["candidate_id"],
                    cand["real_start"],
                    cand["real_end"],
                    tuple(cand.get("fasta_edited_positions", [])),
                    cand["peptide_seq"],
                )
                if key not in seen:
                    unique_candidates.append(cand)
                    seen.add(key)

            if not unique_candidates:
                stats["matrix_rows_first_ids_failed_fasta_substring"] += 1
                continue
            stats["matrix_rows_with_fasta_substring_candidates"] += 1

            intensities = [parse_float(x) for x in row[sample_start_idx:]]
            mrow = {
                "line_no": line_no,
                "peptide_seq": peptide_seq,
                "modified_seq": (
                    row[modified_seq_idx] if modified_seq_idx is not None else ""
                ),
                "precursor_charge": (
                    row[precursor_charge_idx]
                    if precursor_charge_idx is not None
                    else ""
                ),
                "precursor_id": (
                    row[precursor_id_idx] if precursor_id_idx is not None else ""
                ),
                "raw_candidate_ids": ";".join(raw_ids),
                "retained_candidate_ids": ";".join(retained_ids),
                "candidates": unique_candidates,
                "intensities": intensities,
            }

            added_fbpp = set()
            for cand in unique_candidates:
                fbpp = cand["fbpp"]
                if fbpp not in added_fbpp:
                    rows_by_fbpp[fbpp].append(mrow)
                    added_fbpp.add(fbpp)

    return sample_names, rows_by_fbpp, stats


# -----------------------------------------------------------------------------
# Classification / aggregation
# -----------------------------------------------------------------------------


def classify_candidate_at_site(cand, site_group):
    aa_pos = int(site_group["aa_pos"])
    if not (cand["real_start"] <= aa_pos <= cand["real_end"]):
        return None

    local_idx = aa_pos - cand["real_start"]
    pep = cand["peptide_seq"]
    if local_idx < 0 or local_idx >= len(pep):
        return None

    observed = pep[local_idx]
    ref = site_group["refAA1"]
    header_marks_target = aa_pos in cand.get("fasta_edited_positions", [])

    if observed == ref:
        return {
            "base_state": "REF",
            "AA": observed,
            "alt_match_status": "REF",
            "header_marks_target_edit": "yes" if header_marks_target else "no",
        }

    return {
        "base_state": "ALT",
        "AA": observed,
        "alt_match_status": (
            "KNOWN_ALT" if observed in site_group["known_alt_AAs"] else "UNLISTED_ALT"
        ),
        "header_marks_target_edit": "yes" if header_marks_target else "no",
    }


def choose_first_valid_mapping(mrow, fbpp, site_group):
    """Use first candidate mapping for this FBpp that actually covers aa_pos."""
    for cand in mrow["candidates"]:
        if cand["fbpp"] != fbpp:
            continue
        cls = classify_candidate_at_site(cand, site_group)
        if cls is not None:
            return cand, cls
    return None, None


def summarize(site_groups, sample_names, rows_by_fbpp, fbtr_to_fbpp):
    main_rows = []
    detail_rows = []
    stats = {
        "aa_sites_without_fbpp_map": 0,
        "aa_sites_with_fbpp_map": 0,
        "aa_sites_matched": 0,
        "aa_sites_unmatched_after_map": 0,
        "informative_precursor_rows": 0,
    }

    for site_group in site_groups:
        fbtr = normalize_id(site_group["FBtr"])
        mapped_fbpps = fbtr_to_fbpp.get(fbtr, [])
        if not mapped_fbpps:
            stats["aa_sites_without_fbpp_map"] += 1
            continue
        stats["aa_sites_with_fbpp_map"] += 1

        chosen_fbpp = None
        chosen_evidence = []

        # Mapping table is expected to be one-to-one. If more than one FBpp is
        # present for an FBtr, preserve file order and use the first FBpp that
        # yields informative peptide evidence; do not combine proteins.
        for fbpp in mapped_fbpps:
            route_evidence = []
            for mrow in rows_by_fbpp.get(fbpp, []):
                cand, cls = choose_first_valid_mapping(mrow, fbpp, site_group)
                if cand is None:
                    continue
                route_evidence.append((mrow, cand, cls))
            if route_evidence:
                chosen_fbpp = fbpp
                chosen_evidence = route_evidence
                break

        if not chosen_evidence:
            stats["aa_sites_unmatched_after_map"] += 1
            continue

        stats["aa_sites_matched"] += 1
        stats["informative_precursor_rows"] += len(chosen_evidence)

        # Sum precursor rows for the same Stripped.Sequence first.
        # key = (AA, sample_index, stripped_sequence)
        peptide_level = defaultdict(float)
        state_windows = defaultdict(set)
        all_windows = set()
        observed_alt_aas = set()
        detail_indices = []

        for mrow, cand, cls in chosen_evidence:
            aa = cls["AA"]
            if cls["base_state"] == "ALT":
                observed_alt_aas.add(aa)

            win = f"{cand['fasta_start']}-{cand['fasta_end']}"
            state_windows[aa].add(win)
            all_windows.add(win)

            for i, val in enumerate(mrow["intensities"]):
                peptide_level[(aa, i, mrow["peptide_seq"])] += val

            detail_rows.append(
                {
                    "site": site_group["site"],
                    "source_site_ids": ";".join(site_group["source_site_ids"]),
                    "source_cds_pos": ";".join(site_group["source_cds_pos"]),
                    "gene": site_group["gene"],
                    "FBgn": site_group["FBgn"],
                    "FBtr": fbtr,
                    "FBpp": chosen_fbpp,
                    "aa_pos": site_group["aa_pos"],
                    "refAA": site_group["refAA1"],
                    "known_alt_AAs": ";".join(site_group["known_alt_AAs"]),
                    "state": "",
                    "AA": aa,
                    "alt_match_status": cls["alt_match_status"],
                    "header_marks_target_edit": cls.get("header_marks_target_edit", ""),
                    "protein_id_used": cand["candidate_id"],
                    "raw_candidate_ids": mrow["raw_candidate_ids"],
                    "retained_candidate_ids": mrow["retained_candidate_ids"],
                    "peptide_seq": mrow["peptide_seq"],
                    "modified_seq": mrow["modified_seq"],
                    "precursor_charge": mrow["precursor_charge"],
                    "precursor_id": mrow["precursor_id"],
                    "real_start": cand["real_start"],
                    "real_end": cand["real_end"],
                    "window": win,
                    "fasta_edited_positions": ",".join(
                        map(str, cand.get("fasta_edited_positions", []))
                    ),
                    "effective_edited_positions": ",".join(
                        map(str, cand.get("effective_edited_positions", []))
                    ),
                    "raw_fasta_edited_positions": ",".join(
                        map(str, cand.get("raw_fasta_edited_positions", []))
                    ),
                    "substring_offset_0based": cand["substring_offset_0based"],
                    "matrix_line": mrow["line_no"],
                    "fasta_header": cand["fasta_header"],
                    "source_fasta": cand["source_fasta"],
                    "coordinate_system": cand["coordinate_system"],
                    "id_layout": cand["id_layout"],
                    "_intensities": mrow["intensities"],
                }
            )
            detail_indices.append(len(detail_rows) - 1)

        # Stable ALT numbering: known observed AAs first, then unlisted AAs.
        alt_order = []
        for aa in site_group["known_alt_AAs"]:
            if aa in observed_alt_aas and aa not in alt_order:
                alt_order.append(aa)
        for aa in sorted(observed_alt_aas):
            if aa not in alt_order:
                alt_order.append(aa)

        state_by_aa = {site_group["refAA1"]: "REF"}
        for j, aa in enumerate(alt_order, start=1):
            state_by_aa[aa] = f"ALT{j}"

        for di in detail_indices:
            aa = detail_rows[di]["AA"]
            detail_rows[di]["state"] = state_by_aa.get(aa, "ALT")

        output_aas = [site_group["refAA1"]] + alt_order
        for aa in output_aas:
            state = state_by_aa[aa]
            if aa == site_group["refAA1"]:
                alt_status = "REF"
            elif aa in site_group["known_alt_AAs"]:
                alt_status = "KNOWN_ALT"
            else:
                alt_status = "UNLISTED_ALT"

            windows = state_windows.get(aa, set()) or all_windows

            for i, sample in enumerate(sample_names):
                vals = [
                    pep_intensity
                    for (
                        aa_key,
                        sample_idx,
                        _pep,
                    ), pep_intensity in peptide_level.items()
                    if aa_key == aa and sample_idx == i and pep_intensity > 0
                ]

                main_rows.append(
                    {
                        "site": site_group["site"],
                        "source_site_ids": ";".join(site_group["source_site_ids"]),
                        "source_cds_pos": ";".join(site_group["source_cds_pos"]),
                        "gene": site_group["gene"],
                        "FBgn": site_group["FBgn"],
                        "FBtr": fbtr,
                        "FBpp": chosen_fbpp,
                        "aa_pos": site_group["aa_pos"],
                        "refAA": site_group["refAA1"],
                        "known_alt_AAs": ";".join(site_group["known_alt_AAs"]),
                        "window": ";".join(sorted(windows)),
                        "state": state,
                        "AA": aa,
                        "alt_match_status": alt_status,
                        "sample": sample,
                        "sum_intensity": sum(vals) if vals else 0.0,
                        "median_intensity": median(vals) if vals else 0.0,
                        "n_peptides": len(vals),
                    }
                )

    return main_rows, detail_rows, stats


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------


def write_main(path, rows):
    fields = [
        "site",
        "source_site_ids",
        "source_cds_pos",
        "gene",
        "FBgn",
        "FBtr",
        "FBpp",
        "aa_pos",
        "refAA",
        "known_alt_AAs",
        "window",
        "state",
        "AA",
        "alt_match_status",
        "sample",
        "sum_intensity",
        "median_intensity",
        "n_peptides",
    ]
    with open_output(path) as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for row in rows:
            out = dict(row)
            out["sum_intensity"] = format_float(out["sum_intensity"])
            out["median_intensity"] = format_float(out["median_intensity"])
            w.writerow({k: out.get(k, "") for k in fields})


def write_unlisted_alt(path, rows):
    subset = [r for r in rows if r.get("alt_match_status") == "UNLISTED_ALT"]
    write_main(path, subset)
    return len(subset)


def write_detail(path, rows, sample_names):
    fields = [
        "site",
        "source_site_ids",
        "source_cds_pos",
        "gene",
        "FBgn",
        "FBtr",
        "FBpp",
        "aa_pos",
        "refAA",
        "known_alt_AAs",
        "state",
        "AA",
        "alt_match_status",
        "header_marks_target_edit",
        "protein_id_used",
        "raw_candidate_ids",
        "retained_candidate_ids",
        "peptide_seq",
        "modified_seq",
        "precursor_charge",
        "precursor_id",
        "real_start",
        "real_end",
        "window",
        "fasta_edited_positions",
        "effective_edited_positions",
        "raw_fasta_edited_positions",
        "substring_offset_0based",
        "matrix_line",
        "fasta_header",
        "source_fasta",
        "coordinate_system",
        "id_layout",
    ] + sample_names

    with open_output(path) as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for row in rows:
            out = {k: row.get(k, "") for k in fields if k not in sample_names}
            ints = row.get("_intensities", [])
            for i, sample in enumerate(sample_names):
                out[sample] = format_float(ints[i] if i < len(ints) else 0.0)
            w.writerow(out)


def parse_int_list(x):
    return [int(v.strip()) for v in x.split(",") if v.strip()]


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(
        description=(
            "Quantify REF/multiple-ALT amino-acid states from DIA-NN pr_matrix "
            "with FBpp AA-coordinate FASTA IDs and an FBtr<->FBpp conversion table."
        )
    )
    p.add_argument("--sites", required=True, help="Editing/missense site table")
    p.add_argument("--matrix", required=True, help="DIA-NN report.pr_matrix.tsv")
    p.add_argument("--fasta", required=True, nargs="+", help="Search FASTA file(s)")
    p.add_argument(
        "--id-map", required=True, help="FlyBase table containing FBtr and FBpp columns"
    )
    p.add_argument("--out", required=True, help="State-level summary TSV")
    p.add_argument("--detail-out", required=True, help="Precursor-level detail TSV")
    p.add_argument(
        "--unlisted-alt-out",
        default="",
        help="Optional UNLISTED_ALT output; default derived from --out",
    )
    p.add_argument(
        "--candidate-id-cols",
        default="1,2,3,4",
        help="Fallback 1-based ID columns if Protein.Ids is absent/empty",
    )
    p.add_argument(
        "--peptide-seq-col",
        type=int,
        default=7,
        help="1-based Stripped.Sequence column; default 7",
    )
    p.add_argument(
        "--sample-start-col",
        type=int,
        default=11,
        help="1-based first sample intensity column; default 11",
    )
    args = p.parse_args()

    unlisted_out = args.unlisted_alt_out or derive_unlisted_alt_path(args.out)

    sites = load_sites(args.sites)
    site_groups = group_sites_by_amino_acid(sites)
    fbtr_to_fbpp = load_fbtr_fbpp_map(args.id_map)
    fasta_map = load_fasta_coord_map(args.fasta)

    sample_names, rows_by_fbpp, matrix_stats = load_matrix(
        matrix_file=args.matrix,
        fasta_coord_map=fasta_map,
        candidate_id_cols_1based=parse_int_list(args.candidate_id_cols),
        peptide_seq_col_1based=args.peptide_seq_col,
        sample_start_col_1based=args.sample_start_col,
    )

    main_rows, detail_rows, match_stats = summarize(
        site_groups=site_groups,
        sample_names=sample_names,
        rows_by_fbpp=rows_by_fbpp,
        fbtr_to_fbpp=fbtr_to_fbpp,
    )

    write_main(args.out, main_rows)
    write_detail(args.detail_out, detail_rows, sample_names)
    n_unlisted_rows = write_unlisted_alt(unlisted_out, main_rows)

    n_ref_states = len(
        {(r["site"], r["state"]) for r in main_rows if r["state"] == "REF"}
    )
    n_alt_states = len(
        {
            (r["site"], r["state"], r["AA"])
            for r in main_rows
            if str(r["state"]).startswith("ALT")
        }
    )
    n_unlisted_states = len(
        {
            (r["site"], r["state"], r["AA"])
            for r in main_rows
            if r["alt_match_status"] == "UNLISTED_ALT"
        }
    )

    print(f"Original site rows loaded: {len(sites)}")
    print(f"AA-position site groups after same-codon merge: {len(site_groups)}")
    print(f"Samples detected: {len(sample_names)}")
    print(f"FASTA files supplied: {len(args.fasta)}")
    print(f"FASTA ID keys loaded: {len(fasta_map)}")
    print(f"DIA-NN FBpp prefixes indexed: {len(rows_by_fbpp)}")
    print(f"Matrix rows total: {matrix_stats['matrix_rows_total']}")
    print(
        f"Matrix rows with candidate IDs: {matrix_stats['matrix_rows_with_candidate_ids']}"
    )
    print(
        "Matrix rows with retained-ID FASTA substring match: "
        f"{matrix_stats['matrix_rows_with_fasta_substring_candidates']}"
    )
    print(
        "Matrix rows whose first-per-FBpp IDs failed FASTA substring match: "
        f"{matrix_stats['matrix_rows_first_ids_failed_fasta_substring']}"
    )
    print(
        f"AA sites without FBtr->FBpp mapping: {match_stats['aa_sites_without_fbpp_map']}"
    )
    print(f"AA sites with FBtr->FBpp mapping: {match_stats['aa_sites_with_fbpp_map']}")
    print(
        f"AA sites matched to informative FBpp peptides: {match_stats['aa_sites_matched']}"
    )
    print(
        f"AA sites unmatched after FBpp mapping: {match_stats['aa_sites_unmatched_after_map']}"
    )
    print(
        "Informative precursor rows covering an editing AA position: "
        f"{match_stats['informative_precursor_rows']}"
    )
    print(f"REF site states: {n_ref_states}")
    print(f"Observed ALT site states: {n_alt_states}")
    print(f"Observed unlisted ALT site states: {n_unlisted_states}")
    print(f"State-level output: {args.out}")
    print(f"Detail output: {args.detail_out}")
    print(f"Unlisted-ALT output: {unlisted_out} (sample-level rows={n_unlisted_rows})")


if __name__ == "__main__":
    main()
