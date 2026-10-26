#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DIA_NN_proteoform_state_quant_v2.py

Purpose
-------
Summarize DIA-NN report.pr_matrix.tsv peptide intensities into amino-acid
states at RNA-editing-derived protein sites.

Core logic in v2
----------------
1. Coordinate-bearing FASTA IDs:
   - collect all trailing numeric fields;
   - the smallest value (one occurrence) is the window start;
   - the largest value (one occurrence) is the window end;
   - all remaining numeric values are edited positions;
   - for FBtr IDs these values are CDS nucleotide coordinates and are converted
     to 1-based amino-acid coordinates;
   - for non-FBtr IDs the values are treated as amino-acid coordinates.

2. DIA-NN Protein.Ids handling:
   - Protein.Ids is preferred over Protein.Group / other candidate columns;
   - within one matrix row, IDs are de-duplicated by transcript/prefix;
   - for each transcript/prefix, ONLY the first listed FASTA ID is retained;
   - the peptide sequence is still required to be an actual substring of that
     FASTA sequence.

3. Peptide inclusion and state classification:
   - only peptides whose mapped sequence actually covers the target amino-acid
     position are used;
   - the amino acid observed in Stripped.Sequence at the target position is read
     directly;
   - if observed AA == reference AA: REF;
   - otherwise: ALT1, ALT2, ... with the actual AA retained in the AA column;
   - ALT numbering is stable within a site: known ALT AAs from the site table
     come first (site-file order), followed by newly observed/unlisted ALT AAs;
   - an ALT AA not present among the known altAA values for that amino-acid site
     is marked UNLISTED_ALT and is also written to a separate output.

4. Same-codon multiple editing rows:
   - rows sharing FBtr + aa_pos + refAA are treated as ONE amino-acid site;
   - their source site IDs and known altAA values are merged;
   - this prevents duplicated protein intensity for two nucleotide edits in the
     same codon and allows a multi-edit-derived new AA to be recognized.

5. State-level peptide aggregation:
   - DIA-NN pr_matrix is precursor-level; the same Stripped.Sequence may occur
     at multiple charges/modification rows;
   - within site + AA state + sample, precursor intensities belonging to the same
     Stripped.Sequence are first SUMMED into one peptide-level intensity;
   - state sum_intensity = sum of positive peptide-level intensities;
   - state median_intensity = median of positive peptide-level intensities;
   - n_peptides = number of unique Stripped.Sequence peptides with intensity > 0.

Outputs
-------
--out:
  long-format site/state/sample summary with sum, median and n_peptides.

--detail-out:
  wide precursor-level evidence table for auditing peptide-to-FASTA mapping.

--unlisted-alt-out (optional):
  summary rows whose observed ALT AA was not listed in the original site table.
  If omitted, it is automatically derived from --out.
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
# Basic utilities
# -----------------------------------------------------------------------------


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def aa_to_one(x: str) -> str:
    x = str(x).strip()
    if len(x) == 1:
        return x.upper()
    if x in AA3_TO_AA1:
        return AA3_TO_AA1[x]
    raise ValueError(f"Unknown amino acid code: {x}")


def parse_float(x: str) -> float:
    if x is None:
        return 0.0
    x = str(x).strip()
    if x == "" or x.upper() in {"NA", "NAN", "NULL"}:
        return 0.0
    try:
        return float(x)
    except ValueError:
        return 0.0


def format_float(x: float) -> str:
    return f"{x:.10g}"


def normalize_seq(x: str) -> str:
    return str(x).strip().strip('"').strip("'").replace(" ", "").upper()


def split_auto(line: str):
    """Site table can be tab-delimited or whitespace-delimited."""
    line = line.rstrip("\n").rstrip("\r")
    if "\t" in line:
        return line.split("\t")
    return line.split()


def clean_colname(x: str) -> str:
    return x.strip().lstrip("\ufeff").strip('"').strip("'").rstrip("\r")


def normalize_match_key(x: str) -> str:
    """Normalize transcript/gene keys before dictionary lookup."""
    x = str(x).strip().strip('"').strip("'")
    if x.endswith("_DROME"):
        x = x[:-6]
    # Allow FBtr identifiers with a terminal version suffix.
    x = re.sub(r"^(FBtr\d+)\.\d+$", r"\1", x)
    return x


def cds_nt_to_aa_pos(pos_1based: int) -> int:
    """Convert a 1-based CDS nucleotide coordinate to a 1-based AA position."""
    return (int(pos_1based) - 1) // 3 + 1


def derive_unlisted_alt_path(out_path: str) -> str:
    if out_path.endswith(".tsv.gz"):
        return out_path[:-7] + ".unlisted_alt.tsv.gz"
    if out_path.endswith(".tsv"):
        return out_path[:-4] + ".unlisted_alt.tsv"
    return out_path + ".unlisted_alt.tsv"


# -----------------------------------------------------------------------------
# FASTA coordinate parsing
# -----------------------------------------------------------------------------


def parse_coord_id(x: str):
    """
    Parse a coordinate-bearing FASTA/DIA-NN ID.

    v2 coordinate rule:
      - collect all trailing numeric fields;
      - sort numerically;
      - remove ONE smallest value as window start;
      - remove ONE largest value as window end;
      - every remaining value is an edited position.

    This intentionally preserves a boundary edit when duplicate coordinates
    occur. Example:

      FBtrXXXX_100_100_200

    becomes:
      start = 100
      edits = [100]
      end   = 200

    Coordinate systems:
      FBtr prefix    -> numeric fields are CDS nt; converted to AA positions.
      non-FBtr       -> numeric fields are already AA positions.
    """
    original = str(x).strip().strip('"').strip("'")
    clean = original

    if clean.endswith("_DROME"):
        clean = clean[:-6]

    toks = clean.split("_")

    numeric = []
    i = len(toks) - 1
    while i >= 0 and toks[i].isdigit():
        numeric.append(int(toks[i]))
        i -= 1

    numeric = list(reversed(numeric))

    if len(numeric) < 2:
        raise ValueError(f"Cannot parse coordinate ID: {original}")

    prefix = "_".join(toks[: len(toks) - len(numeric)])
    prefix_norm = normalize_match_key(prefix)

    sorted_numeric = sorted(numeric)
    raw_start = sorted_numeric[0]
    raw_end = sorted_numeric[-1]
    # Slice removes exactly one min and exactly one max occurrence.
    raw_edits = sorted_numeric[1:-1]

    if re.match(r"^FBtr\d+$", prefix_norm):
        start = cds_nt_to_aa_pos(raw_start)
        end = cds_nt_to_aa_pos(raw_end)
        edited_positions = [cds_nt_to_aa_pos(p) for p in raw_edits]
        coordinate_system = "FBtr_CDS_nt_to_AA"
    else:
        start = raw_start
        end = raw_end
        edited_positions = list(raw_edits)
        coordinate_system = "protein_AA"

    if start > end:
        raise ValueError(
            f"Parsed start is greater than end for {original}: {start}>{end}"
        )

    return {
        "prefix": prefix_norm,
        "start": start,
        "end": end,
        "edited_positions": edited_positions,
        "raw_numeric": numeric,
        "raw_numeric_sorted": sorted_numeric,
        "raw_start": raw_start,
        "raw_end": raw_end,
        "raw_edited_positions": raw_edits,
        "coordinate_system": coordinate_system,
        "id_layout": "min_start_max_end_remaining_edits",
    }


def find_all_offsets(haystack: str, needle: str):
    """Return all 0-based offsets of needle in haystack."""
    offsets = []
    start = 0
    while True:
        idx = haystack.find(needle, start)
        if idx == -1:
            break
        offsets.append(idx)
        start = idx + 1
    return offsets


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


def extract_fasta_ids(header: str):
    """
    Extract possible DIA-NN IDs from UniProt-like or plain FASTA headers.
    """
    ids = []

    m = re.match(r"^(sp|tr)\|([^|]+)\|([^ ]+)\s*(.*)$", header)
    if m:
        acc = m.group(2).strip()
        entry = m.group(3).strip()
        rest = m.group(4)

        ids.append(acc)
        ids.append(entry)

        if entry.endswith("_DROME"):
            ids.append(entry[:-6])

        gm = re.search(r"\bGN=([^ ]+)", rest)
        if gm:
            ids.append(gm.group(1).strip())
    else:
        first = header.split()[0].strip()
        ids.append(first)

    out = []
    seen = set()
    for x in ids:
        if x and x not in seen:
            out.append(x)
            seen.add(x)
    return out


def extract_gene_prefix_from_header(header: str):
    """Prefer GN=<coordinate-bearing-name> for a gene-level lookup prefix."""
    gm = re.search(r"\bGN=([^ ]+)", header)
    if not gm:
        return None

    gn = gm.group(1).strip()
    try:
        return parse_coord_id(gn)["prefix"]
    except ValueError:
        return None


def load_fasta_coord_map(fasta_files):
    """Merge one or more search FASTAs into an ID -> candidate-record map."""
    if isinstance(fasta_files, (str, bytes)):
        fasta_files = [fasta_files]

    expanded_files = []
    for item in fasta_files:
        for path in str(item).split(","):
            path = path.strip()
            if path:
                expanded_files.append(path)

    if not expanded_files:
        raise ValueError("No FASTA file was supplied")

    id_to_records = defaultdict(list)
    seen_by_id = defaultdict(set)
    total_entries = 0
    total_added_records = 0

    for fasta_file in expanded_files:
        if not os.path.isfile(fasta_file):
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

        file_entries = 0
        file_added_records = 0

        for header, seq in read_fasta(fasta_file):
            file_entries += 1
            total_entries += 1

            ids = extract_fasta_ids(header)
            gene_prefix_from_gn = extract_gene_prefix_from_header(header)

            for one_id in ids:
                try:
                    parsed = parse_coord_id(one_id)
                except ValueError:
                    continue

                gene_prefix = (
                    gene_prefix_from_gn if gene_prefix_from_gn else parsed["prefix"]
                )

                record_key = (
                    one_id,
                    seq,
                    parsed["start"],
                    parsed["end"],
                    tuple(parsed["edited_positions"]),
                    parsed["coordinate_system"],
                    gene_prefix,
                )

                if record_key in seen_by_id[one_id]:
                    continue
                seen_by_id[one_id].add(record_key)

                rec = {
                    "fasta_id": one_id,
                    "fasta_seq": seq,
                    "fasta_start": parsed["start"],
                    "fasta_end": parsed["end"],
                    "fasta_edited_positions": parsed["edited_positions"],
                    "raw_fasta_start": parsed["raw_start"],
                    "raw_fasta_end": parsed["raw_end"],
                    "raw_fasta_edited_positions": parsed["raw_edited_positions"],
                    "coordinate_system": parsed["coordinate_system"],
                    "id_layout": parsed["id_layout"],
                    "gene_prefix": gene_prefix,
                    "fasta_header": header,
                    "source_fasta": fasta_file,
                }

                id_to_records[one_id].append(rec)
                file_added_records += 1
                total_added_records += 1

        print(
            f"Loaded FASTA: {fasta_file} | "
            f"entries={file_entries} | indexed_records={file_added_records}"
        )

    print(
        f"Merged FASTA files: {len(expanded_files)} | "
        f"entries={total_entries} | "
        f"unique ID keys={len(id_to_records)} | "
        f"indexed_records={total_added_records}"
    )

    return id_to_records


# -----------------------------------------------------------------------------
# Site-table loading and same-AA-position grouping
# -----------------------------------------------------------------------------


def load_sites(site_file: str):
    """
    Read editing site table.

    Required columns:
      site_id gene FBgn FBtr refAA altAA trans_pos trans_len
      cds_pos cds_len aa_pos aa_len
    """
    rows = []

    with open_text(site_file) as f:
        header_line = f.readline()
        header = [clean_colname(x) for x in split_auto(header_line)]

        missing = [col for col in SITE_COLUMNS if col not in header]
        if missing:
            raise ValueError(
                "Missing required columns in site file: "
                + ",".join(missing)
                + "\nDetected columns are: "
                + ",".join(header)
            )

        for line_no, line in enumerate(f, start=2):
            if not line.strip():
                continue

            values = split_auto(line)
            if len(values) < len(header):
                values = values + [""] * (len(header) - len(values))

            row_all = dict(zip(header, values))

            row = {}
            for k in SITE_COLUMNS:
                row[k] = row_all.get(k, "").strip().strip('"').rstrip("\r")

            row["aa_pos"] = int(row["aa_pos"])
            row["_ref1"] = aa_to_one(row["refAA"])
            row["_alt1"] = aa_to_one(row["altAA"])
            row["_line_no"] = line_no
            rows.append(row)

    return rows


def group_sites_by_amino_acid(sites):
    """
    Collapse multiple nucleotide editing rows affecting the same amino-acid
    position into one protein-state site.

    Primary grouping key:
      normalized FBtr + aa_pos + refAA

    If FBtr is unavailable, gene is used as the identity key.
    """
    groups = OrderedDict()

    for row in sites:
        fbtr = normalize_match_key(row.get("FBtr", ""))
        gene = normalize_match_key(row.get("gene", ""))
        identity = fbtr if fbtr else gene
        key = (identity, int(row["aa_pos"]), row["_ref1"])

        if key not in groups:
            groups[key] = {
                "site": f"{identity}:{row['aa_pos']}",
                "gene": row.get("gene", ""),
                "FBgn": row.get("FBgn", ""),
                "FBtr": row.get("FBtr", ""),
                "aa_pos": int(row["aa_pos"]),
                "refAA": row.get("refAA", ""),
                "refAA1": row["_ref1"],
                "aa_len": row.get("aa_len", ""),
                "source_site_ids": [],
                "source_cds_pos": [],
                "source_trans_pos": [],
                "known_alt_AAs": [],
                "known_alt_AAs_3letter": [],
                "_source_rows": [],
            }

        g = groups[key]
        g["_source_rows"].append(row)

        sid = row.get("site_id", "")
        if sid and sid not in g["source_site_ids"]:
            g["source_site_ids"].append(sid)

        cds = row.get("cds_pos", "")
        if cds and cds not in g["source_cds_pos"]:
            g["source_cds_pos"].append(cds)

        tpos = row.get("trans_pos", "")
        if tpos and tpos not in g["source_trans_pos"]:
            g["source_trans_pos"].append(tpos)

        alt1 = row["_alt1"]
        if alt1 not in g["known_alt_AAs"]:
            g["known_alt_AAs"].append(alt1)
            g["known_alt_AAs_3letter"].append(row.get("altAA", ""))

    return list(groups.values())


# -----------------------------------------------------------------------------
# DIA-NN matrix parsing
# -----------------------------------------------------------------------------


def parse_col_list(x: str):
    """Split a DIA-NN semicolon-separated ID field."""
    if x is None:
        return []

    x = str(x).strip().strip('"')
    if not x:
        return []

    out = []
    for item in x.split(";"):
        item = item.strip().strip('"').strip("'")
        if item:
            out.append(item)
    return out


def first_id_per_prefix(raw_ids):
    """
    Within one matrix row, retain only the first coordinate ID for each
    transcript/gene prefix.

    Example:
      FBtr0344440_1411_1562_1713;FBtr0344440_1411_1713

    -> retain only:
      FBtr0344440_1411_1562_1713
    """
    kept = OrderedDict()

    for rid in raw_ids:
        variants = [rid]
        if rid.endswith("_DROME"):
            variants.append(rid[:-6])

        parsed = None
        parsed_id = None
        for v in variants:
            try:
                parsed = parse_coord_id(v)
                parsed_id = v
                break
            except ValueError:
                continue

        if parsed is None:
            continue

        prefix = normalize_match_key(parsed["prefix"])
        if prefix not in kept:
            # Keep the original ID spelling first; parsed_id is only a fallback.
            kept[prefix] = rid if rid else parsed_id

    return list(kept.values())


def make_substring_candidates(rec, peptide_seq):
    """
    Map one Stripped.Sequence back to one FASTA record by substring matching.
    """
    fasta_seq = rec["fasta_seq"]
    peptide_seq = normalize_seq(peptide_seq)

    offsets = find_all_offsets(fasta_seq, peptide_seq)

    candidates = []
    for offset in offsets:
        real_start = rec["fasta_start"] + offset
        real_end = real_start + len(peptide_seq) - 1

        effective_edits = [
            p for p in rec["fasta_edited_positions"] if real_start <= p <= real_end
        ]

        candidates.append(
            {
                "candidate_id": rec["fasta_id"],
                "gene_prefix": rec["gene_prefix"],
                "peptide_seq": peptide_seq,
                "real_start": real_start,
                "real_end": real_end,
                "substring_offset_0based": offset,
                "fasta_start": rec["fasta_start"],
                "fasta_end": rec["fasta_end"],
                "fasta_edited_positions": rec["fasta_edited_positions"],
                "effective_edited_positions": effective_edits,
                "raw_fasta_start": rec.get("raw_fasta_start", ""),
                "raw_fasta_end": rec.get("raw_fasta_end", ""),
                "raw_fasta_edited_positions": rec.get("raw_fasta_edited_positions", []),
                "coordinate_system": rec.get("coordinate_system", ""),
                "id_layout": rec.get("id_layout", ""),
                "fasta_header": rec["fasta_header"],
                "source_fasta": rec.get("source_fasta", ""),
            }
        )

    return candidates


def load_matrix(
    matrix_file: str,
    fasta_coord_map,
    candidate_id_cols_1based,
    peptide_seq_col_1based,
    sample_start_col_1based,
):
    """
    Read DIA-NN report.pr_matrix.tsv.

    Important v2 behavior:
      - if a Protein.Ids column exists, it is used FIRST and exclusively when
        non-empty;
      - only if Protein.Ids is empty do we fall back to --candidate-id-cols;
      - within the selected ID list, one first FASTA ID is kept per transcript /
        gene prefix;
      - peptide sequence must be found in that exact retained FASTA entry.
    """
    rows_by_gene = defaultdict(list)

    stats = {
        "matrix_rows_total": 0,
        "matrix_rows_with_candidate_ids": 0,
        "matrix_rows_with_fasta_substring_candidates": 0,
        "matrix_rows_without_candidate_ids": 0,
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

        protein_ids_idx = None
        if "Protein.Ids" in header_clean:
            protein_ids_idx = header_clean.index("Protein.Ids")

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
                row = row + [""] * (len(header) - len(row))

            peptide_seq = normalize_seq(row[seq_idx])
            if not peptide_seq:
                stats["matrix_rows_without_candidate_ids"] += 1
                continue

            # Prefer Protein.Ids because the user-defined FASTA IDs are carried
            # there. Protein.Group can collapse edited + unedited entries.
            raw_ids = []
            if protein_ids_idx is not None:
                raw_ids = parse_col_list(row[protein_ids_idx])

            if not raw_ids:
                for idx in candidate_idxs:
                    if 0 <= idx < len(row):
                        raw_ids.extend(parse_col_list(row[idx]))

            retained_ids = first_id_per_prefix(raw_ids)

            if not retained_ids:
                stats["matrix_rows_without_candidate_ids"] += 1
                continue

            stats["matrix_rows_with_candidate_ids"] += 1

            candidates = []

            for rid in retained_ids:
                lookup_ids = [rid]
                if rid.endswith("_DROME"):
                    lookup_ids.append(rid[:-6])

                # Use the retained ID only. Do not move to a second same-prefix
                # DIA-NN ID if this first ID fails substring matching.
                found_this_rid = False
                for lid in lookup_ids:
                    recs = fasta_coord_map.get(lid, [])
                    if not recs:
                        continue
                    for rec in recs:
                        mapped = make_substring_candidates(rec, peptide_seq)
                        if mapped:
                            candidates.extend(mapped)
                            found_this_rid = True
                    if found_this_rid:
                        break

            # De-duplicate identical mappings generated from duplicated FASTA
            # aliases/files while retaining all distinct substring offsets.
            unique_candidates = []
            seen_key = set()

            for cand in candidates:
                key = (
                    normalize_match_key(cand["gene_prefix"]),
                    cand["candidate_id"],
                    cand["real_start"],
                    cand["real_end"],
                    tuple(cand["fasta_edited_positions"]),
                    cand["peptide_seq"],
                )
                if key not in seen_key:
                    unique_candidates.append(cand)
                    seen_key.add(key)

            if not unique_candidates:
                stats["matrix_rows_first_ids_failed_fasta_substring"] += 1
                continue

            stats["matrix_rows_with_fasta_substring_candidates"] += 1

            intensities = [parse_float(x) for x in row[sample_start_idx:]]

            matrix_row = {
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

            # Add this matrix row once per candidate transcript/gene prefix.
            added_genes = set()
            for cand in unique_candidates:
                gp = normalize_match_key(cand.get("gene_prefix", ""))
                if gp and gp not in added_genes:
                    rows_by_gene[gp].append(matrix_row)
                    added_genes.add(gp)

    return sample_names, rows_by_gene, stats


# -----------------------------------------------------------------------------
# Site-state classification and aggregation
# -----------------------------------------------------------------------------


def classify_candidate_at_site(cand, site_group):
    """
    Classify an actually mapped peptide at one amino-acid site.

    A peptide is informative ONLY when it physically covers site_group['aa_pos'].
    Classification is determined by the actual residue in Stripped.Sequence:
      observed == refAA -> REF
      observed != refAA -> ALT (actual AA retained)
    """
    aa_pos = int(site_group["aa_pos"])

    if not (cand["real_start"] <= aa_pos <= cand["real_end"]):
        return None

    local_idx = aa_pos - cand["real_start"]
    peptide_seq = cand["peptide_seq"]

    if local_idx < 0 or local_idx >= len(peptide_seq):
        return None

    observed = peptide_seq[local_idx]
    ref1 = site_group["refAA1"]

    if observed == ref1:
        base_state = "REF"
        alt_match_status = "REF"
    else:
        base_state = "ALT"
        if observed in site_group["known_alt_AAs"]:
            alt_match_status = "KNOWN_ALT"
        else:
            alt_match_status = "UNLISTED_ALT"

    return {
        "base_state": base_state,
        "AA": observed,
        "alt_match_status": alt_match_status,
        "header_marks_target_edit": (
            "yes" if aa_pos in cand.get("fasta_edited_positions", []) else "no"
        ),
    }


def choose_first_valid_mapping(mrow, lookup_key, site_group):
    """
    For one DIA-NN precursor row and one site, use the first retained FASTA
    mapping that both belongs to lookup_key and covers the site.

    This follows the first-identifier-per-transcript rule and avoids counting
    one precursor intensity more than once for the same site.
    """
    lookup_key = normalize_match_key(lookup_key)

    for cand in mrow["candidates"]:
        candidate_key = normalize_match_key(cand.get("gene_prefix", ""))
        if candidate_key != lookup_key:
            continue

        cls = classify_candidate_at_site(cand, site_group)
        if cls is None:
            continue

        return cand, cls

    return None, None


def summarize(site_groups, sample_names, rows_by_gene):
    """
    Produce state-level summary plus precursor-level audit details.

    Matching order is preserved from the original script:
      1. FBtr
      2. gene fallback only if no informative FBtr peptide exists
    """
    main_rows = []
    detail_rows = []

    match_stats = {
        "aa_sites_matched_by_FBtr": 0,
        "aa_sites_matched_by_gene_fallback": 0,
        "aa_sites_unmatched": 0,
        "informative_precursor_rows": 0,
    }

    for site_group in site_groups:
        fbtr_key = normalize_match_key(site_group.get("FBtr", ""))
        gene_key = normalize_match_key(site_group.get("gene", ""))

        route_candidates = []
        if fbtr_key:
            route_candidates.append(("FBtr", fbtr_key))
        if gene_key and gene_key != fbtr_key:
            route_candidates.append(("gene_fallback", gene_key))

        chosen_mode = None
        chosen_key = None
        chosen_evidence = []

        for match_mode, lookup_key in route_candidates:
            candidate_matrix_rows = rows_by_gene.get(lookup_key, [])
            route_evidence = []

            for mrow in candidate_matrix_rows:
                cand, cls = choose_first_valid_mapping(
                    mrow=mrow,
                    lookup_key=lookup_key,
                    site_group=site_group,
                )

                if cand is None:
                    # Peptides that do not cover the target AA position are not
                    # included in REF/ALT and are not written as evidence.
                    continue

                route_evidence.append((mrow, cand, cls))

            if route_evidence:
                chosen_mode = match_mode
                chosen_key = lookup_key
                chosen_evidence = route_evidence
                break

        if not chosen_evidence:
            match_stats["aa_sites_unmatched"] += 1
            continue

        if chosen_mode == "FBtr":
            match_stats["aa_sites_matched_by_FBtr"] += 1
        else:
            match_stats["aa_sites_matched_by_gene_fallback"] += 1

        match_stats["informative_precursor_rows"] += len(chosen_evidence)

        # precursor -> peptide -> state aggregation
        # Key: (AA, sample_index, stripped_sequence) -> summed precursor intensity
        peptide_level = defaultdict(float)
        state_windows = defaultdict(set)  # AA -> AA-coordinate windows
        state_raw_windows = defaultdict(set)  # AA -> raw coordinate windows
        all_windows = set()
        all_raw_windows = set()
        observed_alt_aas = set()

        # Store detail rows first with state assigned later after ALT numbering.
        site_detail_indices = []

        for mrow, cand, cls in chosen_evidence:
            aa = cls["AA"]

            if cls["base_state"] == "ALT":
                observed_alt_aas.add(aa)

            win = f"{cand['fasta_start']}-{cand['fasta_end']}"
            raw_win = f"{cand['raw_fasta_start']}-{cand['raw_fasta_end']}"
            state_windows[aa].add(win)
            state_raw_windows[aa].add(raw_win)
            all_windows.add(win)
            all_raw_windows.add(raw_win)

            # Sum multiple precursor rows/charges/modifications belonging to the
            # same Stripped.Sequence before state-level median/n calculation.
            for i, val in enumerate(mrow["intensities"]):
                peptide_level[(aa, i, mrow["peptide_seq"])] += val

            detail = {
                "site": site_group["site"],
                "source_site_ids": ";".join(site_group["source_site_ids"]),
                "gene": site_group["gene"],
                "FBgn": site_group["FBgn"],
                "FBtr": site_group["FBtr"],
                "aa_pos": site_group["aa_pos"],
                "refAA": site_group["refAA1"],
                "known_alt_AAs": ";".join(site_group["known_alt_AAs"]),
                "match_mode": chosen_mode,
                "lookup_key": chosen_key,
                "state": "",  # filled after ALT numbering
                "AA": aa,
                "alt_match_status": cls["alt_match_status"],
                "header_marks_target_edit": cls["header_marks_target_edit"],
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
                "raw_window": raw_win,
                "fasta_edited_positions": ",".join(
                    map(str, cand.get("fasta_edited_positions", []))
                ),
                "raw_fasta_edited_positions": ",".join(
                    map(str, cand.get("raw_fasta_edited_positions", []))
                ),
                "substring_offset_0based": cand["substring_offset_0based"],
                "matrix_line": mrow["line_no"],
                "fasta_header": cand["fasta_header"],
                "source_fasta": cand.get("source_fasta", ""),
                "coordinate_system": cand.get("coordinate_system", ""),
                "id_layout": cand.get("id_layout", ""),
                "_intensities": mrow["intensities"],
            }
            detail_rows.append(detail)
            site_detail_indices.append(len(detail_rows) - 1)

        # ALT numbering: known table ALT AAs first in site-file order, then
        # newly observed/unlisted AAs in alphabetical order. Only observed ALT
        # AAs receive an ALT number.
        alt_order = []
        for aa in site_group["known_alt_AAs"]:
            if aa in observed_alt_aas and aa not in alt_order:
                alt_order.append(aa)
        for aa in sorted(observed_alt_aas):
            if aa not in alt_order:
                alt_order.append(aa)

        state_by_aa = {site_group["refAA1"]: "REF"}
        for idx, aa in enumerate(alt_order, start=1):
            state_by_aa[aa] = f"ALT{idx}"

        # Fill final state names into detail rows.
        for di in site_detail_indices:
            aa = detail_rows[di]["AA"]
            detail_rows[di]["state"] = state_by_aa.get(aa, "ALT")

        # Always output the single REF state for a matched AA site, even when
        # its intensity is zero in a particular sample. ALT states are the ALT
        # AAs actually observed in at least one informative precursor row.
        output_aas = [site_group["refAA1"]] + alt_order

        for aa in output_aas:
            state = state_by_aa[aa]
            if aa == site_group["refAA1"]:
                alt_match_status = "REF"
            elif aa in site_group["known_alt_AAs"]:
                alt_match_status = "KNOWN_ALT"
            else:
                alt_match_status = "UNLISTED_ALT"

            windows = state_windows.get(aa, set()) or all_windows
            raw_windows = state_raw_windows.get(aa, set()) or all_raw_windows

            for i, sample in enumerate(sample_names):
                vals = []
                for (
                    aa_key,
                    sample_idx,
                    peptide,
                ), pep_intensity in peptide_level.items():
                    if aa_key == aa and sample_idx == i and pep_intensity > 0:
                        vals.append(pep_intensity)

                sum_intensity = sum(vals) if vals else 0.0
                median_intensity = median(vals) if vals else 0.0
                n_peptides = len(vals)

                main_rows.append(
                    {
                        "site": site_group["site"],
                        "source_site_ids": ";".join(site_group["source_site_ids"]),
                        "source_cds_pos": ";".join(site_group["source_cds_pos"]),
                        "gene": site_group["gene"],
                        "FBgn": site_group["FBgn"],
                        "FBtr": site_group["FBtr"],
                        "aa_pos": site_group["aa_pos"],
                        "refAA": site_group["refAA1"],
                        "known_alt_AAs": ";".join(site_group["known_alt_AAs"]),
                        "window": ";".join(sorted(windows)),
                        "raw_window": ";".join(sorted(raw_windows)),
                        "match_mode": chosen_mode,
                        "state": state,
                        "AA": aa,
                        "alt_match_status": alt_match_status,
                        "sample": sample,
                        "sum_intensity": sum_intensity,
                        "median_intensity": median_intensity,
                        "n_peptides": n_peptides,
                    }
                )

    return main_rows, detail_rows, match_stats


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------


def open_output(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return open(path, "w", encoding="utf-8", newline="")


def write_main(out_file, rows):
    fieldnames = [
        "site",
        "source_site_ids",
        "source_cds_pos",
        "gene",
        "FBgn",
        "FBtr",
        "aa_pos",
        "refAA",
        "known_alt_AAs",
        "window",
        "raw_window",
        "match_mode",
        "state",
        "AA",
        "alt_match_status",
        "sample",
        "sum_intensity",
        "median_intensity",
        "n_peptides",
    ]

    with open_output(out_file) as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            out = dict(row)
            out["sum_intensity"] = format_float(out["sum_intensity"])
            out["median_intensity"] = format_float(out["median_intensity"])
            writer.writerow({k: out.get(k, "") for k in fieldnames})


def write_unlisted_alt(out_file, rows):
    subset = [r for r in rows if r.get("alt_match_status") == "UNLISTED_ALT"]
    write_main(out_file, subset)
    return len(subset)


def write_detail(out_file, rows, sample_names):
    fieldnames = [
        "site",
        "source_site_ids",
        "gene",
        "FBgn",
        "FBtr",
        "aa_pos",
        "refAA",
        "known_alt_AAs",
        "match_mode",
        "lookup_key",
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
        "raw_window",
        "fasta_edited_positions",
        "raw_fasta_edited_positions",
        "substring_offset_0based",
        "matrix_line",
        "fasta_header",
        "source_fasta",
        "coordinate_system",
        "id_layout",
    ] + sample_names

    with open_output(out_file) as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            out = {k: row.get(k, "") for k in fieldnames if k not in sample_names}
            intensities = row.get("_intensities", [])
            for i, sample in enumerate(sample_names):
                val = intensities[i] if i < len(intensities) else 0.0
                out[sample] = format_float(val)
            writer.writerow(out)


def parse_int_list(x: str):
    return [int(i.strip()) for i in x.split(",") if i.strip()]


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Quantify REF and multiple ALT amino-acid states from a DIA-NN "
            "report.pr_matrix.tsv by peptide-to-FASTA substring mapping."
        )
    )

    parser.add_argument("--sites", required=True, help="Editing/missense site table")
    parser.add_argument("--matrix", required=True, help="DIA-NN report.pr_matrix.tsv")
    parser.add_argument(
        "--fasta",
        required=True,
        nargs="+",
        help=(
            "One or more FASTA files used for the DIA-NN search. "
            "Comma-separated paths are also accepted."
        ),
    )
    parser.add_argument("--out", required=True, help="State-level summary TSV")
    parser.add_argument(
        "--detail-out", required=True, help="Precursor-level detail TSV"
    )
    parser.add_argument(
        "--unlisted-alt-out",
        default="",
        help=(
            "Optional output for ALT states whose AA is absent from the original "
            "site table. Default: derive <out>.unlisted_alt.tsv"
        ),
    )

    # Kept for backward compatibility with the original script. Protein.Ids is
    # automatically preferred when that named column exists.
    parser.add_argument(
        "--candidate-id-cols",
        default="1,2,3,4",
        help=(
            "Fallback comma-separated 1-based ID columns if Protein.Ids is empty. "
            "Default: 1,2,3,4"
        ),
    )
    parser.add_argument(
        "--peptide-seq-col",
        type=int,
        default=7,
        help="1-based Stripped.Sequence column. Default: 7",
    )
    parser.add_argument(
        "--sample-start-col",
        type=int,
        default=11,
        help="1-based first sample intensity column. Default: 11",
    )

    args = parser.parse_args()

    unlisted_alt_out = (
        args.unlisted_alt_out
        if args.unlisted_alt_out
        else derive_unlisted_alt_path(args.out)
    )

    sites = load_sites(args.sites)
    site_groups = group_sites_by_amino_acid(sites)
    fasta_coord_map = load_fasta_coord_map(args.fasta)

    sample_names, rows_by_gene, stats = load_matrix(
        matrix_file=args.matrix,
        fasta_coord_map=fasta_coord_map,
        candidate_id_cols_1based=parse_int_list(args.candidate_id_cols),
        peptide_seq_col_1based=args.peptide_seq_col,
        sample_start_col_1based=args.sample_start_col,
    )

    n_fbtr_index_keys = sum(
        1 for key in rows_by_gene if re.match(r"^FBtr\d+$", normalize_match_key(key))
    )
    n_gene_index_keys = len(rows_by_gene) - n_fbtr_index_keys

    if n_fbtr_index_keys > 0 and n_gene_index_keys == 0:
        detected_mapping = "FBtr"
    elif n_gene_index_keys > 0 and n_fbtr_index_keys == 0:
        detected_mapping = "gene"
    elif n_fbtr_index_keys > 0 and n_gene_index_keys > 0:
        detected_mapping = "mixed FBtr/gene"
    else:
        detected_mapping = "undetermined"

    print(
        "Detected matrix mapping: "
        f"{detected_mapping} "
        f"(FBtr keys={n_fbtr_index_keys}, gene keys={n_gene_index_keys})"
    )

    main_rows, detail_rows, match_stats = summarize(
        site_groups=site_groups,
        sample_names=sample_names,
        rows_by_gene=rows_by_gene,
    )

    write_main(args.out, main_rows)
    write_detail(args.detail_out, detail_rows, sample_names)
    n_unlisted_rows = write_unlisted_alt(unlisted_alt_out, main_rows)

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
    print(f"FASTA ID keys loaded after merge: {len(fasta_coord_map)}")
    print(f"Matrix rows total: {stats['matrix_rows_total']}")
    print(f"Matrix rows with candidate IDs: {stats['matrix_rows_with_candidate_ids']}")
    print(
        "Matrix rows with retained-ID FASTA substring match: "
        f"{stats['matrix_rows_with_fasta_substring_candidates']}"
    )
    print(
        "Matrix rows whose first-per-transcript IDs failed FASTA substring match: "
        f"{stats['matrix_rows_first_ids_failed_fasta_substring']}"
    )
    print(f"AA sites matched by FBtr: {match_stats['aa_sites_matched_by_FBtr']}")
    print(
        "AA sites matched by gene fallback: "
        f"{match_stats['aa_sites_matched_by_gene_fallback']}"
    )
    print(f"AA sites unmatched: {match_stats['aa_sites_unmatched']}")
    print(
        "Informative precursor rows covering an editing AA position: "
        f"{match_stats['informative_precursor_rows']}"
    )
    print(f"REF site states: {n_ref_states}")
    print(f"Observed ALT site states: {n_alt_states}")
    print(f"Observed unlisted ALT site states: {n_unlisted_states}")
    print(f"State-level output: {args.out}")
    print(f"Detail output: {args.detail_out}")
    print(
        f"Unlisted-ALT output: {unlisted_alt_out} "
        f"(sample-level rows={n_unlisted_rows})"
    )


if __name__ == "__main__":
    main()
