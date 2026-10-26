import sys
import pysam


def load_txinfo(txinfo_file):
    """
    Read the required transcript metadata columns:
      tx_name, utr5_len, cds_len
    Return:
      txinfo[tx_name] = (utr5_len, cds_len)
    """
    txinfo = {}
    with open(txinfo_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {k: i for i, k in enumerate(header)}

        required = ["tx_name", "utr5_len", "cds_len"]
        for r in required:
            if r not in col:
                raise ValueError(f"Missing required column in txinfo: {r}")

        tx_i = col["tx_name"]
        utr5_i = col["utr5_len"]
        cds_i = col["cds_len"]

        for line in f:
            if not line.strip():
                continue
            arr = line.rstrip("\n").split("\t")
            tx = arr[tx_i]
            utr5 = int(arr[utr5_i])
            cds_len = int(arr[cds_i])
            txinfo[tx] = (utr5, cds_len)
    return txinfo


def fmt_cds_pos(pos, cds_len):
    """
    Format 1-based CDS-relative coordinates:
      Within CDS: 1,2,3...
      5-prime UTR: negative values or zero
      3'UTR: +1,+2,+3...
    """
    if pos < 1:
        return str(pos)
    elif pos > cds_len:
        return f"+{pos - cds_len}"
    else:
        return str(pos)


def build_qalign_to_ref_map(aln, qseq_len):
    """
    Map aligned-query indices to 1-based transcript coordinates.
    Handle insertions, deletions and clipping through aligned pairs.

    Coordinate conventions:
    - query_alignment_sequence excludes soft-clipped bases.
    - get_aligned_pairs() uses full query_sequence indices.
    - Subtract the alignment start to obtain aligned-query indices.
    - Insertions have no reference position; deletions/skips have no query position.
    """
    q2r = [None] * qseq_len

    q_align_start = aln.query_alignment_start
    q_align_end = aln.query_alignment_end

    for qpos, rpos in aln.get_aligned_pairs(matches_only=False, with_seq=False):
        if qpos is None:
            continue

        # Exclude soft-clipped query positions.
        if qpos < q_align_start or qpos >= q_align_end:
            continue

        qidx = qpos - q_align_start
        if 0 <= qidx < qseq_len and rpos is not None:
            q2r[qidx] = rpos + 1  # Convert to 1-based transcript coordinates.

    return q2r


def main():
    if len(sys.argv) < 4 or len(sys.argv) > 6:
        sys.stderr.write(
            f"Usage: {sys.argv[0]} <tagged.bam> <txinfo.tsv> <out.tsv> [bam_threads=4] [flush_n=10000]\n"
        )
        sys.exit(1)

    bam_file = sys.argv[1]
    txinfo_file = sys.argv[2]
    out_file = sys.argv[3]
    bam_threads = int(sys.argv[4]) if len(sys.argv) >= 5 else 4
    flush_n = int(sys.argv[5]) if len(sys.argv) >= 6 else 10000

    txinfo = load_txinfo(txinfo_file)

    buffer = []
    n_total = 0
    n_written = 0
    n_skipped_no_tag = 0
    n_skipped_no_tx = 0
    n_skipped_no_cds = 0
    n_skipped_no_seq = 0
    n_skipped_bad_frame = 0
    n_skipped_oob = 0
    n_skipped_noncontig = 0

    with pysam.AlignmentFile(bam_file, "rb", threads=bam_threads) as bam, open(
        out_file, "w"
    ) as out:
        out.write("transcript_id\tstart\tend\tcodon\tsite_type\tread_name\n")

        for aln in bam:
            n_total += 1

            # Forward is relative to the transcript reference in this BAM.
            if aln.is_unmapped or aln.is_reverse:
                continue
            if not aln.has_tag("PS") or not aln.has_tag("XF"):
                n_skipped_no_tag += 1
                continue

            tx = aln.reference_name
            info = txinfo.get(tx)
            if info is None:
                n_skipped_no_tx += 1
                continue

            utr5_len, cds_len = info
            if cds_len <= 0:
                n_skipped_no_cds += 1
                continue

            qseq = aln.query_alignment_sequence
            if qseq is None:
                n_skipped_no_seq += 1
                continue

            qlen = len(qseq)
            if qlen == 0:
                n_skipped_no_seq += 1
                continue

            ps = int(aln.get_tag("PS"))
            frame = int(aln.get_tag("XF"))

            if frame not in (0, 1, 2):
                n_skipped_bad_frame += 1
                continue

            # Locate the P-site in the aligned query.
            # BAM stores reverse-aligned sequences in reference orientation.
            # Convert indices without reverse-complementing the stored sequence.
            if aln.is_reverse:
                p_idx = qlen - 1 - ps
            else:
                p_idx = ps

            # Map aligned-query indices to transcript coordinates.
            q2r = build_qalign_to_ref_map(aln, qlen)

            for site_type, shift in (("P", 0), ("A", 3), ("E", -3)):
                site_idx = p_idx + shift

                # Adjust to the codon start using the predicted frame.
                if frame == 0:
                    codon_start_idx = site_idx
                elif frame == 1:
                    codon_start_idx = site_idx - 1
                else:  # frame == 2
                    codon_start_idx = site_idx + 1

                codon_end_idx = codon_start_idx + 2

                # Exclude codons outside the read.
                if codon_start_idx < 0 or codon_end_idx >= qlen:
                    n_skipped_oob += 1
                    continue

                codon = qseq[codon_start_idx : codon_start_idx + 3]
                if len(codon) != 3:
                    n_skipped_oob += 1
                    continue

                # Require three consecutive mapped transcript bases.
                r1 = q2r[codon_start_idx]
                r2 = q2r[codon_start_idx + 1]
                r3 = q2r[codon_start_idx + 2]

                if r1 is None or r2 is None or r3 is None:
                    n_skipped_noncontig += 1
                    continue

                if not (r2 == r1 + 1 and r3 == r2 + 1):
                    n_skipped_noncontig += 1
                    continue

                tx_start = r1
                tx_end = r3

                # Convert to CDS-relative coordinates.
                cds_start = tx_start - utr5_len
                cds_end = tx_end - utr5_len

                line = (
                    f"{tx}\t"
                    f"{fmt_cds_pos(cds_start, cds_len)}\t"
                    f"{fmt_cds_pos(cds_end, cds_len)}\t"
                    f"{codon}\t"
                    f"{site_type}\t"
                    f"{aln.query_name}\n"
                )
                buffer.append(line)
                n_written += 1

                if len(buffer) >= flush_n:
                    out.writelines(buffer)
                    buffer.clear()

        if buffer:
            out.writelines(buffer)

    sys.stderr.write(
        f"[DONE] total_reads={n_total}, written_lines={n_written}\n"
        f"[SKIP] no_tag={n_skipped_no_tag}, no_txinfo={n_skipped_no_tx}, no_cds={n_skipped_no_cds}, "
        f"no_seq={n_skipped_no_seq}, bad_frame={n_skipped_bad_frame}, "
        f"out_of_bounds={n_skipped_oob}, noncontig_or_indel={n_skipped_noncontig}\n"
    )


if __name__ == "__main__":
    main()
