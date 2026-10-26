#!/usr/bin/env python3
"""Map CDS sites and calculate focal/flanking codon coverage.

The implementation is designed for large, transcript-grouped files.  Count files
and the longest-transcript annotation are streamed one transcript at a time, so
memory use depends on the longest single CDS rather than on the whole library.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO

ANNOTATION_HEADER = (
    "chrom",
    "genomic_codon_first_pos_1based",
    "transcript",
    "cds_block_index",
    "strand",
    "cds_codon_first_pos_1based",
    "codon_index_1based",
)
COUNT_HEADER = ANNOTATION_HEADER + (
    "rna_count",
    "rpf_count",
    "flank_rna_count",
    "flank_rpf_count",
    "flank_n",
    "library_RNA_sum",
    "library_RPF_sum",
)
PROGRESS_EVERY = 5_000_000


def log(message: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}] {message}", file=sys.stderr, flush=True)


def open_text(path: Path, mode: str = "rt") -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="", buffering=1024 * 1024)


def tmp_path(path: Path) -> Path:
    return path.with_name(f".{path.name}.tmp.{os.getpid()}")


def check_output_paths(paths: list[Path], overwrite: bool) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing and not overwrite:
        joined = "\n  ".join(existing)
        raise FileExistsError(
            f"Output already exists (use --overwrite to replace):\n  {joined}"
        )


def replace_outputs(temp_to_final: list[tuple[Path, Path]]) -> None:
    for temporary, final in temp_to_final:
        os.replace(temporary, final)


def cleanup_temps(paths: list[Path]) -> None:
    for path in paths:
        try:
            path.unlink()
        except FileNotFoundError:
            pass


def load_longest_transcripts(path: Path) -> tuple[set[bytes], dict[bytes, int]]:
    transcripts: set[bytes] = set()
    rank: dict[bytes, int] = {}
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            transcript = raw.strip().split(b"\t", 1)[0]
            if not transcript:
                continue
            if transcript not in transcripts:
                rank[transcript] = len(rank)
                transcripts.add(transcript)
    if not transcripts:
        raise ValueError(f"No transcript IDs found in {path}")
    return transcripts, rank


@dataclass(frozen=True)
class M6ARecord:
    site_key: str
    chrom: str
    genomic_pos: int


def load_m6a_sites(
    path: Path,
) -> tuple[list[M6ARecord], dict[tuple[bytes, int], list[int]]]:
    records: list[M6ARecord] = []
    lookup: dict[tuple[bytes, int], list[int]] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            site_key = line.rstrip("\r\n").split("\t", 1)[0]
            try:
                chrom, position_text = site_key.rsplit("_", 1)
                position = int(position_text)
            except (ValueError, IndexError):
                if line_number == 1:
                    continue
                raise ValueError(
                    f"Invalid m6A site key at {path}:{line_number}: {site_key!r}"
                )
            index = len(records)
            records.append(M6ARecord(site_key, chrom, position))
            lookup.setdefault((chrom.encode(), position), []).append(index)
    if not records:
        raise ValueError(f"No m6A sites found in {path}")
    return records, lookup


@dataclass(frozen=True)
class M6ACandidate:
    transcript_rank: int
    transcript: str
    strand: str
    cds_pos: int
    codon_first_cds_pos: int
    genomic_codon_first_pos: int
    codon_index: int
    codon_offset: int


def prepare_annotation(
    annotation: Path,
    longest_transcripts: Path,
    m6a_sites: Path,
    output_annotation: Path,
    output_m6a_map: Path,
    overwrite: bool = False,
) -> dict[str, int]:
    """Filter annotation and select one longest-transcript mapping per m6A site."""
    check_output_paths([output_annotation, output_m6a_map], overwrite)
    output_annotation.parent.mkdir(parents=True, exist_ok=True)
    output_m6a_map.parent.mkdir(parents=True, exist_ok=True)

    longest, transcript_rank = load_longest_transcripts(longest_transcripts)
    m6a_records, m6a_lookup = load_m6a_sites(m6a_sites)
    selected: list[M6ACandidate | None] = [None] * len(m6a_records)
    candidate_counts = [0] * len(m6a_records)

    annotation_tmp = tmp_path(output_annotation)
    mapping_tmp = tmp_path(output_m6a_map)
    cleanup_temps([annotation_tmp, mapping_tmp])

    input_rows = 0
    kept_rows = 0
    try:
        with annotation.open(
            "rb", buffering=1024 * 1024
        ) as source, annotation_tmp.open("wb", buffering=1024 * 1024) as destination:
            for input_rows, raw in enumerate(source, 1):
                fields = raw.rstrip(b"\r\n").split(b"\t")
                if len(fields) < 7:
                    raise ValueError(
                        f"Expected at least 7 columns at {annotation}:{input_rows}"
                    )
                transcript_bytes = fields[2]
                if transcript_bytes not in longest:
                    if input_rows % PROGRESS_EVERY == 0:
                        log(
                            f"prepare: scanned {input_rows:,} annotation rows; "
                            f"kept {kept_rows:,}"
                        )
                    continue

                destination.write(raw if raw.endswith(b"\n") else raw + b"\n")
                kept_rows += 1
                try:
                    genomic_codon_start = int(fields[1])
                    cds_codon_start = int(fields[5])
                    codon_index = int(fields[6])
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid numeric annotation field at {annotation}:{input_rows}"
                    ) from exc
                strand = fields[4]
                if strand == b"+":
                    direction = 1
                elif strand == b"-":
                    direction = -1
                else:
                    raise ValueError(
                        f"Invalid strand at {annotation}:{input_rows}: {strand!r}"
                    )

                for offset in range(3):
                    indices = m6a_lookup.get(
                        (fields[0], genomic_codon_start + direction * offset)
                    )
                    if not indices:
                        continue
                    candidate = M6ACandidate(
                        transcript_rank=transcript_rank[transcript_bytes],
                        transcript=transcript_bytes.decode(),
                        strand=strand.decode(),
                        cds_pos=cds_codon_start + offset,
                        codon_first_cds_pos=cds_codon_start,
                        genomic_codon_first_pos=genomic_codon_start,
                        codon_index=codon_index,
                        codon_offset=offset,
                    )
                    for index in indices:
                        candidate_counts[index] += 1
                        previous = selected[index]
                        if previous is None or (
                            candidate.transcript_rank,
                            candidate.codon_first_cds_pos,
                        ) < (
                            previous.transcript_rank,
                            previous.codon_first_cds_pos,
                        ):
                            selected[index] = candidate

                if input_rows % PROGRESS_EVERY == 0:
                    log(
                        f"prepare: scanned {input_rows:,} annotation rows; "
                        f"kept {kept_rows:,}"
                    )

        with mapping_tmp.open(
            "wt", encoding="utf-8", newline="", buffering=1024 * 1024
        ) as destination:
            writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
            writer.writerow(
                (
                    "m6a_site_key",
                    "chrom",
                    "genomic_pos_1based",
                    "mapping_status",
                    "selected_transcript",
                    "strand",
                    "cds_pos_1based",
                    "codon_first_cds_pos_1based",
                    "genomic_codon_first_pos_1based",
                    "codon_index_1based",
                    "codon_offset_0based",
                    "candidate_count",
                )
            )
            for record, candidate, count in zip(
                m6a_records, selected, candidate_counts, strict=True
            ):
                if candidate is None:
                    writer.writerow(
                        (
                            record.site_key,
                            record.chrom,
                            record.genomic_pos,
                            "unmapped",
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            0,
                        )
                    )
                else:
                    writer.writerow(
                        (
                            record.site_key,
                            record.chrom,
                            record.genomic_pos,
                            "mapped",
                            candidate.transcript,
                            candidate.strand,
                            candidate.cds_pos,
                            candidate.codon_first_cds_pos,
                            candidate.genomic_codon_first_pos,
                            candidate.codon_index,
                            candidate.codon_offset,
                            count,
                        )
                    )

        replace_outputs(
            [(annotation_tmp, output_annotation), (mapping_tmp, output_m6a_map)]
        )
    except BaseException:
        cleanup_temps([annotation_tmp, mapping_tmp])
        raise

    mapped = sum(candidate is not None for candidate in selected)
    ambiguous = sum(count > 1 for count in candidate_counts)
    stats = {
        "annotation_input_rows": input_rows,
        "annotation_kept_rows": kept_rows,
        "m6a_total": len(m6a_records),
        "m6a_mapped": mapped,
        "m6a_unmapped": len(m6a_records) - mapped,
        "m6a_multiple_candidates": ambiguous,
    }
    log(
        "prepare complete: "
        + ", ".join(f"{key}={value:,}" for key, value in stats.items())
    )
    return stats


def codon_first_position(cds_position: int) -> int:
    if cds_position < 1:
        raise ValueError(f"CDS position must be >= 1, got {cds_position}")
    return 1 + ((cds_position - 1) // 3) * 3


def load_masked_codons(m6a_map: Path, editing_sites: Path) -> set[tuple[str, int]]:
    masked: set[tuple[str, int]] = set()
    with open_text(m6a_map) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "mapping_status",
            "selected_transcript",
            "codon_first_cds_pos_1based",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected m6A mapping header in {m6a_map}")
        for row in reader:
            if row["mapping_status"] != "mapped":
                continue
            masked.add(
                (
                    row["selected_transcript"],
                    int(row["codon_first_cds_pos_1based"]),
                )
            )

    with open_text(editing_sites) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"new_transcript", "new_cds_pos_1based", "target_region"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected editing-site header in {editing_sites}")
        for row in reader:
            if row["target_region"] != "CDS":
                continue
            transcript = row["new_transcript"].strip()
            position_text = row["new_cds_pos_1based"].strip()
            if not transcript or not position_text or position_text == "NA":
                continue
            position = int(position_text)
            if position < 1:
                continue
            masked.add((transcript, codon_first_position(position)))
    return masked


@dataclass(frozen=True)
class AnnotationRow:
    text: str
    cds_codon_start: int


def iter_annotation_groups(path: Path) -> Iterator[tuple[str, list[AnnotationRow]]]:
    current_transcript: str | None = None
    rows: list[AnnotationRow] = []
    previous_cds_position = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            text = line.rstrip("\r\n")
            if not text:
                continue
            fields = text.split("\t")
            if len(fields) < 7:
                raise ValueError(f"Expected at least 7 columns at {path}:{line_number}")
            transcript = fields[2]
            try:
                cds_position = int(fields[5])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid CDS position at {path}:{line_number}"
                ) from exc

            if current_transcript is None:
                current_transcript = transcript
                previous_cds_position = 0
            elif transcript != current_transcript:
                if transcript < current_transcript:
                    raise ValueError(
                        f"Annotation transcript groups are not sorted: "
                        f"{current_transcript!r} before {transcript!r} at line {line_number}"
                    )
                yield current_transcript, rows
                current_transcript = transcript
                rows = []
                previous_cds_position = 0

            if cds_position <= previous_cds_position:
                raise ValueError(
                    f"CDS positions are not strictly increasing for {transcript} "
                    f"at {path}:{line_number}"
                )
            rows.append(AnnotationRow(text, cds_position))
            previous_cds_position = cds_position

    if current_transcript is not None:
        yield current_transcript, rows


def parse_count_line(path: Path, line_number: int, line: str) -> tuple[str, int, int]:
    fields = line.rstrip("\r\n").split("\t", 3)
    if len(fields) != 4:
        raise ValueError(f"Expected 4 columns at {path}:{line_number}")
    try:
        position = int(fields[1])
        count = int(fields[3])
    except ValueError as exc:
        raise ValueError(
            f"Invalid numeric count field at {path}:{line_number}"
        ) from exc
    return fields[0], position, count


def iter_count_groups(path: Path) -> Iterator[tuple[str, dict[int, int]]]:
    current_transcript: str | None = None
    counts: dict[int, int] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            transcript, position, count = parse_count_line(path, line_number, line)
            if current_transcript is None:
                current_transcript = transcript
            elif transcript != current_transcript:
                if transcript < current_transcript:
                    raise ValueError(
                        f"Count transcript groups are not sorted in {path}: "
                        f"{current_transcript!r} before {transcript!r}"
                    )
                yield current_transcript, counts
                current_transcript = transcript
                counts = {}
            counts[position] = counts.get(position, 0) + count
    if current_transcript is not None:
        yield current_transcript, counts


class CountGroupLookup:
    def __init__(self, path: Path):
        self.path = path
        self.groups = iter(iter_count_groups(path))
        self.current: tuple[str, dict[int, int]] | None = next(self.groups, None)

    def get(self, transcript: str) -> dict[int, int]:
        while self.current is not None and self.current[0] < transcript:
            self.current = next(self.groups, None)
        if self.current is None or self.current[0] != transcript:
            return {}
        counts = self.current[1]
        self.current = next(self.groups, None)
        return counts


def sum_count_file(path_text: str) -> tuple[str, int, int]:
    path = Path(path_text)
    total = 0
    rows = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                count = int(line.rsplit("\t", 1)[1].strip())
            except (ValueError, IndexError) as exc:
                raise ValueError(f"Invalid count at {path}:{line_number}") from exc
            total += count
            rows += 1
            if rows % PROGRESS_EVERY == 0:
                log(f"library sum: {path.name}: scanned {rows:,} rows")
    log(f"library sum complete: {path.name}: rows={rows:,}, sum={total:,}")
    return str(path), total, rows


def prefix_sum(values: list[int], valid: list[bool] | None = None) -> list[int]:
    result = [0] * (len(values) + 1)
    running = 0
    if valid is None:
        for index, value in enumerate(values, 1):
            running += value
            result[index] = running
    else:
        for index, (value, keep) in enumerate(zip(values, valid, strict=True), 1):
            if keep:
                running += value
            result[index] = running
    return result


def flank_bounds(index: int, length: int, flank_size: int) -> tuple[int, int]:
    """Return [start, end) around focal index, with edge compensation.

    The returned interval contains at most 2 * flank_size non-focal positions;
    the focal position itself must still be subtracted by the caller.
    """
    left_count = min(flank_size, index)
    right_count = min(flank_size, length - index - 1)
    needed = 2 * flank_size - left_count - right_count
    if needed:
        add_right = min(needed, length - index - 1 - right_count)
        right_count += add_right
        needed -= add_right
    if needed:
        add_left = min(needed, index - left_count)
        left_count += add_left
    return index - left_count, index + right_count + 1


def flank_metrics(
    index: int,
    flank_size: int,
    length: int,
    rna_prefix: list[int],
    rpf_prefix: list[int],
    valid_prefix: list[int],
) -> tuple[int, int, int]:
    start, end = flank_bounds(index, length, flank_size)
    rna = (rna_prefix[index] - rna_prefix[start]) + (
        rna_prefix[end] - rna_prefix[index + 1]
    )
    rpf = (rpf_prefix[index] - rpf_prefix[start]) + (
        rpf_prefix[end] - rpf_prefix[index + 1]
    )
    valid_n = (valid_prefix[index] - valid_prefix[start]) + (
        valid_prefix[end] - valid_prefix[index + 1]
    )
    return rna, rpf, valid_n


@dataclass(frozen=True)
class Sample:
    sample_id: str
    rna_count: Path
    rpf_count: Path


@dataclass(frozen=True)
class CountJob:
    sample: Sample
    annotation: Path
    m6a_map: Path
    editing_sites: Path
    output_dir: Path
    flank_sizes: tuple[int, ...]
    exclude_start: int
    exclude_end: int
    library_rna_sum: int
    library_rpf_sum: int
    overwrite: bool


def parse_sample_sheet(path: Path) -> list[Sample]:
    samples: list[Sample] = []
    with open_text(path) as handle:
        reader = csv.DictReader(
            (
                line
                for line in handle
                if line.strip() and not line.lstrip().startswith("#")
            ),
            delimiter="\t",
        )
        required = {"sample_id", "rna_count", "rpf_count"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Sample sheet must contain columns: {', '.join(sorted(required))}"
            )
        for row in reader:
            sample_id = row["sample_id"].strip()
            if not sample_id or not re.fullmatch(r"[A-Za-z0-9_.-]+", sample_id):
                raise ValueError(
                    f"Invalid sample_id {sample_id!r}; use letters, digits, _, . or -"
                )
            paths: list[Path] = []
            for column in ("rna_count", "rpf_count"):
                value = Path(row[column].strip()).expanduser()
                if not value.is_absolute():
                    value = (path.parent / value).resolve()
                paths.append(value)
            samples.append(Sample(sample_id, paths[0], paths[1]))
    if not samples:
        raise ValueError(f"No samples found in {path}")
    ids = [sample.sample_id for sample in samples]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate sample_id in {path}")
    return samples


def sample_output_paths(
    sample: Sample, output_dir: Path, flank_sizes: tuple[int, ...]
) -> dict[int, Path]:
    return {
        flank_size: output_dir / f"{sample.sample_id}.flank{flank_size}.tsv"
        for flank_size in flank_sizes
    }


def count_one_sample(job: CountJob) -> dict[str, int | str]:
    sample = job.sample
    outputs = sample_output_paths(sample, job.output_dir, job.flank_sizes)
    check_output_paths(list(outputs.values()), job.overwrite)
    temporary = {size: tmp_path(path) for size, path in outputs.items()}
    cleanup_temps(list(temporary.values()))
    masked_codons = load_masked_codons(job.m6a_map, job.editing_sites)
    rna_lookup = CountGroupLookup(sample.rna_count)
    rpf_lookup = CountGroupLookup(sample.rpf_count)

    handles: dict[int, TextIO] = {}
    annotation_rows = 0
    masked_annotation_rows = 0
    try:
        for size, path in temporary.items():
            handle = path.open(
                "wt", encoding="utf-8", newline="", buffering=1024 * 1024
            )
            handles[size] = handle
            handle.write("\t".join(COUNT_HEADER) + "\n")

        for transcript, rows in iter_annotation_groups(job.annotation):
            rna_by_position = rna_lookup.get(transcript)
            rpf_by_position = rpf_lookup.get(transcript)
            rna_counts = [rna_by_position.get(row.cds_codon_start, 0) for row in rows]
            rpf_counts = [rpf_by_position.get(row.cds_codon_start, 0) for row in rows]
            length = len(rows)
            valid = []
            for index, row in enumerate(rows):
                keep = (
                    index >= job.exclude_start
                    and index < length - job.exclude_end
                    and (transcript, row.cds_codon_start) not in masked_codons
                )
                valid.append(keep)
                if not keep:
                    masked_annotation_rows += 1

            rna_prefix = prefix_sum(rna_counts, valid)
            rpf_prefix = prefix_sum(rpf_counts, valid)
            valid_prefix = prefix_sum([1] * length, valid)

            for index, row in enumerate(rows):
                common = f"{row.text}\t{rna_counts[index]}\t{rpf_counts[index]}"
                for flank_size, handle in handles.items():
                    flank_rna, flank_rpf, flank_n = flank_metrics(
                        index,
                        flank_size,
                        length,
                        rna_prefix,
                        rpf_prefix,
                        valid_prefix,
                    )
                    handle.write(
                        f"{common}\t{flank_rna}\t{flank_rpf}\t{flank_n}\t"
                        f"{job.library_rna_sum}\t{job.library_rpf_sum}\n"
                    )
            annotation_rows += length
            if (
                annotation_rows // PROGRESS_EVERY
                != (annotation_rows - length) // PROGRESS_EVERY
            ):
                log(
                    f"{sample.sample_id}: processed {annotation_rows:,} annotation rows"
                )

        for handle in handles.values():
            handle.close()
        handles.clear()
        replace_outputs([(temporary[size], outputs[size]) for size in job.flank_sizes])
    except BaseException:
        for handle in handles.values():
            handle.close()
        cleanup_temps(list(temporary.values()))
        raise

    stats: dict[str, int | str] = {
        "sample_id": sample.sample_id,
        "annotation_rows": annotation_rows,
        "masked_annotation_rows": masked_annotation_rows,
        "library_RNA_sum": job.library_rna_sum,
        "library_RPF_sum": job.library_rpf_sum,
    }
    log(
        f"{sample.sample_id} complete: rows={annotation_rows:,}, "
        f"RNA_sum={job.library_rna_sum:,}, RPF_sum={job.library_rpf_sum:,}"
    )
    return stats


def calculate_library_totals(
    samples: list[Sample], jobs: int
) -> dict[str, tuple[int, int]]:
    paths = sorted(
        {
            str(path)
            for sample in samples
            for path in (sample.rna_count, sample.rpf_count)
        }
    )
    totals: dict[str, tuple[int, int]] = {}
    if jobs == 1 or len(paths) == 1:
        for path in paths:
            path_text, total, rows = sum_count_file(path)
            totals[path_text] = (total, rows)
        return totals
    with ProcessPoolExecutor(max_workers=min(jobs, len(paths))) as executor:
        futures = [executor.submit(sum_count_file, path) for path in paths]
        for future in as_completed(futures):
            path_text, total, rows = future.result()
            totals[path_text] = (total, rows)
    return totals


def write_library_summary(
    path: Path,
    samples: list[Sample],
    totals: dict[str, tuple[int, int]],
    overwrite: bool,
) -> None:
    check_output_paths([path], overwrite)
    temporary = tmp_path(path)
    cleanup_temps([temporary])
    try:
        with temporary.open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                (
                    "sample_id",
                    "library_RNA_sum",
                    "library_RPF_sum",
                    "RNA_input_rows",
                    "RPF_input_rows",
                    "rna_count_file",
                    "rpf_count_file",
                )
            )
            for sample in samples:
                rna_sum, rna_rows = totals[str(sample.rna_count)]
                rpf_sum, rpf_rows = totals[str(sample.rpf_count)]
                writer.writerow(
                    (
                        sample.sample_id,
                        rna_sum,
                        rpf_sum,
                        rna_rows,
                        rpf_rows,
                        sample.rna_count,
                        sample.rpf_count,
                    )
                )
        os.replace(temporary, path)
    except BaseException:
        cleanup_temps([temporary])
        raise


def run_count(args: argparse.Namespace) -> list[dict[str, int | str]]:
    samples = parse_sample_sheet(args.sample_sheet)
    flank_sizes = tuple(sorted(set(args.flanks)))
    if any(size < 1 for size in flank_sizes):
        raise ValueError("All flank sizes must be >= 1")
    if args.exclude_start < 0 or args.exclude_end < 0:
        raise ValueError("CDS exclusion sizes must be >= 0")
    jobs = max(1, args.jobs)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "library_depths.tsv"
    planned_outputs = [summary_path]
    for sample in samples:
        planned_outputs.extend(
            sample_output_paths(sample, args.output_dir, flank_sizes).values()
        )
    check_output_paths(planned_outputs, args.overwrite)

    log(
        f"count: {len(samples)} sample(s), flank sizes={flank_sizes}, jobs={jobs}; "
        "first pass calculates library depths"
    )
    totals = calculate_library_totals(samples, jobs)
    count_jobs = [
        CountJob(
            sample=sample,
            annotation=args.annotation_longest,
            m6a_map=args.m6a_map,
            editing_sites=args.editing_sites,
            output_dir=args.output_dir,
            flank_sizes=flank_sizes,
            exclude_start=args.exclude_start,
            exclude_end=args.exclude_end,
            library_rna_sum=totals[str(sample.rna_count)][0],
            library_rpf_sum=totals[str(sample.rpf_count)][0],
            overwrite=args.overwrite,
        )
        for sample in samples
    ]

    results: list[dict[str, int | str]] = []
    if jobs == 1 or len(count_jobs) == 1:
        for job in count_jobs:
            results.append(count_one_sample(job))
    else:
        with ProcessPoolExecutor(max_workers=min(jobs, len(count_jobs))) as executor:
            futures = [executor.submit(count_one_sample, job) for job in count_jobs]
            for future in as_completed(futures):
                results.append(future.result())

    write_library_summary(summary_path, samples, totals, args.overwrite)
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Map m6A sites and calculate focal/flanking codon coverage."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser(
        "prepare", help="filter annotation and map m6A sites"
    )
    prepare.add_argument("--annotation", type=Path, required=True)
    prepare.add_argument("--longest-transcripts", type=Path, required=True)
    prepare.add_argument("--m6a", type=Path, required=True)
    prepare.add_argument("--output-annotation", type=Path, required=True)
    prepare.add_argument("--output-m6a-map", type=Path, required=True)
    prepare.add_argument("--overwrite", action="store_true")

    count = subparsers.add_parser(
        "count", help="calculate focal and flanking RNA/RPF coverage"
    )
    count.add_argument("--annotation-longest", type=Path, required=True)
    count.add_argument("--m6a-map", type=Path, required=True)
    count.add_argument("--editing-sites", type=Path, required=True)
    count.add_argument("--sample-sheet", type=Path, required=True)
    count.add_argument("--output-dir", type=Path, required=True)
    count.add_argument("--flanks", type=int, nargs="+", default=[20, 30])
    count.add_argument("--exclude-start", type=int, default=15)
    count.add_argument("--exclude-end", type=int, default=5)
    count.add_argument(
        "--jobs",
        type=int,
        default=2,
        help="parallel files/samples (default: 2; low-memory per worker)",
    )
    count.add_argument("--overwrite", action="store_true")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            prepare_annotation(
                annotation=args.annotation,
                longest_transcripts=args.longest_transcripts,
                m6a_sites=args.m6a,
                output_annotation=args.output_annotation,
                output_m6a_map=args.output_m6a_map,
                overwrite=args.overwrite,
            )
        elif args.command == "count":
            run_count(args)
        else:
            parser.error(f"Unknown command: {args.command}")
    except (OSError, ValueError) as exc:
        log(f"ERROR: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
