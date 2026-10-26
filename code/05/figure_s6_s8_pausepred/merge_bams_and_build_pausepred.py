#!/usr/bin/env python3
"""Merge biological-replicate BAMs and build a PausePred sample sheet.

The PS tag is summarized by read length.  PausePred A-site offsets are
calculated as the modal P-site offset plus three nucleotides.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from functools import partial
from pathlib import Path

import pysam

SOURCE_DIR_DEFAULT = Path(
    "/path/to/RNA_editing_DATAprocessing_new/05.psite/psite_2times"
)
OUTPUT_DIR_DEFAULT = (
    Path(__file__).resolve().parents[3] / "results/05_stalling/psite_2times_merged_bam"
)
BAM_SUFFIX = "_RPF.transcriptome.tag.bam"
READ_LENGTHS = tuple(range(25, 36))


@dataclass(frozen=True)
class Sample:
    name: str
    inputs: tuple[Path, ...]
    analysis_bam: Path
    action: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=SOURCE_DIR_DEFAULT)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR_DEFAULT)
    parser.add_argument("--jobs", type=int, default=2, help="Concurrent BAM jobs")
    parser.add_argument(
        "--threads", type=int, default=4, help="Threads per merge/index"
    )
    parser.add_argument(
        "--stage",
        choices=("all", "merge", "offsets"),
        default="all",
        help="Run all steps, merging only, or offset/sample-sheet generation only",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def discover_samples(
    source_dir: Path, output_dir: Path
) -> tuple[list[Sample], list[Path]]:
    bams = sorted(source_dir.glob(f"*{BAM_SUFFIX}"))
    if not bams:
        raise RuntimeError(f"No *{BAM_SUFFIX} files found in {source_dir}")

    groups: dict[str, list[Path]] = defaultdict(list)
    excluded: list[Path] = []
    standard_re = re.compile(
        r"^(?P<sample>.+_[FM])(?P<rep>[0-9]+)" + re.escape(BAM_SUFFIX) + r"$"
    )
    adar_re = re.compile(r"^Dmel_Adar5G1_rep[0-9]+_ribo" + re.escape(BAM_SUFFIX) + r"$")
    male_head_re = re.compile(
        r"^(?P<species>Dere|Dsec|Dsim)_male_head_ribo" + re.escape(BAM_SUFFIX) + r"$"
    )
    female_head_re = re.compile(
        r"^Dsim_female_head_ribo" + re.escape(BAM_SUFFIX) + r"$"
    )

    for bam in bams:
        base = bam.name
        if base.startswith("Dmel_w1118_rep"):
            excluded.append(bam)
            continue

        match = standard_re.match(base)
        if match:
            groups[match.group("sample")].append(bam)
        elif adar_re.match(base):
            groups["Dmel_Adar5G1_M"].append(bam)
        elif match := male_head_re.match(base):
            groups[f"{match.group('species')}_M"].append(bam)
        elif female_head_re.match(base):
            groups["Dsim_F"].append(bam)
        else:
            groups[base[: -len(BAM_SUFFIX)]].append(bam)

    samples: list[Sample] = []
    for name, inputs_list in sorted(groups.items()):
        inputs = tuple(sorted(inputs_list))
        if len(inputs) > 1:
            analysis_bam = output_dir / f"{name}{BAM_SUFFIX}"
            action = "merge"
        else:
            # Keep the BAM data unchanged, but expose it in the output folder
            # via a symlink so its index can be written beside it.
            analysis_bam = output_dir / inputs[0].name
            action = "link_singleton"
        samples.append(Sample(name, inputs, analysis_bam, action))
    return samples, excluded


def atomic_tsv(path: Path, rows: list[list[object]]) -> None:
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with tmp.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)
    os.replace(tmp, path)


def write_manifest(
    output_dir: Path, samples: list[Sample], excluded: list[Path]
) -> None:
    rows: list[list[object]] = [
        ["sample", "action", "n_inputs", "analysis_bam", "input_bams"]
    ]
    for sample in samples:
        rows.append(
            [
                sample.name,
                sample.action,
                len(sample.inputs),
                sample.analysis_bam,
                ";".join(str(path.resolve()) for path in sample.inputs),
            ]
        )
    atomic_tsv(output_dir / "sample_merge_manifest.tsv", rows)
    atomic_tsv(
        output_dir / "excluded_bams.tsv",
        [["reason", "bam"]]
        + [["Dmel_w1118 excluded by request", path.resolve()] for path in excluded],
    )


def quickcheck(path: Path) -> None:
    try:
        result = pysam.quickcheck("-v", str(path))
    except pysam.SamtoolsError as exc:
        raise RuntimeError(f"BAM quickcheck failed for {path}: {exc}") from exc
    if result.strip():
        raise RuntimeError(f"BAM quickcheck failed for {path}: {result.strip()}")


def validate_merge_inputs(sample: Sample) -> None:
    expected_sq: list[tuple[str, int]] | None = None
    for path in sample.inputs:
        quickcheck(path)
        with pysam.AlignmentFile(path, "rb") as bam:
            header = bam.header.to_dict()
            sort_order = header.get("HD", {}).get("SO")
            if sort_order != "coordinate":
                raise RuntimeError(
                    f"Input BAM is not declared coordinate-sorted ({sort_order!r}): {path}"
                )
            sq = [(entry["SN"], int(entry["LN"])) for entry in header.get("SQ", [])]
            if expected_sq is None:
                expected_sq = sq
            elif sq != expected_sq:
                raise RuntimeError(
                    f"Reference sequence headers differ within group {sample.name}"
                )


def merge_one(sample: Sample, threads: int) -> str:
    output = sample.analysis_bam
    output.parent.mkdir(parents=True, exist_ok=True)

    if output.is_file() and output.stat().st_size > 0:
        quickcheck(output)
        bai = Path(f"{output}.bai")
        if not bai.is_file() or bai.stat().st_mtime < output.stat().st_mtime:
            pysam.index("-@", str(threads), str(output))
        return f"SKIP verified existing: {output.name}"

    validate_merge_inputs(sample)
    tmp = output.with_name(f".{output.stem}.tmp.{os.getpid()}.bam")
    try:
        pysam.merge(
            "-f",
            "-@",
            str(threads),
            "-o",
            str(tmp),
            *(str(path) for path in sample.inputs),
            catch_stdout=False,
        )
        quickcheck(tmp)
        os.replace(tmp, output)
        pysam.index("-@", str(threads), str(output))
    finally:
        if tmp.exists():
            tmp.unlink()
    return f"MERGED {len(sample.inputs)} inputs: {output.name}"


def link_singleton(sample: Sample, threads: int) -> str:
    source = sample.inputs[0].resolve()
    output = sample.analysis_bam
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_symlink():
        if output.resolve() != source:
            raise RuntimeError(f"Singleton symlink points to the wrong BAM: {output}")
    elif output.exists():
        raise RuntimeError(f"Singleton target exists but is not a symlink: {output}")
    else:
        output.symlink_to(source)
    quickcheck(output)
    bai = Path(f"{output}.bai")
    if not bai.is_file() or bai.stat().st_mtime < source.stat().st_mtime:
        pysam.index("-@", str(threads), str(output))
    return f"LINKED singleton: {output.name}"


def summarize_one(
    sample: Sample, offsets_dir: Path
) -> tuple[str, str, str, str, list[str]]:
    bam_path = sample.analysis_bam
    quickcheck(bam_path)
    counts: dict[int, Counter[int]] = defaultdict(Counter)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.mapping_quality < 10 or not aln.has_tag("PS"):
                continue
            if aln.query_sequence is None:
                continue
            counts[len(aln.query_sequence)][int(aln.get_tag("PS"))] += 1

    detail_rows: list[list[object]] = [
        [
            "read_length",
            "mode_ps_offset",
            "n_mode",
            "n_total",
            "mode_fraction",
            "all_offsets",
        ]
    ]
    modes: dict[int, int] = {}
    for read_len in sorted(counts):
        total = sum(counts[read_len].values())
        mode_offset, n_mode = counts[read_len].most_common(1)[0]
        modes[read_len] = mode_offset
        all_offsets = ",".join(
            f"{offset}:{n}" for offset, n in counts[read_len].most_common()
        )
        detail_rows.append(
            [read_len, mode_offset, n_mode, total, f"{n_mode / total:.3f}", all_offsets]
        )

    atomic_tsv(offsets_dir / f"{sample.name}.psite_offset_by_length.tsv", detail_rows)
    atomic_tsv(
        offsets_dir / f"{sample.name}.psite_offsets.pausepred.tsv",
        [["length", "psite_offset"]]
        + [[length, modes[length]] for length in READ_LENGTHS if length in modes],
    )

    warnings = [
        f"{sample.name}: missing read length {length}"
        for length in READ_LENGTHS
        if length not in modes
    ]
    # Omit lengths lacking a measured PS-tag mode. Passing a literal NA to
    # PausePred's Perl numeric comparison would incorrectly coerce it to zero.
    available_lengths = [length for length in READ_LENGTHS if length in modes]
    read_lengths = ",".join(str(length) for length in available_lengths)
    asite_offsets = ",".join(str(modes[length] + 3) for length in available_lengths)
    return sample.name, str(bam_path), read_lengths, asite_offsets, warnings


def run_parallel(function, items: list[Sample], jobs: int, label: str) -> list[object]:
    results: list[object] = []
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(function, item): item for item in items}
        for future in as_completed(futures):
            item = futures[future]
            result = future.result()
            print(
                f"[{label}] {item.name}: {result if isinstance(result, str) else 'done'}",
                flush=True,
            )
            results.append(result)
    return results


def main() -> int:
    args = parse_args()
    if args.jobs < 1 or args.threads < 1:
        raise ValueError("--jobs and --threads must be positive integers")
    source_dir = args.source_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    samples, excluded = discover_samples(source_dir, output_dir)
    write_manifest(output_dir, samples, excluded)
    merge_samples = [sample for sample in samples if sample.action == "merge"]
    singleton_samples = [
        sample for sample in samples if sample.action == "link_singleton"
    ]

    print(f"[INFO] source: {source_dir}")
    print(f"[INFO] output: {output_dir}")
    print(
        f"[INFO] samples: {len(samples)} ({len(merge_samples)} merge, "
        f"{len(singleton_samples)} singleton)"
    )
    print(f"[INFO] excluded BAMs: {len(excluded)}")
    for sample in samples:
        print(
            f"[PLAN] {sample.name}\t{sample.action}\t{len(sample.inputs)}\t{sample.analysis_bam}"
        )
    if args.dry_run:
        return 0

    if args.stage in ("all", "merge"):
        run_parallel(
            partial(merge_one, threads=args.threads),
            merge_samples,
            args.jobs,
            "MERGE",
        )
        run_parallel(
            partial(link_singleton, threads=args.threads),
            singleton_samples,
            min(args.jobs, len(singleton_samples)) or 1,
            "SINGLETON",
        )

    if args.stage in ("all", "offsets"):
        missing = [
            sample.analysis_bam
            for sample in samples
            if not sample.analysis_bam.is_file()
        ]
        if missing:
            raise RuntimeError(
                "Required BAMs are missing:\n" + "\n".join(map(str, missing))
            )
        offsets_dir = output_dir / "offsets"
        offsets_dir.mkdir(parents=True, exist_ok=True)
        results = run_parallel(
            partial(summarize_one, offsets_dir=offsets_dir),
            samples,
            args.jobs,
            "OFFSETS",
        )
        parsed_results = sorted(results, key=lambda result: result[0])
        sheet_rows: list[list[object]] = [["sample", "bam", "read_lengths", "offsets"]]
        warnings: list[str] = []
        for name, bam, lengths, offsets, sample_warnings in parsed_results:
            sheet_rows.append([name, bam, lengths, offsets])
            warnings.extend(sample_warnings)
        atomic_tsv(output_dir / "pausepred_sample_sheet.offset_plus3.tsv", sheet_rows)
        warning_path = output_dir / "pausepred_sample_sheet.offset_plus3.warn.log"
        warning_tmp = warning_path.with_name(f".{warning_path.name}.tmp.{os.getpid()}")
        warning_tmp.write_text("".join(f"[WARN] {warning}\n" for warning in warnings))
        os.replace(warning_tmp, warning_path)
        print(f"[DONE] {output_dir / 'pausepred_sample_sheet.offset_plus3.tsv'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
