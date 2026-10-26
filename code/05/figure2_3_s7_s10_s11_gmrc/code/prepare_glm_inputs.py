#!/usr/bin/env python3
"""Expose the eight new single-library tables in the legacy GLM input format."""

import gzip
import hashlib
from pathlib import Path
import recompute as r

SAMPLES = [f"WBg_{sex}{i}" for sex in "FM" for i in range(1, 4)]


def main():
    rows = []
    for sample in SAMPLES:
        src = r.ROOT / "per_library" / f"{sample}.codonTE_AEP.AA.gz"
        dst = r.ROOT / "glm_inputs" / f"{sample}.codonTE_AEP.AA"
        tmp = dst.with_suffix(".partial")
        h = hashlib.sha256()
        with gzip.open(src, "rb") as f, tmp.open("wb") as out:
            while block := f.read(4 * 1024 * 1024):
                h.update(block)
                out.write(block)
        tmp.replace(dst)
        rows.append(
            {
                "sample": sample,
                "source": str(src),
                "GLM_input": str(dst),
                "uncompressed_sha256": h.hexdigest(),
                "bytes": dst.stat().st_size,
            }
        )
        print("GLM INPUT", sample, flush=True)
    r.tsv_write(r.ROOT / "manifests" / "glm_inputs.tsv", rows)


if __name__ == "__main__":
    main()
