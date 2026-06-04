#!/usr/bin/env python3
import functools
import gzip
import logging
import os
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

from fastq_chunk import FastqRecord, get_read_dimensions, calculate_chunk_size, run_parallel

INPUT_PATH = "sample_data/sample_50k_R1.fastq.gz"
OUTPUT_PATH = "output.fastq.gz"
N_WORKERS = 2
MEM_PER_THREAD_MB = 20

logging.basicConfig(level=logging.INFO)


def process_chunk(
    chunk: list[FastqRecord],
    chunk_idx: int,
    *,
    temp_dir: str,
) -> str:
    out_path = os.path.join(temp_dir, f"chunk_{chunk_idx:06d}.fastq.gz")
    print(f"process_chunk: {chunk_idx}; {len(chunk)} sequences to {out_path}", file=sys.stderr)
    with gzip.open(out_path, "wt") as fout:
        for rec in chunk:
            # replace with your processing
            fout.write(f"@{rec.name}\n{rec.sequence}\n+\n{rec.qualities}\n")

    return out_path


def main() -> None:
    dims = get_read_dimensions(INPUT_PATH)
    if dims is None:
        raise SystemExit(f"no reads in {INPUT_PATH}")
    _, _, mem_per_read = dims
    chunk_size = calculate_chunk_size(mem_per_read, MEM_PER_THREAD_MB)

    with tempfile.TemporaryDirectory() as tmp:
        worker = functools.partial(process_chunk, temp_dir=tmp)
        chunk_paths = list(run_parallel(INPUT_PATH, worker,
                                        chunk_size=chunk_size,
                                        n_workers=N_WORKERS,
                                        executor_class = ProcessPoolExecutor))
        with open(OUTPUT_PATH, "wb") as fout:
            for path in chunk_paths:
                with open(path, "rb") as fin:
                    shutil.copyfileobj(fin, fout)

    print(f"wrote {len(chunk_paths)} chunk(s) → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
