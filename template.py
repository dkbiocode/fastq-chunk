#!/usr/bin/env python3
import functools
import gzip
import io
import logging
import os
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import dnaio

from fastq_chunk import run_parallel

INPUT_PATH = "sample_data/sample_50k_R1.fastq.gz"
OUTPUT_PATH = "output.fastq.gz"
N_WORKERS = 2
BUFFER_SIZE = 4 * 1024 * 1024

logging.basicConfig(level=logging.INFO)


def process_chunk(chunk_bytes: bytes, chunk_idx: int, *, temp_dir: str) -> str:
    out_path = os.path.join(temp_dir, f"chunk_{chunk_idx:06d}.fastq.gz")
    print(f"process_chunk: {chunk_idx}; to {out_path}", file=sys.stderr)
    with dnaio.open(io.BytesIO(chunk_bytes)) as fin, gzip.open(out_path, "wt") as fout:
        for rec in fin:
            # replace with your processing
            fout.write(f"@{rec.name}\n{rec.sequence}\n+\n{rec.qualities}\n")
    return out_path


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        worker = functools.partial(process_chunk, temp_dir=tmp)
        chunk_paths = list(run_parallel(INPUT_PATH, worker,
                                        buffer_size=BUFFER_SIZE,
                                        n_workers=N_WORKERS,
                                        executor_class=ProcessPoolExecutor))
        with open(OUTPUT_PATH, "wb") as fout:
            for path in chunk_paths:
                with open(path, "rb") as fin:
                    shutil.copyfileobj(fin, fout)

    print(f"wrote {len(chunk_paths)} chunk(s) → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
