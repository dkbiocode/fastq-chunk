#!/usr/bin/env python3
"""
Benchmark chunked parallel fastq read/split/write using fastq_chunk.run_parallel.

Usage:
    python benchmark_fastq_chunk.py [--threads N] [--mem-per-thread MB]
                                    [--input FILE] [--repeats N]

This script also serves as a usage example for fastq_chunk: define a worker
function with the (chunk, chunk_idx, *, ...) signature, bind extra context with
functools.partial, and pass it to run_parallel.
"""
import argparse
import functools
import gzip
import logging
import os
import shutil
import tempfile
import time

from fastq_chunk import FastqRecord, get_read_dimensions, calculate_chunk_size, run_parallel

# Suppress INFO-level output from fastq_chunk during timed runs; derived values
# are printed explicitly from return values instead.
logging.basicConfig(level=logging.WARNING)


# ---------------------------------------------------------------------------
# Worker — the function passed to run_parallel
# ---------------------------------------------------------------------------

def passthrough_chunk(
    chunk: list[FastqRecord],
    chunk_idx: int,
    *,
    temp_dir: str,
) -> str:
    """Write chunk to a numbered temp gzip file, unchanged.

    temp_dir is keyword-only so functools.partial can bind it, leaving
    (chunk, chunk_idx) as the two-argument worker contract run_parallel expects.
    Replace the body here with any per-chunk processing to adapt this pattern.
    """
    temp_path = os.path.join(temp_dir, f"chunk_{chunk_idx:06d}.fastq.gz")
    with gzip.open(temp_path, 'wt') as fout:
        for rec in chunk:
            fout.write(f"@{rec.name}\n{rec.sequence}\n+\n{rec.qualities}\n")
    return temp_path


# ---------------------------------------------------------------------------
# Benchmark driver
# ---------------------------------------------------------------------------

def run_benchmark(
    input_path: str,
    threads: int,
    mem_per_thread_mb: int,
    repeats: int,
) -> None:
    input_mb = os.path.getsize(input_path) / 1024 / 1024

    dims = get_read_dimensions(input_path)
    if dims is None:
        raise ValueError(f"could not read any records from {input_path}")
    read_len, bytes_per_read, mem_per_read = dims
    chunk_size = calculate_chunk_size(mem_per_read, mem_per_thread_mb)

    print(f"\n=== benchmark_fastq_chunk.py: threads={threads}, mem_per_thread={mem_per_thread_mb} MB ===")
    print(f"Input:       {os.path.basename(input_path)}  ({input_mb:.1f} MB)")
    print(f"read_len:    {read_len} bp  |  mem_per_read: {mem_per_read:,} bytes  |  chunk_size: {chunk_size:,} reads")
    print()

    elapsed_times: list[float] = []

    for run in range(1, repeats + 1):
        working_dir = tempfile.mkdtemp(prefix="bench_")
        try:
            worker = functools.partial(passthrough_chunk, temp_dir=working_dir)
            t0 = time.perf_counter()
            temp_paths = list(
                run_parallel(input_path, worker, chunk_size=chunk_size, n_workers=threads)
            )
            elapsed = time.perf_counter() - t0
        finally:
            shutil.rmtree(working_dir)

        elapsed_times.append(elapsed)
        mbps = input_mb / elapsed
        print(f"  Run {run}/{repeats}:  {len(temp_paths)} chunk(s)   {elapsed:.2f} s  ({mbps:.1f} MB/s)")

    print()
    min_t  = min(elapsed_times)
    mean_t = sum(elapsed_times) / len(elapsed_times)
    max_t  = max(elapsed_times)
    print(f"  min   {min_t:.2f} s   {input_mb / min_t:.1f} MB/s")
    print(f"  mean  {mean_t:.2f} s   {input_mb / mean_t:.1f} MB/s")
    print(f"  max   {max_t:.2f} s   {input_mb / max_t:.1f} MB/s")
    print()


def main() -> None:
    default_input = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "sample_data", "2629LEI-2_S2_L001_R1_001.fastq.gz",
    )

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--threads",        type=int, default=4,            metavar="N")
    parser.add_argument("--mem-per-thread", type=int, default=3840,         metavar="MB")
    parser.add_argument("--input",          type=str, default=default_input, metavar="FILE")
    parser.add_argument("--repeats",        type=int, default=1,            metavar="N")
    args = parser.parse_args()

    run_benchmark(args.input, args.threads, args.mem_per_thread, args.repeats)


if __name__ == "__main__":
    main()
