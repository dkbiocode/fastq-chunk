#!/usr/bin/env python3
"""
Benchmark chunked parallel fastq read/split/write using fastq_chunk.

Single-end (default):
    python benchmark_fastq_chunk.py [--threads N] [--mem-per-thread MB]
                                    [--input FILE] [--repeats N]

Paired-end:
    python benchmark_fastq_chunk.py --paired [--input FILE] [--input-r2 FILE]
                                    [--threads N] [--mem-per-thread MB] [--repeats N]

Note: MB/s figures include gzip recompression in the worker and reflect
overall pipeline throughput, not chunking overhead alone. Run 1 may be
slower than later runs due to OS page-cache warm-up.
"""
import argparse
import functools
import gzip
import logging
import os
import shutil
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

from fastq_chunk import (
    FastqRecord,
    get_read_dimensions,
    calculate_chunk_size,
    run_parallel,
    run_parallel_paired,
)

logging.basicConfig(level=logging.WARNING)


# ---------------------------------------------------------------------------
# Workers
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


def passthrough_chunk_paired(
    r1_chunk: list[FastqRecord],
    r2_chunk: list[FastqRecord],
    chunk_idx: int,
    *,
    temp_dir: str,
) -> tuple[str, str]:
    """Write R1 and R2 chunks to numbered temp gzip files, unchanged.

    temp_dir is keyword-only so functools.partial can bind it, leaving
    (r1_chunk, r2_chunk, chunk_idx) as the worker contract run_parallel_paired expects.
    """
    r1_path = os.path.join(temp_dir, f"chunk_{chunk_idx:06d}_R1.fastq.gz")
    r2_path = os.path.join(temp_dir, f"chunk_{chunk_idx:06d}_R2.fastq.gz")
    for path, chunk in ((r1_path, r1_chunk), (r2_path, r2_chunk)):
        with gzip.open(path, 'wt') as fout:
            for rec in chunk:
                fout.write(f"@{rec.name}\n{rec.sequence}\n+\n{rec.qualities}\n")
    return r1_path, r2_path


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _print_stats(elapsed_times: list[float], input_mb: float) -> None:
    if len(elapsed_times) < 2:
        return
    min_t  = min(elapsed_times)
    mean_t = sum(elapsed_times) / len(elapsed_times)
    max_t  = max(elapsed_times)
    print(f"  min   {min_t:.2f} s   {input_mb / min_t:.1f} MB/s")
    print(f"  mean  {mean_t:.2f} s   {input_mb / mean_t:.1f} MB/s")
    print(f"  max   {max_t:.2f} s   {input_mb / max_t:.1f} MB/s")
    print()


# ---------------------------------------------------------------------------
# Benchmark drivers
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

    print(f"\n=== single-end: threads={threads}, mem_per_thread={mem_per_thread_mb} MB ===")
    print(f"Input:       {os.path.basename(input_path)}  ({input_mb:.1f} MB)")
    print(f"read_len:    {read_len} bp  |  mem_per_read: {mem_per_read:,} bytes  |  chunk_size: {chunk_size:,} reads")
    print(f"(MB/s includes gzip recompression; run 1 may be slower due to OS cache warm-up)")
    print()

    elapsed_times: list[float] = []

    for run in range(1, repeats + 1):
        working_dir = tempfile.mkdtemp(prefix="bench_")
        try:
            worker = functools.partial(passthrough_chunk, temp_dir=working_dir)
            t0 = time.perf_counter()
            temp_paths = list(
                run_parallel(input_path, worker, chunk_size=chunk_size, n_workers=threads,
                             executor_class=ProcessPoolExecutor)
            )
            elapsed = time.perf_counter() - t0
        finally:
            shutil.rmtree(working_dir)

        elapsed_times.append(elapsed)
        print(f"  Run {run}/{repeats}:  {len(temp_paths)} chunk(s)   {elapsed:.2f} s  ({input_mb / elapsed:.1f} MB/s)")

    print()
    _print_stats(elapsed_times, input_mb)


def run_benchmark_paired(
    r1_path: str,
    r2_path: str,
    threads: int,
    mem_per_thread_mb: int,
    repeats: int,
) -> None:
    r1_mb = os.path.getsize(r1_path) / 1024 / 1024
    r2_mb = os.path.getsize(r2_path) / 1024 / 1024
    total_mb = r1_mb + r2_mb

    dims = get_read_dimensions(r1_path)
    if dims is None:
        raise ValueError(f"could not read any records from {r1_path}")
    read_len, bytes_per_read, mem_per_read = dims
    chunk_size = calculate_chunk_size(mem_per_read, mem_per_thread_mb)

    print(f"\n=== paired-end: threads={threads}, mem_per_thread={mem_per_thread_mb} MB ===")
    print(f"R1:          {os.path.basename(r1_path)}  ({r1_mb:.1f} MB)")
    print(f"R2:          {os.path.basename(r2_path)}  ({r2_mb:.1f} MB)")
    print(f"read_len:    {read_len} bp  |  mem_per_read: {mem_per_read:,} bytes  |  chunk_size: {chunk_size:,} reads")
    print(f"(MB/s is R1+R2 combined; includes gzip recompression; run 1 may be slower due to OS cache warm-up)")
    print()

    elapsed_times: list[float] = []

    for run in range(1, repeats + 1):
        working_dir = tempfile.mkdtemp(prefix="bench_paired_")
        try:
            worker = functools.partial(passthrough_chunk_paired, temp_dir=working_dir)
            t0 = time.perf_counter()
            results = list(
                run_parallel_paired(r1_path, r2_path, worker, chunk_size=chunk_size, n_workers=threads,
                                    executor_class=ProcessPoolExecutor)
            )
            elapsed = time.perf_counter() - t0
        finally:
            shutil.rmtree(working_dir)

        elapsed_times.append(elapsed)
        print(f"  Run {run}/{repeats}:  {len(results)} chunk pair(s)   {elapsed:.2f} s  ({total_mb / elapsed:.1f} MB/s combined)")

    print()
    _print_stats(elapsed_times, total_mb)


def main() -> None:
    project_root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    sample_dir = os.path.join(project_root, "sample_data")
    default_r1 = os.path.join(sample_dir, "sample_50k_R1.fastq.gz")
    default_r2 = os.path.join(sample_dir, "sample_50k_R2.fastq.gz")

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--threads",        type=int, default=4,           metavar="N")
    parser.add_argument("--mem-per-thread", type=int, default=40,        metavar="MB")
    parser.add_argument("--input",          type=str, default=default_r1,  metavar="FILE",
                        help="R1 (or single-end) input file")
    parser.add_argument("--input-r2",       type=str, default=None,        metavar="FILE",
                        help="R2 input; triggers paired-end benchmark")
    parser.add_argument("--paired",         action="store_true",
                        help="Run paired-end benchmark using default sample R1/R2")
    parser.add_argument("--repeats",        type=int, default=3,           metavar="N")
    args = parser.parse_args()

    if args.paired and args.input_r2 is None:
        args.input_r2 = default_r2

    if args.input_r2 is not None:
        run_benchmark_paired(args.input, args.input_r2, args.threads, args.mem_per_thread, args.repeats)
    else:
        run_benchmark(args.input, args.threads, args.mem_per_thread, args.repeats)


if __name__ == "__main__":
    main()
