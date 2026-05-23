#!/usr/bin/env python3
"""
Measure actual peak memory per read through degrade_chunk so that
calculate_chunk_size uses a measured value rather than a hand-estimated one.

Paired-end Illumina reads are constant length, so profiling one file gives
an accurate per-read cost for the entire run. Both R1 and R2 are measured
separately since they may have different cycle counts.
"""
import tracemalloc
import os
import dnaio
import numpy as np

from fastq_chunk import FastqRecord, get_read_dimensions
from degrade_fastq import DegradeParams, degrade_chunk

SAMPLE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_data")
FILES = {
    "R1": os.path.join(SAMPLE_DIR, "2629LEI-2_S2_L001_R1_001.fastq.gz"),
    "R2": os.path.join(SAMPLE_DIR, "2629LEI-2_S2_L001_R2_001.fastq.gz"),
}

# Conservative noise params — values don't meaningfully affect allocation
PARAMS = DegradeParams(
    tail_start=50, tail_slope=0.05, global_noise=1.0,
    dropout_rate=0.01, dropout_floor=2.0, seed=42
)

N_READS      = 10_000   # enough reads to average out per-read overhead reliably
MEM_MB       = 3840     # HPC task allocation per thread
SAFETY       = 0.6      # fraction of budget to actually use


def profile_file(label: str, path: str) -> None:
    dims = get_read_dimensions(path)
    if dims is None:
        print(f"{label}: could not read dimensions from {path}")
        return
    read_len, bytes_per_read, get_read_dim_mem = dims

    # Start tracing before loading records so the input buffer is included —
    # that's the dominant cost when sizing chunks for multithreaded use.
    tracemalloc.start()

    records: list[FastqRecord] = []
    with dnaio.open(path) as fin:
        for rec in fin:
            if len(records) >= N_READS:
                break
            records.append(FastqRecord(rec.name, rec.sequence, rec.qualities))

    actual_n = len(records)
    rng = np.random.default_rng(42)

    # degrade_and_write_chunk writes to disk as it goes, so the output is not
    # held in memory — consume the generator without storing results.
    for _ in degrade_chunk(records, PARAMS, rng):
        pass
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    measured_per_read = peak_bytes / actual_n
    usable_bytes      = MEM_MB * 1024 * 1024 * SAFETY
    chunk_size        = int(usable_bytes / measured_per_read)
    chunk_size        = max(100, min(chunk_size, 2_000_000))

    print(f"\n── {label} ──────────────────────────────")
    print(f"  read length             : {read_len} bp")
    print(f"  disk bytes/read         : {bytes_per_read:>8,}")
    print(f"  get_read_dimensions mem : {get_read_dim_mem:>8,.0f}  (load only, no worker)")
    print(f"  measured mem/read       : {measured_per_read:>8,.0f}  (load + degrade_chunk, {actual_n:,} reads)")
    print(f"  peak total ({actual_n:,} reads)  : {peak_bytes/1024/1024:>6.1f} MB")
    print(f"  chunk_size @ {MEM_MB} MB × {SAFETY}  : {chunk_size:>8,} reads")


for label, path in FILES.items():
    profile_file(label, path)
