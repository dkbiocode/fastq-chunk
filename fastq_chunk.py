from __future__ import annotations
import collections
import logging
import os
import tracemalloc
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Callable, Iterator, TypeVar

import dnaio

logger = logging.getLogger(__name__)

T = TypeVar('T')


@dataclass
class FastqRecord:
    name: str
    sequence: str
    qualities: str


def iter_chunks(
    fastq_path: str | os.PathLike,
    chunk_size: int,
) -> Iterator[list[FastqRecord]]:
    """Yield consecutive fixed-size lists of reads from a single fastq file.

    For paired-end coordination, zip two iter_chunks calls at the caller.
    """
    chunk: list[FastqRecord] = []
    with dnaio.open(fastq_path) as fin:
        for rec in fin:
            chunk.append(FastqRecord(rec.name, rec.sequence, rec.qualities))
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
    if chunk:
        yield chunk


def get_read_dimensions(
    fastq_path: str | os.PathLike,
    n_profile: int = 1000,
) -> tuple[int, int, int] | None:
    """Sample up to n_profile reads to measure read length, disk footprint per
    read, and actual in-memory cost per read (chunk buffer overhead).

    n_profile > 1 averages across variable-length reads for heterogeneous data.
    The peak is dominated by holding the chunk buffer, not processing overhead.
    """
    already_tracing = tracemalloc.is_tracing()
    if not already_tracing:
        tracemalloc.start()

    records: list[FastqRecord] = []
    read_len = 0
    bytes_per_read = 0

    with dnaio.open(fastq_path) as fin:
        for rec in fin:
            if not records:
                read_len = len(rec.sequence)
                header_len = len(rec.name)
                bytes_per_read = (
                    header_len + 1 +
                    read_len + 1 +
                    2 +
                    read_len + 1
                )
            records.append(FastqRecord(rec.name, rec.sequence, rec.qualities))
            if len(records) >= n_profile:
                break

    if not records:
        if not already_tracing:
            tracemalloc.stop()
        return None

    _, peak_bytes = tracemalloc.get_traced_memory()
    if not already_tracing:
        tracemalloc.stop()

    mem_per_read = peak_bytes // len(records)
    logger.info(
        "profiled %d reads from %s: read_len=%d bp, disk=%d bytes/read, mem=%d bytes/read",
        len(records), os.path.basename(fastq_path), read_len, bytes_per_read, mem_per_read,
    )
    return read_len, bytes_per_read, mem_per_read


def calculate_chunk_size(
    mem_per_read: int,
    mem_per_thread_mb: int = 3840,
    safety_factor: float = 0.6,
) -> int:
    """Return the number of reads that safely fit within one thread's memory budget."""
    usable_bytes = mem_per_thread_mb * 1024 * 1024 * safety_factor
    raw = int(usable_bytes / mem_per_read)
    chunk_size = max(raw, 100)
    chunk_size = min(chunk_size, 2_000_000)
    logger.info(
        "chunk_size=%d (mem_per_read=%d bytes, budget=%d MB x %.1f = %.0f MB usable, raw=%d%s)",
        chunk_size, mem_per_read, mem_per_thread_mb, safety_factor,
        usable_bytes / 1024 / 1024, raw,
        f", capped at {chunk_size:,}" if raw != chunk_size else "",
    )
    return chunk_size


def run_parallel(
    fastq_path: str | os.PathLike,
    worker: Callable[[list[FastqRecord], int], T],
    *,
    chunk_size: int,
    n_workers: int = 4,
) -> Iterator[T]:
    """Dispatch chunks to worker threads, yielding results in submission order.

    Sliding-window: at most n_workers chunks are in memory at once, so the
    caller's memory budget (used to calculate chunk_size) is respected.

    worker signature: (chunk: list[FastqRecord], chunk_idx: int) -> T
    Bind extra context (params, output paths, etc.) with a lambda or partial.
    """
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        pending: collections.deque = collections.deque()
        for idx, chunk in enumerate(iter_chunks(fastq_path, chunk_size)):
            pending.append(pool.submit(worker, chunk, idx))
            while len(pending) >= n_workers:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()
