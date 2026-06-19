from __future__ import annotations
import collections
import io
import logging
import os
import tracemalloc
from concurrent.futures import Executor, ThreadPoolExecutor, ProcessPoolExecutor
from dataclasses import dataclass
from typing import Callable, Iterator, TypeVar, Tuple, cast

import dnaio
import xopen

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

def iter_chunks_paired(
    fastq_path1: str | os.PathLike,
    fastq_path2: str | os.PathLike,
    chunk_size: int,
) -> Iterator[Tuple[list[FastqRecord], list[FastqRecord]]]:
    """Yield (r1_chunk, r2_chunk) pairs from paired fastq files, keeping reads synchronized.

    Uses dnaio's paired mode, which validates read-name pairing on every record and
    raises immediately if the files are out of sync or have different read counts.
    """
    chunk1: list[FastqRecord] = []
    chunk2: list[FastqRecord] = []
    with dnaio.open(fastq_path1, fastq_path2) as fin:
        for rec1, rec2 in fin:
            chunk1.append(FastqRecord(rec1.name, rec1.sequence, rec1.qualities))
            chunk2.append(FastqRecord(rec2.name, rec2.sequence, rec2.qualities))
            if len(chunk1) >= chunk_size:
                yield chunk1, chunk2
                chunk1, chunk2 = [], []
    if chunk1:
        yield chunk1, chunk2


def iter_byte_chunks_paired(
    fastq_path1: str | os.PathLike,
    fastq_path2: str | os.PathLike,
    buffer_size: int = 4 * 1024 * 1024,
) -> Iterator[tuple[bytes, bytes]]:
    """Yield (r1_bytes, r2_bytes) chunks of raw FASTQ data, paired and synchronized.

    Each chunk is a bytes copy of an internal memoryview, safe to pass to worker
    processes. Workers parse the bytes themselves using dnaio.open(io.BytesIO(...)).
    """
    with xopen.xopen(fastq_path1, "rb", threads=1) as f1, \
         xopen.xopen(fastq_path2, "rb", threads=1) as f2:
        for chunk1_mv, chunk2_mv in dnaio.read_paired_chunks(cast(io.RawIOBase, f1), cast(io.RawIOBase, f2), buffer_size):
            yield bytes(chunk1_mv), bytes(chunk2_mv)


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
    executor_class: Callable[..., Executor] = ThreadPoolExecutor,
) -> Iterator[T]:
    """Dispatch chunks to a pool, yielding results in submission order.

    Sliding-window: at most n_workers chunks are in memory at once, so the
    caller's memory budget (used to calculate chunk_size) is respected.

    worker signature: (chunk: list[FastqRecord], chunk_idx: int) -> T
    Bind extra context (params, output paths, etc.) with functools.partial.

    executor_class: ThreadPoolExecutor (default, I/O-bound) or
        ProcessPoolExecutor (CPU-bound). With ProcessPoolExecutor, worker
        must be picklable — module-level functools.partial works; lambdas
        and closures do not.
    """
    with executor_class(max_workers=n_workers) as pool:
        pending: collections.deque = collections.deque()
        for idx, chunk in enumerate(iter_chunks(fastq_path, chunk_size)):
            pending.append(pool.submit(worker, chunk, idx))
            while len(pending) >= n_workers:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()

def run_parallel_paired(
    fastq_path1: str | os.PathLike,
    fastq_path2: str | os.PathLike,
    worker: Callable[[bytes, bytes, int], T],
    *,
    buffer_size: int = 4 * 1024 * 1024,
    n_workers: int = 4,
    executor_class: Callable[..., Executor] = ProcessPoolExecutor,
) -> Iterator[T]:
    """Dispatch paired-end byte chunks to a pool, yielding results in submission order.

    Sliding-window: at most n_workers chunk-pairs are in memory at once.

    worker signature: (r1_bytes: bytes, r2_bytes: bytes, chunk_idx: int) -> T
    Workers are responsible for parsing the bytes (e.g., with
    dnaio.open(io.BytesIO(r1_bytes), io.BytesIO(r2_bytes))).
    Bind extra context (params, output paths, etc.) with functools.partial.

    executor_class: ProcessPoolExecutor (default, CPU-bound) or
        ThreadPoolExecutor (for I/O-bound user code). With ProcessPoolExecutor,
        worker must be picklable — module-level functools.partial works;
        lambdas and closures do not.
    """
    with executor_class(max_workers=n_workers) as pool:
        pending: collections.deque = collections.deque()
        for idx, (chunk1_bytes, chunk2_bytes) in enumerate(
            iter_byte_chunks_paired(fastq_path1, fastq_path2, buffer_size)
        ):
            pending.append(pool.submit(worker, chunk1_bytes, chunk2_bytes, idx))
            while len(pending) >= n_workers:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()
