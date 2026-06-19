import io

import dnaio
import pytest

from fastq_chunk import run_parallel


def _identity_worker(chunk_bytes, chunk_idx):
    count = sum(1 for _ in dnaio.open(io.BytesIO(chunk_bytes)))
    return (chunk_idx, count)


def test_results_in_order(fastq_5reads):
    results = list(run_parallel(fastq_5reads, _identity_worker, n_workers=2))
    indices = [idx for idx, _ in results]
    assert indices == list(range(len(results)))


def test_total_reads(fastq_5reads):
    results = list(run_parallel(fastq_5reads, _identity_worker, n_workers=2))
    assert sum(count for _, count in results) == 5


def test_single_worker_order(fastq_5reads):
    results = list(run_parallel(fastq_5reads, _identity_worker, n_workers=1))
    indices = [idx for idx, _ in results]
    assert indices == sorted(indices)


def test_large_buffer_total_reads(fastq_5reads):
    results = list(run_parallel(fastq_5reads, _identity_worker,
                                buffer_size=1024 * 1024, n_workers=4))
    assert sum(count for _, count in results) == 5


def test_returns_iterator(fastq_5reads):
    result = run_parallel(fastq_5reads, _identity_worker, n_workers=2)
    assert hasattr(result, "__iter__") and hasattr(result, "__next__")


@pytest.mark.integration
def test_integration_total_reads(sample_data_dir):
    total = sum(
        count for _, count in
        run_parallel(sample_data_dir / "sample_50k_R1.fastq.gz",
                     _identity_worker, n_workers=2)
    )
    assert total == 50_000
