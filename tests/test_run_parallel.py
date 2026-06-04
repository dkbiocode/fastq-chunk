import functools

import pytest
from fastq_chunk import run_parallel


def _identity_worker(chunk, chunk_idx):
    return (chunk_idx, len(chunk))


def _counting_worker(chunk, chunk_idx, *, call_log):
    call_log.append(chunk_idx)
    return chunk_idx


def test_results_in_submission_order(fastq_5reads):
    results = list(run_parallel(fastq_5reads, _identity_worker, chunk_size=2, n_workers=2))
    assert results == [(0, 2), (1, 2), (2, 1)]


def test_total_reads(fastq_5reads):
    results = list(run_parallel(fastq_5reads, _identity_worker, chunk_size=2, n_workers=2))
    assert sum(count for _, count in results) == 5


def test_single_worker_order(fastq_5reads):
    results = list(run_parallel(fastq_5reads, _identity_worker, chunk_size=2, n_workers=1))
    indices = [idx for idx, _ in results]
    assert indices == sorted(indices)


def test_fewer_chunks_than_workers(fastq_5reads):
    results = list(run_parallel(fastq_5reads, _identity_worker, chunk_size=5, n_workers=4))
    assert results == [(0, 5)]


def test_worker_call_count(fastq_5reads):
    call_log = []
    worker = functools.partial(_counting_worker, call_log=call_log)
    list(run_parallel(fastq_5reads, worker, chunk_size=2, n_workers=4))
    assert len(call_log) == 3  # ceil(5/2) == 3 chunks


def test_returns_iterator(fastq_5reads):
    result = run_parallel(fastq_5reads, _identity_worker, chunk_size=5, n_workers=2)
    assert hasattr(result, "__iter__") and hasattr(result, "__next__")


@pytest.mark.integration
def test_integration_total_reads(sample_data_dir):
    total = sum(
        count for _, count in
        run_parallel(sample_data_dir / "sample_50k_R1.fastq.gz",
                     _identity_worker, chunk_size=10_000, n_workers=2)
    )
    assert total == 50_000
