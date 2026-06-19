import io

import dnaio
import pytest

from fastq_chunk import run_parallel_paired


def _paired_identity_worker(r1_bytes, r2_bytes, chunk_idx):
    r1_count = sum(1 for _ in dnaio.open(io.BytesIO(r1_bytes)))
    r2_count = sum(1 for _ in dnaio.open(io.BytesIO(r2_bytes)))
    return (chunk_idx, r1_count, r2_count)


def test_paired_results_in_order(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    results = list(run_parallel_paired(r1, r2, _paired_identity_worker, n_workers=2))
    indices = [idx for idx, _, _ in results]
    assert indices == list(range(len(results)))


def test_paired_total_reads(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    results = list(run_parallel_paired(r1, r2, _paired_identity_worker, n_workers=2))
    assert sum(r1c for _, r1c, _ in results) == 5
    assert sum(r2c for _, _, r2c in results) == 5


def test_paired_r1_r2_equal_per_chunk(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    for _, r1c, r2c in run_parallel_paired(r1, r2, _paired_identity_worker, n_workers=2):
        assert r1c == r2c


def test_paired_single_worker(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    results = list(run_parallel_paired(r1, r2, _paired_identity_worker, n_workers=1))
    indices = [idx for idx, _, _ in results]
    assert indices == list(range(len(results)))


def test_paired_single_chunk(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    results = list(run_parallel_paired(
        r1, r2, _paired_identity_worker, buffer_size=1024 * 1024, n_workers=4,
    ))
    assert results == [(0, 5, 5)]


@pytest.mark.integration
def test_integration_paired_total_reads(sample_data_dir):
    results = list(run_parallel_paired(
        sample_data_dir / "sample_50k_R1.fastq.gz",
        sample_data_dir / "sample_50k_R2.fastq.gz",
        _paired_identity_worker,
        n_workers=2,
    ))
    assert sum(r1c for _, r1c, _ in results) == 50_000
    assert sum(r2c for _, _, r2c in results) == 50_000
