import pytest
import dnaio
from fastq_chunk import iter_chunks_paired


def test_paired_chunk_count(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    assert len(list(iter_chunks_paired(r1, r2, chunk_size=3))) == 2


def test_paired_last_chunk_size(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    pairs = list(iter_chunks_paired(r1, r2, chunk_size=3))
    assert len(pairs[-1][0]) == 2
    assert len(pairs[-1][1]) == 2


def test_paired_total_reads(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    r1_total = sum(len(c1) for c1, _ in iter_chunks_paired(r1, r2, chunk_size=2))
    r2_total = sum(len(c2) for _, c2 in iter_chunks_paired(r1, r2, chunk_size=2))
    assert r1_total == 5
    assert r2_total == 5


def test_paired_r1_r2_same_length_per_chunk(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    for r1_chunk, r2_chunk in iter_chunks_paired(r1, r2, chunk_size=2):
        assert len(r1_chunk) == len(r2_chunk)


def test_paired_synchronization(fastq_paired_5reads):
    r1, r2 = fastq_paired_5reads
    for r1_chunk, r2_chunk in iter_chunks_paired(r1, r2, chunk_size=2):
        for r1_rec, r2_rec in zip(r1_chunk, r2_chunk):
            assert r1_rec.name == r2_rec.name


def test_paired_mismatch_raises(fastq_paired_mismatched):
    r1, r2 = fastq_paired_mismatched
    with pytest.raises(dnaio.FileFormatError):
        for _ in iter_chunks_paired(r1, r2, chunk_size=5):
            pass
