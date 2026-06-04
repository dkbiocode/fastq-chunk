import pytest
from fastq_chunk import iter_chunks


def test_chunk_count_exact_divisor(fastq_5reads):
    assert len(list(iter_chunks(fastq_5reads, chunk_size=5))) == 1


def test_chunk_count_with_remainder(fastq_5reads):
    assert len(list(iter_chunks(fastq_5reads, chunk_size=3))) == 2


def test_last_chunk_smaller_than_chunk_size(fastq_5reads):
    chunks = list(iter_chunks(fastq_5reads, chunk_size=3))
    assert len(chunks[-1]) == 2


def test_total_reads_preserved(fastq_5reads):
    total = sum(len(c) for c in iter_chunks(fastq_5reads, chunk_size=2))
    assert total == 5


@pytest.mark.parametrize("chunk_size", [1, 2, 3])
def test_all_chunks_at_most_chunk_size(fastq_5reads, chunk_size):
    for chunk in iter_chunks(fastq_5reads, chunk_size):
        assert len(chunk) <= chunk_size


def test_record_order_preserved(fastq_5reads):
    chunks = list(iter_chunks(fastq_5reads, chunk_size=3))
    assert chunks[0][0].name == "read1"
    assert chunks[-1][-1].name == "read5"


def test_record_fields_are_strings(fastq_5reads):
    rec = next(iter_chunks(fastq_5reads, chunk_size=5))[0]
    assert isinstance(rec.name, str)
    assert isinstance(rec.sequence, str)
    assert isinstance(rec.qualities, str)


def test_chunk_size_larger_than_file(fastq_5reads):
    chunks = list(iter_chunks(fastq_5reads, chunk_size=999))
    assert len(chunks) == 1
    assert len(chunks[0]) == 5


def test_single_record_chunks(fastq_5reads):
    chunks = list(iter_chunks(fastq_5reads, chunk_size=1))
    assert len(chunks) == 5
    assert all(len(c) == 1 for c in chunks)
