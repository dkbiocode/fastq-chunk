import gzip
import tracemalloc

import pytest
from fastq_chunk import get_read_dimensions


def test_returns_three_tuple(fastq_5reads):
    result = get_read_dimensions(fastq_5reads)
    assert isinstance(result, tuple) and len(result) == 3


def test_read_len_matches_fixture(fastq_5reads):
    read_len, _, _ = get_read_dimensions(fastq_5reads)
    assert read_len == 8  # fixture sequences are 8 bp


def test_bytes_per_read_positive(fastq_5reads):
    _, bytes_per_read, _ = get_read_dimensions(fastq_5reads)
    assert isinstance(bytes_per_read, int) and bytes_per_read > 0


def test_mem_per_read_positive(fastq_5reads):
    _, _, mem_per_read = get_read_dimensions(fastq_5reads)
    assert isinstance(mem_per_read, int) and mem_per_read > 0


def test_empty_file_returns_none(tmp_path):
    empty = tmp_path / "empty.fastq.gz"
    with gzip.open(empty, "wt"):
        pass
    assert get_read_dimensions(empty) is None


def test_n_profile_limits_reads(tmp_path):
    path = tmp_path / "ten.fastq.gz"
    with gzip.open(path, "wt") as f:
        for i in range(10):
            f.write(f"@read{i}\nACGTACGT\n+\nIIIIIIII\n")
    result = get_read_dimensions(path, n_profile=3)
    assert result is not None and len(result) == 3


def test_no_crash_when_tracemalloc_active(fastq_5reads):
    tracemalloc.start()
    try:
        result = get_read_dimensions(fastq_5reads)
        assert result is not None
    finally:
        tracemalloc.stop()


@pytest.mark.integration
def test_integration_sample_read_len(sample_data_dir):
    read_len, _, _ = get_read_dimensions(sample_data_dir / "sample_50k_R1.fastq.gz")
    assert read_len == 251
