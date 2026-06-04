import pytest
from fastq_chunk import calculate_chunk_size


def test_typical_hpc():
    assert calculate_chunk_size(1620, 3840) == 1_491_308


def test_demo_low_budget():
    assert calculate_chunk_size(1620, 40) == 15_534


def test_lower_clamp():
    assert calculate_chunk_size(999_999_999, 3840) == 100


def test_upper_clamp():
    assert calculate_chunk_size(1, 3840) == 2_000_000


def test_zero_budget():
    assert calculate_chunk_size(1620, 0) == 100


def test_return_type_is_int():
    assert isinstance(calculate_chunk_size(1620, 3840), int)


@pytest.mark.parametrize("mem_per_read,budget", [
    (1, 1),
    (1_000_000, 1),
    (500, 100),
    (1620, 3840),
    (999_999_999, 3840),
])
def test_result_always_in_bounds(mem_per_read, budget):
    assert 100 <= calculate_chunk_size(mem_per_read, budget) <= 2_000_000
