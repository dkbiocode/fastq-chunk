import gzip
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent


def _write_fastq_gz(path: Path, records: list[tuple[str, str, str]]) -> None:
    with gzip.open(path, "wt") as fout:
        for name, seq, qual in records:
            fout.write(f"@{name}\n{seq}\n+\n{qual}\n")


@pytest.fixture
def fastq_5reads(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    _write_fastq_gz(path, [(f"read{i}", "ACGTACGT", "IIIIIIII") for i in range(1, 6)])
    return path


@pytest.fixture
def fastq_paired_5reads(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    _write_fastq_gz(r1, [(f"read{i}", "ACGTACGT", "IIIIIIII") for i in range(1, 6)])
    _write_fastq_gz(r2, [(f"read{i}", "TTTTGGGG", "HHHHHHHH") for i in range(1, 6)])
    return r1, r2


@pytest.fixture
def fastq_paired_mismatched(tmp_path):
    r1 = tmp_path / "r1_mismatch.fastq.gz"
    r2 = tmp_path / "r2_mismatch.fastq.gz"
    _write_fastq_gz(r1, [("readA", "ACGTACGT", "IIIIIIII")])
    _write_fastq_gz(r2, [("readB", "TTTTGGGG", "HHHHHHHH")])
    return r1, r2


@pytest.fixture(scope="session")
def sample_data_dir():
    return PROJECT_ROOT / "sample_data"


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "integration: mark test as requiring real sample data files (slow)"
    )
