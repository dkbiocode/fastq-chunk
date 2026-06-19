#!/usr/bin/env python3
import logging
import io
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import dnaio
import pandas as pd

from fastq_chunk import get_read_dimensions, run_parallel_paired

wnv_pattern = re.compile("CT.AC.GT.AC.GT.AC.GC.GC.AC.CT.CT.")

COLUMNS = ["viral_sequence", "cell_sequence", "umi_sequence",
           "viral_quality",  "cell_quality",  "umi_quality",
           "min_cell_quality","min_umi_quality","min_viral_quality",
           "run","flowcell","lane","tile","x","y"]


def min_qual(s):
    return min(ord(c) for c in s) - 33


def extract_codes_from_chunk_pair(chunk1_bytes: bytes, chunk2_bytes: bytes, chunk_ix: int) -> list[dict]:
    data_rows = []

    with dnaio.open(io.BytesIO(chunk1_bytes), io.BytesIO(chunk2_bytes)) as reader:
        for read_ix, (r1, r2) in enumerate(reader):
            assert r1.name == r2.name

            match = wnv_pattern.search(r2.sequence)
            if not match:
                continue
            s, e = match.start(), match.end()

            viral_sequence = r2.sequence[s:e]
            viral_quality  = r2.qualities[s:e]
            cell_sequence  = r1.sequence[:15]
            umi_sequence   = r1.sequence[15:26]
            cell_quality   = r1.qualities[:15]
            umi_quality    = r1.qualities[15:26]

            try:
                position_info = r1.name.split()[1]
                instrument, run, flowcell, lane, tile, x, y = position_info.split(':')
            except ValueError:
                logging.error("Couldn't parse header in chunk_ix=%d read_ix=%d name=%r", chunk_ix, read_ix, r1.name)
                raise

            data_rows.append(dict(zip(COLUMNS, [
                viral_sequence, cell_sequence, umi_sequence,
                viral_quality,  cell_quality,  umi_quality,
                min_qual(cell_quality), min_qual(umi_quality), min_qual(viral_quality),
                run, flowcell, lane, tile, x, y,
            ])))

    return data_rows


def main() -> None:
    INPUT_PATH_R1 = sys.argv[1]
    INPUT_PATH_R2 = sys.argv[2]
    OUTPUT_PATH   = sys.argv[3]
    N_WORKERS     = int(sys.argv[4])

    dims = get_read_dimensions(INPUT_PATH_R1)
    if dims is None:
        raise SystemExit(f"no reads in {INPUT_PATH_R1}")

    all_rows = []
    for rows in run_parallel_paired(INPUT_PATH_R1, INPUT_PATH_R2,
                                    extract_codes_from_chunk_pair,
                                    n_workers=N_WORKERS,
                                    executor_class=ProcessPoolExecutor):
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows, columns=COLUMNS)
    df.to_csv(OUTPUT_PATH, sep="\t", index=False)
    print(f"wrote {len(df)} row(s) → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
