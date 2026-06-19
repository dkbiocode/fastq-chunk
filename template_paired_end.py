#!/usr/bin/env python3
import functools
import logging
import os
import sys
import tempfile
import re
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import dnaio
import io

from fastq_chunk import FastqRecord, get_read_dimensions, calculate_chunk_size, run_parallel_paired

wnv_pattern = re.compile("CT.AC.GT.AC.GT.AC.GC.GC.AC.CT.CT.")

def min_qual(s): # s - string of quality symbols, like "FFFFF:FFF,"
    q = list(map(ord,list(s))) # ord() returns the ASCII value for a character. See ASCII character table for more info.
    return min(q) - 33   # Phred33 is the ascii value minus 33.

def extract_codes_from_chunk_pair(chunk1_bytes: bytes, chunk2_bytes: bytes, chunk_ix:int,
    *, temp_dir: str,
    ) -> str:
    # search for viral barcode
    # if found, add entry for barcode, qualities, and read position info from the header
    # header info: instrument:run:flowcell:lane:tile:x:y

    # the format of the dataframe will contain all of the following as columns
    COLUMNS = ["viral_sequence", "cell_sequence", "umi_sequence",        # sequences
               "viral_quality",  "cell_quality",  "umi_quality",         # qualities
               "min_cell_quality","min_umi_quality","min_viral_quality", # min qualities
               "run","flowcell","lane","tile","x","y"                    # location
               ]
    data_rows = []

    with dnaio.open(io.BytesIO(chunk1_bytes), io.BytesIO(chunk2_bytes)) as reader:
        for read_ix, (r1, r2) in enumerate(reader):
            # ensure pairing
            assert r1.name == r2.name

            # search for virus
            match = wnv_pattern.search(r2.sequence)
            if match:
                s,e = match.start(), match.end()
            else:
                continue

            # get info from the sequences
            viral_sequence = r2.sequence[s:e]
            viral_quality = r2.qualities[s:e]
            # extract barcodes by given positions (0-based; end exclusive)
            cell_sequence = r1.sequence[:15]
            umi_sequence  = r1.sequence[15:26]
            cell_quality = r1.qualities[:15]
            umi_quality  = r1.qualities[15:26]

            # minimum qualities
            min_cell_quality = min_qual(cell_quality)
            min_umi_quality = min_qual(umi_quality)
            min_viral_quality = min_qual(viral_quality)

            # split header info
            try:
                position_info = r1.name.split()[1]
                instrument,run,flowcell,lane,tile,x,y = position_info.split(':')
            except ValueError:
                logging.error(f"Couldn't parse header in {chunk_ix=}; {read_ix=}; {r1.name}")
                raise

            data_rows.append(
                dict(
                    zip(COLUMNS, # list values in the same order as the COLUMNS
                    [
                        viral_sequence, cell_sequence, umi_sequence,        # sequences
                        viral_quality,  cell_quality,  umi_quality,         # qualities
                        min_cell_quality,min_umi_quality,min_viral_quality, # min qualities
                        run,flowcell,lane,tile,x,y                          # location
                    ] 
                    )
                )
            )
    out_path = os.path.join(temp_dir, f"chunk_{chunk_ix:06d}.tsv")
    pd.DataFrame(data_rows, columns=COLUMNS).to_csv(out_path, sep="\t")
    return out_path

def main() -> None:
    INPUT_PATH_R1 = sys.argv[1]
    INPUT_PATH_R2 = sys.argv[2]
    OUTPUT_PATH = sys.argv[3]
    N_WORKERS = int(sys.argv[4])
    MEM_PER_THREAD_MB = 3084

    dims = get_read_dimensions(INPUT_PATH_R1)
    if dims is None: 
        raise SystemExit(f"no reads in {INPUT_PATH_R1}")

    with tempfile.TemporaryDirectory(dir=os.getenv('SLURM_SCRATCH')) as tmp:
        worker = functools.partial(extract_codes_from_chunk_pair, temp_dir=tmp)
        chunk_paths = list(
            run_parallel_paired(INPUT_PATH_R1, INPUT_PATH_R2,
                    worker,
                    n_workers=N_WORKERS,
                    executor_class = ProcessPoolExecutor))


        # merge everything before deleting tmpdir
        chunk_dfs = [pd.read_csv(chunk_path) for chunk_path in chunk_paths]
        df = pd.concat(chunk_dfs, ignore_index=True)
        df.to_csv(OUTPUT_PATH, sep="\t", index=False)

    print(f"wrote {len(chunk_paths)} chunk(s) → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
