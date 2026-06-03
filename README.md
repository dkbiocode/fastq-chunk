## fastq-chunk

Apply a user function to fastq files in a multicore computing environment. 

## Installation 

### pip

You can install with pip, but not all dependencies may be resolved correctly. 

`python -m pip install --no-deps git+https://github.com/dkbiocode/fastq-chunk.git@main`

### As part of a conda environment 

Example environment 

```
name: fq-chunk
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.13
  - numpy
  - dnaio
  - pip:
    - git+https://github.com/dkbiocode/fastq-chunk.git

```

## How to use

This code is inspired by Fastp and other tools that divide input files into "chunks", perform computation on those chunks in separate threads/processes, write the results for the chunk in temporary files, then concatenate the files, preserving order. 

Note: Currently, the file streams all use gzip, so gzipped input and output files are required. 

### How to implement 

There are two primary functions to write:

* Worker function - takes a list of fastq records and writes to a temporary path. The path is the return value of the function. 
* Aggregation function - takes the returned paths and concatenates them to the final output file. This can be in `main()`. 

Return types can be extended to return multiple paths or data, for more sophisticated applications. 

## Examples

### Basic: read-through

```python
import fastq_chunk
import shutil
import os
import gzip

def readthrough(chunk, idx, *, tmpdir):
    temp_path = os.path.join(tmpdir, f"chunk_{chunk_idx:06d}.fastq.gz")
    with gzip.open(temp_path, "wt") as out:
        for read in chunk:
            # work here
            out.write(...)

    return temp_path

def main():
   # workers, memory 
   ... 
   # Input 
   fq_gz = ... 
   # get chunk size estimate from file
   chunksize = estimate_chunk_size(fq_gz) 

   # launch thread or process pool
    worker = functools.partial(readthrough, tempdir=tmpdir) 
   results = list(run_parallel(fq_gz, worker, chunksize, threads)) 

   # concatenate output files
   with open("out.fastq.gz", "wb") as outf:
       for pth in results:
          with open(pth, "rb") as fin:
              shutil.copyfileobj(fin,outf)
```