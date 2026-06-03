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

* Worker function - takes (at least) a list of fastq records and the chunk index, and writes to a temporary path. The path is the return value of the function. 
* Aggregation function - takes the returned paths and concatenates them into the final output file. This can be in `main()`. 

Return types can be extended to return multiple paths or data for more sophisticated applications. 

#### Adding arguments to your worker function

The worker function must take a list of fastq objects (via dnaio) and an index, but you'll want to add your own parameters as well. These are done with keyword arguments and the module functools. 

Define your custom function arguments after the asterisk (*) as below:

```python
def your_func(chunks, idx, *, your_arg1, your_arg2)... 
```

Then use functools to partially call the function with your arguments, returning a new function that only takes the chunk of fastq and chunk index. 

```
worker = functools.partial(your_func, your_arg1 = someValue, your_arg2 = someOtherValue) 
```

The new function `worker` is now ready to be passed into `run_parallel`, where the chunk and chunk index will be set during parallel processing. 

```
results = list(fastq_chunk.run_parallel(infile_fq_gz, worker, chunksize, threads))  
```

## Examples

### Basic: read-through

```python
import shutil
import os
import gzip
import functools
import fastq_chunk 

def readthrough(chunk, idx, *, tmpdir):
    temp_path = os.path.join(tmpdir, f"chunk_{chunk_idx:06d}.fastq.gz")
    with gzip.open(temp_path, "wt") as out:
        for read in chunk:
            # work here
            out.write(...)

    return temp_path

def main():
   # parallelization params, such as hpc resources
    threads = ... # set by SLURM? 
    memory_per_thread # might be fixed on your HPC 
    scratchdir = ... # node-local storage if available 
   ... 
   # Input files(S) 
   fq_gz = ... 
   # get chunk size estimate from file
   chunksize = fastq_chunk.estimate_chunk_size(fq_gz) 

   with tempfile.TemporaryDirectory(dir = scratchdir) as tmpdir:
       # launch thread or process pool
       worker = functools.partial(readthrough, tempdir=tmpdir) 
       results = list(fastq_chunk.run_parallel(fq_gz, worker, chunksize, threads)) 

       # concatenate output files
       with open("out.fastq.gz", "wb") as outf:
           for pth in results:
              with open(pth, "rb") as fin:
                  shutil.copyfileobj(fin,outf)
```