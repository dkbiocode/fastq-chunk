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
