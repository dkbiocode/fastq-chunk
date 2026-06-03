# Benchmark guide

## What is being measured

`benchmark_fastq_chunk.py` exercises the full pipeline: chunked reading →
worker (gzip recompression to temp files) → cleanup. The reported MB/s is
wall-clock throughput including I/O and compression — not chunking overhead
in isolation. Keep this in mind when comparing absolute numbers.

---

## Executor choice: ThreadPoolExecutor vs ProcessPoolExecutor

The benchmark passes `executor_class` to `run_parallel` / `run_parallel_paired`.
The right choice depends on what the worker does:

| Worker type | Better executor | Why |
|---|---|---|
| gzip / zlib compression | `ThreadPoolExecutor` | zlib's C implementation **releases the GIL** during compression, so threads achieve real parallelism with no process overhead |
| Pure-Python CPU work | `ProcessPoolExecutor` | Bypasses the GIL for code that holds it throughout |
| NumPy / SciPy operations | `ThreadPoolExecutor` | NumPy releases the GIL for most array ops |

The passthrough worker used here calls `gzip.open`, which relies on zlib.
Because zlib releases the GIL, `ThreadPoolExecutor` already achieves true
parallelism — adding process-pool overhead only makes things marginally
slower (observed: ~3.3 MB/s with `ProcessPoolExecutor` vs ~3.4 MB/s with
`ThreadPoolExecutor` at 4 workers, 15k-read chunks).

For a **quality-degradation or barcode-analysis worker** written in pure
Python, `ProcessPoolExecutor` is the correct choice and will scale better
across cores.

---

## When ProcessPoolExecutor pays off

Process-pool startup is amortized across the lifetime of the pool (one
`with` block), so per-task overhead is mainly IPC serialization of the
chunk. The practical rule of thumb:

> **ProcessPoolExecutor wins when per-chunk work time ≳ 1–2 seconds.**

At the current benchmark settings (`--mem-per-thread 40`, chunk_size ~15k
reads), each chunk takes ~0.5 s — too small to overcome overhead. The fix
is more reads per chunk (raise `--mem-per-thread`) or more total chunks
(larger input file). Both are addressed below.

---

## Setting `--mem-per-thread` correctly

The benchmark default (40 MB) is artificially low — it exists only to force
multiple chunks out of the small sample files. For real use, pass the actual
per-thread memory budget for your environment:

**HPC:** your scheduler's QOS defines this. In SLURM it is typically
`--mem-per-cpu` (or `--mem / --ntasks`). Read it from your job script or ask
your sysadmin. Pass that number directly:
```bash
python testing/benchmark_fastq_chunk.py --mem-per-thread <QOS_mb> ...
```

**Workstation / laptop:** a reasonable starting point is
`(total_ram_mb × 0.7) / n_workers`. On a 32 GB machine with 4 workers that
is ~5,600 MB. Find total RAM with:
```bash
# macOS
sysctl -n hw.memsize | awk '{printf "%d MB\n", $1/1024/1024}'
# Linux
grep MemTotal /proc/meminfo
```

---

## Recommended input sizes

The benchmark only exercises parallelism when the input spans **multiple
chunks**. The chunk size is set by `--mem-per-thread`; use this formula to
find the minimum file size for meaningful results:

```
chunk_size  = (mem_per_thread_mb × 1024 × 1024 × 0.6) / mem_per_read_bytes
min_reads   = n_workers × chunk_size          # at least one chunk per worker
good_reads  = 4 × n_workers × chunk_size      # four full passes through the pool
```

For 251 bp Illumina reads, `mem_per_read` is ~1,620 bytes (run
`get_read_dimensions` on your data to get the exact figure). Example targets
at common budget levels with 4 workers:

| `--mem-per-thread` | chunk_size | min_reads (1 chunk/worker) | good_reads (4× pool) |
|---|---|---|---|
| 40 MB *(demo only)* | ~15k | ~60k | ~240k |
| 512 MB | ~190k | ~760k | ~3M |
| 3,840 MB *(typical HPC QOS)* | ~1.49M | ~6M | ~24M |

If your input produces fewer chunks than workers, the sliding-window
dispatcher never saturates the pool and the benchmark understates throughput.

The simplest way to build a larger test file is to concatenate the sample
files — see the next section.

---

## Generating a larger test file by concatenation

gzip files can be concatenated at the raw byte level — this is valid per
the gzip spec and is the same technique `degrade_fastq.py` uses to assemble
output without a decompress/recompress round-trip.

```bash
# 10× sample → ~500k reads, ~65 MB R1 + ~79 MB R2
for i in $(seq 10); do cat sample_data/sample_50k_R1.fastq.gz; done \
    > sample_data/sample_500k_R1.fastq.gz
for i in $(seq 10); do cat sample_data/sample_50k_R2.fastq.gz; done \
    > sample_data/sample_500k_R2.fastq.gz

# Benchmark (single-end)
python testing/benchmark_fastq_chunk.py \
    --input sample_data/sample_500k_R1.fastq.gz --repeats 3

# Benchmark (paired-end)
python testing/benchmark_fastq_chunk.py --paired \
    --input sample_data/sample_500k_R1.fastq.gz \
    --input-r2 sample_data/sample_500k_R2.fastq.gz --repeats 3
```

> **Note:** raw concatenation of the same file produces duplicate read
> names. This is fine for throughput benchmarking but not for downstream
> analysis tools that enforce name uniqueness.
