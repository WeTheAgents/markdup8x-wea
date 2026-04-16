# Performance — markdup-wea

This document records measured wall-clock and memory usage on real
production BAM files, captured during Phase D (threading) of the project.

All numbers below come from a single Hetzner cpx62 instance (16 vCPU AMD
EPYC, 30 GB RAM, NVMe-backed `/mnt` volume). Each measurement is a single
run; `/usr/bin/time -v` provided wall, RSS and CPU%.

## K562_REP1 — threading scaling

Input: `K562_REP1.sorted.bam` (8.7 GB, 156,141,236 records, 78,070,618
read pairs, 85.4% duplication rate). Output written, metrics emitted,
duplicate flags compared QNAME-by-QNAME against real Picard 3.4.0 output.

| `-@`  | Wall      | RSS    | CPU%  | Speedup vs `-@ 1` | only_ours / only_picard |
|------:|----------:|-------:|------:|------------------:|:-----------------------:|
|     1 |  9:47.01  | 325 MB | 213 % | 1.00 ×            | 0 / 0                   |
|     4 |  7:46.78  | 330 MB | 341 % | 1.26 ×            | 0 / 0                   |
|     8 |  7:51.84  | 327 MB | 420 % | 1.24 ×            | 0 / 0                   |

Compared to a synchronous single-threaded baseline (commit `38d8622`
before the multithreaded BGZF wiring) the same sample took **16:43.05
wall, 99 % CPU** — so just switching the BAM reader/writer from
`bgzf::Reader/Writer` to `bgzf::MultithreadedReader/Writer` (even with
only one BGZF worker) cuts wall by **1.71×**. With four workers the
wall is **2.16× lower** than the synchronous baseline.

### Interpretation

* **`-@ 4` is the sweet spot.** Going to `-@ 8` does not reduce wall —
  in fact it nudges it up by ~5 s while pulling +20 % more CPU. This is
  expected: BGZF decode/encode parallelism saturates around 3–4 workers
  on this storage, and beyond that the additional worker threads stall
  on the **sequential** scan loop (per-record CIGAR/quality parsing,
  `pending_mates` lookups, group-resolution state — all coordinate-order
  dependent and unsafe to parallelize).
* **`-@ 1` is not really single-threaded.** `MultithreadedReader::with_worker_count(1)`
  still spawns one background worker for BGZF decompression; combined
  with the writer's own worker that's main + 2 background = ~213 % CPU.
  The wall difference between the synchronous reader and `-@ 1` shows
  the cost of synchronously inlining BGZF decode in the scan path.
* **Memory is flat.** RSS stays around **325 MB** regardless of thread
  count. The bgzf workers share fixed-size queues and don't grow with
  thread count.

### vs Picard

Picard 3.4.0 (Java, `-Xmx16g`) on the same input — same machine, same
storage, same day:

* Wall **25:21.73**, RSS **17.6 GB**, CPU **112 %**, exit 0,
  `pair_dups = 66,655,937`.

Side by side:

| Tool                      | Wall    | RSS    | CPU%  | Speedup vs Picard |
|---------------------------|--------:|-------:|------:|------------------:|
| Picard 3.4.0 (Java)       | 25:21.73| 17.6 GB| 112 % |              1.0× |
| markdup-wea `-@ 1`        |  9:47.01|  325 MB| 213 % |          **2.59×**|
| markdup-wea `-@ 4`        |  7:46.78|  330 MB| 341 % |          **3.26×**|
| markdup-wea `-@ 8`        |  7:51.84|  327 MB| 420 % |          **3.22×**|

54× less peak memory across the board.

## Phase C — full-batch correctness on 8 ENCODE samples

Captured during Phase C validation (commit `fa3f49c`, single-threaded
binary — pre-D1 numbers). All eight samples produced QNAME-set
byte-identical output to real Picard 3.4.0 (zero one-sided diffs across
467,145,902 duplicates total).

| Sample        | Input | Ours wall | Picard wall | Ours RSS | Picard RSS |
|---------------|------:|----------:|------------:|---------:|-----------:|
| K562_REP1     |  8.7G | 16:47     | 25:22       | 318 MB   | 17.2 GB    |
| GM12878_REP1  |  9.0G | 18:09     | 26:44       | 207 MB   | 17.5 GB    |
| GM12878_REP2  |  9.1G | 18:24     | 27:06       | 359 MB   | 17.1 GB    |
| H1_REP1       |  12G  | 22:47     | 34:32       | 613 MB   | 16.8 GB    |
| H1_REP2       |  9.2G | 18:60     | 28:12       | 265 MB   | 16.9 GB    |
| K562_REP2     |  12G  | 21:13     | 31:56       | 309 MB   | 16.9 GB    |
| MCF7_REP1     |  13G  | 22:48     | 34:50       | 471 MB   | 16.9 GB    |
| MCF7_REP2     |  12G  | 23:30     | 35:37       | 480 MB   | 17.0 GB    |
| **Mean**      |  —    | **20:20** | **30:32**   | **378 MB**| **17.0 GB**|

These numbers used the synchronous reader/writer (pre-D1). With the new
default (`-@ 1` ≈ 1.7×), full-batch wall would drop to roughly
**1:25 hours** (vs Picard's 4:05 hours).

## Phase D — full-batch with `-@ 4` (post-threading)

Re-run of all 8 samples on the new binary (commit `31399f9`, with D1+D2+D4
landed) using `-@ 4` and `-o /dev/null` (metrics only — output BAM not
needed for this benchmark).

| Sample        | Input | Old wall | -@ 4 wall | Speedup |
|---------------|------:|---------:|----------:|--------:|
| K562_REP1     |  8.7G | 16:47    |  7:47     | 2.16×   |
| GM12878_REP1  |  9.0G | 18:09    |  8:29     | 2.14×   |
| GM12878_REP2  |  9.1G | 18:24    |  8:35     | 2.14×   |
| H1_REP1       |  12G  | 22:47    | 10:54     | 2.09×   |
| H1_REP2       |  9.2G | 18:60    |  8:53     | 2.13×   |
| K562_REP2     |  12G  | 21:13    | 10:36     | 2.00×   |
| MCF7_REP1     |  13G  | 22:48    | 10:46     | 2.12×   |
| MCF7_REP2     |  12G  | 23:30    | 10:52     | 2.16×   |
| **Mean**      |  —    | **20:20**| **9:36**  | **2.12×** |

Per-sample speedup is remarkably consistent (2.00–2.16×) across input
sizes from 8.7 G to 13 G — the BGZF I/O parallelism delivers the same
factor regardless of file size, confirming the bottleneck shifts to the
sequential scan loop equally for each sample.

**Full-batch totals**

| Tool                    | 8-sample wall | vs Picard |
|-------------------------|--------------:|----------:|
| Picard 3.4.0 (Java)     | **4:05 h**    | 1.0×      |
| markdup-wea (single)    | 2:43 h        | 1.50×     |
| markdup-wea (`-@ 4`)    | **1:17 h**    | **3.18×** |

## Phase D — MultiQC parity

All 8 sample-pairs (16 metrics files) parsed by MultiQC 1.33 with no
warnings; the `multiqc_picard_dups` table shows byte-identical values
across all ten columns for every (`*_wea`, `*_picard`) pair:

* `LIBRARY` = `Unknown Library` (D4 fallback matches Picard).
* `READ_PAIRS_EXAMINED`, `READ_PAIR_DUPLICATES`, `PERCENT_DUPLICATION`,
  `ESTIMATED_LIBRARY_SIZE` — exact match (D2 trimmed-float formatting
  matches Picard's `DecimalFormat("0.######")` precision).
* `READ_PAIR_OPTICAL_DUPLICATES = 0` — both (we don't compute optical;
  documented as Deviation #1).


## How to reproduce

The numbers in this document came from `/root/markdup-test/threading_bench.sh`
(threading scaling) and `/root/markdup-test/phasec_batch.sh` /
`/root/markdup-test/picard_batch.sh` (full-batch). Both are simple shell
loops around `/usr/bin/time -v` that capture wall, RSS, CPU and pipe
duplicate-QNAME extraction through `samtools view -f 1024 | sort -u`
before diffing with `comm`.
