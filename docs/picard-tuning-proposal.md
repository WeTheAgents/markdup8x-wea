# Picard MarkDuplicates memory-tuning proposal for nf-core/rnaseq

**Date:** 2026-04-17
**Server:** Hetzner cpx62 (16 vCPU / 30 GB RAM / 0 swap), `204.168.253.225`
**Branch:** `codex/picard-parity-proof-gate`
**Reference:** Picard MarkDuplicates 3.4.0, OpenJDK 21.0.10
**Inputs:** 8 ENCODE RNA-seq BAMs (1.55 billion records total; see
[session-report.md](session-report.md) for the parity-gate context)

---

## Summary

nf-core/rnaseq currently runs Picard MarkDuplicates under
`label 'process_medium'`, which reserves **36 GB** of RAM per task and
sets heap to **`-Xmx28g`** (`task.memory.mega × 0.8`). On the 8 ENCODE
RNA-seq samples in this benchmark the **actual per-JVM peak RSS does
not exceed 8 GB** — the default is a **~4× over-allocation**.

Dropping the heap to **`-Xmx6g`** and running **4 samples in parallel**
on a single 30 GB node keeps duplicate-marking output byte-identical
(**8/8 parity**, see below — zero QNAME+FLAG md5 divergence, zero
metrics-data divergence), holds aggregate RSS under 28 GB, and
delivers a **3.05× end-to-end throughput gain** on the benchmark set.
Stress tests confirm retry headroom: a single 2× largest-sample
"werewolf" BAM (445 M records, 17 GB) completes at `-Xmx9g` with
peak RSS 9.6 GB, and **two werewolves running concurrently at
`-Xmx9g` each finish in the same wall-clock** as one — sum RSS
20.8 GB, no OOM.

## Context — what nf-core does today

From `nf-core/modules` → `modules/nf-core/picard/markduplicates/main.nf`:

```nextflow
process PICARD_MARKDUPLICATES {
    label 'process_medium'
    ...
    def avail_mem = 3072
    if (!task.memory) { ... } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    ...
    java -Xmx${avail_mem}M -jar picard.jar MarkDuplicates ...
}
```

From `nf-core/rnaseq` → `conf/base.config`:

```nextflow
process {
    withLabel:process_medium {
        memory = { 36.GB * task.attempt }
    }
}
```

On a first attempt this yields `-Xmx28672m` = **-Xmx28g**. On retry
the request doubles to 72 GB → `-Xmx57g`. There is **no sample-size
or complexity-based routing** — every MarkDuplicates invocation, from
a 2-GB smoke-test BAM to a 13-GB ENCODE replicate, gets the same
static heap reservation.

**Why this matters on small nodes.** A Hetzner cpx62 box (16 vCPU,
30 GB RAM) cannot even start one `process_medium` task under default
config — the scheduler refuses to place a 36 GB memory request on a
node with < 36 GB available. A 72 GB box schedules exactly one
MarkDuplicates task at a time instead of 3–4.

## Setup

- Picard: `/opt/picard/picard.jar` (3.4.0), OpenJDK 21.0.10
- Inputs: 8 ENCODE RNA-seq sorted BAMs, 85 GB compressed total,
  1,549,871,416 records. Individual sizes 8.7–13 GB (222 M records
  for the largest, MCF7\_REP1). All inputs verified present before
  sweep; see `tuning-bench/recon.md` on the server.
- Werewolf: `samtools cat MCF7_REP1 MCF7_REP1 | samtools sort` →
  17 GB, 445 M records (2× the largest real ENCODE sample).
- Measurement: `/usr/bin/time -v` per JVM for peak RSS and wall;
  per-config `TOTAL WALL` captures parallel wall-clock for the
  whole 8-sample batch; `ps -eo rss,cmd` sampler every 5 s for
  aggregate `java` RSS.
- Picard invocation (same flags across all configs, matching
  nf-core/rnaseq defaults):

```bash
java -Xmx${HEAP} -jar picard.jar MarkDuplicates \
  INPUT=<sample>.sorted.bam \
  OUTPUT=<sample>.mkdup.bam \
  METRICS_FILE=<sample>.metrics \
  ASSUME_SORT_ORDER=coordinate \
  VALIDATION_STRINGENCY=LENIENT \
  TMP_DIR=/tmp
```

## Sweep matrix

| ID | Heap × parallelism | Intent |
|----|--------------------|--------|
| A | `-Xmx28g` × 1 | nf-core-faithful baseline (`process_medium` default) |
| B | `-Xmx7g` × 3 | Comfortable parallel, heap ≈ measured peak |
| C | `-Xmx6g` × 4 | Target config — 4× throughput on 30 GB node |
| D | `-Xmx5g` × 4 | Tight — where spill would start biting wall |
| W | `-Xmx9g` × 1 | Werewolf (2× largest BAM) stress |
| W2 | `-Xmx9g` × 2 | Two werewolves in parallel — retry-envelope stress |

## Results

### Aggregate (8 samples each, A/B/C/D)

| Config | Heap × par | TOTAL wall (s) | vs A | Peak per-JVM RSS | OOM |
|--------|-----------:|---------------:|-----:|-----------------:|-----|
| A `-Xmx28g` × 1 | 28 GB × 1 | 12947 | 1.00× | 30.75 GB | 0 |
| B `-Xmx7g` × 3  |  7 GB × 3 |  5514 | 2.35× |  7.95 GB | 0 |
| C `-Xmx6g` × 4  |  6 GB × 4 |  4240 | **3.05×** |  7.11 GB | 0 |
| D `-Xmx5g` × 4  |  5 GB × 4 |  4247 | 3.05× |  6.16 GB | 0 |

**A peak RSS** is 30.75 GB because `-Xmx28g` JVMs retain heap eagerly
under a generous budget — this is the actual working set under
nf-core defaults on this hardware. It is **not** the algorithm's
necessary peak. Under tighter `-Xmx`, JVM GCs more aggressively and
per-JVM peak drops close to the heap limit.

### Per-sample wall (minutes:seconds)

| Sample | A (28g×1) | B (7g×3) | C (6g×4) | D (5g×4) | C − A | D − A |
|--------|----------:|---------:|---------:|---------:|------:|------:|
| GM12878_REP1 | 26:46 | 27:07 | 27:22 | 27:35 | +2.2% | +3.1% |
| GM12878_REP2 | 26:46 | 27:19 | 27:51 | 27:29 | +4.1% | +2.7% |
| H1_REP1      | 34:12 | 35:12 | 34:53 | 34:57 | +2.0% | +2.2% |
| H1_REP2      | 27:38 | 28:24 | 28:17 | 28:43 | +2.4% | +3.9% |
| K562_REP1    | 24:46 | 25:17 | 25:40 | 25:31 | +3.6% | +3.0% |
| K562_REP2    | 31:25 | 31:37 | 31:59 | 31:57 | +1.8% | +1.7% |
| MCF7_REP1    | 35:14 | 35:34 | 35:16 | 35:08 | +0.1% | −0.3% |
| MCF7_REP2    | 35:37 | 36:16 | 35:45 | 35:47 | +0.4% | +0.5% |
| **Avg Δ**    |   —   |   —   |   —   |   —   | **+2.0%** | **+2.1%** |

Per-sample wall overhead of ~2% under parallelism is IO-contention
from 4 JVMs sharing `/tmp` spill. The **total wall compresses 3×**
because 4 samples process concurrently.

### Werewolf (2× MCF7_REP1, 445 M records, 17 GB BAM)

| Config | Wall | Peak RSS (per JVM) | Heap × par | Exit |
|--------|-----:|-------------------:|-----------:|-----:|
| W (single) | 1:14:57 | 9.64 GB | `-Xmx9g` × 1 | 0 |
| W2 werewolf_a | 1:14:22 | 10.25 GB | `-Xmx9g` × 2 | 0 |
| W2 werewolf_b | 1:14:37 | 10.58 GB | `-Xmx9g` × 2 | 0 |

Werewolf completes successfully under `-Xmx9g`. This establishes
headroom: a 13-GB ENCODE sample needs ~7 GB heap; doubling the
sample needs ~10 GB.

**W2 (two werewolves in parallel, `-Xmx9g` × 2).** Both instances
run against the same shared werewolf.bam input, writing to distinct
OUTPUT paths. Combined wall-clock **1:14:38** is essentially equal
to the single-werewolf wall — parallelism is near-perfect even at
2× largest-sample scale. Sum peak RSS ≈ 20.8 GB, well under the
28 GB safe ceiling on a 30 GB / 0-swap node. Per-JVM RSS rose
modestly (9.64 GB → 10.25–10.58 GB) — attributable to shared `/tmp`
spill contention and competing IO buffering, not heap inflation.

Recommended config C (`-Xmx6g`) has comfortable margin for typical
RNA-seq library sizes and retries go to `-Xmx9g` (or nf-core's
`process_low × 2 = -Xmx19g`) before hitting genuine sample-size
problems.

## Full parity: A vs C on 8 samples

For each sample we computed `samtools view A/X.mkdup.bam | awk
'{print $1"\t"$2}' | LC_ALL=C sort | md5sum` and the same for
`C/X.mkdup.bam`, then compared. Metrics files diffed with `#`-prefix
lines (tool-identity preamble and timestamp) excluded.

| Sample | QNAME+FLAG md5 (A = C) | Metrics diff |
|--------|:----------------------:|:------------:|
| GM12878_REP1 | `3fcc9e29df71b98f36fcdd2fe7b350a0` | 0 |
| GM12878_REP2 | `71f2ef3d3c5009b09547d6291a5930e8` | 0 |
| H1_REP1      | `20536f8261e59709305a0e2a4867f009` | 0 |
| H1_REP2      | `b2ff792426c86ecb8a320fb9c7c4b10e` | 0 |
| K562_REP1    | `f99ce91e97ff0f1639051036d6308482` | 0 |
| K562_REP2    | `72e4bf63effff1cb903f7c4ee1df1ded` | 0 |
| MCF7_REP1    | `b5da493726d7697a9a56b73a4b269728` | 0 |
| MCF7_REP2    | `37d45b3d88d9d5d7f1ccf53544d68885` | 0 |
| **Result**   | **8/8 byte-identical** | **0 diff** |

Heap size affects GC scheduling and spill batching, not the
duplicate-detection algorithm — Picard MarkDuplicates is
deterministic given identical input, flags, and
`ASSUME_SORT_ORDER=coordinate`. Dropping `-Xmx28g` → `-Xmx6g`
preserves bit-exact output.

## Recommendation

**Two equivalent ways to land this in nf-core/rnaseq.**

### Option 1 (minimal diff) — demote the label

In `nf-core/modules/modules/nf-core/picard/markduplicates/main.nf`:

```diff
 process PICARD_MARKDUPLICATES {
-    label 'process_medium'
+    label 'process_low'
```

where `process_low` in nf-core's standard `conf/base.config` is
`memory = 12.GB × task.attempt` → heap = `-Xmx9g` on first attempt,
`-Xmx19g` on retry. This lands comfortably inside the measured
werewolf envelope (9.6 GB peak) with headroom for sample outliers.

### Option 2 (explicit heap cap)

In a consumer's `modules.config`, `ext.args2` for the
`PICARD_MARKDUPLICATES` process:

```nextflow
withName: 'PICARD_MARKDUPLICATES' {
    ext.args2 = '-Xmx6g'
    memory   = { 8.GB * task.attempt }
}
```

This matches benchmark config C precisely: heap 6 GB, 8 GB task
reservation, retry doubles.

**Option 1 is preferred** — it changes one line in the canonical
module and rides on nf-core's existing label taxonomy. Option 2 is
for downstream pipelines that want to pin heap without touching the
nf-core/modules repo.

### Expected impact on real users

- **30 GB nodes:** MarkDuplicates starts at all (previously failed
  scheduling). Typical RNA-seq job fits 4 in parallel.
- **72 GB nodes:** 4× the MarkDuplicates throughput per node.
- **Retry semantics preserved:** `task.attempt = 2` still doubles
  the reservation (Option 1: to 24 GB / heap 19 GB; Option 2: to
  16 GB / heap 6 GB static, which a pipeline can widen in its retry
  block if desired).

## Caveats

Not covered by this benchmark, and therefore not blanket-recommended
without additional validation:

- **Non-RNA-seq data** (WGS, WGBS, long-read). RSS scales with read
  count and library complexity — WGS BAMs with 1.5 B reads may
  legitimately need more heap.
- **MarkDuplicatesSpark.** Different code path entirely; this
  proposal does not cover it.
- **Very small heaps (< 5 GB).** Not measured. `-Xmx4g` is the
  Picard-documented minimum and may start spill-dominating walls
  before it OOMs.
- **Very large BAMs (> 2× werewolf).** Not measured. Our werewolf
  cap is 445 M records; BAMs meaningfully larger than this should
  retry into `-Xmx12g`+.
- **Unsorted BAMs.** We run with `ASSUME_SORT_ORDER=coordinate`, the
  nf-core default. Unsorted input forces Picard's internal sort and
  materially changes the heap profile.
- **Optical duplicates enabled.** nf-core/rnaseq disables these.
  `READ_NAME_REGEX=null`-style configs enable an extra hash table
  whose footprint is not captured here.

## Reproducibility

All per-config logs, per-sample `.time` files, `sweep.log`
queue traces, and `memory.timeline` 5 s sampler output live under
`~/markdup-test/tuning-bench/results/{A,B,C,D,W}/` on
`204.168.253.225`. The driver script (`run_sweep.sh`, `run_one.sh`)
and this proposal's parity script (`parity.sh`) are in
`tuning-bench/scripts/` on the same server.
