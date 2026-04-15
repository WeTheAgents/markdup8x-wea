# Proposal: markdup-wea MVP — Rust BAM Duplicate Marker

## Why

Picard MarkDuplicates is the 3rd most expensive step in nf-core/rnaseq (~9% of total pipeline runtime). It is Java-based, single-threaded, allocates 3.8-25GB JVM heap, and OOM'd on 5 of 8 human genome samples at 32GB RAM. The alternative, samtools markdup, requires a 4-step pipeline (sort-n, fixmate, sort, markdup) and is slower than Picard at t=1 on human data (5,070s vs 1,791s). sambamba is fast but has file descriptor limits and memory tuning requirements that are unacceptable for a pipeline with 1200+ users on diverse infrastructure (Docker, Singularity, HPC, cloud).

The nf-core maintainer has stated that any replacement must be a straight drop-in, not another option — the pipeline already has too many tool-choice parameters. We need a tool that is faster, uses less memory, and requires zero configuration.

## What Changes

- New standalone Rust binary `markdup-wea` that marks duplicate reads in coordinate-sorted BAM files
- Two-pass architecture: scan (build duplicate groups) then write (apply flags)
- Picard-compatible duplicate detection: same grouping key, same scoring, same FLAG 0x400 semantics
- Picard-compatible metrics output parseable by MultiQC without changes
- Multi-threaded BAM I/O via htslib thread pool
- Static binary with zero runtime dependencies (no JVM, no Python, no shared libs)
- **NOT in MVP**: optical duplicate detection (reports 0), UMI support, queryname-sorted input

## Capabilities

### New Capabilities

- `duplicate-detection`: Mark duplicate reads in coordinate-sorted BAM using Picard-compatible algorithm
- `paired-end-grouping`: Group paired-end reads by unclipped 5' positions of both mates with library-aware keying
- `single-end-grouping`: Group single-end and mate-unmapped reads by unclipped 5' position and orientation
- `quality-scoring`: Score reads by sum of base qualities >= 15 (Picard-compatible)
- `picard-metrics`: Output Picard-format DuplicationMetrics file with histogram, parseable by MultiQC
- `sort-order-validation`: Validate and enforce coordinate sort order at runtime
- `stdin-buffering`: Accept piped input by transparent temp-file buffering for two-pass processing
- `dupset-strategies`: Two dupset implementations (BitVec, FxHashSet) behind a trait for benchmarking

### Modified Capabilities

(none — this is a new tool)

### Impact

- New repository: `peachgabba22/markdup-wea`
- New binary: `markdup-wea` (static, <10MB)
- Integration target: nf-core/modules `BAM_MARKDUPLICATES` subworkflow (future PR)
- Test data dependency: nf-core/rnaseq test profile (yeast, 5 samples) and ENCODE human RNA-seq (8 samples)
