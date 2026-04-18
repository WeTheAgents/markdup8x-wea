# Session Report: Picard 3.4.0 Parity Gate — 8/8 ENCODE Samples Passed

**Date:** 2026-04-17  
**Server:** Hetzner cpx62 (16 vCPU / 30 GB RAM), `204.168.253.225`  
**Branch:** `codex/picard-parity-proof-gate`  
**Reference:** Picard MarkDuplicates 3.4.0 (`e76128c283...`, OpenJDK 21.0.10)

---

## Summary

markdup-wea produces **byte-identical duplicate flags and metrics** to
Picard MarkDuplicates 3.4.0 across all 8 ENCODE RNA-seq samples
(1.55 billion records, 934 million duplicates). Zero per-record flag
mismatches. Zero metrics data divergence.

## Parity Results

| Sample | Records | Duplicates | Flag mismatches | Metrics match |
|--------|--------:|-----------:|:---------------:|:-------------:|
| K562_REP1 | 156,141,236 | 133,311,874 | 0 | YES |
| GM12878_REP1 | 166,920,844 | 125,959,040 | 0 | YES |
| GM12878_REP2 | 167,447,264 | 126,547,312 | 0 | YES |
| H1_REP1 | 223,421,774 | 111,561,954 | 0 | YES |
| H1_REP2 | 186,588,962 | 97,362,862 | 0 | YES |
| K562_REP2 | 195,837,680 | 160,964,674 | 0 | YES |
| MCF7_REP1 | 222,752,794 | 80,424,932 | 0 | YES |
| MCF7_REP2 | 230,760,862 | 98,159,156 | 0 | YES |
| **Total** | **1,549,871,416** | **934,291,804** | **0** | **8/8** |

## Performance (K562_REP1 head-to-head, same box, `-@ 4`)

| | Picard 3.4.0 | markdup-wea | Ratio |
|---|---:|---:|---:|
| Wall time | 25:00 | 7:54 | **3.2x faster** |
| Peak RSS | 8,568 MB | 323 MB | **26x less memory** |

## What was verified

- **Per-record flags:** every QNAME + FLAG pair compared across both
  BAM outputs, position-by-position. Zero mismatches on all 8 samples.
- **Metrics data:** all numeric columns (LIBRARY, READ_PAIRS_EXAMINED,
  READ_PAIR_DUPLICATES, PERCENT_DUPLICATION, ESTIMATED_LIBRARY_SIZE)
  and the full ROI histogram are byte-identical after excluding the
  tool-identity preamble (line 2: command-line, line 4: timestamp).
- **BAM header:** `@HD VN:1.6`, `SO:coordinate`, `@SQ`, `@RG` lines
  identical. `@PG ID:MarkDuplicates` present in both (our CL text
  differs — honest tool identity).
- **Per-record tags:** `PG:Z:MarkDuplicates` written on every record
  (matching Picard `ADD_PG_TAG_TO_READS=true` default). `DT` tags
  cleared (matching Picard `CLEAR_DT=true` default).

## What intentionally differs

These are inherent tool-identity differences, not correctness gaps:

1. `@PG CL` field: we write `markdup-wea 0.1.0 INPUT=... OUTPUT=...`,
   Picard writes its full 60-parameter command line.
2. Metrics line 2: tool name and arguments.
3. Metrics line 4: run timestamp (same format, different time).

## Picard invocation for reproducibility

```bash
java -Xmx4g -jar /opt/picard/picard.jar MarkDuplicates \
  INPUT=<sample>.sorted.bam \
  OUTPUT=<sample>.picard340.bam \
  METRICS_FILE=<sample>.picard340.metrics.txt \
  VALIDATION_STRINGENCY=LENIENT \
  ASSUME_SORTED=true \
  CLEAR_DT=true \
  ADD_PG_TAG_TO_READS=true \
  REMOVE_DUPLICATES=false \
  COMPRESSION_LEVEL=5
```

Note: nf-core/rnaseq additionally passes `REFERENCE_SEQUENCE=genome.fa`
and uses a larger heap (`-Xmx17g`). `REFERENCE_SEQUENCE` in Picard
MarkDuplicates is used only for CRAM input support — it does not affect
the duplicate-detection algorithm or flag output on BAM inputs.

## Fixes applied in this session

| Fix | File | Description |
|-----|------|-------------|
| `@HD VN:1.6` | `src/markdup.rs` | Bump SAM spec version to match Picard |
| `@PG MarkDuplicates` | `src/markdup.rs` | Add program header line, chained to last PG |
| `PG:Z` per-record tag | `src/markdup.rs` | Write `PG:Z:MarkDuplicates` on every record |
| Metrics blank line | `src/metrics.rs` | Add blank line after "Started on" preamble |
| Timestamp `UTC` | `src/metrics.rs` | Change zone label from GMT to UTC |
| Trailing blank line | `src/metrics.rs` | Add trailing newline at metrics EOF |
| Histogram cap | `src/metrics.rs` | CoverageMult = 0 for BIN > 100 (Picard behavior) |

## Known out-of-scope gaps (documented in `docs/deviations.md`)

- Optical duplicate detection (column always 0 — nf-core disables it)
- Multi-library metrics rows (single-library in all test data)
- `@PG CL` field text parity (inherent tool identity)

None of these affect duplicate flag output.
