# Tasks

## 1. Project Skeleton (Phase 1)

- [ ] 1.1 Initialize Cargo.toml with rust-htslib (static), clap 4, anyhow, log, env_logger, rustc-hash, bitvec, tempfile
- [ ] 1.2 Configure release profile: lto=true, codegen-units=1, strip=true, opt-level=3
- [ ] 1.3 Create src/main.rs with clap CLI: INPUT, -o, -M, -@, --remove-duplicates, --assume-sort-order
- [ ] 1.4 Create src/lib.rs with public module re-exports
- [ ] 1.5 Create LICENSE (MIT)
- [ ] 1.6 Create README.md (project description, usage, build instructions)
- [ ] 1.7 Create .github/workflows/ci.yml (build, test, clippy, fmt on Linux/macOS/Windows)

### GATE: `cargo build --release` compiles. CI green.

## 2. Core Primitives (Phase 1)

- [ ] 2.1 Implement src/position.rs: `unclipped_5prime(record) -> i64` with forward/reverse/clip handling
- [ ] 2.2 Unit tests for position.rs: forward no-clip, forward soft-clip, reverse soft-clip, hard-clip, mixed, empty CIGAR, spliced alignment
- [ ] 2.3 Implement src/scoring.rs: `quality_sum(record) -> u32` (sum bases with Q >= 15)
- [ ] 2.4 Unit tests for scoring.rs: all-high-Q, all-low-Q, mixed, empty quals, edge at Q=15
- [ ] 2.5 Implement src/io.rs: open_reader (file + stdin detection), open_writer, validate_sort_order
- [ ] 2.6 Implement runtime sort-order check in io.rs: track (prev_ref_id, prev_pos), error on violation
- [ ] 2.7 Implement stdin temp-file buffering using tempfile crate
- [ ] 2.8 Unit tests for io.rs: valid SO, missing SO, wrong SO, runtime out-of-order, stdin detection

### GATE: All unit tests pass. `markdup-wea input.bam -o output.bam` copies BAM unchanged. Clippy clean.

## 3. Paired-End Duplicate Detection (Phase 2)

- [ ] 3.1 Implement src/pending_mates.rs: PendingMate struct (#[repr(C, packed)], 52 bytes), FxHashMap<u64, PendingMate>
- [ ] 3.2 Implement QNAME hashing: 64-bit FxHash primary, 32-bit check_hash (first 4 bytes of QNAME)
- [ ] 3.3 Implement hash collision detection: check_hash mismatch → treat as new pending mate
- [ ] 3.4 Implement src/groups.rs: PairedEndKey, SingleEndKey, ScoredPair, GroupTracker
- [ ] 3.5 Implement PairedEndKey ordering: lo/hi by (ref_id, unclipped_5prime) ascending
- [ ] 3.6 Implement incremental group resolution: BTreeMap<i64, Vec<GroupId>>, drain on stream advance
- [ ] 3.7 Implement chromosome boundary resolution: resolve all groups for prev chromosome
- [ ] 3.8 Implement EOF resolution: resolve all remaining groups (cross-chromosome, stragglers)
- [ ] 3.9 Implement src/scan.rs: Pass 1 orchestrator — classify records, route to grouping, return dup_bits
- [ ] 3.10 Implement record classification: unmapped, secondary, supplementary, paired-mate-unmapped, paired-both-mapped, single-end
- [ ] 3.11 Implement library detection from @RG headers (ID → LB mapping)
- [ ] 3.12 Implement single-end inline grouping (consecutive reads with same SingleEndKey)
- [ ] 3.13 Implement mate-unmapped routing to single-end grouping
- [ ] 3.14 Implement orphan read handling: warn to stderr at EOF, pass through unflagged
- [ ] 3.15 Implement src/markdup.rs: two-pass orchestrator (scan → dup_bits → reset → write)
- [ ] 3.16 Implement pre-existing flag clearing: clear FLAG 0x400 on input before grouping
- [ ] 3.17 Implement --remove-duplicates mode: skip records in dup_bits during Pass 2

### GATE: All unit tests pass.

## 4. DupSet Strategies (Phase 2)

- [ ] 4.1 Define DupSet trait: insert(record_id), contains(record_id) -> bool
- [ ] 4.2 Implement BitVecDupSet: bitvec::BitVec, grow-on-demand
- [ ] 4.3 Implement HashDupSet: FxHashSet<u64>
- [ ] 4.4 Wire strategy selection (feature flag or runtime flag)
- [ ] 4.5 Test: both strategies produce identical output on synthetic BAM

### GATE: Both strategies compile and pass identical-output test.

## 5. Synthetic BAM Integration Tests (Phase 2)

- [ ] 5.1 Create test fixture: two identical pairs at same position → one flagged
- [ ] 5.2 Create test fixture: pair with higher quality wins over lower
- [ ] 5.3 Create test fixture: single-end reads mixed with paired
- [ ] 5.4 Create test fixture: mate on different chromosome → both reads flagged
- [ ] 5.5 Create test fixture: mate unmapped → single-end treatment
- [ ] 5.6 Create test fixture: supplementary read of dup primary → NOT flagged
- [ ] 5.7 Create test fixture: pre-existing FLAG 0x400 → cleared and re-evaluated
- [ ] 5.8 Create test fixture: orphan read (no mate in BAM) → pass through with warning
- [ ] 5.9 Create test fixture: multiple libraries → separate grouping
- [ ] 5.10 Create test fixture: runtime sort-order violation → error exit
- [ ] 5.11 Wire all fixtures into tests/integration.rs with assertions on FLAG values

### GATE: All 10 synthetic test scenarios pass on both DupSet strategies.

## 6. Picard-Compatible Metrics (Phase 3)

- [ ] 6.1 Implement src/metrics.rs: MetricsCounters struct (per-library counters)
- [ ] 6.2 Implement counter accumulation during Pass 1 scan
- [ ] 6.3 Implement PERCENT_DUPLICATION calculation
- [ ] 6.4 Implement ESTIMATED_LIBRARY_SIZE via Lander-Waterman / Newton's method
- [ ] 6.5 Implement edge cases: zero dups, all dups, zero examined
- [ ] 6.6 Implement duplicate group size histogram tracking
- [ ] 6.7 Implement metrics file writer with exact Picard header format (TAB-separated, htsjdk headers)
- [ ] 6.8 Implement per-library rows when multiple libraries present
- [ ] 6.9 Unit tests: PERCENT_DUPLICATION formula, library size estimation, edge cases
- [ ] 6.10 Integration test: metrics file format matches Picard template character-by-character

### GATE: Metrics file is parseable by MultiQC (manual test with `multiqc .`).

## 7. Threading + Optimization (Phase 4)

- [ ] 7.1 Wire -@ flag to hts_set_threads() for BGZF decompression/compression
- [ ] 7.2 Implement record reuse: reader.read(&mut record) pattern in both passes
- [ ] 7.3 Pre-size FxHashMap from BAM file size estimate
- [ ] 7.4 Add logging: record count, pending_mates peak, groups resolved, wall time per pass

### GATE: `markdup-wea -@ 4` produces correct output. Performance log shows thread utilization.

## 8. Yeast Validation (Phase 5a)

- [ ] 8.1 Acquire yeast test BAMs: run `nf-core/rnaseq -profile test,docker --skip_markduplicates`
- [ ] 8.2 Document BAM acquisition in tests/README.md
- [ ] 8.3 Run markdup-wea on all 5 yeast samples, compare dup counts to ground truth:
  - WT_REP1: 36,810 / 180,342
  - WT_REP2: 11,688 / 90,962
  - RAP1_UNINDUCED_REP2: 78,929 / 98,201
  - RAP1_UNINDUCED_REP1: 36,294 / 48,977
  - RAP1_IAA_30M_REP1: 11,094 / 48,347
- [ ] 8.4 If any count differs: diff flagged reads (`samtools view -f 1024`), isolate divergence
- [ ] 8.5 Run MultiQC on markdup-wea metrics, verify duplication plot
- [ ] 8.6 Compare metrics values to Picard's metrics on same samples
- [ ] 8.7 Benchmark both DupSet strategies on all 5 samples: wall time (hyperfine), peak RSS (/usr/bin/time -v)
- [ ] 8.8 Pick DupSet winner. If negligible difference → keep BitVec.
- [ ] 8.9 Test stdin piping: `samtools view -h yeast.bam | markdup-wea - -o out.bam`
- [ ] 8.10 Test --remove-duplicates produces correct count
- [ ] 8.11 Verify static binary size < 10MB

### GATE: All 5 yeast samples match Picard. Metrics parseable by MultiQC. DupSet strategy decided. This gate MUST pass before proceeding to Phase 5b.

## 9. Human Genome Validation (Phase 5b)

- [ ] 9.1 Deploy markdup-wea to cloud benchmark server (existing infra in domains/rnaseq/benchmark/)
- [ ] 9.2 Run on all 8 ENCODE human RNA-seq samples (166M-230M reads each)
- [ ] 9.3 Compare dup counts to Picard ground truth for all 8 samples
- [ ] 9.4 If any count differs: investigate and document
- [ ] 9.5 Benchmark: Picard vs samtools (t=1,4,8) vs markdup-wea (t=1,4,8), 3 runs each
- [ ] 9.6 Measure peak RSS with /usr/bin/time -v
- [ ] 9.7 Log pending_mates peak count for memory analysis
- [ ] 9.8 Create benchmarks/ directory with scripts and CSV results
- [ ] 9.9 Publish performance comparison table in README

### GATE: All 8 human samples match Picard. Peak RSS < target. Performance competitive with Picard.

## 10. Release Preparation

- [ ] 10.1 Update README with benchmark table, usage examples, deviation table
- [ ] 10.2 Add CHANGELOG.md with v0.1.0 entry
- [ ] 10.3 Configure GitHub Actions release workflow: build static binaries for Linux x86_64, macOS arm64, Windows x86_64
- [ ] 10.4 Tag v0.1.0 and create GitHub Release with binaries
- [ ] 10.5 Publish to crates.io (cargo publish)

### GATE: v0.1.0 released. Binary downloadable. `cargo install markdup-wea` works.
