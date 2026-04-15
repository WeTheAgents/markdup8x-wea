# Design: markdup-wea MVP

## Context

We are building a Rust-based BAM duplicate marker as a drop-in replacement for Picard MarkDuplicates in the nf-core/rnaseq pipeline. The tool must produce identical duplicate counts to Picard on all test data, use <200MB RSS (vs Picard's 21.5GB), and require zero configuration.

Benchmark data from our RNA-seq pipeline profiling (5 yeast samples, 8 human genome samples) provides ground truth for correctness validation and performance targets.

## Goals / Non-Goals

**Goals:**
- Identical dup counts to Picard on all 13 test samples
- Peak RSS < 200MB (BitVec) or < 400MB (FxHashSet) on 200M-read human BAMs
- Wall time < 1,500s at t=1 and < 900s at t=4 on human data
- Static binary < 10MB, zero runtime dependencies
- Picard-format metrics parseable by MultiQC
- Stdin piping support for pipeline integration

**Non-Goals:**
- Optical duplicate detection (MVP outputs 0)
- UMI/molecular barcode support
- Queryname-sorted input support
- Supplementary read flagging (differs from Picard, matches samtools)
- Parallel group resolution (core stays single-threaded)
- Disk spill for pending_mates (out of scope; add if >500MB in practice)

## Decisions

### Decision 1: Two-pass architecture
Using two-pass (scan then write) because:
- In coordinate-sorted BAM, read1 and read2 appear at distant positions
- We need both mates' CIGARs to compute unclipped 5' (BAM MPOS is aligned, not unclipped)
- Single-pass would require buffering all BAM records per chromosome (~1-2GB for chr1)
- Two-pass uses only metadata (~40-52 bytes per pending mate) — 100x less memory
- Considered: single-pass with chromosome buffering (too much memory), write-then-seek-back (breaks stdout), fixmate preprocessing (defeats simplicity goal)

### Decision 2: rust-htslib (not noodles)
Using rust-htslib 0.47 with static feature because:
- C bindings to htslib are battle-tested (samtools, bcftools, millions of users)
- BGZF threading via `hts_set_threads()` is mature
- Static linking produces self-contained binary
- Considered: noodles (pure Rust, but less mature BAM write path, no built-in BGZF threading)

### Decision 3: QNAME hashing with collision detection
Using 64-bit FxHash + 32-bit check hash because:
- Full QNAME storage for 200M reads = ~8GB (40 bytes avg QNAME)
- 64-bit hash: collision probability = ~2.7 * 10^-12 per pair at 200M reads (birthday paradox: n^2 / 2^65)
- 32-bit check hash catches the rare collision (combined false match: ~6 * 10^-22)
- If check_hash mismatch detected: treat as new pending mate (safe fallback)
- Considered: storing full QNAMEs (too much memory), 128-bit hash (overkill, slower)

### Decision 4: Incremental group resolution
Using BTreeMap<position, Vec<GroupId>> to resolve groups as stream advances because:
- Holding all groups until EOF wastes memory (100M pairs = ~8GB of group tracker)
- RNA-seq inserts are tight (~150-500bp) so groups resolve quickly
- BTreeMap allows efficient "drain all groups with max_pos <= current_pos"
- Chromosome boundaries force full resolution of prior chromosome
- Cross-chromosome pairs stay pending until mate's chromosome (rare, <1%)
- Considered: resolve only at chromosome boundary (holds too many groups), resolve at EOF (worst case memory)

### Decision 5: DupSet dual strategy
Implementing both BitVec and FxHashSet behind a trait because:
- BitVec: 25MB for 200M reads (best memory)
- FxHashSet: ~230MB for 20M dups (simpler, no record_id tracking needed)
- Benchmarking on yeast will reveal if the memory difference matters at small scale
- Human genome benchmark will stress-test the winner
- Considered: only BitVec (premature optimization), only FxHashSet (may bust memory budget)

### Decision 6: PendingMate struct packing
Using `#[repr(C, packed)]` for 52-byte PendingMate because:
- Without packing, Rust padding inflates to 64-72 bytes
- At 2M pending mates: 52 bytes = 104MB vs 72 bytes = 144MB (30% waste)
- Packed struct access is slightly slower on some architectures (unaligned loads)
- RNA-seq pending set is small enough that the speed penalty is negligible
- Considered: natural alignment (simpler, but 30% more memory), arena allocation (complex)

## Data Flow

```
INPUT (file or stdin)
    │
    ├─ stdin? ──► temp file ──┐
    │                         │
    ▼                         ▼
Pass 1: SCAN ◄── seekable file
    │
    ├── For each record:
    │   ├── Assign record_id (sequential counter)
    │   ├── Classify (unmapped/secondary/supplementary/paired/single-end)
    │   ├── Single-end → inline group resolution
    │   └── Paired-end:
    │       ├── Mate in pending? → form PairedEndKey → add to GroupTracker
    │       └── Mate not found? → insert into pending_mates
    │
    ├── Incremental resolution: when stream advances past group's max_pos
    │   └── Best score wins → losers' record_ids → dup_bits
    │
    ├── Chromosome boundary → resolve all groups for prev chrom
    │
    └── EOF → resolve all remaining groups (cross-chrom, stragglers)
         └── Result: dup_bits (BitVec or FxHashSet)

Pass 2: WRITE ◄── same seekable file (reset to start)
    │
    ├── For each record:
    │   ├── Increment record_id counter
    │   ├── If record_id in dup_bits → FLAG |= 0x400
    │   ├── If --remove-duplicates and dup → skip
    │   └── Write to output
    │
    └── Output BAM (file or stdout)

METRICS (if -M specified)
    │
    └── Write Picard-format DuplicationMetrics
        ├── Header (htsjdk format)
        ├── Data rows (per library)
        ├── ESTIMATED_LIBRARY_SIZE (Lander-Waterman / Newton's method)
        └── Histogram (group sizes)
```

## Memory Budget (200M reads, human RNA-seq)

| Component | BitVec mode | FxHashSet mode |
|-----------|-------------|----------------|
| dup_bits / dupset | 25 MB | 230 MB |
| pending_mates (peak ~2M, RNA-seq) | 130 MB | 130 MB |
| GroupTracker (active ~100K) | 8 MB | 8 MB |
| BAM reader buffers | 4 MB | 4 MB |
| Library table, counters | < 1 MB | < 1 MB |
| **Total** | **~168 MB** | **~373 MB** |

## Risks / Trade-offs

**Risk: Dup counts don't match Picard**
Impact: fatal — the tool is useless if correctness fails.
Mitigation: yeast validation gate before human genome. Diff individual flagged reads to isolate divergence. Test each edge case (cross-chrom, mate-unmapped, soft-clips) in isolation.

**Risk: Two-pass is slower than Picard's single pass**
Impact: medium — defeats the performance argument.
Mitigation: BGZF caching means second read is partially in OS page cache. htslib I/O threads offset the double-read overhead. If still slow on human data, investigate memory-mapped read for Pass 2.

**Risk: Pending mates exceed memory on WGS data**
Impact: low for MVP (RNA-seq only) but blocks future WGS support.
Mitigation: log peak pending_mates count. Add `--max-pending` with disk spill in post-MVP. For RNA-seq, insert sizes are tight (150-500bp) so pending set stays small.

**Risk: rust-htslib build fails on CI or user machines**
Impact: medium — blocks adoption.
Mitigation: static feature bundles htslib. CI builds release binaries for Linux/macOS/Windows. Users download prebuilt binaries from GitHub Releases.

**Risk: MultiQC doesn't parse our metrics**
Impact: medium — breaks pipeline integration.
Mitigation: test exact format against MultiQC 1.19+. The critical line is `## METRICS CLASS\tpicard.sam.markduplicates.DuplicationMetrics` — must be character-identical.

## Open Questions

1. **Tie-breaking:** When two pairs have identical combined_score, which wins? Picard uses "first encountered" (index-based). We should match this — use the pair with the lower record_id_1. Verify on test data.

2. **ESTIMATED_LIBRARY_SIZE without optical adjustment:** Picard subtracts optical dups before computing library size. We report 0 optical dups → our library size estimate will be slightly lower than Picard's when optical dups exist. Document this as known deviation. For standard Illumina flowcells (our test data), optical = 0, so estimates should match.

3. **Pre-existing duplicate flags:** Should we clear FLAG 0x400 on input reads before grouping? Picard does (it re-marks from scratch). We should too — adds one bitwise AND per record in Pass 2.
