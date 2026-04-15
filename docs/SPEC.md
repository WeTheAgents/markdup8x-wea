# markdup-wea — MVP Specification

**Status: APPROVED — Spec is law.**

## Context

Picard MarkDuplicates is the 3rd most expensive step in nf-core/rnaseq (~9% runtime). Java, single-threaded, 3.8-25GB JVM heap, OOM'd on 5/8 human samples at 32GB. samtools markdup requires 4 passes and is slower than Picard at t=1 (5,070s vs 1,791s). sambamba has ulimit/memory tuning issues.

**Goal:** Single static Rust binary. Two passes over coordinate-sorted BAM. Drop-in Picard replacement for nf-core. No JVM, no config, no surprises.

**MVP scope:** Paired-end + single-end duplicate detection, Picard-compatible metrics, multi-threaded I/O. **Not in MVP:** optical duplicates (output 0, document), UMI support, queryname-sorted input.

## Targets

| Metric | Picard | samtools t=1 | samtools t=4 | **markdup-wea** |
|--------|--------|-------------|-------------|-----------------|
| Wall time (human) | 1,791s | 5,070s | 1,220s | **<1,500s (t=1), <900s (t=4)** |
| Peak RSS | 21.5GB | 3.1GB | 3.6GB | **<200MB** |
| Passes | 1 (but buffers internally) | 4 | 4 | **2** |

Correctness: identical dup counts to Picard on 5 yeast + 8 human samples.

---

## Architecture: Two-Pass with Dupset

The dupset (set of record IDs to mark as duplicates) uses one of two strategies, selected at build time or via feature flag:
- **BitVec**: 0.125 bytes/read -> 25MB for 200M reads. Requires knowing max record count upfront (grow on demand). Best for memory.
- **FxHashSet<u64>**: 8+ bytes/entry but only stores actual dups -> 160-230MB for 20M dups. Simpler code, no index tracking. Best for simplicity.

Both implementations behind a `DupSet` trait. Benchmark both on yeast, pick winner for human genome runs.

### Why two passes?

In coordinate-sorted BAM, a read's mate appears at a distant position. We need the mate's CIGAR to compute its unclipped 5' position (BAM's MPOS field only stores the aligned position, not unclipped). Without seeing the mate, we can't form the complete grouping key. Therefore we must scan the entire BAM before we can determine which reads are duplicates.

Single-pass alternatives (buffering all records per chromosome, or writing then seeking back to patch flags) either blow the memory budget or require seekable output (breaks stdout piping).

### Pass 1 -- Scan

Read every record sequentially. Assign each record a sequential `record_id: u64` (counter from 0). Classify each record:

| Record type | How to detect | Action |
|-------------|--------------|--------|
| Unmapped | FLAG & 0x4 | Count for metrics, skip grouping |
| Secondary | FLAG & 0x100 | Count for metrics, skip grouping |
| Supplementary | FLAG & 0x800 | Count for metrics, skip grouping |
| Paired, mate unmapped | FLAG & 0x8 | Treat as **single-end** for grouping |
| Paired, both mapped | !(FLAG & 0x4) && !(FLAG & 0x8) | Paired-end grouping |
| Unpaired (single-end) | !(FLAG & 0x1) | Single-end grouping |

**Single-end grouping:** Group consecutive reads by `SingleEndKey { library_idx, ref_id, unclipped_5prime, is_reverse }`. Since BAM is coordinate-sorted, reads with the same key are adjacent. Resolve immediately: highest `quality_sum` wins, losers' record_ids -> set bit in `dup_bits`.

**Paired-end grouping:**
1. Compute own `unclipped_5prime` from CIGAR and `quality_sum` (bases with Q >= 15)
2. Look up QNAME hash in `pending_mates: FxHashMap<u64, PendingMate>`
3. If mate found -> remove from pending, form `PairedEndKey`, add `ScoredPair` to `GroupTracker`
4. If mate not found -> insert into `pending_mates` with own data

**PairedEndKey construction:** When both mates are known, order them by (ref_id, unclipped_5prime) ascending. The "lower" mate becomes pos_1, the "higher" becomes pos_2. This ensures the same key regardless of which mate we encounter first.

**Incremental group resolution:** Track each active group's `max_pos` (the higher of the two positions). Maintain a `BTreeMap<i64, Vec<GroupId>>` of groups keyed by their max_pos. When the stream advances past a position, resolve all groups with max_pos <= that position. Resolution: sort pairs by combined_score descending, first pair is primary, rest -> set both record_ids in `dup_bits`.

**At chromosome boundary:** Resolve all groups for the previous chromosome. At EOF: resolve all remaining groups (catches cross-chromosome pairs and any stragglers). Pending mates with no matching mate at EOF -> warn and skip (orphan reads, pass through unflagged).

### Pass 2 -- Write

Re-read the BAM from the beginning, maintaining the same record_id counter. For each record: if `dup_bits[record_id]` is set -> `FLAG |= 0x400`. If `--remove-duplicates` -> skip the record entirely. Write to output.

### Memory Budget

| Component | Bytes/entry | Count (200M reads) | BitVec mode | FxHashSet mode |
|-----------|------------|---------------------|-------------|----------------|
| dupset | 0.125 / 16 | 200M / 20M dups | **25MB** | **230MB** |
| `pending_mates` | ~52 + overhead | Peak ~2M (RNA-seq) | **~130MB** | **~130MB** |
| `GroupTracker` active | ~80 | ~100K concurrent | **~8MB** | **~8MB** |
| **Total** | | | **~163MB** | **~368MB** |

Both variants fit under 500MB. BitVec is clearly better for memory; FxHashSet is simpler (no record_id counter, no index tracking). We implement both behind a trait and benchmark on yeast to decide.

Worst case (WGS, wide inserts): pending_mates peak ~20M -> ~1.3GB. Out of MVP scope; add `--max-pending` with disk spill later.

### Stdin Handling

Two-pass requires seekable input. When input is `-` or `/dev/stdin`:
1. Copy stdin to a temp file (auto-deleted on drop via `tempfile` crate)
2. Two-pass over the temp file
3. Log warning: "Streaming from stdin; buffering to temp file at {path}"

---

## Phase 1: Skeleton + Core Primitives (~3 days)

**Goal:** Binary that reads BAM, validates sort order, passes through unchanged. All building blocks tested.

**Create:**

`Cargo.toml`:
```toml
[dependencies]
rust-htslib = { version = "0.47", features = ["static"] }
clap = { version = "4", features = ["derive"] }
anyhow = "1"
log = "0.4"
env_logger = "0.11"
rustc-hash = "2"       # FxHashMap
bitvec = "1"            # BitVec for dupset
tempfile = "3"          # stdin buffering

[profile.release]
lto = true
codegen-units = 1
strip = true
opt-level = 3
```

`src/main.rs` -- CLI:
```
markdup-wea [OPTIONS] <INPUT>

Arguments:
  <INPUT>   Input BAM file (coordinate-sorted). Use - for stdin.

Options:
  -o, --output <FILE>          Output BAM [default: stdout]
  -M, --metrics <FILE>         Picard-format metrics file
  -@, --threads <N>            I/O threads for BAM reading/writing [default: 1]
      --remove-duplicates      Exclude duplicates from output
      --assume-sort-order <SO> Override header sort order check
```

`src/position.rs` -- `unclipped_5prime(record: &Record) -> i64`:
- Forward read: `pos - left_clips` (soft + hard at 5' end)
- Reverse read: `cigar.end_pos() + right_clips` (soft + hard at 3' end, which is 5' for reverse)
- Edge case: empty CIGAR (unmapped reads) -> return pos as-is

`src/scoring.rs` -- `quality_sum(record: &Record) -> u32`:
- Sum of all base qualities >= 15
- For paired scoring: `pair_score = quality_sum(read1) + quality_sum(read2)`

`src/io.rs`:
- `open_reader(path, threads)` -> `bam::Reader` (detect stdin -> temp file)
- `open_writer(path, threads, header)` -> `bam::Writer`
- `validate_sort_order(header)` -> Result (check @HD SO:coordinate)
- Runtime sort-order enforcement: track `(prev_ref_id, prev_pos)`, error if current < prev

**Tests:**
- `position.rs`: forward no-clip, forward soft-clip, reverse soft-clip, hard-clip, mixed clip, empty CIGAR
- `scoring.rs`: all Q>=15, all Q<15, mixed, empty quals
- `io.rs`: valid SO:coordinate, missing SO, wrong SO, runtime out-of-order detection

**Done when:** `markdup-wea input.bam -o output.bam` produces byte-identical copy. CI green (build + test + clippy + fmt).

---

## Phase 2: Paired-End Duplicate Detection (~1 week)

**Goal:** Correct duplicate marking. This is the core of the tool.

**Create:**

`src/pending_mates.rs`:
```rust
#[repr(C, packed)]  // force no padding
struct PendingMate {
    name_hash: u64,        // 8  -- FxHash of QNAME
    check_hash: u32,       // 4  -- secondary hash (first 4 bytes of QNAME)
    ref_id: i32,           // 4
    unclipped_5prime: i64,  // 8
    mate_ref_id: i32,      // 4
    mate_pos: i64,         // 8  -- MPOS from BAM (aligned, not unclipped)
    quality_sum: u32,      // 4
    record_id: u64,        // 8
    library_idx: u8,       // 1  -- up to 255 libraries
    is_reverse: bool,      // 1
    _pad: [u8; 2],         // 2  -- explicit padding to align to 8
}                          // = 52 bytes total
// Stored in FxHashMap<u64, PendingMate> keyed by name_hash
```

When mate arrives, collision check: compare `check_hash`. If mismatch (hash collision), fall back to full QNAME comparison by re-reading the record (extremely rare with 64-bit primary hash).

`src/groups.rs`:
```rust
struct PairedEndKey {
    library_idx: u8,
    ref_id_lo: i32,        // lower-coordinate mate's ref
    pos_lo: i64,           // lower-coordinate mate's unclipped 5'
    is_reverse_lo: bool,
    ref_id_hi: i32,        // higher-coordinate mate's ref
    pos_hi: i64,           // higher-coordinate mate's unclipped 5'
    is_reverse_hi: bool,
}

struct ScoredPair {
    combined_score: u32,
    record_id_1: u64,      // read1's record_id
    record_id_2: u64,      // read2's record_id
}

struct GroupTracker {
    // Active groups: key -> vec of scored pairs
    groups: FxHashMap<PairedEndKey, Vec<ScoredPair>>,
    // Index: max_pos -> list of keys to resolve at that position
    resolve_at: BTreeMap<i64, Vec<PairedEndKey>>,
}
```

Group resolution: for each group, sort by combined_score desc. Index 0 = primary (not flagged). All others -> set `dup_bits[record_id_1]` and `dup_bits[record_id_2]`.

`src/scan.rs` -- Pass 1 orchestrator:
- Read all records, classify, route to single-end or paired-end grouping
- Trigger incremental resolution as stream advances
- Return `dup_bits: BitVec`

`src/markdup.rs` -- Two-pass orchestrator:
- Pass 1: `scan(reader) -> dup_bits`
- Reset reader to beginning
- Pass 2: iterate records, apply flags from dup_bits, write

**Library handling:** Parse `@RG` headers -> map `RG:ID -> LB`. If no `@RG`, single library. Each record's `RG` aux tag -> library_idx.

**Edge cases explicitly handled:**
- Mate unmapped (FLAG 0x8): treat as single-end
- Both unmapped: skip grouping entirely
- MAPQ=0: handled normally (quality_sum may be low -> likely flagged as dup)
- Orphan reads (mate never found): warn to stderr, pass through unflagged
- Cross-chromosome pairs: stay pending until mate's chromosome; resolved at EOF if not before

**Supplementary/secondary reads: NOT flagged as duplicates in MVP.** Only primary alignments get FLAG 0x400. This matches samtools markdup behavior. Document in README: "Unlike Picard, supplementary alignments of duplicate reads are not flagged. This does not affect downstream analysis tools (featureCounts, HTSeq, RSEM) which only use primary alignments."

**Tests:**
1. Synthetic BAM: two identical pairs at same position -> one flagged
2. Synthetic BAM: pair with higher quality wins over lower
3. Synthetic BAM: single-end reads mixed with paired
4. Synthetic BAM: mate on different chromosome -> both flagged
5. Synthetic BAM: mate unmapped -> single-end treatment
6. **Real data: match Picard dup counts on 5 yeast samples exactly**
   - WT_REP1: 36,810 dups / 180,342 total
   - WT_REP2: 11,688 / 90,962
   - RAP1_UNINDUCED_REP2: 78,929 / 98,201
   - RAP1_UNINDUCED_REP1: 36,294 / 48,977
   - RAP1_IAA_30M_REP1: 11,094 / 48,347

**Test BAM acquisition:** Run `nf-core/rnaseq -profile test,docker --outdir test_out --skip_markduplicates` to get coordinate-sorted BAMs without Picard flags. Document in `tests/README.md`.

**Done when:** All synthetic tests pass. Dup counts match Picard on all 5 yeast samples.

---

## Phase 3: Picard-Compatible Metrics (~3 days)

**Goal:** Metrics file that MultiQC parses as Picard output.

**Create:** `src/metrics.rs`

Counters accumulated during Pass 1:
```rust
struct MetricsCounters {
    unpaired_reads_examined: u64,
    read_pairs_examined: u64,
    secondary_or_supplementary: u64,
    unmapped_reads: u64,
    unpaired_read_duplicates: u64,
    read_pair_duplicates: u64,
    // optical = 0 in MVP
    dup_group_sizes: Vec<u64>,  // histogram: index = group_size, value = count
}
```

**Output format (exact -- MultiQC parses by header line):**
```
## htsjdk.samtools.metrics.StringHeader
# markdup-wea {version} INPUT={input} OUTPUT={output}
## METRICS CLASS	picard.sam.markduplicates.DuplicationMetrics
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
{lib}	{u}	{p}	{ss}	{um}	{ud}	{pd}	0	{pct}	{est}

## HISTOGRAM	java.lang.Double
BIN	VALUE
1.0	{n}
2.0	{n}
...
```

**PERCENT_DUPLICATION** = `(unpaired_dups + read_pair_dups * 2) / (unpaired_examined + read_pairs_examined * 2)`

**ESTIMATED_LIBRARY_SIZE** -- Lander-Waterman equation:
- n = read_pairs_examined, c = read_pairs_examined - read_pair_duplicates
- Solve `c/x = 1 - exp(-n/x)` for x by Newton's method (~10 iterations)
- Special cases: c == 0 -> return 0, c == n -> return n (no dups)

**READ_PAIR_OPTICAL_DUPLICATES** = 0 in MVP. ESTIMATED_LIBRARY_SIZE uses c = pairs - dups (no optical adjustment).

**Test:** Run `multiqc .` on output metrics -> verify duplication rate plot appears with correct values. Compare our metrics values to Picard's metrics on same yeast samples.

**Done when:** MultiQC generates identical duplication plot from our metrics.

---

## Phase 4: Threading + Optimization (~3 days)

**I/O threads:** Pass `-@ N` to `hts_set_threads()` for BGZF decompression/compression. This is ~40% of wall time. Core algorithm stays single-threaded.

**Optimizations:**
- Record reuse: `reader.read(&mut record)` (amortize allocation)
- Pre-size FxHashMap: estimate read count from BAM file size / avg record size
- Use `BufWriter` for metrics output

**NOT in MVP:** Parallel group resolution, SIMD scoring, memory-mapped I/O, disk spill for pending_mates.

**Done when:** `markdup-wea -@ 4` on yeast data is faster than Picard at t=1.

---

## Phase 5a: Yeast Validation + Dupset Benchmark (~2 days)

**Gate: everything works on yeast before touching human data.**

**Test data:** 5 yeast samples from `nf-core/rnaseq -profile test,docker --skip_markduplicates`

**Correctness:** Dup counts must match Picard exactly on all 5 samples.

**Dupset comparison:** Run both BitVec and FxHashSet modes on all 5 samples. Compare:
- Peak RSS (`/usr/bin/time -v`)
- Wall time (`hyperfine --runs 5`)
- Code complexity (subjective)

Pick winner for human genome runs. If negligible difference -> keep BitVec (lower memory ceiling).

**Done when:** All 5 yeast samples pass. Dupset strategy decided.

## Phase 5b: Human Genome Validation + Full Benchmark (~1 week)

**Only after Phase 5a is green.**

**Test data:** 8 ENCODE human RNA-seq samples (existing cloud benchmark infra)

**Correctness gate:** All 8 samples produce identical dup counts to Picard. If any differ -> investigate and document.

**Benchmark protocol:**
```bash
hyperfine --warmup 1 --runs 3 \
  'markdup-wea input.bam -o /dev/null -M metrics.txt -@ {threads}' \
  --export-csv results.csv

/usr/bin/time -v markdup-wea input.bam -o out.bam -M metrics.txt -@ 4
```

**Deliverables:**
- Performance comparison table (Picard vs samtools vs markdup-wea) in README
- Scripts in `benchmarks/`
- Memory profile confirming <200MB RSS (BitVec) or <400MB (FxHashSet)

---

## Repo Structure

```
markdup-wea/
├── Cargo.toml
├── LICENSE (MIT)
├── README.md
├── docs/
│   └── SPEC.md              # this file
├── src/
│   ├── main.rs              # CLI
│   ├── lib.rs               # public API
│   ├── markdup.rs           # two-pass orchestrator
│   ├── scan.rs              # Pass 1: scan + classify + group
│   ├── pending_mates.rs     # PendingMate + FxHashMap buffer
│   ├── groups.rs            # PairedEndKey, SingleEndKey, GroupTracker
│   ├── position.rs          # unclipped_5prime()
│   ├── scoring.rs           # quality_sum()
│   ├── metrics.rs           # Picard-format metrics + Lander-Waterman
│   └── io.rs                # BAM reader/writer, sort validation, stdin->temp
├── tests/
│   ├── integration.rs       # end-to-end on synthetic BAMs
│   ├── fixtures/            # committed small BAMs
│   └── README.md            # how to get real test data
├── benchmarks/              # scripts + results
└── .github/workflows/ci.yml
```

## Documented Deviations from Picard

| Behavior | Picard | markdup-wea (MVP) | Impact |
|----------|--------|-------------------|--------|
| Supplementary flagging | Flags supplementaries of dup reads | Does NOT flag supplementaries | None for RNA-seq (downstream tools ignore supplementaries) |
| Optical duplicates | Detects and counts | Reports 0 | ESTIMATED_LIBRARY_SIZE slightly different on NovaSeq data |
| Queryname-sorted input | Supported | Not supported (error) | None for nf-core (always coordinate-sorted) |
| Cross-chromosome dups | Detected | Detected (resolved at EOF) | Identical behavior |

## Verification Checklist

### Stage 1 -- Yeast (gate for human)
- [ ] `cargo test` -- unit + integration tests pass
- [ ] `cargo clippy` -- zero warnings
- [ ] Dup counts match Picard on all 5 yeast samples
- [ ] `multiqc .` -> correct duplication plot from our metrics
- [ ] Dupset strategy benchmarked (BitVec vs FxHashSet), winner picked
- [ ] `markdup-wea -@ 4` faster than Picard on yeast
- [ ] Runtime sort-order check catches unsorted input
- [ ] Stdin piping works (`samtools sort in.bam | markdup-wea - -o out.bam`)
- [ ] `--remove-duplicates` produces clean BAM
- [ ] Orphan reads (missing mate) produce warning, not crash
- [ ] Static binary < 10MB

### Stage 2 -- Human genome
- [ ] Dup counts match Picard on all 8 human samples
- [ ] Peak RSS < target on human data
- [ ] Wall time competitive with Picard at t=1, faster at t=4
- [ ] Performance table published in README
