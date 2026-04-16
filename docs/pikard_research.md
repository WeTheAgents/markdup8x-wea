# Picard MarkDuplicates: complete internals reference for byte-compatible reimplementation

**Picard's MarkDuplicates groups reads by unclipped 5ʹ positions and orientation, scores them by sum of base qualities ≥ Q15, and breaks ties by file order — but its behavior diverges sharply between coordinate-sorted and query-sorted input modes.** This report documents every known behavioral nuance, from source-code-level algorithm details to community-discovered quirks, that could cause divergence in a Rust replacement. The findings draw from direct analysis of the Picard and htsjdk source, 30+ GitHub issues, GATK forum threads, Biostars discussions, and comparisons with samtools/sambamba/biobambam implementations.

---

## 1. How unclipped 5ʹ positions are computed

The 5ʹ coordinate calculation is the single most critical piece of the algorithm. Picard delegates to htsjdk's `SAMRecord.getUnclippedStart()` and `SAMRecord.getUnclippedEnd()`, which in turn call `SAMUtils` methods.

**Forward-strand reads** use `getUnclippedStart()`:
```
read1Coordinate = rec.getUnclippedStart()  // when !getReadNegativeStrandFlag()
```
The algorithm walks the CIGAR from the **left**, subtracting the length of every leading `S` (soft-clip) or `H` (hard-clip) element from `alignmentStart`. It stops at the first non-clip operator. Example: `alignmentStart=100`, CIGAR `3H5S50M` → unclippedStart = 100 − 5 − 3 = **92**.

**Reverse-strand reads** use `getUnclippedEnd()`:
```
read1Coordinate = rec.getUnclippedEnd()  // when getReadNegativeStrandFlag()
```
This walks the CIGAR from the **right**, adding the length of every trailing `S`/`H` element to `alignmentEnd`. The key subtlety is that **`alignmentEnd` itself is computed as `alignmentStart + referenceLength − 1`**, where `referenceLength` sums all reference-consuming operators: **M, D, N, =, X**. This means the **N (intron-skip) operator inflates the unclipped end** for reverse-strand reads.

For a reverse-strand read with CIGAR `50M5000N50M3S`, the calculation is: `alignmentEnd = start + 50 + 5000 + 50 − 1`, then `unclippedEnd = alignmentEnd + 3`. Two reverse-strand reads at the same start position but with different intron sizes will have **different 5ʹ coordinates** — a correct behavior for RNA-seq (different splice isoforms should not be duplicates) but a common source of confusion.

Both soft and hard clips are treated **identically** in the unclipped position calculation. All coordinates are **1-based** (SAM convention). `getUnclippedStart()` can theoretically return values < 1 if clips extend before position 1, though this is extremely rare in practice. Unmapped reads (`getReadUnmappedFlag()`) are skipped entirely and never enter the duplicate detection pipeline.

---

## 2. Scoring, tie-breaking, and the short overflow trap

**Default scoring strategy is `SUM_OF_BASE_QUALITIES`**, set in the `MarkDuplicates` constructor. The score is computed by `DuplicateScoringStrategy.getSumOfBaseQualities()`:

```java
short score = 0;
for (final byte b : rec.getBaseQualities()) {
    if (b >= 15) score += b;  // threshold is >= 15, inclusive
}
```

Critical details for byte-compatibility:

- The threshold is **Q ≥ 15** (inclusive), not Q > 15. All bases in the full quality array are evaluated.
- The return type is **`short` (Java `i16`)**, which overflows silently at 32,767. For reads longer than ~400 bp with uniformly high qualities, this can wrap negative — a latent bug that matters for long-read data.
- When QUAL is `"*"` (missing), `getBaseQualities()` returns an **empty array**, yielding score **0**. These reads always lose to reads with quality scores.
- **Pair score = read1_score + read2_score**, computed by adding the second mate's score when it is encountered during the coordinate-sorted traversal.

**Tie-breaking in MarkDuplicates** (not MarkDuplicatesWithMateCigar) uses a simple two-step chain:

1. **Highest `score` wins** (the read/pair kept as non-duplicate).
2. **Lowest `read1IndexInFile` wins** — meaning the first record encountered in file order survives.

There is **no QNAME lexicographic comparison**, **no read-name hash**, and **no mapping quality** in MarkDuplicates tie-breaking. This is a crucial distinction from MarkDuplicatesWithMateCigar, which uses a richer comparison chain: score → MAPQ (sum of read + mate MAPQ) → read name (lexicographic) → first-of-pair flag. The elPrep paper confirms that this file-order dependency makes MarkDuplicates inherently **non-deterministic** when reads with identical scores arrive in different orders.

**MAPQ=0 reads are fully included** in duplicate detection with no special handling. Picard issue #1285 (filed by developer @yfarjoun) documents a known consequence: chimeric pairs where one mate has MQ=0 get a random placement by the aligner, causing duplicates to be missed. This bug is intentionally unremediated in current Picard.

---

## 3. The grouping key and coordinate-ordering rules

### Paired-end grouping key
For pairs, the duplicate group key is the tuple: **(libraryId, read1ReferenceIndex, read1Coordinate, orientation, read2ReferenceIndex, read2Coordinate)**, plus optional barcodes. Two `ReadEndsForMarkDuplicates` entries must match on all these fields to be considered duplicates (via `areComparableForDuplicates`).

### Fragment grouping key
For fragments (single-end or mate-unmapped): **(libraryId, read1ReferenceIndex, read1Coordinate, orientation)**. No read2 fields are compared.

### Coordinate ordering of read1 / read2
The read with the **lower reference index** becomes read1. If reference indices are equal, the read with the **lower unclipped coordinate** becomes read1. When the mate arrives during traversal:

- If `matesRefIndex > read1RefIndex` OR (`==` AND `matesCoord >= read1Coord`): mate becomes read2. Note the **`>=`** — when coordinates are equal, the first-encountered read stays as read1.
- Otherwise: the reads are **flipped** — the mate becomes read1, the original read1 becomes read2.

### The same-position RF → FR normalization
After coordinate ordering, if both reads map to the **same reference and same coordinate** and the orientation is `RF`, it is **forced to `FR`**. This ensures that RF and FR pairs at identical positions are grouped together rather than forming separate duplicate sets.

### Orientation encoding
Constants in `ReadEnds`: **F=0, R=1, FF=2, FR=3, RR=4, RF=5**. The `getOrientationByte(read1NegStrand, read2NegStrand)` method uses the coordinate-ordered read1/read2 strands. A separate field, `orientationForOpticalDuplicates`, preserves the **sequencing-order** orientation (first-of-pair always comes first), computed before coordinate reordering. This distinction matters: PCR duplicate detection groups FR and RF together at same positions, but optical duplicate detection keeps them separate because true optical duplicates preserve read-order orientation.

### Sort comparator (ReadEndsMDComparator)
The full sort order is: `libraryId → [barcodes] → read1ReferenceIndex → read1Coordinate → orientation → read2ReferenceIndex → read2Coordinate → read1IndexInFile → read2IndexInFile`.

---

## 4. Fragment-vs-pair priority and the "fragments always lose" rule

When a genomic position has **both paired and unpaired** reads, all unpaired fragments at that position are **unconditionally marked as duplicates**, regardless of their scores. Pairs always dominate fragments. This is implemented in `markDuplicateFragments()`:

```java
if (containsPairs) {
    for (ReadEndsForMarkDuplicates end : list) {
        if (!end.isPaired()) addIndexAsDuplicate(end.read1IndexInFile);
    }
} else {
    // only fragments: highest score wins
}
```

This means a single-end read with a perfect quality score will still be marked as a duplicate if any pair maps to the same position. This is an intentional design decision — pairs provide stronger evidence of unique molecular origin.

---

## 5. The coordinate-vs-queryname behavioral split

This is the **most consequential architectural difference** for reimplementation. The same BAM produces different outputs depending on sort order:

| Behavior | Coordinate-sorted | Query-sorted/grouped |
|---|---|---|
| Supplementary alignments marked as dup | **No** — skipped entirely | **Yes** — inherit primary's dup flag |
| Secondary alignments marked as dup | **No** | **Yes** |
| Unmapped mates of dup reads | **Not marked** | **Marked as dup** |
| Mechanism | Only `read1IndexInFile`/`read2IndexInFile` flagged | `duplicateQueryName` tracking — all records with matching QNAME get flagged |

For query-sorted input, the duplicate name propagation code checks:
```java
boolean isDuplicate = recordIndex == nextDuplicateIndex ||
    (sortOrder == queryname && recordIndex > nextDuplicateIndex && 
     rec.getReadName().equals(duplicateQueryName));
```

**Important quirk**: Picard expects **lexicographic** queryname sorting, but `samtools sort -n` uses **numeric** sorting. This mismatch can cause `IllegalArgumentException: Alignments added out of order` (GitHub issue #1194). The user must use `ASSUME_SORT_ORDER=queryname` carefully.

---

## 6. Library detection and the @RG fallback chain

**Correction (Phase D, 2026-04-16):** the original write-up below incorrectly
described step 2 as falling back to the @RG ID. Verified against Picard 3.4.0
source (`picard.sam.markduplicates.util.LibraryIdGenerator.getLibraryName`)
and against real-data byte-diff on 8 ENCODE samples (Phase C):

1. Read's `RG` tag → look up `@RG` header record → use `LB` field
2. If `LB` is missing: fall back to the literal string **"Unknown Library"**
   (NOT the @RG ID). All reads with absent LB therefore collapse into a
   single library bucket — which means duplicates ARE detected across @RG
   entries that all lack LB.
3. If the read has **no RG tag at all**: in theory returns a constant
   "unknown library" string, but in practice **throws NullPointerException**
   at `buildSortedReadEndLists` (confirmed by GATK community reports).

Multiple `@RG` lines sharing the same `LB` value receive the **same `libraryId` (short)**, meaning reads from different read groups with the same library are potential duplicates of each other. Library ID is the **first field checked** in `areComparableForDuplicates` — reads from different libraries are never compared.

The `libraryId` is a `short` (max 32,767 unique libraries), assigned incrementally by `LibraryIdGenerator.getLibraryId()`.

---

## 7. Supplementary, secondary, unmapped, and orphan read handling

**Supplementary (0x800) and secondary (0x100)** alignments are excluded from `buildReadEnds()` with an explicit `!rec.isSecondaryOrSupplementary()` check. They never participate in duplicate scoring or grouping. Their flagging depends on sort order (see section 5).

**Mate-unmapped reads (FLAG 0x8 set)**: When a read is paired (`0x1`) but its mate is unmapped, the `read2ReferenceIndex` stays at −1, causing `isPaired()` to return false. The read is added to `fragSort` only, participating in **fragment-level** duplicate detection. It does not enter `pairSort`.

**Chimeric/inter-chromosomal pairs**: Treated as normal pairs. The read on the lower reference index becomes read1. Both mates are flagged if the pair is a duplicate. Cross-chromosome pairs in **split BAMs** (separate files per chromosome) will not be identified because the mate is never encountered.

**Orphan reads** (paired flag set but mate absent from file): Added to the temporary pair-matching map. When the mate never appears, an INFO-level log message reports the count: `"Read N records. M pairs never matched."` No error is thrown. The orphan read remains in `fragSort` as a fragment.

**Pre-existing FLAG 0x400**: Picard **always clears and recomputes** duplicate flags. Non-duplicates explicitly get `setDuplicateReadFlag(false)`. The algorithm ignores any pre-existing duplicate marking.

---

## 8. Optical duplicate detection internals

### READ_NAME_REGEX optimization
The "default regex" is **not a regex at all**. Picard's `ReadNameParser` uses an optimized colon-split algorithm: it splits the read name on `':'` and extracts the last three numeric fields as tile, x, y. For 5-element names (e.g., `HWUSI:6:73:941:1973`), fields 3/4/5 are used. For 7-element names (standard Illumina, e.g., `HISEQ:1:1101:1234:5678:ACGT:0`), fields 5/6/7 are used. The header prints this as `--READ_NAME_REGEX <optimized capture of last three ':' separated fields as numeric values>`.

When read names don't match (e.g., SRA-style names like `SRR5183158.19928893`), a warning is logged and optical detection is silently disabled for those reads.

### PhysicalLocation short overflow bug
In older Picard versions, tile x/y coordinates were stored as **`short` (max 32,767)**. HiSeq X and NovaSeq flowcells produce coordinates exceeding this range, causing silent overflow. GitHub issues #1216, #921, and #1441 document this: coordinates ~65K apart could erroneously appear within pixel distance, creating **false-positive optical duplicates**. Developer @yfarjoun acknowledged this as a known issue. A fix (PR #2034) switched to `PhysicalLocationInt` in newer versions. **For byte-compatibility, you must match the target Picard version's coordinate storage type.**

### Default OPTICAL_DUPLICATE_PIXEL_DISTANCE
The default is **100 pixels**, appropriate for unpatterned flowcells. For patterned flowcells (NovaSeq 6000/X, HiSeq X, HiSeq 4000), **2500** is recommended (issue #1252). The O(n²) comparison algorithm can hang on extremely large duplicate sets (issue #472 reports problems at n=500K).

### Optical duplicates are pair-only
The `READ_PAIR_OPTICAL_DUPLICATES` metric tracks only **pairs**, never fragments. Single-end reads cannot be classified as optical duplicates.

---

## 9. Metrics file format and library size estimation

### File structure
```
## htsjdk.samtools.metrics.StringHeader
# <full command line>
## htsjdk.samtools.metrics.StringHeader
# Started on: <timestamp>

## METRICS CLASS	picard.sam.DuplicationMetrics
<TAB-separated column headers>
<TAB-separated data rows, one per library>

## HISTOGRAM	java.lang.Double
<histogram data>
```

### Column order (10 columns in current master)
`LIBRARY`, `UNPAIRED_READS_EXAMINED`, `READ_PAIRS_EXAMINED`, `SECONDARY_OR_SUPPLEMENTARY_RDS`, `UNMAPPED_READS`, `UNPAIRED_READ_DUPLICATES`, `READ_PAIR_DUPLICATES`, `READ_PAIR_OPTICAL_DUPLICATES`, `PERCENT_DUPLICATION`, `ESTIMATED_LIBRARY_SIZE`

The `SECONDARY_OR_SUPPLEMENTARY_RDS` column was **added in a later version** — older Picard had only 9 columns. `PERCENT_DUPLICATION` is a **fraction** (0.05, not 5%), despite its name.

### Counting quirk: per-read then divide by 2
`READ_PAIRS_EXAMINED` and `READ_PAIR_DUPLICATES` are incremented **once per read** during traversal (the code comments say "will need to be divided by 2 at the end"). The finalization step divides these counts by 2 to convert from per-read to per-pair counts.

### ESTIMATED_LIBRARY_SIZE (Lander-Waterman bisection)
```
readPairs = READ_PAIRS_EXAMINED − READ_PAIR_OPTICAL_DUPLICATES
uniqueReadPairs = READ_PAIRS_EXAMINED − READ_PAIR_DUPLICATES
f(x, c, n) = c/x − 1 + exp(−n/x)
```
The algorithm uses **bisection** (not Newton's method): initial bracket `m=1.0, M=100.0` (as multipliers of `uniqueReadPairs`), `M` grows by 10× until `f(M × uniqueReadPairs) > 0`, then 40 bisection iterations. Returns `(long)(uniqueReadPairs × (m + M) / 2.0)`.

Edge cases:
- **Zero duplicates** (`readPairDuplicates == 0`): returns `null` (empty field in output)
- **All duplicates** (`uniqueReadPairs >= readPairs`): throws `IllegalStateException`
- **Zero examined** (denominator = 0): `PERCENT_DUPLICATION` is set to `0.0` (explicit check, not NaN)
- **Single-end only data**: `ESTIMATED_LIBRARY_SIZE` is `null`; histogram is omitted

### MultiQC compatibility
MultiQC identifies Picard files by searching for strings like `picard.sam.DuplicationMetrics` or `picard.sam.markduplicates.MarkDuplicates`. It dynamically reads column headers from the `UNPAIRED_READ_DUPLICATES` line. When multiple libraries exist, it **sums numeric values** across rows and recomputes derived fields. The library size re-implementation in MultiQC mirrors Picard's bisection method exactly (40 iterations, same bracket logic). MultiQC skips samples where `READ_PAIRS_EXAMINED == 0 AND UNPAIRED_READS_EXAMINED == 0`.

---

## 10. RNA-seq: how intron-spanning reads interact with duplicate detection

Picard has **no special RNA-seq logic**. Spliced alignments are handled by the generic CIGAR processing rules. The key consequence:

- **Forward-strand spliced reads**: The N operator does not affect `getUnclippedStart()` because N is an internal operator, not a leading clip. Two forward-strand reads at the same start position with different intron structures will have the **same** 5ʹ coordinate and may be grouped as duplicates.
- **Reverse-strand spliced reads**: The N operator **does** affect `getUnclippedEnd()` because `alignmentEnd = alignmentStart + Σ(M,D,N,=,X) − 1`. Different intron sizes produce different 5ʹ coordinates, so different splice isoforms are correctly separated.

**Practical RNA-seq recommendations** from Picard documentation and community:
- Use MarkDuplicates, **not** MarkDuplicatesWithMateCigar (which exhausts memory on large intron skips)
- Set `READ_NAME_REGEX=null` to disable optical duplicate detection (extremely large duplicate sets at highly-expressed loci cause O(n²) slowdowns)
- Note that **64% of single-end** and **19% of paired-end** "duplicates" in RNA-seq may be valid biological transcripts (Salzberg et al.), not PCR artifacts
- nf-core/rnaseq uses `ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT` but does **not** disable optical detection by default

---

## 11. Key GitHub issues documenting intentional quirks

| Issue | Status | Behavioral quirk | Reimplementation impact |
|---|---|---|---|
| **#166** | Closed | Developer Tim Fennell explains the complete algorithm: unclipped 5ʹ, Q≥15, per-library grouping | Canonical reference |
| **#1119** | Closed | Coordinate vs query-sorted mode differences for supplementary/unmapped flagging | Must implement both paths |
| **#451** | Closed | Unmapped mates NOT marked in coordinate mode (intentional) | Must reproduce |
| **#1138** | Open | Developer @yfarjoun argues unmapped mates should be marked; current fix only for queryname mode | Known inconsistency |
| **#1285** | Open | MQ=0 reads cause missed duplicates (intentional — not treated as unmapped) | Must reproduce |
| **#1216** | Closed | Optical dup coordinate short overflow — acknowledged, won't fix (until PR #2034) | Version-dependent |
| **#1194** | Open | Queryname sort: Picard expects lexicographic, samtools uses numeric ordering | Sort validation |
| **#1604** | Info | Read names not matching regex → optical detection silently disabled | Must reproduce |
| **#902** | Closed | `ADD_PG_TAG_TO_READS=true` (default) adds PG tag to every read, increasing file size | Output format |
| **#1944** | Open | Single-end metrics omit `ESTIMATED_LIBRARY_SIZE` and histogram | Metrics edge case |
| **#472** | Open | O(n²) optical dup detection hangs on large sets (>500K) | Performance |

---

## 12. How samtools, sambamba, and GATK Spark diverge from Picard

### samtools markdup
Requires pre-processing with `samtools fixmate -m` to add `MC` and `ms` tags. Uses the `ms` tag for scoring (sum of quality values, with a different threshold formula). Supports two modes: **template** (default, measures template start/end with R1/R2 distinction) and **sequence** (measures 5ʹ ends, closer to Picard). The `-S` flag optionally marks supplementary alignments of duplicates. samtools was "written to match Picard 2.10.3" according to community reports but is not byte-compatible.

### sambamba markdup
Explicitly implements "the same criteria as Picard" for grouping (unclipped 5ʹ + orientation). ~99% concordant with Picard in practice. Key difference: **no optical duplicate detection**. The Bismark issue #170 revealed a name-sort bug: sambamba sorts only by read name string (no SAM flag tiebreaker), so when Bismark assigns identical names to R1/R2, sambamba's name sort can swap R1/R2 order, corrupting methylation extraction. samtools avoids this by using SAM flags as a secondary sort key.

### GATK MarkDuplicatesSpark
Internally queryname-sorts all input, so it **always** marks supplementary/secondary/unmapped mates — even for coordinate-sorted input. This produces **different duplicate counts** than Picard on coordinate-sorted BAMs (one user reported 10,571 vs 10,622 duplicates). Results converge when input is queryname-sorted.

### biobambam2 bammarkduplicates
Similar algorithm to Picard with much better memory efficiency. Results are "slightly different but similar" — the biobambam paper focuses on runtime rather than exact concordance.

---

## 13. The CLEAR_DT tag, ASSUME_SORT_ORDER, and other parameter quirks

**CLEAR_DT** (default `true`): Removes existing `DT` (Duplicate Type) tags from all input records before processing. The `DT` tag values are `LB` (library/PCR duplicate) and `SQ` (sequencing/optical duplicate). Writing new DT tags is controlled by `TAGGING_POLICY`: `DontTag` (default), `OpticalOnly`, or `All`.

**ASSUME_SORT_ORDER** vs **ASSUME_SORTED**: These are **mutually exclusive** (`mutex` annotation). `ASSUME_SORTED` (deprecated) is a boolean that assumes coordinate sort. `ASSUME_SORT_ORDER` accepts `coordinate`, `queryname`, or `unsorted`. If neither is specified and the header's `SO` field doesn't indicate coordinate or queryname, Picard throws: `"This program requires input that are either coordinate or query sorted."` When `ASSUME_SORT_ORDER=queryname`, the output header is set to `SO:unknown GO:query`.

**ADD_PG_TAG_TO_READS** (default `true`): Adds `PG:Z:MarkDuplicates` to every output read (not just duplicates), increasing file size by ~16% according to benchmark testing.

**REMOVE_DUPLICATES quirk**: When `true`, `REMOVE_DUPLICATES` removes more reads than simply filtering by FLAG 0x400 — in queryname mode, it also removes associated supplementary, secondary, and unmapped mate records of duplicate primaries.

---

## 14. ReadEnds data structure: the complete field inventory

```
PhysicalLocationInt → ReadEnds → ReadEndsForMarkDuplicates → ReadEndsForMarkDuplicatesWithBarcodes

ReadEnds fields:
  libraryId         : short    — library index (grouping)
  orientation       : byte     — F/R/FF/FR/RR/RF (grouping)
  read1ReferenceIndex: int     — 0-based ref index of read1 (grouping)
  read1Coordinate   : int      — 1-based unclipped 5ʹ of read1 (grouping)
  read2ReferenceIndex: int     — -1 if unpaired (grouping for pairs)
  read2Coordinate   : int      — (grouping for pairs)
  readGroup         : short    — read group index (optical only)
  orientationForOpticalDuplicates: byte — sequencing-order (optical only)

ReadEndsForMarkDuplicates adds:
  score             : short    — quality score (scoring)
  read1IndexInFile  : long     — file position (tie-breaking + flagging)
  read2IndexInFile  : long     — file position (flagging)

PhysicalLocation fields:
  tile              : int      — flowcell tile (optical)
  x                 : short/int — x coordinate (optical, type is version-dependent!)
  y                 : short/int — y coordinate (optical)

ReadEndsForMarkDuplicatesWithBarcodes adds:
  barcode           : int      — Objects.hash(normalizedUMI)
  readOneBarcode    : int
  readTwoBarcode    : int
```

The `isPaired()` method simply checks `read2ReferenceIndex != -1`. The `readGroup` field is the 0-based index of the RG in the header's read group list, used exclusively for optical duplicate detection (reads on different tiles from different read groups are never optical duplicates of each other).

---

## 15. Specific technical edge cases for completeness

**Empty CIGAR strings**: Unmapped reads may have `*` as CIGAR. These reads are filtered out by the `getReadUnmappedFlag()` check before `buildReadEnds()` is called, so they never reach CIGAR parsing.

**Reads at position 0**: In SAM, position 0 means unmapped. `getAlignmentStart()` returns 0 for unmapped reads. These are caught by the unmapped-flag check.

**Invalid FLAG combinations**: Picard does not validate FLAG consistency exhaustively. With `VALIDATION_STRINGENCY=LENIENT` (commonly used), malformed flags may produce surprising results — e.g., a read with both 0x4 (unmapped) and a valid POS could bypass the unmapped check if the flag is not set but the position is 0.

**Reference index**: `read1ReferenceIndex` stores the **0-based** reference sequence index from the BAM header (not the reference name), while `read1Coordinate` stores the **1-based** alignment position.

**Mate lookup key**: The temporary map for pairing uses `readGroupId + readName` as the lookup key, keyed by the **mate's reference index**. This means reads with identical names in different read groups will not be paired with each other.

---

## Conclusion: the eight hardest things to get right

Building a byte-compatible Picard MarkDuplicates replacement requires getting these items exactly right, in priority order:

1. **The coordinate-vs-queryname behavioral split** — two fundamentally different code paths for supplementary/secondary/unmapped mate flagging, which many reimplementations miss entirely.

2. **Tie-breaking by file order** — non-deterministic by nature; the Rust replacement must process records in exactly the same traversal order as Picard's `SortingCollection` to match output.

3. **The same-position RF → FR normalization** — a single conditional that silently merges two orientation groups, easily overlooked but essential for correct pair grouping.

4. **N operator inflating reverse-strand 5ʹ positions** — a side effect of how `getAlignmentEnd()` sums all reference-consuming operators, with major implications for RNA-seq data.

5. **Fragment-vs-pair priority** — unpaired fragments unconditionally lose to pairs at the same position, regardless of quality scores.

6. **Score type overflow** — the `short` return type silently wraps for long reads, and the optical duplicate coordinate `short` overflow in older versions creates reproducible-but-wrong results.

7. **Metrics counter division** — `READ_PAIRS_EXAMINED` and `READ_PAIR_DUPLICATES` are counted per-read during traversal and divided by 2 at finalization, a detail that affects every metric downstream.

8. **READ_NAME_REGEX colon-split optimization** — the default optical detection does not use regex at all, and failure to match silently disables optical detection rather than erroring.