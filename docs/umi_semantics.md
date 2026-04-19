# UMI / BARCODE_TAG semantics — design & parity spec

Source-of-truth design doc for Track A (UMI feature work). Every semantic
claim is backed by a Picard 3.4.0 source citation. Implementers of A.1–A.3
should not need to re-read Picard source — this doc is the contract.

Status: **A.0 done** (2026-04-19). Phase 4 Hetzner probe ran GREEN — all 5
fixtures confirmed §3 findings empirically. See §9.1 for results.

---

## 1. Scope

**In scope for Track A:**
- `BARCODE_TAG` (single UMI tag shared by both mates; e.g. `RX`)
- `READ_ONE_BARCODE_TAG` / `READ_TWO_BARCODE_TAG` (per-mate tags; e.g. `BX`)
- `MOLECULAR_IDENTIFIER_TAG` (MI tag emission on output records)

**Deferred to A.7 (not in A.1–A.6):**
- `DUPLEX_UMI` — see §7 for the source-cited reason.

**Out of scope for markdup-wea entirely (unchanged):**
- Optical duplicates (`READ_NAME_REGEX` / `OPTICAL_DUPLICATE_PIXEL_DISTANCE`)
- `UmiAwareMarkDuplicates` UMI consensus / assignment (separate Picard tool)
- `EstimateLibraryComplexityWithBarcodes` (separate Picard tool; we reuse
  the same `getReadBarcodeValue` hash function for parity with
  `MarkDuplicates BARCODE_TAG=…`).

---

## 2. Picard sources (3.4.0)

| File | Role | Permalink |
|---|---|---|
| `MarkDuplicates.java` | CLI args, `buildReadEnds`, pair-merge, `areComparableForDuplicates`, MI emission | [picard@3.4.0:MarkDuplicates.java](https://github.com/broadinstitute/picard/blob/3.4.0/src/main/java/picard/sam/markduplicates/MarkDuplicates.java) |
| `util/ReadEndsForMarkDuplicatesWithBarcodes.java` | 3-int field layout | [picard@3.4.0:util/ReadEndsForMarkDuplicatesWithBarcodes.java](https://github.com/broadinstitute/picard/blob/3.4.0/src/main/java/picard/sam/markduplicates/util/ReadEndsForMarkDuplicatesWithBarcodes.java) |
| `util/ReadEndsForMarkDuplicatesWithBarcodesCodec.java` | On-disk field write order | [picard@3.4.0:util/ReadEndsForMarkDuplicatesWithBarcodesCodec.java](https://github.com/broadinstitute/picard/blob/3.4.0/src/main/java/picard/sam/markduplicates/util/ReadEndsForMarkDuplicatesWithBarcodesCodec.java) |
| `UmiUtil.java` | `getTopStrandNormalizedUmi`, `setMolecularIdentifier`, `getStrand`, regex | [picard@3.4.0:UmiUtil.java](https://github.com/broadinstitute/picard/blob/3.4.0/src/main/java/picard/sam/markduplicates/UmiUtil.java) |
| `EstimateLibraryComplexity.java` | `getReadBarcodeValue` (shared hash for BC/RX/BX/read1/read2) | [picard@3.4.0:EstimateLibraryComplexity.java](https://github.com/broadinstitute/picard/blob/3.4.0/src/main/java/picard/sam/markduplicates/EstimateLibraryComplexity.java) |

Anchor line numbers below reference these 3.4.0 snapshots.

---

## 3. Semantics — the 9 questions

### Q1. Hash function

**Claim:** Picard stores three barcodes as `int` values computed by two distinct
functions:
- `barcode` ← `Objects.hash(topStrandNormalizedUmi)` where `topStrandNormalizedUmi`
  is either the raw tag string or `null`.
- `readOneBarcode` / `readTwoBarcode` ← `attr.hashCode()` directly (no
  `Objects.hash` wrapping).

**Primary citation:**
- `MarkDuplicates.java:695-696` — `final String topStrandNormalizedUmi =
  UmiUtil.getTopStrandNormalizedUmi(rec, BARCODE_TAG, DUPLEX_UMI);
  endsWithBarcode.barcode = Objects.hash(topStrandNormalizedUmi);`
- `EstimateLibraryComplexity.java:406-411` — `public static int
  getReadBarcodeValue(final SAMRecord record, final String tag) {
  if (null == tag) return 0; final String attr =
  record.getStringAttribute(tag); if (null == attr) return 0; else return
  attr.hashCode(); }`

**Cross-check:** `ReadEndsForMarkDuplicatesWithBarcodesCodec.java:51-53`
writes exactly these three int fields in order (`barcode`,
`readOneBarcode`, `readTwoBarcode`) — confirms the typing.

**Java details that matter for the Rust port:**
- `Objects.hash(x)` for a single argument is **NOT** `x.hashCode()`. Per
  JDK javadoc: `Objects.hash(x)` = `Arrays.hashCode(new Object[]{x})` =
  `31 * 1 + (x == null ? 0 : x.hashCode())` = `31 + stringHash`.
- `Objects.hash(null)` = `31`. `Objects.hash("")` = `31 + 0 = 31`. So
  `barcode=31` for both null and empty — an **asymmetry** with
  `readOneBarcode=0` / `readTwoBarcode=0` on missing tag.
- `String.hashCode()` is the classic polynomial-31 over 16-bit chars:
  `h = 31*h + charAt(i)`. Since `ALLOWED_UMI` restricts UMI bytes to ASCII
  (see Q2), the Rust port can treat the tag as bytes and use the same
  formula on `u8 as i32`.

**Rust port (canonical):**
```rust
fn java_string_hashcode(s: &[u8]) -> i32 {
    let mut h: i32 = 0;
    for &c in s { h = h.wrapping_mul(31).wrapping_add(c as i32); }
    h
}

/// BARCODE_TAG path (MarkDuplicates.java:695-696).
fn picard_barcode_hash(normalized_umi: Option<&[u8]>) -> i32 {
    let inner = normalized_umi.map_or(0, java_string_hashcode);
    31i32.wrapping_add(inner)
}

/// READ_ONE/READ_TWO_BARCODE_TAG path (EstimateLibraryComplexity.java:406-411).
fn read_barcode_value(tag_attr: Option<&[u8]>) -> i32 {
    tag_attr.map_or(0, java_string_hashcode)
}
```

**Pitfall:** Using FxHash, xxHash, or any non-Java hash guarantees silent
parity divergence — no panic, just wrong dup counts. Must replicate Java
`String.hashCode` byte-for-byte. Must preserve the `Objects.hash` `31+`
seed for the `barcode` field specifically.

---

### Q2. Normalization

**Claim:** The UMI is passed **verbatim** (case-preserving, no trimming, no
dash-collapse). Only characters `[ATCGNatcgn-]` are allowed; anything else
throws `PicardException`. `AGCT` and `agct` hash to DIFFERENT values.

**Primary citation:**
- `UmiUtil.java:48` — `static final Pattern ALLOWED_UMI =
  Pattern.compile("^[ATCGNatcgn-]*$");`
- `UmiUtil.java:72-74` — `if (!ALLOWED_UMI.matcher(umi).matches()) {
  throw new PicardException("UMI found with illegal characters…"); }`
- `UmiUtil.java:76-78` — `if (!duplexUmi) { return umi; }` — with
  DUPLEX_UMI=false (our MVP), the string is returned unmodified.

**Cross-check:** `EstimateLibraryComplexity.java:406-411` hashes the tag
string directly with no pre-processing (`attr.hashCode()`), matching the
"verbatim" claim for READ_ONE/READ_TWO_BARCODE_TAG.

**Dual-UMI form (`AGCT-AGCT`):** With DUPLEX_UMI=false, returned unchanged
(the `-` is hashed as a regular character). With DUPLEX_UMI=true, the UMI
is split on `-` into two parts and, for bottom-strand reads, swapped; see
Q7. **DUPLEX_UMI is deferred**, so for A.1–A.6 treat the UMI verbatim.

**Pitfall:** Do NOT uppercase, trim, or strip `-`. A helpful pre-normalizer
silently breaks parity.

**A.0 decision:** markdup-wea will **match Picard exactly** — verbatim
value, validate against `^[ATCGNatcgn-]*$`, error on mismatch. No
deviations.md entry needed.

---

### Q3. Missing-tag policy

**Claim:** "Missing tag" has **two different fallbacks** in Picard depending
on which field:
- `BARCODE_TAG` missing (or value is null) → `topStrandNormalizedUmi = null`
  → `Objects.hash(null) = 31`. Stored as `barcode = 31`.
- `READ_ONE_BARCODE_TAG` / `READ_TWO_BARCODE_TAG` missing or tag
  unspecified (CLI arg null) → `getReadBarcodeValue` returns `0`.

Missing ≠ empty-string for the barcode field specifically: `Objects.hash("")
= 31 + 0 = 31`, coincidentally equal to `Objects.hash(null) = 31`. So empty
and missing collide for `barcode`, which is benign. But **missing barcode
(=31) does NOT equal missing readOne/readTwo (=0)**, which matters if we
ever cross-compare fields.

**Primary citation:**
- `UmiUtil.java:62-64` — `if (umiTag == null) { return null; }`
- `UmiUtil.java:68-70` — `if (umi == null) { return null; }` (record
  lacks the tag).
- `MarkDuplicates.java:696` — `Objects.hash(null)` = 31 (JDK javadoc).
- `EstimateLibraryComplexity.java:406-411` — two early returns of `0`.

**Cross-check:** `ReadEndsForMarkDuplicatesWithBarcodes.java:28-30`
initializes all three fields to `0`, so "no useBarcodes mode" matches
"useBarcodes + all READ_ONE/TWO missing" for read_one/two but NOT for
`barcode` (which would be 31 when useBarcodes and BARCODE_TAG is set but
the record lacks the tag). In practice `useBarcodes` can only be true
when at least one of BARCODE_TAG / READ_ONE_BARCODE_TAG / READ_TWO_BARCODE_TAG
is non-null (`MarkDuplicates.java:261`), so the `0` defaults are only
observed when the user didn't ask for the barcode dimension at all.

**Pitfall:** Treating missing BARCODE_TAG as "same group as empty" via
`hash("")=0` will be wrong — Picard's `Objects.hash("") = 31` not 0.
Similarly, substituting `0` for missing BARCODE_TAG breaks parity when
the user's data has a mix of records with and without RX.

---

### Q4. Key composition (equality)

**Claim:** When `useBarcodes=true`, a duplicate group's equality is
extended additively: the usual (lib, ref1, pos1, orient, ref2, pos2) are
AND-ed with `barcode == barcode && readOneBarcode == readOneBarcode &&
readTwoBarcode == readTwoBarcode`. All three barcode ints participate in
equality.

**Primary citation:**
- `MarkDuplicates.java:810-819` — `areComparableForDuplicates(...)` — `if
  (useBarcodes && areComparable) { ... areComparable = lhsWithBarcodes.barcode
  == rhsWithBarcodes.barcode && lhsWithBarcodes.readOneBarcode ==
  rhsWithBarcodes.readOneBarcode && lhsWithBarcodes.readTwoBarcode ==
  rhsWithBarcodes.readTwoBarcode; }`

**Cross-check (sort comparator, different purpose):**
`MarkDuplicates.java:984-1008` — `ReadEndsMDComparator` also compares
barcodes in the same order (`barcode`, then `readOneBarcode`, then
`readTwoBarcode`). This is for Picard's disk-backed `SortingCollection`
traversal order; markdup-wea uses equality-based grouping
(`FxHashMap<PairedEndKey, Vec<ScoredPair>>` at
[src/groups.rs:57](src/groups.rs:57)), so **the sort comparator does not
translate to a required code change** — but the equality extension does.

**Codec confirms field inclusion independently:**
`ReadEndsForMarkDuplicatesWithBarcodesCodec.java:51-53` writes all three
`int`s after the parent ReadEnds fields. If any field weren't part of the
logical record, the codec wouldn't persist it.

**A.1 action:** Extend [src/groups.rs:13-21](src/groups.rs:13)
`PairedEndKey` with `barcode: i32`, `read_one_barcode: i32`,
`read_two_barcode: i32` (all default `0`). Extend
[src/groups.rs:26-30](src/groups.rs:26) `SingleEndKey` with `barcode: i32`
only (single-end reads have no read1/2 concept — Picard's
`buildReadEnds` sets `readOneBarcode` for unpaired reads too, see Q5, so
we may need to reconsider — actually unpaired-as-fragment should only use
the `barcode` dimension because that's all MarkDuplicates consults for
fragment comparisons in `areComparableForDuplicates`; readOne/Two exist
on the object but are only compared once paired equality is reached. For
safety, carry `read_one_barcode` on `SingleEndKey` too and set it per
Q5).

**Clarification on SingleEndKey scope:** `areComparableForDuplicates` is
called for both fragments and pairs (different `compareRead2` arg). The
barcode equality block runs regardless of fragment-vs-pair — so yes,
unpaired fragments compare `readOneBarcode` and `readTwoBarcode` too.
Per Q5, only `readOneBarcode` will be set on unpaired fragments (from the
if-branch at `MarkDuplicates.java:698`); `readTwoBarcode` stays 0. So
SingleEndKey must also carry `read_one_barcode: i32` (default 0) and
`read_two_barcode: i32` (always 0 for unpaired) to match
`areComparableForDuplicates` exactly. In practice `readTwoBarcode` adds
nothing to equality for unpaired reads (constant 0 on both sides), but
we include it for correctness if two different single-end inputs from
the same library somehow had distinct readTwo via `READ_TWO_BARCODE_TAG`
(they can't — `buildReadEnds` only writes `readTwoBarcode` for
second-of-pair, line 700-701).

**Simplification for A.1 (accepted):** Omit `read_two_barcode` from
`SingleEndKey` since Picard never writes it for non-paired reads.
Document this as an implementation note, not a deviation — both tools
have 0 on both sides.

**Pitfall:** Forgetting to include barcodes in `PartialEq` / `Hash` of
`PairedEndKey` — the compiler won't catch this because the existing
derives need manual extension to cover new fields, and missing a field
silently degrades the key (all barcodes alias to 0).

---

### Q5. READ_ONE vs READ_TWO assignment

**Claim:** `readOneBarcode` is always set from the record whose
`firstOfPair` flag is TRUE; `readTwoBarcode` is from the `secondOfPair`
record. Assignment is **coupled to the BAM FLAG, not to lo/hi coordinate
ordering**.

**Primary citation:**
- `MarkDuplicates.java:698-702` (inside `buildReadEnds`, per-record
  processing): `if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
  endsWithBarcode.readOneBarcode = getReadOneBarcodeValue(rec); } else {
  endsWithBarcode.readTwoBarcode = getReadTwoBarcodeValue(rec); }`
- `MarkDuplicates.java:574-584` (pair-merge when the mate arrives): `if
  (rec.getFirstOfPairFlag()) { ... ((...)pairedEnds).readOneBarcode =
  getReadOneBarcodeValue(rec); } else { ... readTwoBarcode =
  getReadTwoBarcodeValue(rec); }`

Both assignment sites gate on `firstOfPair`, not on which read is
lo-coordinate vs hi-coordinate.

**Cross-check:** The same pair-merge block at
`MarkDuplicates.java:586-613` *separately* computes lo/hi ordering (via
`matesRefIndex` / `matesCoordinate` comparisons) and may swap `read1*`
↔ `read2*` fields. But nowhere in that block does it touch
`readOneBarcode` / `readTwoBarcode`. The two concerns are orthogonal in
Picard.

**Implication for markdup-wea:** The same-position RF→FR normalization at
[src/scan.rs:270-279](src/scan.rs:270) (which currently mirrors Picard
lo/hi flip) does **not** need to swap barcodes. Our `PendingMate` must
carry the pending record's barcode values tagged by `firstOfPair`-ness,
not by lo/hi-ness.

**Pitfall:** Assigning `readOneBarcode` to the lo-coordinate mate breaks
parity on reversed pairs (same-locus FR/RF swap) where the firstOfPair
read happens to land at the hi coordinate.

**A.1 implementation note:** Store `is_first_of_pair: bool` alongside
`barcode` values when constructing `PendingMate`. When the mate arrives,
assign readOneBarcode/readTwoBarcode using the firstOfPair flag, not
lo/hi ordering.

---

### Q6. MI tag (MOLECULAR_IDENTIFIER_TAG)

**Claim:** MI is **emission-only** — it is set on output records via
`rec.setAttribute(...)` and is **not** part of any duplicate grouping key
or comparator. The emitted value is `"{contig}:{pos}/{assignedUmi}"` with
optional `/A` or `/B` duplex suffix. In plain `MarkDuplicates`,
`assignedUmi` is always the empty string, so the emitted MI tag is
`"{contig}:{pos}/"` literally — this is the "empty-UMI quirk".

**Primary citation:**
- `MarkDuplicates.java:404-407` — sole call site:
  `if (BARCODE_TAG != null) { UmiUtil.setMolecularIdentifier(rec, "",
  MOLECULAR_IDENTIFIER_TAG, DUPLEX_UMI); }`
- `UmiUtil.java:137-169` — function body constructing the MI string.
- `UmiUtil.java:139-141` — `if (molecularIdentifierTag == null) { return; }`
  — MI only set when user passed `MOLECULAR_IDENTIFIER_TAG=...`.
- `UmiUtil.java:144-148` — `molecularIdentifier.append(rec.getContig())`;
  then `:`; then `rec.getReadNegativeStrandFlag() ? rec.getAlignmentStart()
  : rec.getMateAlignmentStart()`; then `/`; then `assignedUmi`.

**Cross-check (no grouping):** `MarkDuplicates.java:810-819`
(`areComparableForDuplicates`) compares only the three `int` barcode
fields — never touches `MOLECULAR_IDENTIFIER_TAG` or its value. The sort
comparator at `MarkDuplicates.java:984-1008` likewise does not reference
MI. Confirmed: MI is orthogonal to duplicate detection.

**Position chosen:** `reverseStrand ? alignmentStart : mateAlignmentStart`.
For reverse-strand reads, MI records the read's OWN aligned start; for
forward-strand reads, MI records the MATE's aligned start. This is
Picard's canonical "pair's representative position" convention.

**Gating:** MI is only emitted when `BARCODE_TAG != null` (so without
`BARCODE_TAG`, setting `MOLECULAR_IDENTIFIER_TAG` alone is a no-op —
confirmed at `MarkDuplicates.java:405`).

**DUPLEX_UMI suffix:** When `DUPLEX_UMI=true` and `UmiUtil.getStrand`
returns TOP/BOTTOM, `/A` or `/B` is appended. UNKNOWN → no suffix
(`UmiUtil.java:150-166`). Since DUPLEX_UMI is deferred, A.1–A.6 emit no
suffix.

**A.3 implementation:** Emit MI as a side-effect in Pass 2 (output
writer), not in Pass 1 (scan). Requires contig name lookup from the
header (we have it) and the mate's aligned position for forward-strand
reads (we have `mate_pos` in BAM's MPOS field, aligned — no unclipped
adjustment needed since Picard uses `getMateAlignmentStart`, not
`unclipped`).

**Pitfall:** Setting MI during Pass 1 (scan) conflates two passes; MI
belongs in Pass 2 (write). Also, reading `getMateUnclippedStart` instead
of `getMateAlignmentStart` would require MC tag parsing — do NOT do
that; use MPOS directly.

**Format clarification — 1-based vs 0-based:** Picard's
`rec.getAlignmentStart()` returns **1-based** (SAM spec). Our internal
pos is 0-based (noodles convention). When emitting MI, convert back to
1-based (`pos + 1`).

---

### Q7. DUPLEX_UMI verdict

**DUPLEX_UMI: DEFER to A.7.**

**Source-cited reason:** Duplex-strand assignment in
`UmiUtil.getStrand(rec)` at `UmiUtil.java:104-128` depends on the
**mate's unclipped 5' coordinate**:
```java
final int mate5PrimeStart = (rec.getMateNegativeStrandFlag())
    ? SAMUtils.getMateUnclippedEnd(rec)
    : SAMUtils.getMateUnclippedStart(rec);
```
`SAMUtils.getMateUnclippedStart/End` requires the **MC tag** (mate CIGAR)
on the record to correctly adjust for soft/hard clipping. Without MC, it
silently returns the aligned start — producing wrong strand assignments
for clipped mates.

markdup-wea's scan does NOT currently parse the MC tag. Adding DUPLEX_UMI
requires one of:
1. Parse MC tag and compute mate unclipped 5' per-record (new code in
   `src/scan.rs`), OR
2. Defer `barcode` computation for paired-end reads until both mates are
   seen (we already know each mate's own unclipped 5' by then), breaking
   Picard's per-record `buildReadEnds` structure, OR
3. Require inputs to carry MC tags (nf-core rnaseq may not guarantee this).

None is a one-hour change. Pilot datasets for Track A (GSE75823, GSE134031,
PRJNA416930, PRJNA515063) are all **non-duplex** (single-strand UMIs);
DUPLEX_UMI parity is not blocking.

**A.7 work (when scheduled):** Decide path 1 vs path 2, extend
`PendingMate` to carry `mate_unclipped_5prime` (path 2) or add
MC-tag-based unclipped computation to `src/scan.rs` (path 1), add
`getStrand` helper, extend `picard_barcode_hash` to invoke the
top-strand-normalization swap for bottom-strand reads.

**A.0 action:** Add DUPLEX_UMI to `docs/deviations.md` under
"Out-of-scope" with the MC-tag reason.

---

### Q8. RF→FR flip + barcodes

**Claim:** Our same-locus RF→FR normalization at
[src/scan.rs:270-279](src/scan.rs:270) does **not** need to swap barcodes.
Barcodes are assigned by `firstOfPair` flag (Q5), independent of lo/hi
ordering and orientation.

**Primary citation:** Same as Q5 — `MarkDuplicates.java:574-584` and
`MarkDuplicates.java:698-702`. Neither site references the lo/hi flip at
`MarkDuplicates.java:586-613`.

**Picard's own RF→FR:** `MarkDuplicates.java:599-603` — sets
`pairedEnds.orientation = ReadEnds.FR` when same-pos RF — but does not
touch `readOneBarcode` / `readTwoBarcode`.

**A.1 action:** Leave the flip at [src/scan.rs:270-279](src/scan.rs:270)
unchanged when extending keys with barcodes; only add barcode fields to
the output `PairedEndKey` construction, sourced from `firstOfPair`-tagged
pending/current-record values.

**Pitfall:** Assuming barcodes follow lo/hi — would silently diverge on
the subset of same-pos RF pairs.

---

### Q9. Metrics impact

**Claim:** Enabling `BARCODE_TAG` does NOT change the shape of the
`DuplicationMetrics` output — same columns, same histogram. What changes
is the **content** (more groups, higher PERCENT_DUPLICATION on non-UMI
data, different ESTIMATED_LIBRARY_SIZE). Picard emits `LIBRARY` rows
unchanged; there is no barcode-partitioned metrics view in
`MarkDuplicates`. For that, users run `EstimateLibraryComplexityWithBarcodes`
separately.

**Primary citation:** `MarkDuplicates.java:404-407` shows MI emission
gated on BARCODE_TAG but no metrics gating. The metrics output in
`MarkDuplicates.finalizeMetrics` / `writeMetrics` paths (not read in detail
here but referenced for completeness) does not branch on `useBarcodes`.

**Cross-check:** Picard has a separate class
`EstimateLibraryComplexityWithBarcodes` for per-barcode library size —
existence of a separate tool is itself evidence that `MarkDuplicates`
itself does not do barcode-partitioned metrics.

**A.1–A.3 action:** No changes to `src/metrics.rs`. Histogram
accumulation naturally reflects the new group structure because keys
changed.

**Pitfall:** Accidentally splitting LIBRARY rows by barcode — diverges
from Picard.

---

## 4. CLI flags (verbatim from Picard)

From `MarkDuplicates.java:178-218`:

| Picard option | Default | Our CLI flag (proposed A.2) | Notes |
|---|---|---|---|
| `BARCODE_TAG` | null | `--barcode-tag <TAG>` | 2-char SAM tag name (e.g. `RX`, `BC`). When set, `useBarcodes=true`. |
| `READ_ONE_BARCODE_TAG` | null | `--read-one-barcode-tag <TAG>` | Per-mate tag; read from firstOfPair only. |
| `READ_TWO_BARCODE_TAG` | null | `--read-two-barcode-tag <TAG>` | Per-mate tag; read from secondOfPair only. |
| `MOLECULAR_IDENTIFIER_TAG` | null | `--molecular-identifier-tag <TAG>` | Output MI tag. No-op unless `--barcode-tag` also set (matches Picard:405). |
| `DUPLEX_UMI` | false | `--duplex-umi` | **Deferred to A.7**; CLI flag not added in A.2. If added prematurely, it must error unless A.7 is landed. |

All four flags default to "feature off" — markdup-wea without any of
them is byte-identical to today's behavior (Track A should not regress
the 8/8 ENCODE parity gate).

---

## 5. Implementation hooks (for A.1–A.3)

| File | Lines | Hook |
|---|---|---|
| [src/groups.rs:13-21](src/groups.rs:13) | `PairedEndKey` | append `barcode: i32, read_one_barcode: i32, read_two_barcode: i32` (all default 0 when feature off). Add to `PartialEq`/`Eq`/`Hash` derives. |
| [src/groups.rs:26-30](src/groups.rs:26) | `SingleEndKey` | append `barcode: i32, read_one_barcode: i32` (no read_two for unpaired; see Q4 clarification). |
| [src/pending_mates.rs:10-21](src/pending_mates.rs:10) | `PendingMate` | add `barcode: i32`, `read_one_barcode: i32`, `read_two_barcode: i32`, `is_first_of_pair: bool`. Barcode values here are the PENDING record's own computed hashes (so when the mate arrives, we combine). |
| [src/scan.rs:79-97](src/scan.rs:79) | clone `extract_read_group` pattern | new helpers: `extract_barcode_tag(record, tag_name: &[u8;2]) -> Option<Vec<u8>>` + pure-fn `picard_barcode_hash(bytes) -> i32` + `read_barcode_value(bytes) -> i32`. |
| [src/scan.rs:218-260](src/scan.rs:218) | single-end key construction | compute `barcode` from BARCODE_TAG, `read_one_barcode` from READ_ONE_BARCODE_TAG (see Q5 — unpaired reads get readOne only, readTwo stays 0). |
| [src/scan.rs:262-292](src/scan.rs:262) | paired-end mate-merge | compute own barcode values; combine with pending's via firstOfPair flag (Q5); barcode swap on lo/hi flip is NOT required (Q8). |
| [src/scan.rs:270-279](src/scan.rs:270) | RF→FR flip | unchanged (Q8). |
| [src/main.rs:11-34](src/main.rs:11) | `Cli` struct | add 4 flags per §4 above. |
| [src/markdup.rs](src/markdup.rs) | Pass 2 writer | A.3: emit MI tag via `record.data_mut().insert(...)` when `--molecular-identifier-tag` set AND `--barcode-tag` set (gated per Q6). Format: `{contig}:{pos_1based}/`. |
| [src/scoring.rs](src/scoring.rs) | confirm-only | no changes; Picard scoring is barcode-agnostic (memory confirmed, Picard source at `MarkDuplicates.java:646-648` only uses `DuplicateScoringStrategy.computeDuplicateScore(rec, strategy)`). |
| [src/io.rs](src/io.rs) | confirm-only | no changes; B.0's AlignmentReader enum already present. |

---

## 6. CLI flag names (verbatim)

Per `MarkDuplicates.java:179-218`:
- `BARCODE_TAG` (doc: "Barcode SAM tag (ex. BC for 10X Genomics)")
- `READ_ONE_BARCODE_TAG` (doc: "Read one barcode SAM tag (ex. BX for 10X Genomics)")
- `READ_TWO_BARCODE_TAG` (doc: "Read two barcode SAM tag (ex. BX for 10X Genomics)")
- `MOLECULAR_IDENTIFIER_TAG` (doc: "SAM tag to uniquely identify the molecule from which a read was derived. Use of this option requires that the BARCODE_TAG option be set to a non null value. Default null.")

Our CLI should mirror the semantics in hyphenated form. Lift doc strings
into clap `help` attributes verbatim.

---

## 7. Test-vector strategy (informational; A.4 owns execution)

- **Unit tests (A.1):** 3–4 synthetic BAMs constructed in-code via
  `noodles`. Records at identical coords, controlled RX values. Cases:
  (a) same RX → one group; (b) distinct RX → two groups; (c) one record
  missing RX → `barcode=31` group distinct from records with RX=""
  (impossible via real pipelines but exercises the asymmetry); (d)
  RX=AGCT vs rx=agct → two groups (case-sensitive).
- **Smoke fixture (A.4):** WARP `rna_with_umis` test BAMs at
  `gs://broad-gotc-test-storage/rna_with_umis/` — public, <5GB, already
  RX-tagged. Compare markdup-wea `--barcode-tag RX` output vs Picard
  `BARCODE_TAG=RX` output byte-for-byte on flagged-QNAME set + metrics.
- **Pilot dataset (A.5):** GSE75823 (Parekh 2016) — matched UMI-seq +
  TruSeq on UHRR. UMI arm tests `--barcode-tag RX` parity; TruSeq arm
  regression-tests non-UMI parity (no-op gate).

A.0 does not run any of these; A.0 only confirms the Picard jar used for
parity gate is still at 3.4.0 on Hetzner (see §9) so A.4 onward can use it.

---

## 8. DUPLEX_UMI verdict (summary)

**Decision:** DEFER to A.7.

**Reason (one line, source-cited):** `UmiUtil.getStrand` at
`UmiUtil.java:104-128` requires mate unclipped 5' position via
`SAMUtils.getMateUnclippedStart/End`, which depends on the MC tag;
markdup-wea does not currently parse MC and folding this into A.1
requires either MC parsing in scan.rs or restructuring barcode
computation to wait for the mate.

`docs/deviations.md` will gain a new "Out-of-scope" entry for
DUPLEX_UMI.

---

## 9. Phase 4 probe (not yet run — optional A.0 cross-check)

A Hetzner probe would run a synthetic ~10-record BAM with injected RX
tags through Picard 3.4.0 (same jar at
`/root/markdup-test/picard.jar` per `project_parity_gate_result.md`
memory) and confirm:

1. Two records, same coords, distinct RX → separate dup groups (confirms Q4).
2. Same coords, one missing RX vs one with RX=anything → separate groups
   (confirms Q3: missing hashes to `31` ≠ any realistic UMI hash).
3. Same coords, `RX=AGCT` vs `rx=agct` → two groups (confirms Q2 case).
4. Same coords, both `RX=AGCT-TGCA` → one group; one `RX=AGCT-TGCA`, one
   `RX=AGCT-AGCT` → two groups (confirms Q2 dash).
5. Metrics output shape unchanged vs plain run (confirms Q9).

Probe is a confidence check on top of source-reading. If it contradicts
§3, stop and re-read source. Verdict will be appended to this doc once
run.

### 9.1 Phase 4 probe results — 2026-04-19

**Verdict: GREEN. All 5 fixtures match §3 predictions. A.1 may proceed on source-reading alone.**

Setup: Hetzner cpx62 (`204.168.253.225`), Picard `/opt/picard/picard.jar` v3.4.0, samtools 1.19.2, workspace `/mnt/HC_Volume_105344878/tmp/umi-probe/`. Each fixture is a hand-rolled SAM (heredoc) with 2–3 paired records at identical coords (chr1:1000/1100, CIGAR=100M, mapq=60), differing only in RX tag. Picard run with `BARCODE_TAG=RX VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true`. Probe script kept at `/mnt/HC_Volume_105344878/tmp/umi-probe/probe.sh` for reuse by A.3 (MI) and A.7 (DUPLEX_UMI).

**Fixture A — distinct RX (Q4: barcode in equality key) — ✓ confirmed.**

| QNAME | RX | observed dup | predicted dup |
|---|---|:---:|:---:|
| aa | AAAAAAAA | 0 | 0 (representative) |
| ab | AAAAAAAA | **1** | 1 (dup of aa) |
| bb | TTTTTTTT | 0 | 0 (separate group) |

Picard kept `aa` as representative, flagged `ab` as duplicate, left `bb` unique. Same-coord reads with different RX → not grouped. Confirms barcode participates in the equality key.

**Fixture B — missing RX (Q3: missing-tag deterministic, hashes to 31) — ✓ confirmed.**

| QNAME | RX | observed dup | predicted dup |
|---|---|:---:|:---:|
| aa | AAAAAAAA | 0 | 0 (separate group from nn/nm) |
| nn | _absent_ | 0 | 0 (representative of missing-tag group) |
| nm | _absent_ | **1** | 1 (dup of nn — same missing-tag hash) |

`nn` and `nm` share a dup group despite both lacking the tag — confirms missing-tag is deterministic (`Objects.hash(null) = 31`), NOT skipped (would have grouped with `aa`) and NOT erroring. `aa` stays in its own group because hash(31 + "AAAAAAAA".hashCode()) ≠ 31. **This is the highest-pitfall finding empirically locked in.**

**Fixture C — case sensitivity (Q2: verbatim, no uppercase-fold) — ✓ confirmed.**

| QNAME | RX | observed dup | predicted dup |
|---|---|:---:|:---:|
| uc | AGCTAGCT | 0 | 0 (separate from lc) |
| lc | agctagct | 0 | 0 (separate from uc) |

Two separate dup groups — Picard preserves case in the hash input. A future Rust port MUST NOT call `.to_ascii_uppercase()` on the UMI before hashing.

**Fixture D — dual-UMI dash (Q2: dash survives in hash verbatim) — ✓ confirmed.**

| QNAME | RX | observed dup | predicted dup |
|---|---|:---:|:---:|
| d1 | AGCT-TGCA | 0 | 0 (representative) |
| d2 | AGCT-TGCA | **1** | 1 (dup of d1) |
| d3 | AGCT-AGCT | 0 | 0 (separate group) |

`d3` is NOT collapsed with `d1`/`d2` despite having one half identical (`AGCT`). Confirms the dash is part of the literal hashed string — no split-and-merge logic in the non-DUPLEX path.

**Fixture E — metrics shape invariance (Q9) — ✓ confirmed.**

Comparing `fixture_E.picard.metrics.txt` (BARCODE_TAG=RX) vs `fixture_E.plain.metrics.txt` (no BARCODE_TAG) on the same 3-pair input:

- **METRICS CLASS row:** identical (`picard.sam.DuplicationMetrics`).
- **Column header:** identical (`LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED SECONDARY_OR_SUPPLEMENTARY_RDS UNMAPPED_READS UNPAIRED_READ_DUPLICATES READ_PAIR_DUPLICATES READ_PAIR_OPTICAL_DUPLICATES PERCENT_DUPLICATION ESTIMATED_LIBRARY_SIZE`).
- **HISTOGRAM header:** identical (`BIN CoverageMult all_sets non_optical_sets`).
- **Values differ** as expected: BARCODE_TAG run shows `READ_PAIR_DUPLICATES=1, PERCENT_DUPLICATION=0.333, ESTIMATED_LIBRARY_SIZE=3`; plain run shows `READ_PAIR_DUPLICATES=2, PERCENT_DUPLICATION=0.667, ESTIMATED_LIBRARY_SIZE=1`. The value delta is itself a side-confirmation of Fixture A.

`MarkDuplicates` does not gain or lose columns when `BARCODE_TAG` is set. The separate `EstimateLibraryComplexityWithBarcodes` tool (with the per-barcode rows) is not invoked from `MarkDuplicates`.

**Implications for A.1:**

1. The `picard_barcode_hash` / `read_barcode_value` Rust snippets in §3 Q1 are correct as written — no modifications needed.
2. The asymmetric missing-tag fallback (BARCODE_TAG missing → 31, READ_ONE/TWO missing → 0) MUST be implemented exactly as documented; Fixture B proves the BARCODE_TAG side is non-zero in practice.
3. The case-preserving + dash-preserving normalization rule MUST hold. Any Rust normalization helper that uppercases or splits on `-` will diverge silently.
4. Metrics writer in `src/metrics.rs` (or wherever the writer lives) needs no schema change for A.1 — column shape is invariant under BARCODE_TAG.

---

## 10. Open questions for A.1

None blocking. The following are expected to land via implementation
details rather than semantic decisions:

1. How exactly does noodles expose `record.data().get(&Tag::OTHER(b"RX"))`
   as `Option<Value<String>>`? (Trivial; same pattern as `extract_read_group`
   at [src/scan.rs:79-97](src/scan.rs:79).)
2. When `BARCODE_TAG` is set but a record's RX fails the `^[ATCGNatcgn-]*$`
   regex, do we match Picard's throw-PicardException or downgrade to a
   warning? A.0 default: **match Picard — error out**. User can explicitly
   override via a future `--lenient-barcodes` flag if needed.
3. Do we validate that the specified tag is exactly 2 characters (SAM
   spec)? A.0 default: **yes, validate at CLI parse time**.

---

## 11. Verification checklist for A.0 Done

- [x] `docs/umi_semantics.md` created with all 8 TOC sections (this file).
- [x] 9 questions answered, each with ≥1 `picard@3.4.0:` permalink and ≥1
      cross-check (codec, sibling file, or CLI doc-string).
- [x] Q7 (DUPLEX_UMI) has explicit defer-verdict with source-cited reason.
- [ ] `docs/SPEC.md` updated (§12 pending — next step).
- [ ] `docs/deviations.md` updated (§12 pending — next step).
- [x] Phase 4 probe run on Hetzner — GREEN, all 5 fixtures match §3
      predictions (see §9.1, 2026-04-19).
- [ ] Linear RNA-2 closed, RNA-3 opened with this doc linked.
- [ ] Memory `project_track_a_decision.md` status line appended.

---

## Appendix — Picard source excerpts (for offline reference)

### A. Barcode field layout
```java
// ReadEndsForMarkDuplicatesWithBarcodes.java:28-30
public int barcode = 0;
public int readOneBarcode = 0;
public int readTwoBarcode = 0;
```

### B. Equality check
```java
// MarkDuplicates.java:810-819
protected boolean areComparableForDuplicates(..., useBarcodes) {
    boolean areComparable = /* lib, ref, pos, orient */;
    if (useBarcodes && areComparable) {
        areComparable = lhs.barcode == rhs.barcode
                     && lhs.readOneBarcode == rhs.readOneBarcode
                     && lhs.readTwoBarcode == rhs.readTwoBarcode;
    }
    return areComparable;
}
```

### C. barcode hash (BARCODE_TAG path)
```java
// MarkDuplicates.java:695-696, UmiUtil.java:61-91
String umi = rec.getStringAttribute(BARCODE_TAG);
// validate ^[ATCGNatcgn-]*$ or throw
String topStrandNormalizedUmi = DUPLEX_UMI ? swapOnBottom(umi) : umi;
int barcode = Objects.hash(topStrandNormalizedUmi);  // 31 + umi.hashCode() or 31 if null
```

### D. readOne/readTwo hash (per-mate tag path)
```java
// EstimateLibraryComplexity.java:406-411
public static int getReadBarcodeValue(SAMRecord record, String tag) {
    if (null == tag) return 0;
    String attr = record.getStringAttribute(tag);
    if (null == attr) return 0;
    return attr.hashCode();
}
```

### E. First-of-pair barcode assignment
```java
// MarkDuplicates.java:698-702 (per-record buildReadEnds)
if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
    endsWithBarcode.readOneBarcode = getReadOneBarcodeValue(rec);
} else {
    endsWithBarcode.readTwoBarcode = getReadTwoBarcodeValue(rec);
}

// MarkDuplicates.java:574-584 (pair-merge)
if (rec.getFirstOfPairFlag()) {
    pairedEnds.readOneBarcode = getReadOneBarcodeValue(rec);
} else {
    pairedEnds.readTwoBarcode = getReadTwoBarcodeValue(rec);
}
```

### F. MI tag emission
```java
// MarkDuplicates.java:404-407
if (BARCODE_TAG != null) {
    UmiUtil.setMolecularIdentifier(rec, "", MOLECULAR_IDENTIFIER_TAG, DUPLEX_UMI);
}

// UmiUtil.java:137-169 (body)
if (molecularIdentifierTag == null) return;
StringBuilder mi = new StringBuilder();
mi.append(rec.getContig()).append(":");
mi.append(rec.getReadNegativeStrandFlag() ? rec.getAlignmentStart() : rec.getMateAlignmentStart());
mi.append("/").append(assignedUmi);  // assignedUmi="" from plain MarkDuplicates
if (duplexUmis) {
    switch (getStrand(rec)) {
        case TOP:    mi.append("/A"); break;
        case BOTTOM: mi.append("/B"); break;
        default:     break;  // UNKNOWN: no suffix
    }
}
rec.setAttribute(molecularIdentifierTag, mi.toString());
```
