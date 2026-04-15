# Deviations from Picard MarkDuplicates

This document enumerates every known place where **markdup-wea** behaves
differently from `picard.jar MarkDuplicates` (v3.x reference). Each entry
describes the scenario, Picard's behavior, our behavior, the rationale, and
the downstream-impact assessment. Maintainers adopting this tool should read
this document end-to-end before deciding whether it meets their needs.

The guiding principle: **we reproduce Picard's output, not its bugs — except
when reproducing a bug is required for byte-identical duplicate flags.**

---

## Intentional deviations

### 1. `READ_PAIR_OPTICAL_DUPLICATES` always reports `0`

- **Picard:** Enabled by default. Splits QNAME on colons, treats `tile:x:y`
  coordinates as optical-cluster location, flags intra-cluster duplicates
  separately from PCR duplicates.
- **Ours:** Not implemented. Column is always `0`.
- **Why:** nf-core/rnaseq disables optical-duplicate detection by default
  (no `--READ_NAME_REGEX`). Implementing it would ~double the scope of
  Pass 1, add a QNAME parser, and contribute zero information on the
  workflows this tool targets.
- **Impact:** Metrics column is conservative — zero optical dups never
  causes a downstream false negative. MultiQC renders the column
  without error.
- **If needed later:** straightforward to add, localized to
  `PairedGroupTracker::resolve_group`.

### 2. Queryname-sorted input is rejected

- **Picard:** Supports `queryname` sort input with a different Pass-1
  traversal; also has different semantics around supplementary-record
  flagging in this mode.
- **Ours:** `validate_sort_order` rejects anything except `coordinate`
  with a clear error, pointing the user to `samtools sort`.
- **Why:** nf-core/rnaseq delivers coordinate-sorted BAM. Supporting
  both paths would double our test matrix for a fraction of downstream
  users, and Picard's queryname path has its own set of edge cases
  (different 0x400 behavior on supplementaries).
- **Impact:** An upstream `samtools sort -o ...` is required for
  non-coordinate-sorted inputs. No silent wrong answer — we fail fast.

### 3. `ADD_PG_TAG_TO_READS` not implemented (no `PG:Z` tag on records)

- **Picard:** Defaults to `true`; writes `PG:Z:MarkDuplicates` on every
  output record.
- **Ours:** Not written.
- **Why:** Adds ~16% to BAM file size on typical rnaseq data. nf-core
  downstream consumers do not inspect record-level PG tags; the
  `@PG` header line is what tooling keys on.
- **Impact:** `samtools view -H` still shows the @PG header. Per-record
  lineage is lost. No known downstream consumer affected.

### 4. `CLEAR_DT` not implemented (pre-existing DT tags preserved)

- **Picard:** Defaults to `true`; strips `DT:Z:...` tags from input before
  writing output.
- **Ours:** Preserves whatever was in the input.
- **Why:** We don't write DT tags ourselves, and stripping them would
  require walking every record's data field. The only way input has DT
  is a prior Picard run — in which case the tags are stale anyway, but
  harmless.
- **Impact:** None unless the input was already Picard-processed. In
  that rare case, DT may be misleading; prefer running from a clean
  aligner output.

### 5. Missing `RG` tag → `library_idx = 0`, not NPE

- **Picard:** NullPointerException / `ClassCastException` on records
  lacking the `RG` field.
- **Ours:** Assigns library index `0` (the default library), warns via
  the orphan/telemetry path.
- **Why:** Strict improvement; no reason to crash on malformed input
  when a sensible default exists.
- **Impact:** For single-library workflows (the overwhelming common
  case), zero divergence. For multi-library workflows with
  mis-annotated reads, those reads get pooled into the default
  library's dup-group space; documented here so users know to
  validate their RG annotation.

### 6. Score type is `u32`, not Picard's `short` with silent overflow

- **Picard:** Accumulates `SUM_OF_BASE_QUALITIES` in a `short`, wrapping
  at 32,767 silently — so a read of length > ~440bp with uniformly high
  quality produces a score that wraps negative and loses against
  shorter reads.
- **Ours:** `u32`, saturating at ~4.3 billion (unreachable for any
  realistic read).
- **Why:** Picard's overflow is a documented bug, not an algorithmic
  choice. Reproducing it would require deliberately truncating.
- **Impact:** For Illumina 100–150bp reads, scores are always < 32,767
  and results are identical. For long reads (PacBio/ONT passed through
  an aligner), our tie-breaking may differ from Picard on reads whose
  true score exceeds 32,767. If you are processing long-read data
  through markdup-wea and want byte-identical Picard output, this is
  the one place to expect differences.

### 7. MQ=0 chimeric-pair missed-duplicate quirk (Picard issue #1285) — reproduced

- **Picard:** A documented bug where chimeric pairs with one mate at
  MAPQ=0 can be missed during duplicate detection, because Picard does
  not treat MQ=0 as unmapped.
- **Ours:** Matches. We also do not special-case MQ=0.
- **Why:** Intentional match — reproducing the quirk is required for
  byte-identical flagged-QNAME output against a Picard reference.
- **Impact:** Same as Picard. If you want MQ=0 treated as unmapped,
  neither tool does it — file a shared ticket upstream.

---

## Approximations (small divergences from Picard's exact algorithm)

### 8. Orphan-vs-fragment same-locus competition

- **Picard (research §7):** An orphan paired read (its mate never
  arrived; unresolved at EOF) is pushed into `fragSort` and competes
  against same-locus fragments by quality score, just like any other
  fragment. The orphan may win, lose, or tie depending on its score.
- **Ours:** Our Pass 1 eagerly inserts a paired-presence marker at
  every paired-primary read's single-end locus (see A2 in the plan).
  When the read's mate never arrives, that marker stays in place and
  unconditionally flags any co-located fragment — the orphan never
  competes by score.
- **Why:** In a streaming coordinate-sorted pass we cannot
  retroactively modify group resolution. Eager insertion lets us
  honor Picard §4 (pairs always beat fragments at a locus) without a
  global fragSort pass. The trade-off is orphans behave as if they
  had "won" fragment competition at their locus.
- **Impact:** Bounded by the orphan rate. In well-behaved
  nf-core/rnaseq BAMs, orphan rate is <0.1% and
  orphan-co-located-with-fragment rate is vanishingly small. On
  datasets with heavy chimerism or aggressive upstream filtering,
  this could flag slightly more fragments than Picard does.
- **Detection:** if Phase C byte-diff against Picard shows
  fragment-only false positives, inspect whether those positions
  harbor orphans in the input.

### 9. `estimate_library_size` bail-out when unique ≥ total pairs

- **Picard:** Throws `IllegalStateException`, bubbling up as a fatal
  error.
- **Ours:** Returns `None`, writes an empty `ESTIMATED_LIBRARY_SIZE`
  field (same as Picard's null-serialization path for zero dups).
- **Why:** Graceful handling of an arithmetically impossible input is
  strictly better than crashing. The empty-field serialization is
  indistinguishable to MultiQC from Picard's zero-dup empty field.
- **Impact:** Users running a pathological input see an empty field
  instead of an exception. No silent wrong answer.

---

## Out-of-scope (not implemented, not planned for MVP)

- Multi-threaded Pass 1 / Pass 2 — Phase 4 in tasks.md, deferred.
- `TAGGING_POLICY` / `DT` tag emission — not used by nf-core/rnaseq.
- `BARCODE_TAG`, `READ_ONE_BARCODE_TAG`, `READ_TWO_BARCODE_TAG` — UMI
  handling belongs in `umi_tools`, not MarkDuplicates.
- `MAX_FILE_HANDLES_FOR_READ_ENDS_MAP` and other tuning knobs — we
  use an in-memory BTreeMap with position-keyed incremental
  resolution; no disk spill.
- `DUPLICATE_SCORING_STRATEGY` other than `SUM_OF_BASE_QUALITIES`.

---

## How to read this document

If you are evaluating whether markdup-wea is safe to substitute for
Picard in your pipeline, the acceptance test is:

1. Entries 1–7 (intentional) — read and confirm the rationale applies
   to your data.
2. Entries 8–9 (approximations) — confirm the impact bounds are
   acceptable for your workflow.
3. Phase C byte-diff against Picard on a representative sample of
   your own data — the ultimate ground-truth test. See `BENCHMARK.md`
   for the procedure.
