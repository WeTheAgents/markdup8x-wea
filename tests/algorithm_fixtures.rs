//! Phase B — synthetic-BAM fixture tests verifying Phase A algorithmic fixes.
//!
//! Each test constructs a minimal coordinate-sorted BAM via `BamBuilder`, runs
//! the full markdup-wea pipeline, and asserts either FLAG-0x400 QNAME sets or
//! metrics fields. One test per research finding (A1..A7) plus deviations and
//! regressions. See docs/deviations.md for divergences tested by `b03_*`.
//!
//! Naming: `bNN_<snake_case>` — run `cargo test b0` for the whole phase,
//! `cargo test b03` for a single test.

#![allow(clippy::bool_assert_comparison)]

mod common;

use common::{
    dup_qnames_set, pair_count_by_qname, parse_metrics, read_flags_by_qname_all,
    read_string_tag_by_qname, run_markdup, run_markdup_with_metrics, BamBuilder, ReadSpec,
};
use tempfile::tempdir;

// =============================================================================
// B1 — A1: RF→FR normalization when both reads of a pair sit at same unclipped_5'
// =============================================================================
//
// Verifies: src/scan.rs:225-231 (MarkDuplicates §3 same-position RF→FR fix).
//
// Construction trick: a pair with R1 at pos=100 reverse (1M) and R2 at pos=100
// forward (1M) — both reads have unclipped_5'=100 on chr1. Coordinate sort keeps
// them at alignment_start=100; write order decides which is scanned first,
// which decides whether `<=` lo/hi tie-break produces FR or RF after orientation.
//
// Pair "fr" (FR): scanner sees rev first, fwd second → current=fwd=lo,
//                 rev_lo=false, rev_hi=true → natural FR (no normalization).
// Pair "rf" (RF): scanner sees fwd first, rev second → current=rev=lo,
//                 rev_lo=true, rev_hi=false → raw RF → A1 forces to FR.
//
// Both pairs end up with the same PairedEndKey → same duplicate group →
// exactly one pair flagged. If A1 is broken, RF stays RF, two groups, zero dups.
#[test]
fn b01_rf_at_same_position_merges_with_fr() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    // Pair "fr": write rev first, fwd second → scanner sees rev first.
    let fr_rev = ReadSpec {
        qname: "fr",
        flags: 0x10, // reverse (add_pair will OR 0x41)
        ref_name: Some("chr1"),
        pos: Some(100),
        cigar: "1M",
        ..Default::default()
    };
    let fr_fwd = ReadSpec {
        qname: "fr",
        flags: 0,
        ref_name: Some("chr1"),
        pos: Some(100),
        cigar: "1M",
        ..Default::default()
    };
    // Pair "rf": write fwd first, rev second → scanner sees fwd first.
    let rf_fwd = ReadSpec {
        qname: "rf",
        flags: 0,
        ref_name: Some("chr1"),
        pos: Some(100),
        cigar: "1M",
        ..Default::default()
    };
    let rf_rev = ReadSpec {
        qname: "rf",
        flags: 0x10,
        ref_name: Some("chr1"),
        pos: Some(100),
        cigar: "1M",
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_pair(fr_rev, fr_fwd) // R1=rev, R2=fwd in add_pair order
        .add_pair(rf_fwd, rf_rev) // R1=fwd, R2=rev in add_pair order
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let dups = dup_qnames_set(&output);
    // Exactly one of the two pairs must be flagged (tie-break determines which).
    assert_eq!(
        dups.len(),
        1,
        "expected exactly 1 duplicate QNAME (RF/FR merged by A1), got {:?}",
        dups
    );
    assert!(
        dups.contains("fr") || dups.contains("rf"),
        "expected either fr or rf flagged, got {:?}",
        dups
    );
}

// =============================================================================
// B2 — A2: fragment-vs-pair priority via eager paired-marker
// =============================================================================
//
// Verifies: src/groups.rs:42-47 + src/scan.rs:198-214. When a paired read
// is observed, a zero-score "paired marker" is eagerly inserted into the
// SingleEndTracker at its locus. Any true fragment at that locus is then
// unconditionally flagged, regardless of its quality score.
//
// Setup: SE fragment "frag" at chr1:100 with Q60 (high score). Pair "pair"
// R1 at chr1:100 fwd Q20 (low score), R2 off at chr1:300 rev. The fragment's
// score exceeds the pair's — without A2 the fragment would win. With A2 the
// eager marker at chr1:100:fwd wins regardless.
#[test]
fn b02_fragment_and_pair_same_position_fragment_flagged() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    let frag = ReadSpec {
        qname: "frag",
        flags: 0, // single-end, forward
        ref_name: Some("chr1"),
        pos: Some(100),
        mapq: 60,
        cigar: "100M",
        qual: vec![60u8; 100],
        ..Default::default()
    };
    let pair_r1 = ReadSpec {
        qname: "pair",
        flags: 0, // add_pair will OR 0x41
        ref_name: Some("chr1"),
        pos: Some(100),
        mapq: 20,
        cigar: "100M",
        qual: vec![20u8; 100],
        ..Default::default()
    };
    let pair_r2 = ReadSpec {
        qname: "pair",
        flags: 0x10, // reverse
        ref_name: Some("chr1"),
        pos: Some(300),
        mapq: 20,
        cigar: "100M",
        qual: vec![20u8; 100],
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(frag)
        .add_pair(pair_r1, pair_r2)
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let dups = dup_qnames_set(&output);
    assert!(
        dups.contains("frag"),
        "fragment must be flagged by eager paired-marker (A2), got {:?}",
        dups
    );
    assert!(
        !dups.contains("pair"),
        "paired read itself must NOT be flagged (no other pair to compete with), got {:?}",
        dups
    );
}

// =============================================================================
// B3 — Deviation #8: orphan paired read vs co-located fragment
// =============================================================================
//
// NOT a Picard-parity test. Tests the documented divergence in
// docs/deviations.md §8: orphan paired reads (mate never arrives) leave an
// eager paired-marker that unconditionally beats any co-located fragment —
// score comparison never happens, unlike Picard.
//
// Verifies: src/scan.rs:271 (orphan passes through unflagged) +
//           src/scan.rs:198-214 (eager marker stays in single_tracker).
//
// Construction: SE fragment "frag" Q60 at chr1:100. Orphan "orph" added as a
// single read with paired-flag set (0x1 | 0x40) but NO R2 ever written —
// scanner inserts eager marker + pending_mate, drains pending_mate at EOF
// as orphan. The marker at chr1:100:fwd remains and flags "frag".
#[test]
fn b03_orphan_paired_read_treated_as_fragment() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    let frag = ReadSpec {
        qname: "frag",
        flags: 0,
        ref_name: Some("chr1"),
        pos: Some(100),
        mapq: 60,
        cigar: "100M",
        qual: vec![60u8; 100],
        ..Default::default()
    };
    let orphan = ReadSpec {
        qname: "orph",
        flags: 0x1 | 0x40, // paired + first-in-pair, but R2 never added
        ref_name: Some("chr1"),
        pos: Some(100),
        mapq: 20,
        cigar: "100M",
        mate_ref_name: Some("chr1"),
        mate_pos: Some(999),
        qual: vec![20u8; 100],
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(frag)
        .add_read(orphan) // single read; no add_pair → no R2 written
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let dups = dup_qnames_set(&output);
    assert!(
        dups.contains("frag"),
        "fragment loses to orphan's paired-marker (deviation #8), got {:?}",
        dups
    );
    assert!(
        !dups.contains("orph"),
        "orphan passes through unflagged (src/scan.rs:271), got {:?}",
        dups
    );
}

// =============================================================================
// B4 — Regression: chimeric pair flagging (cross-chromosome)
// =============================================================================
//
// Two identical chimeric pairs (R1 on chr1, R2 on chr2) must still group by
// PairedEndKey and produce one duplicate. Exercises that A1 same-locus
// normalization doesn't accidentally interfere with cross-chromosome pairs.
//
// Implementation note: R2 must be REVERSE so its unclipped_5' extends beyond
// alignment_start, otherwise pair_A would be resolved immediately upon its R2
// arrival (before pair_B's R2) and the two would never group. This is the
// same timing invariant normal non-chimeric paired dedup relies on.
#[test]
fn b04_chimeric_pair_both_flagged() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    // Manual 4-read construction so coordinate order is correct across chromosomes.
    let pa_r1 = ReadSpec {
        qname: "pA",
        flags: 0x1 | 0x40 | 0x20, // paired, first, mate-reverse
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr2"),
        mate_pos: Some(500),
        ..Default::default()
    };
    let pb_r1 = ReadSpec {
        qname: "pB",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr2"),
        mate_pos: Some(500),
        ..Default::default()
    };
    let pa_r2 = ReadSpec {
        qname: "pA",
        flags: 0x1 | 0x80 | 0x10, // paired, second, reverse
        ref_name: Some("chr2"),
        pos: Some(500),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        ..Default::default()
    };
    let pb_r2 = ReadSpec {
        qname: "pB",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr2"),
        pos: Some(500),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .reference("chr2", 100_000)
        .read_group("rg1", "lib1")
        .add_read(pa_r1)
        .add_read(pb_r1)
        .add_read(pa_r2)
        .add_read(pb_r2)
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let dups = dup_qnames_set(&output);
    assert_eq!(
        dups.len(),
        1,
        "expected exactly 1 duplicate chimeric pair, got {:?}",
        dups
    );
    assert!(dups.contains("pA") || dups.contains("pB"));
}

// =============================================================================
// B5 — Regression: splice CIGAR with N, reverse-strand → different 5' → not dup
// =============================================================================
//
// Two paired reads with R1 reverse on chr1:100 but different N-sizes. For
// reverse reads, unclipped_5' = alignment_start + ref_consumed - 1, which
// DOES include the N-skip in ref_consumed. Different intron sizes → different
// unclipped_5' → different PairedEndKey → NOT duplicates.
//
// Calculations:
//   "a" R1 rev pos=100 CIGAR=50M100N50M: ref_consumed=200, 5'=100+200-1=299
//   "b" R1 rev pos=100 CIGAR=50M200N50M: ref_consumed=300, 5'=100+300-1=399
//   R2 fwd pos=1000: 5'=1000 for both
// "a" key: (chr1, 299, rev) ↔ (chr1, 1000, fwd)
// "b" key: (chr1, 399, rev) ↔ (chr1, 1000, fwd) — different lo → different group.
#[test]
fn b05_spliced_reverse_different_introns_not_dup() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    let a_r1 = ReadSpec {
        qname: "a",
        flags: 0x1 | 0x40 | 0x10, // paired, first, reverse (mate is fwd)
        ref_name: Some("chr1"),
        pos: Some(100),
        cigar: "50M100N50M",
        mate_ref_name: Some("chr1"),
        mate_pos: Some(1000),
        ..Default::default()
    };
    let b_r1 = ReadSpec {
        qname: "b",
        flags: 0x1 | 0x40 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(100),
        cigar: "50M200N50M",
        mate_ref_name: Some("chr1"),
        mate_pos: Some(1000),
        ..Default::default()
    };
    let a_r2 = ReadSpec {
        qname: "a",
        flags: 0x1 | 0x80 | 0x20, // paired, second, mate-reverse (R1 rev)
        ref_name: Some("chr1"),
        pos: Some(1000),
        cigar: "100M",
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        ..Default::default()
    };
    let b_r2 = ReadSpec {
        qname: "b",
        flags: 0x1 | 0x80 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(1000),
        cigar: "100M",
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(a_r1)
        .add_read(b_r1)
        .add_read(a_r2)
        .add_read(b_r2)
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let dups = dup_qnames_set(&output);
    assert!(
        dups.is_empty(),
        "different intron sizes on reverse R1 must not group (unclipped_5' differs), got {:?}",
        dups
    );
}

// =============================================================================
// B6 — Regression: splice CIGAR with N, forward-strand → same 5' → are dup
// =============================================================================
//
// Two paired reads with R1 forward on chr1:100 with different N-sizes. For
// forward reads, unclipped_5' = alignment_start (plus soft-clip, none here).
// N-size doesn't affect forward unclipped_5' → same PairedEndKey → ARE dups.
#[test]
fn b06_spliced_forward_same_start_are_dup() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    let a_r1 = ReadSpec {
        qname: "a",
        flags: 0x1 | 0x40 | 0x20, // paired, first, fwd, mate-reverse
        ref_name: Some("chr1"),
        pos: Some(100),
        cigar: "50M100N50M",
        mate_ref_name: Some("chr1"),
        mate_pos: Some(2000),
        ..Default::default()
    };
    let b_r1 = ReadSpec {
        qname: "b",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        cigar: "50M200N50M",
        mate_ref_name: Some("chr1"),
        mate_pos: Some(2000),
        ..Default::default()
    };
    let a_r2 = ReadSpec {
        qname: "a",
        flags: 0x1 | 0x80 | 0x10, // paired, second, reverse
        ref_name: Some("chr1"),
        pos: Some(2000),
        cigar: "100M",
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        ..Default::default()
    };
    let b_r2 = ReadSpec {
        qname: "b",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(2000),
        cigar: "100M",
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(a_r1)
        .add_read(b_r1)
        .add_read(a_r2)
        .add_read(b_r2)
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let dups = dup_qnames_set(&output);
    assert_eq!(
        dups.len(),
        1,
        "forward-strand same alignment_start + differing N → same unclipped_5' → 1 dup, got {:?}",
        dups
    );
}

// =============================================================================
// B7 — Scoring: explicit Q0 bases must score 0 (not auto-filled to Q30)
// =============================================================================
//
// Two identical single-end reads at chr1:100. "high" has Q30 qualities, "low"
// has explicit Q0 qualities. Q30 must win; "low" is the duplicate. The trick
// is that common::ReadSpec auto-fills an empty `qual` vec with Q30 (see
// tests/common/mod.rs:252) — we must pass `vec![0u8;100]` explicitly to make
// the low-quality read actually score zero.
#[test]
fn b07_missing_qual_scores_as_zero() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    let high = ReadSpec {
        qname: "high",
        flags: 0,
        ref_name: Some("chr1"),
        pos: Some(100),
        mapq: 60,
        cigar: "100M",
        qual: vec![30u8; 100],
        ..Default::default()
    };
    let low = ReadSpec {
        qname: "low",
        flags: 0,
        ref_name: Some("chr1"),
        pos: Some(100),
        mapq: 60,
        cigar: "100M",
        qual: vec![0u8; 100],
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(high)
        .add_read(low)
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let dups = dup_qnames_set(&output);
    assert_eq!(dups.len(), 1);
    assert!(
        dups.contains("low"),
        "Q0 read must lose to Q30 read, got {:?}",
        dups
    );
}

// =============================================================================
// B8 — A3: two RGs with same LB group as one library
// =============================================================================
//
// Two @RG entries both tagged LB:libA. Two identical pairs, one per RG. The
// library-idx resolution (src/scan.rs:26-36) must map both RGs to the same
// library, so the pairs land in the same PairedEndKey group and one is
// flagged as a duplicate.
#[test]
fn b08_multiple_rg_same_lb_same_library() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("out.metrics.txt");

    let pa_r1 = ReadSpec {
        qname: "pA",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(500),
        rg: Some("rg1"),
        ..Default::default()
    };
    let pb_r1 = ReadSpec {
        qname: "pB",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(500),
        rg: Some("rg2"),
        ..Default::default()
    };
    let pa_r2 = ReadSpec {
        qname: "pA",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(500),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        rg: Some("rg1"),
        ..Default::default()
    };
    let pb_r2 = ReadSpec {
        qname: "pB",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(500),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        rg: Some("rg2"),
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "libA")
        .read_group("rg2", "libA")
        .add_read(pa_r1)
        .add_read(pb_r1)
        .add_read(pa_r2)
        .add_read(pb_r2)
        .write(&input)
        .unwrap();

    run_markdup_with_metrics(&input, &output, &metrics).unwrap();

    let dups = dup_qnames_set(&output);
    assert_eq!(
        dups.len(),
        1,
        "same LB across RGs must group as one library, got {:?}",
        dups
    );

    let recs = parse_metrics(&metrics).unwrap();
    assert_eq!(recs.len(), 1, "writer emits a single aggregated row");
    assert_eq!(recs[0].read_pair_duplicates, 1);
    assert_eq!(recs[0].read_pairs_examined, 2);
}

// =============================================================================
// B9 — A3 fallback: two RGs with NO LB tag → both collapse to "Unknown Library"
// =============================================================================
//
// Two @RG entries without LB. Verifies that scanner's library fallback
// (src/scan.rs::build_library_map) matches Picard's
// `LibraryIdGenerator.getLibraryName` exactly — when LB is absent the library
// name resolves to the literal "Unknown Library", NOT the @RG ID. Both RGs
// therefore share a single library bucket, and two identical pairs (one per
// RG) collapse into one duplicate group → exactly 1 pair flagged.
#[test]
fn b09_missing_lb_fallback_unknown_library() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("out.metrics.txt");

    let pa_r1 = ReadSpec {
        qname: "pA",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(500),
        rg: Some("rg1"),
        ..Default::default()
    };
    let pb_r1 = ReadSpec {
        qname: "pB",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(500),
        rg: Some("rg2"),
        ..Default::default()
    };
    let pa_r2 = ReadSpec {
        qname: "pA",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(500),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        rg: Some("rg1"),
        ..Default::default()
    };
    let pb_r2 = ReadSpec {
        qname: "pB",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(500),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        rg: Some("rg2"),
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group_no_lb("rg1")
        .read_group_no_lb("rg2")
        .add_read(pa_r1)
        .add_read(pb_r1)
        .add_read(pa_r2)
        .add_read(pb_r2)
        .write(&input)
        .unwrap();

    run_markdup_with_metrics(&input, &output, &metrics).unwrap();

    let dups = dup_qnames_set(&output);
    assert_eq!(
        dups.len(),
        1,
        "RGs without LB must collapse to single \"Unknown Library\" bucket → \
         exactly 1 pair flagged, got {:?}",
        dups
    );
    assert!(dups.contains("pA") || dups.contains("pB"));

    let recs = parse_metrics(&metrics).unwrap();
    assert_eq!(recs[0].read_pair_duplicates, 1);
    assert_eq!(recs[0].library, "Unknown Library");
    // Sanity: no pair lost.
    let counts = pair_count_by_qname(&output);
    assert_eq!(counts["pA"], 2);
    assert_eq!(counts["pB"], 2);
}

// =============================================================================
// B10 — Intentional deviation match: MAPQ=0 on one mate does not suppress dup
// =============================================================================
//
// Picard issue #1285: MQ=0 on a mate does not exclude the pair from
// duplicate detection. We reproduce this behavior by design (documented
// deviation). Two identical chimeric pairs, R1 at chr1:100 with MAPQ=0;
// expected: 1 duplicate. If MQ=0 filtering sneaks in, neither pair is
// flagged and this test catches it.
#[test]
fn b10_mq_zero_chimeric_pair() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    let make_r1 = |qname: &'static str| ReadSpec {
        qname,
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mapq: 0, // <-- Picard #1285: does not suppress dedup
        mate_ref_name: Some("chr2"),
        mate_pos: Some(500),
        ..Default::default()
    };
    let make_r2 = |qname: &'static str| ReadSpec {
        qname,
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr2"),
        pos: Some(500),
        mapq: 60,
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .reference("chr2", 100_000)
        .read_group("rg1", "lib1")
        .add_read(make_r1("pA"))
        .add_read(make_r1("pB"))
        .add_read(make_r2("pA"))
        .add_read(make_r2("pB"))
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let dups = dup_qnames_set(&output);
    assert_eq!(
        dups.len(),
        1,
        "MAPQ=0 must not suppress dedup (Picard issue #1285 match), got {:?}",
        dups
    );
}

// =============================================================================
// B11 — A4: cross-RG same-QNAME mate lookup isolation
// =============================================================================
//
// Two reads in different RGs with the IDENTICAL QNAME "shared", at DIFFERENT
// positions. The A4 fix (src/pending_mates.rs:23-36) keys pending-mate
// lookup by (RG, QNAME), so the two same-named pairs do not collide in the
// pending buffer and both pairs resolve correctly.
//
// Assertion: no duplicates (positions differ) AND all 4 records are in
// output (pair_count_by_qname["shared"] == 4). If A4 is broken, one pair's
// R1 would be matched against the other pair's R2 and produce garbage or
// lose records.
#[test]
fn b11_cross_rg_same_qname_not_matched() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("out.metrics.txt");

    // rg1 pair at chr1:100 / chr1:500
    let a_r1 = ReadSpec {
        qname: "shared",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(500),
        rg: Some("rg1"),
        ..Default::default()
    };
    let a_r2 = ReadSpec {
        qname: "shared",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(500),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        rg: Some("rg1"),
        ..Default::default()
    };
    // rg2 pair at chr1:200 / chr1:700 (different positions)
    let b_r1 = ReadSpec {
        qname: "shared",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(200),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(700),
        rg: Some("rg2"),
        ..Default::default()
    };
    let b_r2 = ReadSpec {
        qname: "shared",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(700),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(200),
        rg: Some("rg2"),
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .read_group("rg2", "lib2")
        .add_read(a_r1)
        .add_read(b_r1)
        .add_read(a_r2)
        .add_read(b_r2)
        .write(&input)
        .unwrap();

    run_markdup_with_metrics(&input, &output, &metrics).unwrap();

    let dups = dup_qnames_set(&output);
    assert!(
        dups.is_empty(),
        "distinct libraries + distinct positions → no dups, got {:?}",
        dups
    );
    let counts = pair_count_by_qname(&output);
    assert_eq!(
        counts.get("shared").copied().unwrap_or(0),
        4,
        "all 4 records must be preserved (no collision via QNAME hash)"
    );
    let recs = parse_metrics(&metrics).unwrap();
    assert_eq!(recs[0].read_pair_duplicates, 0);
    assert_eq!(recs[0].read_pairs_examined, 2);
}

// =============================================================================
// B11b — A4: cross-RG same-QNAME isolation still holds when LB is shared
// =============================================================================
//
// This is the bug-shaped case: rg1 and rg2 share the same LB, so dedup grouping
// is intentionally shared, but mate lookup must still stay isolated per RG.
//
// Construction:
// - Two pairs named "shared" live in different RGs but the SAME library.
// - Their R1s arrive first (chr1:100 and chr1:200).
// - The rg2 R2 arrives before the rg1 R2, so a buggy lookup keyed only by
//   library+QNAME will match it to the wrong pending mate.
// - A third control pair lives at the exact mixed key that the buggy cross-pair
//   would synthesize.
//
// Correct behavior: all three true pairs are distinct → zero duplicates.
// Buggy behavior: one of the cross-paired "shared" combinations collides with
// "control" and produces a false duplicate.
#[test]
fn b11b_cross_rg_same_qname_same_lb_not_cross_matched() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("out.metrics.txt");

    let a_r1 = ReadSpec {
        qname: "shared",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(900),
        rg: Some("rg1"),
        ..Default::default()
    };
    let a_r2 = ReadSpec {
        qname: "shared",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(900),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        rg: Some("rg1"),
        ..Default::default()
    };
    let b_r1 = ReadSpec {
        qname: "shared",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(200),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(800),
        rg: Some("rg2"),
        ..Default::default()
    };
    let b_r2 = ReadSpec {
        qname: "shared",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(800),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(200),
        rg: Some("rg2"),
        ..Default::default()
    };
    let c_r1 = ReadSpec {
        qname: "control",
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(800),
        rg: Some("rg1"),
        ..Default::default()
    };
    let c_r2 = ReadSpec {
        qname: "control",
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(800),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        rg: Some("rg1"),
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .read_group("rg2", "lib1")
        .add_read(a_r1)
        .add_read(c_r1)
        .add_read(b_r1)
        .add_read(b_r2)
        .add_read(c_r2)
        .add_read(a_r2)
        .write(&input)
        .unwrap();

    run_markdup_with_metrics(&input, &output, &metrics).unwrap();

    let dups = dup_qnames_set(&output);
    assert!(
        dups.is_empty(),
        "shared-QNAME mates in different RGs must not cross-match even when LB is shared, got {:?}",
        dups
    );
    let counts = pair_count_by_qname(&output);
    assert_eq!(counts.get("shared").copied().unwrap_or(0), 4);
    assert_eq!(counts.get("control").copied().unwrap_or(0), 2);

    let recs = parse_metrics(&metrics).unwrap();
    assert_eq!(recs[0].read_pair_duplicates, 0);
    assert_eq!(recs[0].read_pairs_examined, 3);
}

// =============================================================================
// B12 — A6: bisection does not panic on near-degenerate input
// =============================================================================
//
// Four identical pairs → 3 duplicates, 1 unique ("winner"). This is the
// smallest-unique-ratio input realizable with standard grouping and exercises
// the Lander-Waterman bisection's growing upper bracket (src/metrics.rs:53-117).
// The bisection must not panic and must produce a finite `Some(_)` estimate.
// True `None` (undefined) is tested by B13 where pair_duplicates == 0.
#[test]
fn b12_all_duplicates_metric_not_panicking() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("out.metrics.txt");

    let make_r1 = |qname: &'static str| ReadSpec {
        qname,
        flags: 0x1 | 0x40 | 0x20,
        ref_name: Some("chr1"),
        pos: Some(100),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(500),
        ..Default::default()
    };
    let make_r2 = |qname: &'static str| ReadSpec {
        qname,
        flags: 0x1 | 0x80 | 0x10,
        ref_name: Some("chr1"),
        pos: Some(500),
        mate_ref_name: Some("chr1"),
        mate_pos: Some(100),
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(make_r1("p1"))
        .add_read(make_r1("p2"))
        .add_read(make_r1("p3"))
        .add_read(make_r1("p4"))
        .add_read(make_r2("p1"))
        .add_read(make_r2("p2"))
        .add_read(make_r2("p3"))
        .add_read(make_r2("p4"))
        .write(&input)
        .unwrap();

    // Must not panic on bisection edge (unique=0).
    run_markdup_with_metrics(&input, &output, &metrics).unwrap();

    let dups = dup_qnames_set(&output);
    assert_eq!(
        dups.len(),
        3,
        "4 identical pairs → 3 duplicates, got {:?}",
        dups
    );

    let recs = parse_metrics(&metrics).unwrap();
    assert_eq!(recs[0].read_pair_duplicates, 3);
    assert_eq!(recs[0].read_pairs_examined, 4);
    assert!(
        recs[0].estimated_library_size.is_some(),
        "bisection must yield a finite estimate, not panic, got {:?}",
        recs[0].estimated_library_size
    );
}

// =============================================================================
// B13 — A5 + A7: single-end-only input produces valid metrics with empty estimate
// =============================================================================
//
// Three single-end reads at distinct positions (no duplicates among them).
// Zero paired reads → zero READ_PAIRS_EXAMINED → zero READ_PAIR_DUPLICATES.
// A7 must serialize ESTIMATED_LIBRARY_SIZE as an empty field (None in our
// parser), and UNPAIRED_READS_EXAMINED (A5) must reflect all three fragments.
#[test]
fn b13_zero_pairs_metric_produces_valid_file() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("out.metrics.txt");

    let mk = |qname: &'static str, pos: usize| ReadSpec {
        qname,
        flags: 0,
        ref_name: Some("chr1"),
        pos: Some(pos),
        cigar: "100M",
        ..Default::default()
    };

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(mk("s1", 100))
        .add_read(mk("s2", 200))
        .add_read(mk("s3", 300))
        .write(&input)
        .unwrap();

    run_markdup_with_metrics(&input, &output, &metrics).unwrap();

    let dups = dup_qnames_set(&output);
    assert!(
        dups.is_empty(),
        "distinct positions → no dups, got {:?}",
        dups
    );

    let recs = parse_metrics(&metrics).unwrap();
    assert_eq!(
        recs[0].unpaired_reads_examined, 3,
        "A5: orphan counter = 3 SE reads"
    );
    assert_eq!(recs[0].read_pairs_examined, 0);
    assert_eq!(recs[0].read_pair_duplicates, 0);
    assert_eq!(recs[0].percent_duplication, 0.0);
    assert_eq!(
        recs[0].estimated_library_size, None,
        "A7: undefined estimate serializes as empty field"
    );
}

// =============================================================================
// B14 — Output-side fidelity: pre-existing duplicate flags are cleared/recomputed
// =============================================================================
//
// One unique single-end read arrives already marked duplicate (0x400). Picard
// clears stale duplicate flags and recomputes from scratch; the output must
// therefore NOT retain 0x400.
#[test]
fn b14_stale_duplicate_flag_is_cleared() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(ReadSpec {
            qname: "stale",
            flags: 0x400,
            ref_name: Some("chr1"),
            pos: Some(100),
            cigar: "100M",
            ..Default::default()
        })
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let flags = read_flags_by_qname_all(&output);
    assert_eq!(flags["stale"].len(), 1);
    assert_eq!(
        flags["stale"][0] & 0x400,
        0,
        "stale input duplicate bit must be cleared when the record is not a duplicate"
    );
}

// =============================================================================
// B15 — Output-side fidelity: CLEAR_DT strips stale DT tags before writing
// =============================================================================
//
// Picard defaults CLEAR_DT=true. A pre-existing DT:Z:LB tag on input must not
// survive to output when we are not explicitly re-tagging duplicates.
#[test]
fn b15_stale_dt_tag_is_cleared() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(ReadSpec {
            qname: "dt",
            ref_name: Some("chr1"),
            pos: Some(100),
            cigar: "100M",
            dt: Some("LB"),
            ..Default::default()
        })
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let tags = read_string_tag_by_qname(&output, *b"DT");
    assert_eq!(
        tags["dt"],
        vec![None],
        "stale DT tag must be removed from output"
    );
}

// =============================================================================
// B16 — Failure mode fidelity: chromosome regressions must hard-fail
// =============================================================================
//
// coordinate sort means references are globally nondecreasing. A stream
// chr1 -> chr2 -> chr1 must error rather than quietly flushing the earlier chr1
// duplicate groups and continuing with a wrong answer.
#[test]
fn b16_chromosome_regression_hard_fails() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");

    BamBuilder::new()
        .reference("chr1", 100_000)
        .reference("chr2", 100_000)
        .read_group("rg1", "lib1")
        .add_read(ReadSpec {
            qname: "a",
            ref_name: Some("chr1"),
            pos: Some(100),
            cigar: "100M",
            ..Default::default()
        })
        .add_read(ReadSpec {
            qname: "b",
            ref_name: Some("chr2"),
            pos: Some(100),
            cigar: "100M",
            ..Default::default()
        })
        .add_read(ReadSpec {
            qname: "c",
            ref_name: Some("chr1"),
            pos: Some(200),
            cigar: "100M",
            ..Default::default()
        })
        .write(&input)
        .unwrap();

    let err = run_markdup(&input, &output).expect_err("cross-chromosome regression must fail");
    let msg = err.to_string();
    assert!(
        msg.contains("Sort order violation") || msg.contains("Not coordinate-sorted"),
        "expected coordinate-sort failure, got: {}",
        msg
    );
}
