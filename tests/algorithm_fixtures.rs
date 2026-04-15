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

use common::{dup_qnames_set, run_markdup, BamBuilder, ReadSpec};
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
