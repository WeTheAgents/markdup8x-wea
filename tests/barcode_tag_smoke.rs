//! Smoke tests for Track A.2 BARCODE_TAG support.
//!
//! Validates that `--barcode-tag RX` splits reads at the same coordinates
//! into distinct duplicate groups by RX value. When RX differs, pairs are
//! NOT duplicates. When RX matches, pairs ARE duplicates. This is the
//! defining Picard semantic (see docs/umi_semantics.md Q4).

mod common;

use common::{
    read_flags_by_qname_all, run_markdup, run_markdup_with_barcode_tag, BamBuilder, ReadSpec,
};
use tempfile::tempdir;

fn is_dup(flag: u16) -> bool {
    flag & 0x400 != 0
}

/// Two paired-end reads at the same coordinates with DIFFERENT RX values
/// must NOT be marked as duplicates when `--barcode-tag RX` is active.
#[test]
fn different_rx_splits_duplicate_group() {
    let tmp = tempdir().unwrap();
    let input = tmp.path().join("in.bam");
    let output = tmp.path().join("out.bam");

    BamBuilder::new()
        .coord_sort()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_pair(
            ReadSpec {
                qname: "pair_a",
                ref_name: Some("chr1"),
                pos: Some(1000),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1200),
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
            ReadSpec {
                qname: "pair_a",
                flags: 0x10, // R2 reverse
                ref_name: Some("chr1"),
                pos: Some(1200),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1000),
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
        )
        .add_pair(
            ReadSpec {
                qname: "pair_b",
                ref_name: Some("chr1"),
                pos: Some(1000),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1200),
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"CCCCCC")],
                ..Default::default()
            },
            ReadSpec {
                qname: "pair_b",
                flags: 0x10,
                ref_name: Some("chr1"),
                pos: Some(1200),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1000),
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"CCCCCC")],
                ..Default::default()
            },
        )
        .write(&input)
        .unwrap();

    run_markdup_with_barcode_tag(&input, &output, b"RX").unwrap();

    let flags = read_flags_by_qname_all(&output);
    // No reads should be flagged — different RX means different dedup groups.
    for (qname, fs) in &flags {
        for f in fs {
            assert!(!is_dup(*f), "{} flagged as dup with flag {:#x}", qname, f);
        }
    }
}

/// Two paired-end reads at the same coordinates with MATCHING RX values
/// MUST be marked as duplicates (one pair kept, the other flagged).
#[test]
fn matching_rx_marks_duplicate() {
    let tmp = tempdir().unwrap();
    let input = tmp.path().join("in.bam");
    let output = tmp.path().join("out.bam");

    BamBuilder::new()
        .coord_sort()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_pair(
            ReadSpec {
                qname: "pair_a",
                ref_name: Some("chr1"),
                pos: Some(1000),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1200),
                qual: vec![40; 100], // higher quality → keep this pair
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
            ReadSpec {
                qname: "pair_a",
                flags: 0x10,
                ref_name: Some("chr1"),
                pos: Some(1200),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1000),
                qual: vec![40; 100],
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
        )
        .add_pair(
            ReadSpec {
                qname: "pair_b",
                ref_name: Some("chr1"),
                pos: Some(1000),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1200),
                qual: vec![20; 100], // lower quality → flag this pair
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
            ReadSpec {
                qname: "pair_b",
                flags: 0x10,
                ref_name: Some("chr1"),
                pos: Some(1200),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1000),
                qual: vec![20; 100],
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
        )
        .write(&input)
        .unwrap();

    run_markdup_with_barcode_tag(&input, &output, b"RX").unwrap();

    let flags = read_flags_by_qname_all(&output);
    // pair_a (higher qual) — both reads kept.
    for f in &flags["pair_a"] {
        assert!(!is_dup(*f), "pair_a unexpectedly flagged: {:#x}", f);
    }
    // pair_b (lower qual) — both reads flagged as dups.
    for f in &flags["pair_b"] {
        assert!(is_dup(*f), "pair_b not flagged as dup: {:#x}", f);
    }
}

/// Default-off (no `--barcode-tag` passed) must ignore RX tags and mark
/// the two same-coord pairs as duplicates just like pre-Track-A behavior.
/// This is the A.1 byte-parity regression gate expressed as a test.
#[test]
fn default_off_ignores_rx_tag() {
    let tmp = tempdir().unwrap();
    let input = tmp.path().join("in.bam");
    let output = tmp.path().join("out.bam");

    BamBuilder::new()
        .coord_sort()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_pair(
            ReadSpec {
                qname: "pair_a",
                ref_name: Some("chr1"),
                pos: Some(1000),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1200),
                qual: vec![40; 100],
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
            ReadSpec {
                qname: "pair_a",
                flags: 0x10,
                ref_name: Some("chr1"),
                pos: Some(1200),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1000),
                qual: vec![40; 100],
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
        )
        .add_pair(
            ReadSpec {
                qname: "pair_b",
                ref_name: Some("chr1"),
                pos: Some(1000),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1200),
                qual: vec![20; 100],
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"CCCCCC")], // different RX, but ignored
                ..Default::default()
            },
            ReadSpec {
                qname: "pair_b",
                flags: 0x10,
                ref_name: Some("chr1"),
                pos: Some(1200),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1000),
                qual: vec![20; 100],
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"CCCCCC")],
                ..Default::default()
            },
        )
        .write(&input)
        .unwrap();

    run_markdup(&input, &output).unwrap();

    let flags = read_flags_by_qname_all(&output);
    // pair_b flagged as dup despite different RX — feature was off.
    for f in &flags["pair_b"] {
        assert!(is_dup(*f), "pair_b not flagged (default-off): {:#x}", f);
    }
}

/// With `--barcode-tag RX` on, a record missing RX hashes to Picard's
/// null-fallback (31) while a record with `RX=AAAAAA` hashes to
/// `31 + java_string_hashcode("AAAAAA")`. Those differ → the two pairs
/// end up in distinct duplicate groups and neither is flagged.
#[test]
fn missing_rx_still_distinct_from_present_rx() {
    let tmp = tempdir().unwrap();
    let input = tmp.path().join("in.bam");
    let output = tmp.path().join("out.bam");

    BamBuilder::new()
        .coord_sort()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_pair(
            ReadSpec {
                qname: "pair_a",
                ref_name: Some("chr1"),
                pos: Some(1000),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1200),
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
            ReadSpec {
                qname: "pair_a",
                flags: 0x10,
                ref_name: Some("chr1"),
                pos: Some(1200),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1000),
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"AAAAAA")],
                ..Default::default()
            },
        )
        .add_pair(
            ReadSpec {
                qname: "pair_b",
                ref_name: Some("chr1"),
                pos: Some(1000),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1200),
                rg: Some("rg1"),
                // no RX tag
                ..Default::default()
            },
            ReadSpec {
                qname: "pair_b",
                flags: 0x10,
                ref_name: Some("chr1"),
                pos: Some(1200),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1000),
                rg: Some("rg1"),
                ..Default::default()
            },
        )
        .write(&input)
        .unwrap();

    run_markdup_with_barcode_tag(&input, &output, b"RX").unwrap();

    let flags = read_flags_by_qname_all(&output);
    for (qname, fs) in &flags {
        for f in fs {
            assert!(!is_dup(*f), "{} flagged: {:#x}", qname, f);
        }
    }
}

/// Invalid UMI character → markdup-wea must error out (Picard throws).
#[test]
fn invalid_umi_character_errors() {
    let tmp = tempdir().unwrap();
    let input = tmp.path().join("in.bam");
    let output = tmp.path().join("out.bam");

    BamBuilder::new()
        .coord_sort()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_pair(
            ReadSpec {
                qname: "pair_a",
                ref_name: Some("chr1"),
                pos: Some(1000),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1200),
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"XYZBAD")], // X/Y/Z/B not in [ATCGNatcgn-]
                ..Default::default()
            },
            ReadSpec {
                qname: "pair_a",
                flags: 0x10,
                ref_name: Some("chr1"),
                pos: Some(1200),
                cigar: "100M",
                mate_ref_name: Some("chr1"),
                mate_pos: Some(1000),
                rg: Some("rg1"),
                aux_tags: vec![(*b"RX", b"XYZBAD")],
                ..Default::default()
            },
        )
        .write(&input)
        .unwrap();

    let err = run_markdup_with_barcode_tag(&input, &output, b"RX").unwrap_err();
    let msg = format!("{:#}", err);
    assert!(
        msg.contains("Invalid UMI") || msg.contains("illegal character"),
        "unexpected error message: {}",
        msg
    );
}
