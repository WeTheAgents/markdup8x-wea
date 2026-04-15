//! Smoke test: validate the fixture builder + end-to-end markdup pipeline.
//!
//! Scenario: two identical single-end reads at the same unclipped 5' position.
//! Expected: the lower-quality read is flagged with FLAG 0x400, the higher-quality one is not.
//!
//! This is INFRASTRUCTURE validation, not a Picard-correctness claim. Real edge-case
//! fixtures (Phase 6) arrive after the edge-case research doc.

mod common;

use common::{read_flags_by_qname, run_markdup, BamBuilder, ReadSpec};

#[test]
fn two_identical_single_end_reads_one_flagged() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("in.bam");
    let output = tmp.path().join("out.bam");

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        // Two identical SE reads at pos 100, 100M CIGAR, forward strand.
        // read_hi has higher quality (Q40 vs Q20) → should keep, mark read_lo as dup.
        .add_read(ReadSpec {
            qname: "read_hi",
            flags: 0, // not paired, mapped, forward
            ref_name: Some("chr1"),
            pos: Some(100),
            cigar: "100M",
            qual: vec![40; 100],
            rg: Some("rg1"),
            ..Default::default()
        })
        .add_read(ReadSpec {
            qname: "read_lo",
            flags: 0,
            ref_name: Some("chr1"),
            pos: Some(100),
            cigar: "100M",
            qual: vec![20; 100],
            rg: Some("rg1"),
            ..Default::default()
        })
        .write(&input)
        .expect("write input BAM");

    run_markdup(&input, &output).expect("markdup run");

    let flags = read_flags_by_qname(&output);
    assert_eq!(flags.len(), 2, "expected both reads in output");

    let hi = flags["read_hi"];
    let lo = flags["read_lo"];

    assert_eq!(
        hi & 0x400,
        0,
        "read_hi (Q40) should NOT be flagged; FLAG=0x{:x}",
        hi
    );
    assert_eq!(
        lo & 0x400,
        0x400,
        "read_lo (Q20) SHOULD be flagged as duplicate; FLAG=0x{:x}",
        lo
    );
}

#[test]
fn single_read_never_flagged() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("in.bam");
    let output = tmp.path().join("out.bam");

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(ReadSpec {
            qname: "solo",
            flags: 0,
            ref_name: Some("chr1"),
            pos: Some(500),
            cigar: "100M",
            rg: Some("rg1"),
            ..Default::default()
        })
        .write(&input)
        .expect("write input BAM");

    run_markdup(&input, &output).expect("markdup run");

    let flags = read_flags_by_qname(&output);
    assert_eq!(flags["solo"] & 0x400, 0, "solo read must not be flagged");
}
