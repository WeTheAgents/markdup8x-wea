//! Smoke tests for Track A.3 MOLECULAR_IDENTIFIER_TAG emission.
//!
//! Validates that `--molecular-identifier-tag MI` (plus the required
//! `--barcode-tag RX`) writes a per-record MI aux tag whose value matches
//! Picard's `UmiUtil.setMolecularIdentifier` format: `{contig}:{pos}/`
//! where pos is 1-based and is the record's own alignment_start on reverse
//! reads, otherwise the mate's alignment_start. `assignedUmi` is always
//! empty in plain MarkDuplicates — hence the trailing slash with nothing
//! after it. See `docs/umi_semantics.md` Q6.

mod common;

use common::{read_string_tag_by_qname, run_markdup_with_mi_tag, BamBuilder, ReadSpec};
use tempfile::tempdir;

#[test]
fn mi_tag_emitted_with_picard_format() {
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
        .write(&input)
        .unwrap();

    run_markdup_with_mi_tag(&input, &output, b"RX", b"MI").unwrap();

    let mi_values = read_string_tag_by_qname(&output, *b"MI");
    let pair_a_mi = mi_values.get("pair_a").expect("pair_a must have MI");
    // Two mates for pair_a. Forward mate (no 0x10) uses mate_alignment_start
    // (=1200 1-based); reverse mate (0x10) uses its own alignment_start
    // (=1200 1-based). Both end up "chr1:1200/".
    for v in pair_a_mi {
        assert_eq!(
            v.as_deref(),
            Some("chr1:1200/"),
            "pair_a MI must be chr1:1200/, got {:?}",
            v
        );
    }
}

#[test]
fn mi_not_emitted_without_barcode_tag() {
    // CLI layer rejects --molecular-identifier-tag without --barcode-tag,
    // and the runtime gate is defense-in-depth: a BarcodeTags with mi=Some
    // and barcode=None must not emit MI. We exercise the library path
    // directly (bypassing CLI validation) to lock the runtime gate.
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
                ..Default::default()
            },
        )
        .write(&input)
        .unwrap();

    let mi = *b"MI";
    markdup_wea::markdup::run(
        input.to_str().unwrap(),
        Some(output.to_str().unwrap()),
        None,
        1,
        false,
        None,
        markdup_wea::barcode_tags::BarcodeTags {
            mi: Some(&mi),
            ..Default::default()
        },
    )
    .unwrap();

    let mi_values = read_string_tag_by_qname(&output, *b"MI");
    for (qname, vs) in &mi_values {
        for v in vs {
            assert!(
                v.is_none(),
                "{} unexpectedly has MI tag {:?} when barcode_tag was not set",
                qname,
                v
            );
        }
    }
}
