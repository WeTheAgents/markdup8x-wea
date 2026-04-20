mod common;

use common::{run_markdup_with_metrics, BamBuilder, ReadSpec};
use serde_json::Value;
use std::process::Command;
use tempfile::tempdir;

#[test]
fn validation_compare_passes_on_identical_outputs() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("out.metrics.txt");
    let report = dir.path().join("compare.json");

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(ReadSpec {
            qname: "read1",
            ref_name: Some("chr1"),
            pos: Some(100),
            cigar: "100M",
            ..Default::default()
        })
        .write(&input)
        .unwrap();

    run_markdup_with_metrics(&input, &output, &metrics).unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_validation_compare"))
        .args([
            "--sample-id",
            "smoke",
            "--expected-bam",
            output.to_str().unwrap(),
            "--actual-bam",
            output.to_str().unwrap(),
            "--expected-metrics",
            metrics.to_str().unwrap(),
            "--actual-metrics",
            metrics.to_str().unwrap(),
            "--output-json",
            report.to_str().unwrap(),
        ])
        .status()
        .unwrap();

    assert!(
        status.success(),
        "comparator should succeed on identical files"
    );

    let json: Value = serde_json::from_str(&std::fs::read_to_string(&report).unwrap()).unwrap();
    assert_eq!(json["passed"], true);
    assert_eq!(json["sample_id"], "smoke");
}

#[test]
fn validation_compare_reports_metrics_mismatch() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("out.metrics.txt");
    let mutated_metrics = dir.path().join("mutated.metrics.txt");
    let report = dir.path().join("compare.json");

    BamBuilder::new()
        .reference("chr1", 100_000)
        .read_group("rg1", "lib1")
        .add_read(ReadSpec {
            qname: "read1",
            ref_name: Some("chr1"),
            pos: Some(100),
            cigar: "100M",
            ..Default::default()
        })
        .write(&input)
        .unwrap();

    run_markdup_with_metrics(&input, &output, &metrics).unwrap();
    let original = std::fs::read_to_string(&metrics).unwrap();
    std::fs::write(&mutated_metrics, original.replace("lib1", "MutatedLibrary")).unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_validation_compare"))
        .args([
            "--sample-id",
            "mismatch",
            "--expected-bam",
            output.to_str().unwrap(),
            "--actual-bam",
            output.to_str().unwrap(),
            "--expected-metrics",
            metrics.to_str().unwrap(),
            "--actual-metrics",
            mutated_metrics.to_str().unwrap(),
            "--output-json",
            report.to_str().unwrap(),
        ])
        .status()
        .unwrap();

    assert!(
        !status.success(),
        "comparator should fail on mismatched metrics"
    );

    let json: Value = serde_json::from_str(&std::fs::read_to_string(&report).unwrap()).unwrap();
    assert_eq!(json["passed"], false);
    assert_eq!(json["mismatch"]["class"], "metrics_text");
}
