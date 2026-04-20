use anyhow::{Context, Result};
use clap::Parser;
use noodles::bam;
use noodles::bgzf;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use serde::Serialize;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "validation_compare")]
#[command(about = "Structured Picard-vs-markdup-wea BAM and metrics comparator")]
struct Cli {
    #[arg(long)]
    sample_id: String,

    #[arg(long)]
    expected_bam: PathBuf,

    #[arg(long)]
    actual_bam: PathBuf,

    #[arg(long)]
    expected_metrics: Option<PathBuf>,

    #[arg(long)]
    actual_metrics: Option<PathBuf>,

    #[arg(long)]
    output_json: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
struct ComparisonResult {
    sample_id: String,
    passed: bool,
    comparator: String,
    expected_bam: String,
    actual_bam: String,
    expected_metrics: Option<String>,
    actual_metrics: Option<String>,
    records_compared: u64,
    mismatch: Option<Mismatch>,
}

#[derive(Debug, Serialize)]
struct Mismatch {
    class: String,
    record_index: Option<u64>,
    field: Option<String>,
    expected: Option<String>,
    actual: Option<String>,
    details: Option<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let result = compare(&cli)?;
    let json = serde_json::to_string_pretty(&result)?;

    if let Some(path) = &cli.output_json {
        std::fs::write(path, json.as_bytes())
            .with_context(|| format!("Failed to write comparator JSON to {}", path.display()))?;
    }

    println!("{json}");
    if result.passed {
        Ok(())
    } else {
        std::process::exit(1);
    }
}

fn compare(cli: &Cli) -> Result<ComparisonResult> {
    let comparator = format!("validation_compare {}", env!("CARGO_PKG_VERSION"));
    let mut result = ComparisonResult {
        sample_id: cli.sample_id.clone(),
        passed: true,
        comparator,
        expected_bam: cli.expected_bam.display().to_string(),
        actual_bam: cli.actual_bam.display().to_string(),
        expected_metrics: cli
            .expected_metrics
            .as_ref()
            .map(|path| path.display().to_string()),
        actual_metrics: cli
            .actual_metrics
            .as_ref()
            .map(|path| path.display().to_string()),
        records_compared: 0,
        mismatch: None,
    };

    let bam_mismatch = compare_bams(
        &cli.expected_bam,
        &cli.actual_bam,
        &mut result.records_compared,
    )?;
    if let Some(mismatch) = bam_mismatch {
        result.passed = false;
        result.mismatch = Some(mismatch);
        return Ok(result);
    }

    if let (Some(expected_metrics), Some(actual_metrics)) =
        (&cli.expected_metrics, &cli.actual_metrics)
    {
        if let Some(mismatch) = compare_metrics(expected_metrics, actual_metrics)? {
            result.passed = false;
            result.mismatch = Some(mismatch);
            return Ok(result);
        }
    }

    Ok(result)
}

fn compare_bams(
    expected_bam: &PathBuf,
    actual_bam: &PathBuf,
    records_compared: &mut u64,
) -> Result<Option<Mismatch>> {
    let (mut expected_reader, expected_header) = open_bam(expected_bam)?;
    let (mut actual_reader, actual_header) = open_bam(actual_bam)?;

    let expected_header_text = render_header(&expected_header)?;
    let actual_header_text = render_header(&actual_header)?;
    if expected_header_text != actual_header_text {
        return Ok(Some(Mismatch {
            class: "header_text".to_string(),
            record_index: None,
            field: None,
            expected: Some(expected_header_text),
            actual: Some(actual_header_text),
            details: Some("SAM header text differs".to_string()),
        }));
    }

    let mut expected_record = bam::Record::default();
    let mut actual_record = bam::Record::default();
    let mut index = 0u64;

    loop {
        let expected_len = expected_reader.read_record(&mut expected_record)?;
        let actual_len = actual_reader.read_record(&mut actual_record)?;

        if expected_len == 0 && actual_len == 0 {
            *records_compared = index;
            return Ok(None);
        }

        if expected_len == 0 || actual_len == 0 {
            return Ok(Some(Mismatch {
                class: "record_count".to_string(),
                record_index: Some(index),
                field: None,
                expected: Some(if expected_len == 0 {
                    "EOF".to_string()
                } else {
                    "record".to_string()
                }),
                actual: Some(if actual_len == 0 {
                    "EOF".to_string()
                } else {
                    "record".to_string()
                }),
                details: Some("One BAM ended before the other".to_string()),
            }));
        }

        let expected_buf = RecordBuf::try_from_alignment_record(&expected_header, &expected_record)
            .with_context(|| format!("Failed to materialize expected BAM record {index}"))?;
        let actual_buf = RecordBuf::try_from_alignment_record(&actual_header, &actual_record)
            .with_context(|| format!("Failed to materialize actual BAM record {index}"))?;

        if let Some(mismatch) = compare_record(index, &expected_buf, &actual_buf) {
            return Ok(Some(mismatch));
        }

        index += 1;
    }
}

fn compare_record(index: u64, expected: &RecordBuf, actual: &RecordBuf) -> Option<Mismatch> {
    compare_field(index, "qname", expected.name(), actual.name())
        .or_else(|| compare_field(index, "flags", expected.flags(), actual.flags()))
        .or_else(|| {
            compare_field(
                index,
                "reference_sequence_id",
                expected.reference_sequence_id(),
                actual.reference_sequence_id(),
            )
        })
        .or_else(|| {
            compare_field(
                index,
                "alignment_start",
                expected.alignment_start(),
                actual.alignment_start(),
            )
        })
        .or_else(|| {
            compare_field(
                index,
                "mapping_quality",
                expected.mapping_quality(),
                actual.mapping_quality(),
            )
        })
        .or_else(|| compare_field(index, "cigar", expected.cigar(), actual.cigar()))
        .or_else(|| {
            compare_field(
                index,
                "mate_reference_sequence_id",
                expected.mate_reference_sequence_id(),
                actual.mate_reference_sequence_id(),
            )
        })
        .or_else(|| {
            compare_field(
                index,
                "mate_alignment_start",
                expected.mate_alignment_start(),
                actual.mate_alignment_start(),
            )
        })
        .or_else(|| {
            compare_field(
                index,
                "template_length",
                expected.template_length(),
                actual.template_length(),
            )
        })
        .or_else(|| compare_field(index, "sequence", expected.sequence(), actual.sequence()))
        .or_else(|| {
            compare_field(
                index,
                "quality_scores",
                expected.quality_scores(),
                actual.quality_scores(),
            )
        })
        .or_else(|| compare_field(index, "aux", expected.data(), actual.data()))
}

fn compare_field<T>(index: u64, field: &str, expected: T, actual: T) -> Option<Mismatch>
where
    T: PartialEq + std::fmt::Debug,
{
    if expected == actual {
        None
    } else {
        Some(Mismatch {
            class: field.to_string(),
            record_index: Some(index),
            field: Some(field.to_string()),
            expected: Some(format!("{expected:?}")),
            actual: Some(format!("{actual:?}")),
            details: Some("Alignment record field differs".to_string()),
        })
    }
}

fn compare_metrics(expected_path: &PathBuf, actual_path: &PathBuf) -> Result<Option<Mismatch>> {
    let expected = read_text(expected_path)?;
    let actual = read_text(actual_path)?;

    if expected == actual {
        return Ok(None);
    }

    let expected_lines: Vec<&str> = expected.lines().collect();
    let actual_lines: Vec<&str> = actual.lines().collect();
    let max_len = expected_lines.len().max(actual_lines.len());
    for idx in 0..max_len {
        let expected_line = expected_lines.get(idx).copied();
        let actual_line = actual_lines.get(idx).copied();
        if expected_line != actual_line {
            return Ok(Some(Mismatch {
                class: "metrics_text".to_string(),
                record_index: None,
                field: Some(format!("line {}", idx + 1)),
                expected: expected_line.map(ToOwned::to_owned),
                actual: actual_line.map(ToOwned::to_owned),
                details: Some("Metrics text differs".to_string()),
            }));
        }
    }

    Ok(Some(Mismatch {
        class: "metrics_text".to_string(),
        record_index: None,
        field: None,
        expected: Some(expected),
        actual: Some(actual),
        details: Some("Metrics text differs".to_string()),
    }))
}

fn open_bam(path: &PathBuf) -> Result<(bam::io::Reader<bgzf::Reader<File>>, sam::Header)> {
    let file =
        File::open(path).with_context(|| format!("Failed to open BAM {}", path.display()))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .with_context(|| format!("Failed to read BAM header from {}", path.display()))?;
    Ok((reader, header))
}

fn render_header(header: &sam::Header) -> Result<String> {
    let mut writer = sam::io::Writer::new(Vec::new());
    writer.write_header(header)?;
    let bytes = writer.get_ref();
    String::from_utf8(bytes.clone()).context("Header serialization was not UTF-8")
}

fn read_text(path: &PathBuf) -> Result<String> {
    let mut buf = String::new();
    File::open(path)
        .with_context(|| format!("Failed to open text file {}", path.display()))?
        .read_to_string(&mut buf)
        .with_context(|| format!("Failed to read text file {}", path.display()))?;
    Ok(buf)
}
