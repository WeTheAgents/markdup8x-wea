//! Two-pass orchestrator: scan (Pass 1) → dup_bits → write (Pass 2).

use crate::dupset::DupSet;
use crate::io;
use crate::metrics;
use crate::scan;
use anyhow::{Context, Result};
use log::info;
use noodles::bam;
use noodles::bgzf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::RecordBuf;
use std::num::NonZeroUsize;

pub fn run(
    input: &str,
    output: Option<&str>,
    metrics_path: Option<&str>,
    threads: u32,
    remove_duplicates: bool,
    assume_sort_order: Option<&str>,
) -> Result<()> {
    let (actual_path, _tmp) = io::resolve_input(input)?;
    let n_threads = (threads as usize).max(1);

    // === Pass 1 ===
    let (mut reader, header) = io::open_bam(&actual_path, n_threads)?;
    io::validate_sort_order(&header, assume_sort_order)?;

    info!("Pass 1: scanning for duplicates...");
    let scan_result = scan::scan_pass(&mut reader, &header)?;
    info!(
        "Pass 1 done: {} records, {} dups",
        scan_result.total_records, scan_result.dup_bits.len()
    );
    drop(reader);

    // === Metrics ===
    if let Some(mp) = metrics_path {
        let out_str = output.unwrap_or("stdout");
        metrics::write_metrics(mp, &scan_result.counters, input, out_str)?;
        info!("Metrics written to {}", mp);
    }

    // === Pass 2: read records, convert to RecordBuf, modify flags, write ===
    info!("Pass 2: writing output...");
    let (mut reader2, header2) = io::open_bam(&actual_path, n_threads)?;

    // Pass-2 writer wraps the destination in MultithreadedWriter for parallel
    // BGZF compression. We then construct bam::io::Writer::from(bgzf_writer)
    // — `from` accepts a pre-wrapped inner writer (vs `new` which would wrap
    // a second time in single-threaded BGZF).
    let workers = NonZeroUsize::new(n_threads).unwrap();
    type BoxedRaw = Box<dyn std::io::Write + Send + 'static>;
    let raw: BoxedRaw = if let Some(p) = output {
        let f = std::fs::File::create(p).with_context(|| format!("Failed to create: {}", p))?;
        Box::new(std::io::BufWriter::new(f))
    } else {
        Box::new(std::io::BufWriter::new(std::io::stdout()))
    };
    let bgzf_w = bgzf::MultithreadedWriter::with_worker_count(workers, raw);
    let mut writer: bam::io::Writer<bgzf::MultithreadedWriter<BoxedRaw>> =
        bam::io::Writer::from(bgzf_w);

    writer.write_header(&header2)?;

    let mut record = bam::Record::default();
    let mut record_id: u64 = 0;
    let mut dups_written: u64 = 0;
    let mut records_written: u64 = 0;

    while reader2.read_record(&mut record)? > 0 {
        let is_dup = scan_result.dup_bits.contains(record_id);

        // Convert to RecordBuf to modify flags
        let mut rec_buf = RecordBuf::try_from_alignment_record(&header2, &record)?;

        let current_flags = u16::from(rec_buf.flags());
        let new_flags = if is_dup {
            dups_written += 1;
            current_flags | 0x400
        } else {
            current_flags & !0x400
        };
        *rec_buf.flags_mut() = noodles::sam::alignment::record::Flags::from(new_flags);

        if !(remove_duplicates && is_dup) {
            writer.write_alignment_record(&header2, &rec_buf)?;
            records_written += 1;
        }

        record_id += 1;
    }

    // MultithreadedWriter::finish flushes pending blocks, joins worker threads,
    // and writes the BGZF EOF marker. Single-threaded bgzf::Writer used `try_finish`;
    // the multithreaded variant exposes `finish` returning the inner W.
    writer.get_mut().finish()?;

    info!(
        "Pass 2 done: {} written, {} flagged{}",
        records_written,
        dups_written,
        if remove_duplicates { format!(", {} removed", dups_written) } else { String::new() }
    );

    Ok(())
}
