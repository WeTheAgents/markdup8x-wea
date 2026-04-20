//! Two-pass orchestrator: scan (Pass 1) → dup_bits → write (Pass 2).

use crate::barcode_tags::{format_mi_value, BarcodeTags};
use crate::dupset::DupSet;
use crate::io;
use crate::metrics;
use crate::scan;
use anyhow::{bail, Context, Result};
use bstr::BString;
use log::info;
use noodles::bam;
use noodles::bgzf;
use noodles::sam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::map::header::Version;
use noodles::sam::header::record::value::map::{self, Program};
use noodles::sam::header::record::value::Map;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

fn normalized_path(path: &str) -> Result<PathBuf> {
    let path = Path::new(path);
    if path.is_absolute() {
        Ok(path.to_path_buf())
    } else {
        Ok(std::env::current_dir()?.join(path))
    }
}

pub fn run(
    input: &str,
    output: Option<&str>,
    metrics_path: Option<&str>,
    threads: u32,
    remove_duplicates: bool,
    assume_sort_order: Option<&str>,
    barcode_tags: BarcodeTags<'_>,
) -> Result<()> {
    let (actual_path, _tmp) = io::resolve_input(input)?;
    let n_threads = (threads as usize).max(1);

    if let Some(output_path) = output {
        if normalized_path(output_path)? == normalized_path(&actual_path)? {
            bail!("Output path must differ from input path to avoid truncating the source BAM");
        }
    }

    // === Pass 1 ===
    let (mut reader, header) = io::open_bam(&actual_path, n_threads)?;
    io::validate_sort_order(&header, assume_sort_order)?;

    info!("Pass 1: scanning for duplicates...");
    let scan_result = scan::scan_pass(&mut reader, &header, barcode_tags)?;
    info!(
        "Pass 1 done: {} records, {} dups",
        scan_result.total_records,
        scan_result.dup_bits.len()
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
    let (mut reader2, mut header2) = io::open_bam(&actual_path, n_threads)?;

    // Picard parity: bump @HD VN to 1.6, preserving SO:coordinate.
    {
        use sam::header::record::value::map::header::tag::SORT_ORDER;
        let mut hd = Map::<map::Header>::new(Version::new(1, 6));
        if let Some(existing) = header2.header() {
            if let Some(so) = existing.other_fields().get(&SORT_ORDER) {
                hd.other_fields_mut().insert(SORT_ORDER, so.clone());
            }
        }
        *header2.header_mut() = Some(hd);
    }

    // Picard parity: add @PG ID:MarkDuplicates, chained to the last existing PG.
    {
        use sam::header::record::value::map::program::tag::{
            COMMAND_LINE, NAME, PREVIOUS_PROGRAM_ID, VERSION,
        };
        let mut pg = Map::<Program>::default();
        pg.other_fields_mut()
            .insert(VERSION, BString::from(env!("CARGO_PKG_VERSION")));
        pg.other_fields_mut()
            .insert(NAME, BString::from("MarkDuplicates"));
        // PP: chain to the last existing @PG ID (leaf of the PG chain).
        if let Ok(leaves) = header2.programs().leaves() {
            if let Some((leaf_id, _)) = leaves.last() {
                let pp_bytes: &[u8] = leaf_id.as_ref();
                pg.other_fields_mut()
                    .insert(PREVIOUS_PROGRAM_ID, BString::from(pp_bytes));
            }
        }
        // CL: Picard writes its invocation; we write ours (honest provenance).
        let out_str = output.unwrap_or("-");
        let cl = format!(
            "markdup-wea {} INPUT={} OUTPUT={}",
            env!("CARGO_PKG_VERSION"),
            input,
            out_str
        );
        pg.other_fields_mut()
            .insert(COMMAND_LINE, BString::from(cl));
        header2
            .programs_mut()
            .add("MarkDuplicates", pg)
            .context("Failed to add @PG MarkDuplicates to header")?;
    }

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

    // MI emission is gated on BOTH --molecular-identifier-tag AND --barcode-tag
    // being set (Picard MarkDuplicates.java:404-407). The CLI layer already
    // rejects MI-without-barcode; this `&&` is defense-in-depth.
    let mi_tag: Option<[u8; 2]> = match (barcode_tags.mi, barcode_tags.barcode) {
        (Some(t), Some(_)) => Some(*t),
        _ => None,
    };

    while reader2.read_record(&mut record)? > 0 {
        let is_dup = scan_result.dup_bits.contains(record_id);

        // Convert to RecordBuf to modify flags
        let mut rec_buf = RecordBuf::try_from_alignment_record(&header2, &record)?;
        rec_buf.data_mut().remove(&Tag::new(b'D', b'T'));
        // Picard parity: ADD_PG_TAG_TO_READS=true (Picard default).
        rec_buf.data_mut().insert(
            Tag::new(b'P', b'G'),
            Value::String(BString::from("MarkDuplicates")),
        );

        if let Some(tag_bytes) = mi_tag {
            let flags_u16 = u16::from(rec_buf.flags());
            let is_reverse = flags_u16 & 0x10 != 0;
            let contig_name: Option<String> = rec_buf.reference_sequence_id().and_then(|idx| {
                header2
                    .reference_sequences()
                    .get_index(idx)
                    .map(|(name, _)| name.to_string())
            });
            // Picard's pair-representative position: reverse → this record's
            // alignment_start; forward → its mate's alignment_start. 1-based.
            let pos_1based: i64 = if is_reverse {
                rec_buf
                    .alignment_start()
                    .map(|p| usize::from(p) as i64)
                    .unwrap_or(0)
            } else {
                rec_buf
                    .mate_alignment_start()
                    .map(|p| usize::from(p) as i64)
                    .unwrap_or(0)
            };
            let mi_value = format_mi_value(contig_name.as_deref(), is_reverse, pos_1based);
            rec_buf.data_mut().insert(
                Tag::new(tag_bytes[0], tag_bytes[1]),
                Value::String(BString::from(mi_value)),
            );
        }

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
        if remove_duplicates {
            format!(", {} removed", dups_written)
        } else {
            String::new()
        }
    );

    Ok(())
}
