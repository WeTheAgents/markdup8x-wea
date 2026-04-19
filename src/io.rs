//! BAM I/O helpers.

use anyhow::{bail, Context, Result};
use log::warn;
use noodles::bam;
use noodles::bgzf;
use noodles::sam;
use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use std::fs::File;
use std::io::{self, Write};
use std::num::NonZeroUsize;
use tempfile::NamedTempFile;

/// The concrete BAM reader type (over a multithreaded-bgzf-wrapped File).
/// `MultithreadedReader` with worker_count=1 behaves like the plain reader,
/// so this type alias works for both single- and multi-threaded I/O.
pub type BamReader = bam::io::Reader<bgzf::MultithreadedReader<File>>;
/// The concrete BAM writer type (multithreaded BGZF wrapping inner W).
pub type BamWriter<W> = bam::io::Writer<bgzf::MultithreadedWriter<W>>;

/// Format-agnostic alignment reader. Hot-path consumers depend on this enum
/// rather than a concrete reader type, so additional formats (CRAM, etc.) can
/// be added later without touching scan/markdup signatures.
pub enum AlignmentReader {
    Bam(BamReader),
}

impl AlignmentReader {
    #[inline]
    pub fn read_record(&mut self, record: &mut bam::Record) -> io::Result<usize> {
        match self {
            AlignmentReader::Bam(r) => r.read_record(record),
        }
    }
}

/// Open a BAM reader from a file path with `threads` BGZF decode workers
/// (clamped to ≥1). Worker threads run BGZF decompression in the background;
/// the scan loop itself remains sequential by algorithm design.
pub fn open_bam(path: &str, threads: usize) -> Result<(AlignmentReader, sam::Header)> {
    let file = File::open(path).with_context(|| format!("Failed to open BAM: {}", path))?;
    let workers = NonZeroUsize::new(threads.max(1)).unwrap();
    let bgzf_reader = bgzf::MultithreadedReader::with_worker_count(workers, file);
    let mut reader = bam::io::Reader::from(bgzf_reader);
    let header = reader.read_header().context("Failed to read BAM header")?;
    Ok((AlignmentReader::Bam(reader), header))
}

/// If input is "-", copy stdin to temp file. Returns (actual_path, Option<temp_file>).
pub fn resolve_input(path: &str) -> Result<(String, Option<NamedTempFile>)> {
    if path == "-" {
        warn!("Reading from stdin; buffering to temp file for two-pass processing");
        let mut tmp = NamedTempFile::new().context("Failed to create temp file")?;
        let mut stdin = io::stdin().lock();
        io::copy(&mut stdin, &mut tmp).context("Failed to buffer stdin")?;
        tmp.as_file_mut().flush()?;
        let p = tmp.path().to_str().unwrap().to_string();
        Ok((p, Some(tmp)))
    } else {
        Ok((path.to_string(), None))
    }
}

/// Validate sort order from the SAM header.
pub fn validate_sort_order(header: &sam::Header, assume_sort_order: Option<&str>) -> Result<()> {
    if let Some(assumed) = assume_sort_order {
        if assumed != "coordinate" {
            bail!("Only coordinate sort order is supported, got --assume-sort-order {}", assumed);
        }
        return Ok(());
    }

    if let Some(hdr) = header.header() {
        if let Some(so) = hdr.other_fields().get(&SORT_ORDER) {
            let so_str = std::str::from_utf8(so.as_ref()).unwrap_or("");
            match so_str {
                "coordinate" => return Ok(()),
                "queryname" => bail!("Input BAM is sorted by queryname. markdup-wea requires coordinate-sorted input."),
                "unsorted" => bail!("Input BAM is unsorted. Sort with: samtools sort input.bam -o sorted.bam"),
                other => bail!("Unknown sort order: {}. Use --assume-sort-order coordinate.", other),
            }
        }
    }
    bail!("No @HD SO field. Use --assume-sort-order coordinate to override.")
}

/// Runtime sort-order checker.
pub struct SortOrderEnforcer {
    seen_any: bool,
    seen_unmapped_tail: bool,
    prev_tid: Option<usize>,
    prev_pos: i64,
}

impl SortOrderEnforcer {
    pub fn new() -> Self {
        Self { seen_any: false, seen_unmapped_tail: false, prev_tid: None, prev_pos: -1 }
    }

    pub fn check(&mut self, tid: Option<usize>, pos: i64) -> Result<()> {
        match tid {
            Some(cur_tid) => {
                if self.seen_unmapped_tail {
                    bail!("Not coordinate-sorted: mapped record on reference {} follows unmapped tail", cur_tid);
                }

                if let Some(prev_tid) = self.prev_tid {
                    if cur_tid < prev_tid {
                        bail!(
                            "Not coordinate-sorted: reference {} follows later reference {}",
                            cur_tid,
                            prev_tid
                        );
                    }

                    if cur_tid == prev_tid && pos < self.prev_pos {
                        bail!(
                            "Not coordinate-sorted: pos {} follows {} on same ref",
                            pos,
                            self.prev_pos
                        );
                    }
                }

                self.prev_tid = Some(cur_tid);
                self.prev_pos = pos;
            }
            None => {
                if self.seen_any {
                    self.seen_unmapped_tail = true;
                }
                self.prev_tid = None;
                self.prev_pos = -1;
            }
        }
        self.seen_any = true;
        Ok(())
    }
}

impl Default for SortOrderEnforcer {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sort_enforcer_valid() {
        let mut e = SortOrderEnforcer::new();
        assert!(e.check(Some(0), 100).is_ok());
        assert!(e.check(Some(0), 200).is_ok());
    }

    #[test]
    fn sort_enforcer_chrom_boundary() {
        let mut e = SortOrderEnforcer::new();
        assert!(e.check(Some(0), 50000).is_ok());
        assert!(e.check(Some(1), 100).is_ok());
    }

    #[test]
    fn sort_enforcer_rejects_chromosome_regression() {
        let mut e = SortOrderEnforcer::new();
        assert!(e.check(Some(0), 100).is_ok());
        assert!(e.check(Some(1), 100).is_ok());
        assert!(e.check(Some(0), 200).is_err());
    }

    #[test]
    fn sort_enforcer_rejects_mapped_after_unmapped_tail() {
        let mut e = SortOrderEnforcer::new();
        assert!(e.check(Some(0), 100).is_ok());
        assert!(e.check(None, -1).is_ok());
        assert!(e.check(Some(1), 100).is_err());
    }

    #[test]
    fn sort_enforcer_invalid() {
        let mut e = SortOrderEnforcer::new();
        assert!(e.check(Some(0), 5000).is_ok());
        assert!(e.check(Some(0), 3000).is_err());
    }
}
