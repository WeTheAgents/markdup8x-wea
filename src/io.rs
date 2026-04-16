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

/// Open a BAM reader from a file path with `threads` BGZF decode workers
/// (clamped to ≥1). Worker threads run BGZF decompression in the background;
/// the scan loop itself remains sequential by algorithm design.
pub fn open_bam(path: &str, threads: usize) -> Result<(BamReader, sam::Header)> {
    let file = File::open(path).with_context(|| format!("Failed to open BAM: {}", path))?;
    let workers = NonZeroUsize::new(threads.max(1)).unwrap();
    let bgzf_reader = bgzf::MultithreadedReader::with_worker_count(workers, file);
    let mut reader = bam::io::Reader::from(bgzf_reader);
    let header = reader.read_header().context("Failed to read BAM header")?;
    Ok((reader, header))
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
    prev_tid: Option<usize>,
    prev_pos: i64,
}

impl SortOrderEnforcer {
    pub fn new() -> Self {
        Self { prev_tid: None, prev_pos: -1 }
    }

    pub fn check(&mut self, tid: Option<usize>, pos: i64) -> Result<()> {
        if let (Some(cur), Some(prev)) = (tid, self.prev_tid) {
            if cur == prev && pos < self.prev_pos {
                bail!("Not coordinate-sorted: pos {} follows {} on same ref", pos, self.prev_pos);
            }
        }
        if tid != self.prev_tid { self.prev_pos = -1; }
        self.prev_tid = tid;
        self.prev_pos = pos;
        Ok(())
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
    fn sort_enforcer_invalid() {
        let mut e = SortOrderEnforcer::new();
        assert!(e.check(Some(0), 5000).is_ok());
        assert!(e.check(Some(0), 3000).is_err());
    }
}
