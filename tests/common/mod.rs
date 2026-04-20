//! Shared test infrastructure for markdup-wea integration tests.
//!
//! Provides:
//! - `BamBuilder` — fluent API to construct synthetic coordinate-sorted BAMs
//! - `read_flags_by_qname` — assertion helper to check FLAG 0x400 on output
//!
//! Use pattern:
//!
//! ```ignore
//! mod common;
//! use common::{BamBuilder, read_flags_by_qname};
//!
//! #[test]
//! fn my_fixture() {
//!     let tmp = tempfile::tempdir().unwrap();
//!     let input = tmp.path().join("in.bam");
//!     let output = tmp.path().join("out.bam");
//!
//!     BamBuilder::new()
//!         .reference("chr1", 100_000)
//!         .read_group("rg1", "lib1")
//!         .add_read(ReadSpec { /* ... */ })
//!         .write(&input).unwrap();
//!
//!     markdup_wea::markdup::run(
//!         input.to_str().unwrap(), Some(output.to_str().unwrap()),
//!         None, 1, false, None,
//!     ).unwrap();
//!
//!     let flags = read_flags_by_qname(&output);
//!     assert!(flags["read1"] & 0x400 != 0);
//! }
//! ```

#![allow(dead_code)] // helpers used selectively across test files

use bstr::BString;
use noodles::bam;
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record::{Flags, MappingQuality};
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use noodles::sam::header::record::value::map::read_group::tag::LIBRARY;
use noodles::sam::header::record::value::map::{self, Map, ReadGroup, ReferenceSequence};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufWriter;
use std::num::NonZeroUsize;
use std::path::Path;

// ---------- Read specification ----------

/// Declarative description of one alignment record.
///
/// Coordinates are 1-based (SAM convention). `pos` / `mate_pos` of `None` means unmapped.
#[derive(Clone, Debug)]
pub struct ReadSpec {
    pub qname: &'static str,
    pub flags: u16,
    pub ref_name: Option<&'static str>,
    pub pos: Option<usize>,
    pub mapq: u8,
    pub cigar: &'static str,
    pub mate_ref_name: Option<&'static str>,
    pub mate_pos: Option<usize>,
    pub tlen: i32,
    /// Sequence as ASCII bytes (A/C/G/T/N). If empty, defaults to all-'A' of CIGAR query length.
    pub seq: &'static [u8],
    /// Per-base Phred qualities. If empty, defaults to Q30 for every base.
    pub qual: Vec<u8>,
    /// Optional RG tag value. If set, a matching @RG must have been added to the builder.
    pub rg: Option<&'static str>,
    /// Optional DT tag value for output-side fidelity tests.
    pub dt: Option<&'static str>,
    /// Arbitrary 2-byte SAM aux string tags (e.g. RX for barcode tests).
    pub aux_tags: Vec<([u8; 2], &'static [u8])>,
}

impl Default for ReadSpec {
    fn default() -> Self {
        Self {
            qname: "read",
            flags: 0,
            ref_name: None,
            pos: None,
            mapq: 60,
            cigar: "100M",
            mate_ref_name: None,
            mate_pos: None,
            tlen: 0,
            seq: b"",
            qual: Vec::new(),
            rg: None,
            dt: None,
            aux_tags: Vec::new(),
        }
    }
}

// ---------- BAM builder ----------

pub struct BamBuilder {
    refs: Vec<(String, usize)>,
    read_groups: Vec<(String, String)>, // (id, library_name)
    reads: Vec<ReadSpec>,
    override_sort_order: Option<String>,
    coord_sort: bool,
}

impl BamBuilder {
    pub fn new() -> Self {
        Self {
            refs: Vec::new(),
            read_groups: Vec::new(),
            reads: Vec::new(),
            override_sort_order: None,
            coord_sort: false,
        }
    }

    /// Opt-in: stably sort records by (ref_id, pos) before writing. Use when
    /// adding multiple pairs at interleaving coordinates where manual ordering
    /// would be tedious. Tests that deliberately emit an unsorted BAM (B16)
    /// must NOT call this.
    pub fn coord_sort(mut self) -> Self {
        self.coord_sort = true;
        self
    }

    /// Add a reference sequence (@SQ).
    pub fn reference(mut self, name: &str, length: usize) -> Self {
        self.refs.push((name.to_string(), length));
        self
    }

    /// Add a read group (@RG) with an associated library (LB).
    pub fn read_group(mut self, id: &str, library: &str) -> Self {
        self.read_groups.push((id.to_string(), library.to_string()));
        self
    }

    /// Add a read group (@RG) **without an LB tag** — used to test the
    /// Picard-exact LB fallback where the library name defaults to the literal
    /// string "Unknown Library" (NOT the @RG ID, despite older docs that
    /// claimed otherwise — verified against Picard 3.4.0 source). Encoded here
    /// as a sentinel empty LB string which `build_header` skips when emitting LB.
    pub fn read_group_no_lb(mut self, id: &str) -> Self {
        self.read_groups.push((id.to_string(), String::new()));
        self
    }

    /// Override the @HD SO field (default: "coordinate"). Used for negative tests.
    pub fn sort_order(mut self, so: &str) -> Self {
        self.override_sort_order = Some(so.to_string());
        self
    }

    /// Add one alignment record.
    pub fn add_read(mut self, spec: ReadSpec) -> Self {
        self.reads.push(spec);
        self
    }

    /// Add a pair of reads (R1 + R2). Fills cross-mate fields automatically.
    /// R1 gets FLAG 0x41 (paired + first-in-pair), R2 gets 0x81 (paired + second-in-pair),
    /// plus any extra flags passed. Caller must still set 0x10/0x20 for reverse orientation.
    pub fn add_pair(mut self, mut r1: ReadSpec, mut r2: ReadSpec) -> Self {
        r1.flags |= 0x1 | 0x40; // paired, first-in-pair
        r2.flags |= 0x1 | 0x80; // paired, second-in-pair
                                // Cross-link mate positions if not already set.
        if r1.mate_ref_name.is_none() {
            r1.mate_ref_name = r2.ref_name;
        }
        if r1.mate_pos.is_none() {
            r1.mate_pos = r2.pos;
        }
        if r2.mate_ref_name.is_none() {
            r2.mate_ref_name = r1.ref_name;
        }
        if r2.mate_pos.is_none() {
            r2.mate_pos = r1.pos;
        }
        // Mate-reverse flag derived from partner's reverse flag.
        if r2.flags & 0x10 != 0 {
            r1.flags |= 0x20;
        }
        if r1.flags & 0x10 != 0 {
            r2.flags |= 0x20;
        }
        self.reads.push(r1);
        self.reads.push(r2);
        self
    }

    /// Write the assembled BAM to `path`. Records are stably sorted by
    /// (ref_id, pos) so callers can add pairs in any order and still produce
    /// a coordinate-sorted BAM. Unmapped reads (no pos) sort to the end.
    pub fn write(mut self, path: &Path) -> anyhow::Result<()> {
        let header = self.build_header()?;
        let ref_name_to_id: HashMap<String, usize> = self
            .refs
            .iter()
            .enumerate()
            .map(|(i, (n, _))| (n.clone(), i))
            .collect();

        if self.coord_sort {
            let sort_key = |spec: &ReadSpec| -> (i64, i64) {
                let tid = spec
                    .ref_name
                    .and_then(|n| ref_name_to_id.get(n).map(|&i| i as i64))
                    .unwrap_or(i64::MAX);
                let pos = spec.pos.map(|p| p as i64).unwrap_or(i64::MAX);
                (tid, pos)
            };
            self.reads.sort_by_key(sort_key);
        }

        let file = File::create(path)?;
        // bam::io::Writer::new wraps in BGZF; Writer::from does NOT (important!).
        let mut writer = bam::io::Writer::new(BufWriter::new(file));
        writer.write_header(&header)?;

        for spec in &self.reads {
            let rec = build_record(spec, &ref_name_to_id)?;
            writer.write_alignment_record(&header, &rec)?;
        }
        writer.try_finish()?; // flush BGZF EOF block
        Ok(())
    }

    fn build_header(&self) -> anyhow::Result<sam::Header> {
        use sam::header::record::value::map::header::Version;

        let mut hdr_map = Map::<map::Header>::new(Version::new(1, 6));
        let so = self.override_sort_order.as_deref().unwrap_or("coordinate");
        hdr_map
            .other_fields_mut()
            .insert(SORT_ORDER, BString::from(so));

        let mut builder = sam::Header::builder().set_header(hdr_map);

        for (name, length) in &self.refs {
            let len = NonZeroUsize::new(*length)
                .ok_or_else(|| anyhow::anyhow!("reference length must be > 0"))?;
            builder =
                builder.add_reference_sequence(name.clone(), Map::<ReferenceSequence>::new(len));
        }

        for (id, lib) in &self.read_groups {
            let mut rg = Map::<ReadGroup>::default();
            // Empty library string = sentinel from `read_group_no_lb`: omit the LB tag.
            if !lib.is_empty() {
                rg.other_fields_mut()
                    .insert(LIBRARY, BString::from(lib.as_str()));
            }
            builder = builder.add_read_group(id.clone(), rg);
        }

        Ok(builder.build())
    }
}

// ---------- Record construction ----------

fn build_record(spec: &ReadSpec, ref_map: &HashMap<String, usize>) -> anyhow::Result<RecordBuf> {
    let mut b = RecordBuf::builder()
        .set_name(spec.qname)
        .set_flags(Flags::from(spec.flags))
        .set_template_length(spec.tlen);

    if let Some(ref_name) = spec.ref_name {
        let id = *ref_map
            .get(ref_name)
            .ok_or_else(|| anyhow::anyhow!("unknown ref: {}", ref_name))?;
        b = b.set_reference_sequence_id(id);
    }
    if let Some(pos) = spec.pos {
        b = b.set_alignment_start(Position::try_from(pos)?);
    }
    if let Some(mref) = spec.mate_ref_name {
        let id = *ref_map
            .get(mref)
            .ok_or_else(|| anyhow::anyhow!("unknown mate ref: {}", mref))?;
        b = b.set_mate_reference_sequence_id(id);
    }
    if let Some(mpos) = spec.mate_pos {
        b = b.set_mate_alignment_start(Position::try_from(mpos)?);
    }
    if let Some(mq) = MappingQuality::new(spec.mapq) {
        b = b.set_mapping_quality(mq);
    }

    let cigar = parse_cigar(spec.cigar)?;
    let query_len = cigar_query_len(&cigar);
    b = b.set_cigar(cigar);

    let seq: Vec<u8> = if spec.seq.is_empty() {
        vec![b'A'; query_len]
    } else {
        spec.seq.to_vec()
    };
    b = b.set_sequence(Sequence::from(seq));

    let qual: Vec<u8> = if spec.qual.is_empty() {
        vec![30; query_len]
    } else {
        spec.qual.clone()
    };
    b = b.set_quality_scores(QualityScores::from(qual));

    if let Some(rg) = spec.rg {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::Data;
        let data: Data = [(Tag::READ_GROUP, Value::String(BString::from(rg)))]
            .into_iter()
            .collect();
        b = b.set_data(data);
    }

    let mut record = b.build();
    if let Some(dt) = spec.dt {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        record
            .data_mut()
            .insert(Tag::new(b'D', b'T'), Value::String(BString::from(dt)));
    }

    for (tag, value) in &spec.aux_tags {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        record.data_mut().insert(
            Tag::new(tag[0], tag[1]),
            Value::String(BString::from(*value)),
        );
    }

    Ok(record)
}

fn parse_cigar(s: &str) -> anyhow::Result<Cigar> {
    let mut ops: Vec<Op> = Vec::new();
    let mut num = 0usize;
    for c in s.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as usize - '0' as usize);
        } else {
            let kind = match c {
                'M' => Kind::Match,
                'I' => Kind::Insertion,
                'D' => Kind::Deletion,
                'N' => Kind::Skip,
                'S' => Kind::SoftClip,
                'H' => Kind::HardClip,
                'P' => Kind::Pad,
                '=' => Kind::SequenceMatch,
                'X' => Kind::SequenceMismatch,
                other => anyhow::bail!("unknown CIGAR op: {}", other),
            };
            ops.push(Op::new(kind, num));
            num = 0;
        }
    }
    Ok(ops.into_iter().collect())
}

/// Number of query bases consumed by the CIGAR (for SEQ/QUAL length).
fn cigar_query_len(cigar: &Cigar) -> usize {
    cigar
        .as_ref()
        .iter()
        .filter(|op| {
            matches!(
                op.kind(),
                Kind::Match
                    | Kind::Insertion
                    | Kind::SoftClip
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
            )
        })
        .map(|op| op.len())
        .sum()
}

// ---------- Output inspection ----------

/// Read an output BAM and return a map from QNAME → FLAG for the FIRST occurrence
/// of each QNAME. For paired-end fixtures, prefer `read_flags_by_qname_all`.
pub fn read_flags_by_qname(path: &Path) -> HashMap<String, u16> {
    let mut map = HashMap::new();
    let file = File::open(path).expect("open output BAM");
    let mut reader = bam::io::Reader::new(file);
    let _header = reader.read_header().expect("read header");
    let mut rec = bam::Record::default();
    while reader.read_record(&mut rec).expect("read record") > 0 {
        let name = rec.name().map(|n| {
            let b: &[u8] = n.as_ref();
            String::from_utf8_lossy(b).to_string()
        });
        if let Some(n) = name {
            map.entry(n).or_insert_with(|| u16::from(rec.flags()));
        }
    }
    map
}

/// Read all occurrences of each QNAME (R1 + R2 + supplementary etc.).
pub fn read_flags_by_qname_all(path: &Path) -> HashMap<String, Vec<u16>> {
    let mut map: HashMap<String, Vec<u16>> = HashMap::new();
    let file = File::open(path).expect("open output BAM");
    let mut reader = bam::io::Reader::new(file);
    let _header = reader.read_header().expect("read header");
    let mut rec = bam::Record::default();
    while reader.read_record(&mut rec).expect("read record") > 0 {
        let name = rec.name().map(|n| {
            let b: &[u8] = n.as_ref();
            String::from_utf8_lossy(b).to_string()
        });
        if let Some(n) = name {
            map.entry(n).or_default().push(u16::from(rec.flags()));
        }
    }
    map
}

/// Read a string-valued SAM tag from every record, grouped by QNAME.
pub fn read_string_tag_by_qname(path: &Path, tag: [u8; 2]) -> HashMap<String, Vec<Option<String>>> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record::data::field::Value;

    let mut map: HashMap<String, Vec<Option<String>>> = HashMap::new();
    let file = File::open(path).expect("open output BAM");
    let mut reader = bam::io::Reader::new(file);
    let _header = reader.read_header().expect("read header");
    let mut rec = bam::Record::default();
    let tag = Tag::new(tag[0], tag[1]);

    while reader.read_record(&mut rec).expect("read record") > 0 {
        let name = rec.name().map(|n| {
            let b: &[u8] = n.as_ref();
            String::from_utf8_lossy(b).to_string()
        });
        if let Some(n) = name {
            let value = rec
                .data()
                .get(&tag)
                .transpose()
                .expect("decode tag")
                .and_then(|v| match v {
                    Value::String(s) => Some(String::from_utf8_lossy(s.as_ref()).to_string()),
                    _ => None,
                });
            map.entry(n).or_default().push(value);
        }
    }

    map
}

/// Convenience: run markdup-wea on `input`, write output to `output`.
pub fn run_markdup(input: &Path, output: &Path) -> anyhow::Result<()> {
    markdup_wea::markdup::run(
        input.to_str().unwrap(),
        Some(output.to_str().unwrap()),
        None,
        1,
        false,
        None,
        markdup_wea::barcode_tags::BarcodeTags::default(),
    )
}

/// Like `run_markdup` but also writes a Picard-format metrics file.
pub fn run_markdup_with_metrics(input: &Path, output: &Path, metrics: &Path) -> anyhow::Result<()> {
    markdup_wea::markdup::run(
        input.to_str().unwrap(),
        Some(output.to_str().unwrap()),
        Some(metrics.to_str().unwrap()),
        1,
        false,
        None,
        markdup_wea::barcode_tags::BarcodeTags::default(),
    )
}

/// Run markdup-wea with a `--barcode-tag` (BARCODE_TAG) set.
pub fn run_markdup_with_barcode_tag(
    input: &Path,
    output: &Path,
    barcode_tag: &[u8; 2],
) -> anyhow::Result<()> {
    markdup_wea::markdup::run(
        input.to_str().unwrap(),
        Some(output.to_str().unwrap()),
        None,
        1,
        false,
        None,
        markdup_wea::barcode_tags::BarcodeTags {
            barcode: Some(barcode_tag),
            ..Default::default()
        },
    )
}

/// Run markdup-wea with both `--barcode-tag` and `--molecular-identifier-tag`.
pub fn run_markdup_with_mi_tag(
    input: &Path,
    output: &Path,
    barcode_tag: &[u8; 2],
    mi_tag: &[u8; 2],
) -> anyhow::Result<()> {
    markdup_wea::markdup::run(
        input.to_str().unwrap(),
        Some(output.to_str().unwrap()),
        None,
        1,
        false,
        None,
        markdup_wea::barcode_tags::BarcodeTags {
            barcode: Some(barcode_tag),
            mi: Some(mi_tag),
            ..Default::default()
        },
    )
}

/// Picard-format DuplicationMetrics row. `estimated_library_size` is `None`
/// when the underlying field is empty (undefined estimate — A7 serialization).
#[derive(Debug, Default, Clone)]
pub struct MetricsRecord {
    pub library: String,
    pub unpaired_reads_examined: u64,
    pub read_pairs_examined: u64,
    pub secondary_or_supplementary_rds: u64,
    pub unmapped_reads: u64,
    pub unpaired_read_duplicates: u64,
    pub read_pair_duplicates: u64,
    pub read_pair_optical_duplicates: u64,
    pub percent_duplication: f64,
    pub estimated_library_size: Option<u64>,
}

/// Parse a Picard-format metrics file. Binds columns by header name (robust to
/// column-order changes across Picard versions). Returns one `MetricsRecord`
/// per data row — our writer currently emits a single "default" row, but Phase
/// C will read real Picard output which may have multiple libraries.
pub fn parse_metrics(path: &Path) -> anyhow::Result<Vec<MetricsRecord>> {
    use std::io::Read;
    let mut s = String::new();
    File::open(path)?.read_to_string(&mut s)?;

    let mut lines = s.lines();
    // Scan for "## METRICS CLASS" marker.
    loop {
        match lines.next() {
            Some(l) if l.starts_with("## METRICS CLASS") => break,
            Some(_) => continue,
            None => anyhow::bail!("no METRICS CLASS section found in {:?}", path),
        }
    }
    let header_line = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("missing column header line after METRICS CLASS"))?;
    let cols: Vec<&str> = header_line.split('\t').collect();
    let idx = |name: &str| -> anyhow::Result<usize> {
        cols.iter()
            .position(|c| *c == name)
            .ok_or_else(|| anyhow::anyhow!("column {} missing in metrics header", name))
    };
    let i_lib = idx("LIBRARY")?;
    let i_uex = idx("UNPAIRED_READS_EXAMINED")?;
    let i_rpx = idx("READ_PAIRS_EXAMINED")?;
    let i_sos = idx("SECONDARY_OR_SUPPLEMENTARY_RDS")?;
    let i_unm = idx("UNMAPPED_READS")?;
    let i_urd = idx("UNPAIRED_READ_DUPLICATES")?;
    let i_rpd = idx("READ_PAIR_DUPLICATES")?;
    let i_rpo = idx("READ_PAIR_OPTICAL_DUPLICATES")?;
    let i_pct = idx("PERCENT_DUPLICATION")?;
    let i_els = idx("ESTIMATED_LIBRARY_SIZE")?;

    let mut out = Vec::new();
    for line in lines {
        if line.is_empty() || line.starts_with("## ") {
            break;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let els = f.get(i_els).copied().unwrap_or("");
        out.push(MetricsRecord {
            library: f.get(i_lib).copied().unwrap_or("").to_string(),
            unpaired_reads_examined: f.get(i_uex).and_then(|v| v.parse().ok()).unwrap_or(0),
            read_pairs_examined: f.get(i_rpx).and_then(|v| v.parse().ok()).unwrap_or(0),
            secondary_or_supplementary_rds: f.get(i_sos).and_then(|v| v.parse().ok()).unwrap_or(0),
            unmapped_reads: f.get(i_unm).and_then(|v| v.parse().ok()).unwrap_or(0),
            unpaired_read_duplicates: f.get(i_urd).and_then(|v| v.parse().ok()).unwrap_or(0),
            read_pair_duplicates: f.get(i_rpd).and_then(|v| v.parse().ok()).unwrap_or(0),
            read_pair_optical_duplicates: f.get(i_rpo).and_then(|v| v.parse().ok()).unwrap_or(0),
            percent_duplication: f.get(i_pct).and_then(|v| v.parse().ok()).unwrap_or(0.0),
            estimated_library_size: if els.is_empty() {
                None
            } else {
                els.parse().ok()
            },
        });
    }
    Ok(out)
}

/// Count of records per QNAME in `path`. Used to assert that no pair was lost
/// to a cross-RG QNAME collision (see B11).
pub fn pair_count_by_qname(path: &Path) -> HashMap<String, usize> {
    let mut map: HashMap<String, usize> = HashMap::new();
    let file = File::open(path).expect("open output BAM");
    let mut reader = bam::io::Reader::new(file);
    let _header = reader.read_header().expect("read header");
    let mut rec = bam::Record::default();
    while reader.read_record(&mut rec).expect("read record") > 0 {
        if let Some(n) = rec.name() {
            let b: &[u8] = n.as_ref();
            *map.entry(String::from_utf8_lossy(b).to_string())
                .or_insert(0) += 1;
        }
    }
    map
}

/// Set of QNAMEs for which at least one record in `path` has FLAG & 0x400 != 0.
/// This matches the Picard-diff QNAME-set semantics (a pair is "marked" if either
/// half carries the duplicate flag). Phase C3 will consume the same helper to diff
/// our output against real Picard output.
pub fn dup_qnames_set(path: &Path) -> HashSet<String> {
    let mut set = HashSet::new();
    let file = File::open(path).expect("open output BAM");
    let mut reader = bam::io::Reader::new(file);
    let _header = reader.read_header().expect("read header");
    let mut rec = bam::Record::default();
    while reader.read_record(&mut rec).expect("read record") > 0 {
        if u16::from(rec.flags()) & 0x400 != 0 {
            if let Some(n) = rec.name() {
                let b: &[u8] = n.as_ref();
                set.insert(String::from_utf8_lossy(b).to_string());
            }
        }
    }
    set
}
