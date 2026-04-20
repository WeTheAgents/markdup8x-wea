//! Pass 1: Scan BAM, build duplicate groups, resolve, return dupset.

use crate::barcode_tags::BarcodeTags;
use crate::barcodes;
use crate::dupset::BitVecDupSet;
use crate::groups::{
    PairedEndKey, PairedGroupTracker, ScoredPair, ScoredSingle, SingleEndKey, SingleEndTracker,
};
use crate::io::SortOrderEnforcer;
use crate::metrics::MetricsCounters;
use crate::pending_mates::{check_hash, qname_hash, PendingMate, PendingMateBuffer};
use crate::position::unclipped_5prime;
use crate::scoring::quality_sum;
use anyhow::{bail, Context, Result};
use log::{info, warn};
use noodles::bam;
use noodles::sam;
use noodles::sam::header::record::value::map::read_group::tag::LIBRARY;
use rustc_hash::FxHashMap;

/// Build library mapping from @RG headers.
/// Returns (rg_to_lib, lib_names) — the second is used to emit a Picard-matching
/// LIBRARY name in the metrics file (see `metrics_library_name`).
fn build_library_map(header: &sam::Header) -> (FxHashMap<Vec<u8>, u8>, Vec<String>) {
    let mut lib_names: Vec<String> = Vec::new();
    let mut rg_to_lib: FxHashMap<Vec<u8>, u8> = FxHashMap::default();

    for (id, rg) in header.read_groups().iter() {
        // Picard library fallback (verified in real Picard 3.4.0 source —
        // `picard.sam.markduplicates.util.LibraryIdGenerator.getLibraryName`):
        // when @RG LB tag is absent, the library name falls back to the literal
        // string "Unknown Library", NOT the @RG ID. All reads with missing LB
        // therefore collapse into a single library bucket. Confirmed by real-data
        // diff against Picard on 8 ENCODE samples (Phase C).
        let _ = id; // @RG ID is unused for library naming; kept for diagnostics if needed.
        let lib_name = rg
            .other_fields()
            .get(&LIBRARY)
            .map(|v| String::from_utf8_lossy(v.as_ref()).to_string())
            .unwrap_or_else(|| "Unknown Library".to_string());

        let lib_idx = if let Some(idx) = lib_names.iter().position(|n| n == &lib_name) {
            idx as u8
        } else {
            let idx = lib_names.len() as u8;
            lib_names.push(lib_name);
            idx
        };
        let id_bytes: &[u8] = id.as_ref();
        rg_to_lib.insert(id_bytes.to_vec(), lib_idx);
    }

    (rg_to_lib, lib_names)
}

/// Choose the LIBRARY value to emit in Picard metrics. Picard's
/// `LibraryIdGenerator.getLibraryName` returns the @RG LB tag verbatim or
/// "Unknown Library" when LB is absent — it does NOT fall back to the @RG ID
/// (that's only our internal grouping convention for duplicate detection).
fn metrics_library_name(header: &sam::Header) -> String {
    for (_id, rg) in header.read_groups().iter() {
        if let Some(lb) = rg.other_fields().get(&LIBRARY) {
            return String::from_utf8_lossy(lb.as_ref()).to_string();
        }
    }
    "Unknown Library".to_string()
}

fn get_library_idx(read_group: Option<&[u8]>, rg_to_lib: &FxHashMap<Vec<u8>, u8>) -> u8 {
    if rg_to_lib.is_empty() {
        return 0;
    }
    if let Some(read_group) = read_group {
        if let Some(&idx) = rg_to_lib.get(read_group) {
            return idx;
        }
    }
    0
}

fn extract_read_group(record: &bam::Record) -> Result<Option<Vec<u8>>> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record::data::field::Value;

    let data = record.data();
    let value = data
        .get(&Tag::READ_GROUP)
        .transpose()
        .context("Failed to decode RG tag")?;

    match value {
        Some(Value::String(s)) => {
            let bytes: &[u8] = s.as_ref();
            Ok(Some(bytes.to_vec()))
        }
        Some(_) => bail!("RG tag must be a string"),
        None => Ok(None),
    }
}

/// Look up an arbitrary 2-byte SAM aux tag and return its value as an owned
/// byte vector. Returns `None` if the tag is absent; errors if present but
/// not a SAM String (Z/H). Used for BARCODE_TAG / READ_ONE_BARCODE_TAG /
/// READ_TWO_BARCODE_TAG lookup. Mirrors `extract_read_group`'s copy-out
/// pattern to sidestep noodles' record-scoped borrow.
fn lookup_string_tag(record: &bam::Record, tag: &[u8; 2]) -> Result<Option<Vec<u8>>> {
    use noodles::sam::alignment::record::data::field::{Tag, Value};

    let data = record.data();
    let value = data
        .get(&Tag::new(tag[0], tag[1]))
        .transpose()
        .with_context(|| {
            format!(
                "Failed to decode {} tag",
                std::str::from_utf8(tag).unwrap_or("??")
            )
        })?;

    match value {
        None => Ok(None),
        Some(Value::String(s)) => {
            let bytes: &[u8] = s.as_ref();
            Ok(Some(bytes.to_vec()))
        }
        Some(_) => bail!(
            "{} tag must be a string",
            std::str::from_utf8(tag).unwrap_or("??")
        ),
    }
}

/// Compute (barcode_hash, read_barcode_hash) for a single record under the
/// given `BarcodeTags` config.
///
/// **Feature off** (`tags.*` = None) → returns 0. This is how A.1's default-off
/// byte-parity contract is preserved: never invoke `picard_barcode_hash`
/// unless the CLI flag is explicitly set (31 ≠ 0 and would diverge).
///
/// **Feature on, tag missing on this record** → returns Picard's asymmetric
/// fallback: 31 for BARCODE_TAG (`Objects.hash(null)` = 31), 0 for READ_ONE/TWO
/// (bare `String.hashCode` on null → 0). See docs/umi_semantics.md Q1.
///
/// **read-one/read-two routing**: `is_first_of_pair` selects which of
/// `tags.read_one` / `tags.read_two` applies to this record. Picard keys that
/// routing off the BAM firstOfPair flag, NOT lo/hi coord ordering — see
/// docs/umi_semantics.md Q5.
fn extract_barcode_hashes(
    record: &bam::Record,
    is_first_of_pair: bool,
    tags: BarcodeTags<'_>,
) -> Result<(i32, i32)> {
    let barcode = if let Some(t) = tags.barcode {
        let umi = lookup_string_tag(record, t)?;
        if let Some(u) = umi.as_deref() {
            barcodes::validate_umi(u).with_context(|| {
                format!(
                    "Invalid UMI in {} tag",
                    std::str::from_utf8(t).unwrap_or("??")
                )
            })?;
        }
        barcodes::picard_barcode_hash(umi.as_deref())
    } else {
        0
    };

    let which = if is_first_of_pair {
        tags.read_one
    } else {
        tags.read_two
    };
    let read_barcode = if let Some(t) = which {
        let v = lookup_string_tag(record, t)?;
        barcodes::read_barcode_value(v.as_deref())
    } else {
        0
    };

    Ok((barcode, read_barcode))
}

pub struct ScanResult {
    pub dup_bits: BitVecDupSet,
    pub counters: MetricsCounters,
    pub total_records: u64,
    pub pending_peak: usize,
}

/// Extract CIGAR ops from BAM record.
fn extract_cigar_ops(
    record: &bam::Record,
) -> std::io::Result<Vec<(noodles::sam::alignment::record::cigar::op::Kind, usize)>> {
    record
        .cigar()
        .iter()
        .map(|r| r.map(|op| (op.kind(), op.len())))
        .collect()
}

/// Extract quality scores from BAM record as owned Vec.
fn extract_quals(record: &bam::Record) -> Vec<u8> {
    let qs = record.quality_scores();
    let bytes: &[u8] = qs.as_ref();
    bytes.to_vec()
}

fn extract_reference_sequence_id(record: &bam::Record) -> Result<Option<usize>> {
    record
        .reference_sequence_id()
        .transpose()
        .context("Failed to decode reference sequence ID")
}

fn extract_alignment_start(record: &bam::Record) -> Result<Option<i64>> {
    record
        .alignment_start()
        .transpose()
        .context("Failed to decode alignment start")
        .map(|pos| pos.map(|p| (usize::from(p) as i64) - 1))
}

fn extract_mate_reference_sequence_id(record: &bam::Record) -> Result<Option<i32>> {
    record
        .mate_reference_sequence_id()
        .transpose()
        .context("Failed to decode mate reference sequence ID")
        .map(|tid| tid.map(|t| t as i32))
}

fn extract_mate_alignment_start(record: &bam::Record) -> Result<Option<i64>> {
    record
        .mate_alignment_start()
        .transpose()
        .context("Failed to decode mate alignment start")
        .map(|pos| pos.map(|p| (usize::from(p) as i64) - 1))
}

/// Run Pass 1.
pub fn scan_pass(
    reader: &mut crate::io::AlignmentReader,
    header: &sam::Header,
    barcode_tags: BarcodeTags<'_>,
) -> Result<ScanResult> {
    let (rg_to_lib, _lib_names) = build_library_map(header);

    let mut dup_bits = BitVecDupSet::new();
    let mut pending = PendingMateBuffer::new();
    let mut paired_tracker = PairedGroupTracker::new();
    let mut single_tracker = SingleEndTracker::new();
    let mut enforcer = SortOrderEnforcer::new();
    let mut counters = MetricsCounters::new();
    counters.library_name = metrics_library_name(header);

    let mut record_id: u64 = 0;
    let mut prev_tid: Option<usize> = None;
    let mut record = bam::Record::default();

    while reader.read_record(&mut record)? > 0 {
        let flags = record.flags();
        let tid = extract_reference_sequence_id(&record)?;
        let pos = extract_alignment_start(&record)?.unwrap_or(-1);

        enforcer.check(tid, pos).context("Sort order violation")?;

        if tid != prev_tid && prev_tid.is_some() {
            single_tracker.flush(&mut dup_bits);
            if let Some(pt) = prev_tid {
                paired_tracker.resolve_chromosome(pt as i32, &mut dup_bits);
            }
        }
        prev_tid = tid;

        let current_id = record_id;
        record_id += 1;

        // Classification (priority order per spec)
        if flags.is_unmapped() {
            counters.unmapped_reads += 1;
            continue;
        }
        if flags.is_secondary() {
            counters.secondary_or_supplementary += 1;
            continue;
        }
        if flags.is_supplementary() {
            counters.secondary_or_supplementary += 1;
            continue;
        }

        let is_paired = flags.is_segmented();
        let mate_unmapped = flags.is_mate_unmapped();
        let is_reverse = flags.is_reverse_complemented();
        let is_first_of_pair = flags.is_first_segment();

        let cigar_ops = extract_cigar_ops(&record)?;
        let uc5 = unclipped_5prime(pos, &cigar_ops, is_reverse);
        let quals = extract_quals(&record);
        let qsum = quality_sum(&quals);
        let read_group = extract_read_group(&record)?;
        let lib_idx = get_library_idx(read_group.as_deref(), &rg_to_lib);
        let tid_i32 = tid.map(|t| t as i32).unwrap_or(-1);

        let (barcode_hash, read_barcode_hash) =
            extract_barcode_hashes(&record, is_first_of_pair, barcode_tags)?;

        if !is_paired || mate_unmapped {
            single_tracker.add_read(
                SingleEndKey {
                    library_idx: lib_idx,
                    barcode_hash,
                    ref_id: tid_i32,
                    unclipped_5prime: uc5,
                    is_reverse,
                },
                ScoredSingle {
                    score: qsum,
                    record_id: current_id,
                    is_paired_marker: false,
                },
                &mut dup_bits,
            );
        } else {
            let qname_bytes: &[u8] = record
                .name()
                .map(|n| {
                    let b: &[u8] = n.as_ref();
                    b
                })
                .unwrap_or(b"");
            let nh = qname_hash(qname_bytes, read_group.as_deref());
            let ch = check_hash(qname_bytes, read_group.as_deref());
            let mate_tid = extract_mate_reference_sequence_id(&record)?.unwrap_or(-1);
            let mate_pos = extract_mate_alignment_start(&record)?.unwrap_or(-1);

            // Insert a fragment-group marker at this read's single-end key.
            // Picard (§4 of research): paired reads always beat fragments at
            // the same locus, so any fragment landing in this group must be
            // flagged unconditionally. Insert inline at observation time —
            // that way retroactive insertion is never needed, because
            // coordinate-sort guarantees the SingleEndTracker is at this key
            // right now. Marker score=0, record_id=current_id (unused unless
            // something dereferences the Vec entry; the resolve code skips
            // markers before flagging).
            single_tracker.add_read(
                SingleEndKey {
                    library_idx: lib_idx,
                    barcode_hash,
                    ref_id: tid_i32,
                    unclipped_5prime: uc5,
                    is_reverse,
                },
                ScoredSingle {
                    score: 0,
                    record_id: current_id,
                    is_paired_marker: true,
                },
                &mut dup_bits,
            );

            if let Some(mate) = pending.remove(nh, ch) {
                let (ref_lo, pos_lo, rev_lo_raw, ref_hi, pos_hi, rev_hi_raw) =
                    if (tid_i32, uc5) <= (mate.ref_id, mate.unclipped_5prime) {
                        (
                            tid_i32,
                            uc5,
                            is_reverse,
                            mate.ref_id,
                            mate.unclipped_5prime,
                            mate.is_reverse,
                        )
                    } else {
                        (
                            mate.ref_id,
                            mate.unclipped_5prime,
                            mate.is_reverse,
                            tid_i32,
                            uc5,
                            is_reverse,
                        )
                    };

                // Picard same-position RF→FR normalization (MarkDuplicates §3 of research).
                // When both reads land at identical (ref, unclipped_5') and the orientation
                // after lo/hi ordering is RF (rev_lo=true, rev_hi=false), Picard forces it
                // to FR so RF and FR at the same locus group as one duplicate set.
                let (rev_lo, rev_hi) =
                    if ref_lo == ref_hi && pos_lo == pos_hi && rev_lo_raw && !rev_hi_raw {
                        (false, true)
                    } else {
                        (rev_lo_raw, rev_hi_raw)
                    };

                // BARCODE_TAG: both mates carry the same value per Picard
                // convention, so either side's hash is equivalent; we use the
                // current record's. READ_ONE/TWO: assignment is firstOfPair-
                // flag-driven, NOT lo/hi-coord driven (umi_semantics.md Q5).
                let (read1_barcode_hash, read2_barcode_hash) = if is_first_of_pair {
                    (read_barcode_hash, mate.read_barcode_hash)
                } else {
                    (mate.read_barcode_hash, read_barcode_hash)
                };
                paired_tracker.add_pair(
                    PairedEndKey {
                        library_idx: lib_idx,
                        barcode_hash,
                        read1_barcode_hash,
                        read2_barcode_hash,
                        ref_id_lo: ref_lo,
                        pos_lo,
                        is_reverse_lo: rev_lo,
                        ref_id_hi: ref_hi,
                        pos_hi,
                        is_reverse_hi: rev_hi,
                    },
                    ScoredPair {
                        combined_score: qsum + mate.quality_sum,
                        record_id_1: mate.record_id,
                        record_id_2: current_id,
                    },
                );
                paired_tracker.resolve_up_to(tid_i32, pos, &mut dup_bits);
            } else {
                pending.insert(PendingMate {
                    name_hash: nh,
                    check_hash: ch,
                    ref_id: tid_i32,
                    unclipped_5prime: uc5,
                    is_reverse,
                    mate_ref_id: mate_tid,
                    mate_pos,
                    quality_sum: qsum,
                    record_id: current_id,
                    library_idx: lib_idx,
                    barcode_hash,
                    read_barcode_hash,
                    is_first_of_pair,
                });
            }
        }
    }

    single_tracker.flush(&mut dup_bits);
    paired_tracker.resolve_all(&mut dup_bits);

    let orphans = pending.drain_all();
    if !orphans.is_empty() {
        warn!("{} orphan reads passed through unflagged", orphans.len());
    }

    counters.read_pairs_examined = paired_tracker.pairs_examined;
    counters.read_pair_duplicates = paired_tracker.pair_duplicates;
    // Picard (research §7): an orphan paired read (mate never arrives)
    // is counted as one `UNPAIRED_READS_EXAMINED`. We also leave the
    // paired-marker inserted at its locus — see deviations.md for the
    // known divergence where an orphan cannot compete with a co-located
    // fragment by score; markers always win at that locus.
    counters.unpaired_reads_examined = single_tracker.reads_examined + orphans.len() as u64;
    counters.unpaired_read_duplicates = single_tracker.read_duplicates;
    counters.merge_histogram(&paired_tracker.group_sizes);
    counters.merge_histogram(&single_tracker.group_sizes);

    info!(
        "Pass 1: {} records, {} pair dups, {} single dups, peak pending {}",
        record_id,
        counters.read_pair_duplicates,
        counters.unpaired_read_duplicates,
        pending.peak()
    );

    Ok(ScanResult {
        dup_bits,
        counters,
        total_records: record_id,
        pending_peak: pending.peak(),
    })
}
