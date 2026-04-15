//! Pass 1: Scan BAM, build duplicate groups, resolve, return dupset.

use crate::dupset::BitVecDupSet;
use crate::groups::{
    PairedEndKey, PairedGroupTracker, ScoredPair, ScoredSingle, SingleEndKey, SingleEndTracker,
};
use crate::io::SortOrderEnforcer;
use crate::metrics::MetricsCounters;
use crate::pending_mates::{check_hash, qname_hash, PendingMate, PendingMateBuffer};
use crate::position::unclipped_5prime;
use crate::scoring::quality_sum;
use anyhow::{Context, Result};
use log::{info, warn};
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::record::Flags;
use noodles::sam::header::record::value::map::read_group::tag::LIBRARY;
use rustc_hash::FxHashMap;

/// Build library mapping from @RG headers.
fn build_library_map(header: &sam::Header) -> FxHashMap<Vec<u8>, u8> {
    let mut lib_names: Vec<String> = Vec::new();
    let mut rg_to_lib: FxHashMap<Vec<u8>, u8> = FxHashMap::default();

    for (id, rg) in header.read_groups().iter() {
        let lib_name = rg
            .other_fields()
            .get(&LIBRARY)
            .map(|v| String::from_utf8_lossy(v.as_ref()).to_string())
            .unwrap_or_else(|| "unknown".to_string());

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

    rg_to_lib
}

fn get_library_idx(record: &bam::Record, rg_to_lib: &FxHashMap<Vec<u8>, u8>) -> u8 {
    if rg_to_lib.is_empty() {
        return 0;
    }
    use noodles::sam::alignment::record::data::field::Tag;
    if let Some(Ok(value)) = record.data().get(&Tag::READ_GROUP) {
        use noodles::sam::alignment::record::data::field::Value;
        if let Value::String(s) = value {
            let s_bytes: &[u8] = s.as_ref();
            if let Some(&idx) = rg_to_lib.get(s_bytes) {
                return idx;
            }
        }
    }
    0
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

/// Run Pass 1.
pub fn scan_pass(
    reader: &mut crate::io::BamReader,
    header: &sam::Header,
) -> Result<ScanResult> {
    let rg_to_lib = build_library_map(header);

    let mut dup_bits = BitVecDupSet::new();
    let mut pending = PendingMateBuffer::new();
    let mut paired_tracker = PairedGroupTracker::new();
    let mut single_tracker = SingleEndTracker::new();
    let mut enforcer = SortOrderEnforcer::new();
    let mut counters = MetricsCounters::new();

    let mut record_id: u64 = 0;
    let mut prev_tid: Option<usize> = None;
    let mut record = bam::Record::default();

    while reader.read_record(&mut record)? > 0 {
        let flags = record.flags();
        let tid = record.reference_sequence_id().transpose().ok().flatten();
        let pos = record
            .alignment_start()
            .transpose()
            .ok()
            .flatten()
            .map(|p| (usize::from(p) as i64) - 1)
            .unwrap_or(-1);

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

        let cigar_ops = extract_cigar_ops(&record)?;
        let uc5 = unclipped_5prime(pos, &cigar_ops, is_reverse);
        let quals = extract_quals(&record);
        let qsum = quality_sum(&quals);
        let lib_idx = get_library_idx(&record, &rg_to_lib);
        let tid_i32 = tid.map(|t| t as i32).unwrap_or(-1);

        if !is_paired || mate_unmapped {
            single_tracker.add_read(
                SingleEndKey {
                    library_idx: lib_idx,
                    ref_id: tid_i32,
                    unclipped_5prime: uc5,
                    is_reverse,
                },
                ScoredSingle { score: qsum, record_id: current_id },
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
            let nh = qname_hash(qname_bytes);
            let ch = check_hash(qname_bytes);
            let mate_tid = record
                .mate_reference_sequence_id()
                .transpose()
                .ok()
                .flatten()
                .map(|t| t as i32)
                .unwrap_or(-1);
            let mate_pos = record
                .mate_alignment_start()
                .transpose()
                .ok()
                .flatten()
                .map(|p| (usize::from(p) as i64) - 1)
                .unwrap_or(-1);

            if let Some(mate) = pending.remove(nh, ch) {
                let (ref_lo, pos_lo, rev_lo, ref_hi, pos_hi, rev_hi) =
                    if (tid_i32, uc5) <= (mate.ref_id, mate.unclipped_5prime) {
                        (tid_i32, uc5, is_reverse, mate.ref_id, mate.unclipped_5prime, mate.is_reverse)
                    } else {
                        (mate.ref_id, mate.unclipped_5prime, mate.is_reverse, tid_i32, uc5, is_reverse)
                    };

                paired_tracker.add_pair(
                    PairedEndKey {
                        library_idx: lib_idx,
                        ref_id_lo: ref_lo, pos_lo, is_reverse_lo: rev_lo,
                        ref_id_hi: ref_hi, pos_hi, is_reverse_hi: rev_hi,
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
    counters.unpaired_reads_examined = single_tracker.reads_examined;
    counters.unpaired_read_duplicates = single_tracker.read_duplicates;
    counters.merge_histogram(&paired_tracker.group_sizes);
    counters.merge_histogram(&single_tracker.group_sizes);

    info!(
        "Pass 1: {} records, {} pair dups, {} single dups, peak pending {}",
        record_id, counters.read_pair_duplicates, counters.unpaired_read_duplicates, pending.peak()
    );

    Ok(ScanResult {
        dup_bits, counters,
        total_records: record_id,
        pending_peak: pending.peak(),
    })
}
