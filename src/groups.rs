//! Duplicate group tracking and resolution.
//!
//! Groups reads by PairedEndKey or SingleEndKey. Resolves groups by score:
//! highest combined_score wins, rest are flagged as duplicates.

use crate::dupset::DupSet;
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;

/// Grouping key for paired-end reads.
/// Ordered: lo = lower (ref_id, position), hi = higher.
///
/// Field order mirrors Picard's `ReadEndsMDComparator` (library → barcode →
/// read1/read2 → coords). Barcode hash defaults to 0 when no `BARCODE_TAG` /
/// `READ_ONE_BARCODE_TAG` / `READ_TWO_BARCODE_TAG` is configured — keys then
/// behave identically to the pre-Track-A layout. See `docs/umi_semantics.md`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PairedEndKey {
    pub library_idx: u8,
    pub barcode_hash: i32,
    pub read1_barcode_hash: i32,
    pub read2_barcode_hash: i32,
    pub ref_id_lo: i32,
    pub pos_lo: i64,
    pub is_reverse_lo: bool,
    pub ref_id_hi: i32,
    pub pos_hi: i64,
    pub is_reverse_hi: bool,
}

/// Grouping key for single-end reads.
///
/// `barcode_hash` defaults to 0 when no `BARCODE_TAG` is configured.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct SingleEndKey {
    pub library_idx: u8,
    pub barcode_hash: i32,
    pub ref_id: i32,
    pub unclipped_5prime: i64,
    pub is_reverse: bool,
}

/// A scored pair in a duplicate group.
#[derive(Clone, Debug)]
pub struct ScoredPair {
    pub combined_score: u32,
    pub record_id_1: u64,
    pub record_id_2: u64,
}

/// A scored single-end read in a duplicate group.
///
/// `is_paired_marker=true` denotes a placeholder inserted when a read of a
/// properly-paired pair is encountered at this position. Per Picard (research
/// §4), when any such marker sits in a fragment group, every non-marker
/// (i.e. true single-end/orphan fragment) in that group is unconditionally
/// flagged as a duplicate regardless of score — paired reads always beat
/// fragments at the same locus.
#[derive(Clone, Debug)]
pub struct ScoredSingle {
    pub score: u32,
    pub record_id: u64,
    pub is_paired_marker: bool,
}

/// Tracks active paired-end duplicate groups and resolves them incrementally.
pub struct PairedGroupTracker {
    groups: FxHashMap<PairedEndKey, Vec<ScoredPair>>,
    /// Maps max_pos -> list of keys to resolve when stream passes that position.
    /// Uses (ref_id, pos) encoded as i64 for same-chromosome groups.
    resolve_at: BTreeMap<(i32, i64), Vec<PairedEndKey>>,
    pub pairs_examined: u64,
    pub pair_duplicates: u64,
    pub group_sizes: Vec<u64>, // histogram: index = group_size, value = count
}

impl PairedGroupTracker {
    pub fn new() -> Self {
        Self {
            groups: FxHashMap::default(),
            resolve_at: BTreeMap::new(),
            pairs_examined: 0,
            pair_duplicates: 0,
            group_sizes: Vec::new(),
        }
    }

    /// Add a scored pair to its group.
    pub fn add_pair(&mut self, key: PairedEndKey, pair: ScoredPair) {
        self.pairs_examined += 1;
        let max_pos = (key.ref_id_hi, key.pos_hi);
        self.resolve_at
            .entry(max_pos)
            .or_default()
            .push(key.clone());
        self.groups.entry(key).or_default().push(pair);
    }

    /// Resolve all groups whose max_pos is <= (ref_id, pos).
    /// Returns number of groups resolved.
    pub fn resolve_up_to(
        &mut self,
        ref_id: i32,
        pos: i64,
        dup_bits: &mut dyn DupSet,
    ) -> usize {
        let cutoff = (ref_id, pos);
        let keys_to_resolve: Vec<Vec<PairedEndKey>> = self
            .resolve_at
            .range(..=cutoff)
            .map(|(_, keys)| keys.clone())
            .collect();

        // Remove resolved entries from BTreeMap
        let positions: Vec<(i32, i64)> = self
            .resolve_at
            .range(..=cutoff)
            .map(|(k, _)| *k)
            .collect();
        for p in &positions {
            self.resolve_at.remove(p);
        }

        let mut resolved = 0;
        for keys in keys_to_resolve {
            for key in keys {
                if let Some(pairs) = self.groups.remove(&key) {
                    self.resolve_group(pairs, dup_bits);
                    resolved += 1;
                }
            }
        }
        resolved
    }

    /// Resolve all groups for a specific chromosome (at chromosome boundary).
    pub fn resolve_chromosome(&mut self, ref_id: i32, dup_bits: &mut dyn DupSet) -> usize {
        self.resolve_up_to(ref_id, i64::MAX, dup_bits)
    }

    /// Resolve all remaining groups (at EOF).
    pub fn resolve_all(&mut self, dup_bits: &mut dyn DupSet) -> usize {
        let all_keys: Vec<PairedEndKey> = self.groups.keys().cloned().collect();
        let mut resolved = 0;
        for key in all_keys {
            if let Some(pairs) = self.groups.remove(&key) {
                self.resolve_group(pairs, dup_bits);
                resolved += 1;
            }
        }
        self.resolve_at.clear();
        resolved
    }

    fn resolve_group(&mut self, mut pairs: Vec<ScoredPair>, dup_bits: &mut dyn DupSet) {
        // Track histogram
        let size = pairs.len();
        if size >= self.group_sizes.len() {
            self.group_sizes.resize(size + 1, 0);
        }
        self.group_sizes[size] += 1;

        if pairs.len() <= 1 {
            return; // No duplicates in a group of 1
        }

        // Sort: highest score first, then lowest record_id (tie-breaking: first encountered wins)
        pairs.sort_by(|a, b| {
            b.combined_score
                .cmp(&a.combined_score)
                .then(a.record_id_1.cmp(&b.record_id_1))
        });

        // First pair is primary; all others are duplicates
        for pair in pairs.iter().skip(1) {
            dup_bits.insert(pair.record_id_1);
            dup_bits.insert(pair.record_id_2);
            self.pair_duplicates += 1;
        }
    }

    pub fn active_groups(&self) -> usize {
        self.groups.len()
    }
}

impl Default for PairedGroupTracker {
    fn default() -> Self {
        Self::new()
    }
}

/// Inline single-end group tracker. Resolves consecutive groups immediately.
pub struct SingleEndTracker {
    current_key: Option<SingleEndKey>,
    current_group: Vec<ScoredSingle>,
    pub reads_examined: u64,
    pub read_duplicates: u64,
    pub group_sizes: Vec<u64>,
}

impl SingleEndTracker {
    pub fn new() -> Self {
        Self {
            current_key: None,
            current_group: Vec::new(),
            reads_examined: 0,
            read_duplicates: 0,
            group_sizes: Vec::new(),
        }
    }

    /// Add a single-end read. If key differs from current group, resolve the current group first.
    ///
    /// `reads_examined` counts only real fragments (not paired markers), to
    /// match Picard's `UNPAIRED_READS_EXAMINED` semantics.
    pub fn add_read(
        &mut self,
        key: SingleEndKey,
        scored: ScoredSingle,
        dup_bits: &mut dyn DupSet,
    ) {
        if !scored.is_paired_marker {
            self.reads_examined += 1;
        }

        if self.current_key.as_ref() != Some(&key) {
            // New group — resolve previous
            self.resolve_current(dup_bits);
            self.current_key = Some(key);
        }
        self.current_group.push(scored);
    }

    /// Flush the current group (call at chromosome boundary and EOF).
    pub fn flush(&mut self, dup_bits: &mut dyn DupSet) {
        self.resolve_current(dup_bits);
        self.current_key = None;
    }

    fn resolve_current(&mut self, dup_bits: &mut dyn DupSet) {
        let size = self.current_group.len();
        if size == 0 {
            return;
        }

        // Histogram counts only real fragments, not paired markers.
        let fragment_count = self
            .current_group
            .iter()
            .filter(|s| !s.is_paired_marker)
            .count();
        if fragment_count >= self.group_sizes.len() {
            self.group_sizes.resize(fragment_count + 1, 0);
        }
        if fragment_count > 0 {
            self.group_sizes[fragment_count] += 1;
        }

        let has_paired_marker = self.current_group.iter().any(|s| s.is_paired_marker);

        if has_paired_marker {
            // Picard §4: every real fragment at a locus with any paired-read
            // presence is unconditionally flagged regardless of score.
            for read in self.current_group.iter() {
                if !read.is_paired_marker {
                    dup_bits.insert(read.record_id);
                    self.read_duplicates += 1;
                }
            }
            self.current_group.clear();
            return;
        }

        if fragment_count <= 1 {
            self.current_group.clear();
            return;
        }

        // Sort: highest score first, then lowest record_id
        self.current_group.sort_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then(a.record_id.cmp(&b.record_id))
        });

        for read in self.current_group.iter().skip(1) {
            dup_bits.insert(read.record_id);
            self.read_duplicates += 1;
        }

        self.current_group.clear();
    }
}

impl Default for SingleEndTracker {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dupset::BitVecDupSet;

    #[test]
    fn paired_group_single_pair_no_dup() {
        let mut tracker = PairedGroupTracker::new();
        let mut dup_bits = BitVecDupSet::new();

        let key = PairedEndKey {
            library_idx: 0,
            barcode_hash: 0,
            read1_barcode_hash: 0,
            read2_barcode_hash: 0,
            ref_id_lo: 0,
            pos_lo: 1000,
            is_reverse_lo: false,
            ref_id_hi: 0,
            pos_hi: 1500,
            is_reverse_hi: true,
        };
        tracker.add_pair(
            key,
            ScoredPair {
                combined_score: 800,
                record_id_1: 0,
                record_id_2: 5,
            },
        );

        tracker.resolve_all(&mut dup_bits);
        assert_eq!(dup_bits.len(), 0); // No duplicates
        assert_eq!(tracker.pair_duplicates, 0);
    }

    #[test]
    fn paired_group_two_pairs_one_dup() {
        let mut tracker = PairedGroupTracker::new();
        let mut dup_bits = BitVecDupSet::new();

        let key = PairedEndKey {
            library_idx: 0,
            barcode_hash: 0,
            read1_barcode_hash: 0,
            read2_barcode_hash: 0,
            ref_id_lo: 0,
            pos_lo: 1000,
            is_reverse_lo: false,
            ref_id_hi: 0,
            pos_hi: 1500,
            is_reverse_hi: true,
        };
        tracker.add_pair(
            key.clone(),
            ScoredPair {
                combined_score: 800,
                record_id_1: 0,
                record_id_2: 5,
            },
        );
        tracker.add_pair(
            key,
            ScoredPair {
                combined_score: 600,
                record_id_1: 1,
                record_id_2: 6,
            },
        );

        tracker.resolve_all(&mut dup_bits);
        assert_eq!(dup_bits.len(), 2); // Both reads of losing pair
        assert!(!dup_bits.contains(0)); // Winner read1
        assert!(!dup_bits.contains(5)); // Winner read2
        assert!(dup_bits.contains(1)); // Loser read1
        assert!(dup_bits.contains(6)); // Loser read2
        assert_eq!(tracker.pair_duplicates, 1);
    }

    #[test]
    fn paired_group_tie_breaking_by_record_id() {
        let mut tracker = PairedGroupTracker::new();
        let mut dup_bits = BitVecDupSet::new();

        let key = PairedEndKey {
            library_idx: 0,
            barcode_hash: 0,
            read1_barcode_hash: 0,
            read2_barcode_hash: 0,
            ref_id_lo: 0,
            pos_lo: 1000,
            is_reverse_lo: false,
            ref_id_hi: 0,
            pos_hi: 1500,
            is_reverse_hi: true,
        };
        // Same score, pair with lower record_id_1 wins
        tracker.add_pair(
            key.clone(),
            ScoredPair {
                combined_score: 800,
                record_id_1: 10,
                record_id_2: 15,
            },
        );
        tracker.add_pair(
            key,
            ScoredPair {
                combined_score: 800,
                record_id_1: 2,
                record_id_2: 7,
            },
        );

        tracker.resolve_all(&mut dup_bits);
        // record_id_1=2 has lower id → wins
        assert!(!dup_bits.contains(2));
        assert!(!dup_bits.contains(7));
        assert!(dup_bits.contains(10));
        assert!(dup_bits.contains(15));
    }

    #[test]
    fn incremental_resolution() {
        let mut tracker = PairedGroupTracker::new();
        let mut dup_bits = BitVecDupSet::new();

        // Group at pos 1000/1500
        let key1 = PairedEndKey {
            library_idx: 0,
            barcode_hash: 0,
            read1_barcode_hash: 0,
            read2_barcode_hash: 0,
            ref_id_lo: 0,
            pos_lo: 1000,
            is_reverse_lo: false,
            ref_id_hi: 0,
            pos_hi: 1500,
            is_reverse_hi: true,
        };
        tracker.add_pair(
            key1.clone(),
            ScoredPair {
                combined_score: 800,
                record_id_1: 0,
                record_id_2: 1,
            },
        );
        tracker.add_pair(
            key1,
            ScoredPair {
                combined_score: 600,
                record_id_1: 2,
                record_id_2: 3,
            },
        );

        // Group at pos 2000/3000
        let key2 = PairedEndKey {
            library_idx: 0,
            barcode_hash: 0,
            read1_barcode_hash: 0,
            read2_barcode_hash: 0,
            ref_id_lo: 0,
            pos_lo: 2000,
            is_reverse_lo: false,
            ref_id_hi: 0,
            pos_hi: 3000,
            is_reverse_hi: true,
        };
        tracker.add_pair(
            key2,
            ScoredPair {
                combined_score: 900,
                record_id_1: 4,
                record_id_2: 5,
            },
        );

        // Resolve up to pos 2000 — should resolve key1 but not key2
        let resolved = tracker.resolve_up_to(0, 2000, &mut dup_bits);
        assert_eq!(resolved, 1);
        assert!(dup_bits.contains(2)); // key1 loser
        assert!(!dup_bits.contains(4)); // key2 not resolved yet
        assert_eq!(tracker.active_groups(), 1);
    }

    #[test]
    fn single_end_inline_resolution() {
        let mut tracker = SingleEndTracker::new();
        let mut dup_bits = BitVecDupSet::new();

        let key = SingleEndKey {
            library_idx: 0,
            barcode_hash: 0,
            ref_id: 0,
            unclipped_5prime: 1000,
            is_reverse: false,
        };

        tracker.add_read(
            key.clone(),
            ScoredSingle {
                score: 500,
                record_id: 0,
                is_paired_marker: false,
            },
            &mut dup_bits,
        );
        tracker.add_read(
            key.clone(),
            ScoredSingle {
                score: 300,
                record_id: 1,
                is_paired_marker: false,
            },
            &mut dup_bits,
        );

        // Different key triggers resolution of previous group
        let key2 = SingleEndKey {
            library_idx: 0,
            barcode_hash: 0,
            ref_id: 0,
            unclipped_5prime: 2000,
            is_reverse: false,
        };
        tracker.add_read(
            key2,
            ScoredSingle {
                score: 400,
                record_id: 2,
                is_paired_marker: false,
            },
            &mut dup_bits,
        );

        assert!(dup_bits.contains(1)); // Loser from first group
        assert!(!dup_bits.contains(0)); // Winner
        assert!(!dup_bits.contains(2)); // Only member of second group

        tracker.flush(&mut dup_bits);
        assert_eq!(tracker.read_duplicates, 1);
    }

    #[test]
    fn different_library_different_group() {
        let mut tracker = PairedGroupTracker::new();
        let mut dup_bits = BitVecDupSet::new();

        let key_lib0 = PairedEndKey {
            library_idx: 0,
            barcode_hash: 0,
            read1_barcode_hash: 0,
            read2_barcode_hash: 0,
            ref_id_lo: 0,
            pos_lo: 1000,
            is_reverse_lo: false,
            ref_id_hi: 0,
            pos_hi: 1500,
            is_reverse_hi: true,
        };
        let key_lib1 = PairedEndKey {
            library_idx: 1,
            barcode_hash: 0,
            read1_barcode_hash: 0,
            read2_barcode_hash: 0,
            ..key_lib0.clone()
        };

        tracker.add_pair(
            key_lib0,
            ScoredPair {
                combined_score: 800,
                record_id_1: 0,
                record_id_2: 1,
            },
        );
        tracker.add_pair(
            key_lib1,
            ScoredPair {
                combined_score: 600,
                record_id_1: 2,
                record_id_2: 3,
            },
        );

        tracker.resolve_all(&mut dup_bits);
        // Different libraries → different groups → no dups
        assert_eq!(dup_bits.len(), 0);
    }

    #[test]
    fn fragment_loses_to_paired_marker_regardless_of_score() {
        // Picard §4: a paired marker in the group flags every fragment,
        // even if the fragment's quality score is higher.
        let mut tracker = SingleEndTracker::new();
        let mut dup_bits = BitVecDupSet::new();

        let key = SingleEndKey {
            library_idx: 0,
            barcode_hash: 0,
            ref_id: 0,
            unclipped_5prime: 1000,
            is_reverse: false,
        };

        // High-score fragment first
        tracker.add_read(
            key.clone(),
            ScoredSingle { score: 999, record_id: 7, is_paired_marker: false },
            &mut dup_bits,
        );
        // Then a paired marker (zero score) — marker must still win
        tracker.add_read(
            key,
            ScoredSingle { score: 0, record_id: 42, is_paired_marker: true },
            &mut dup_bits,
        );
        tracker.flush(&mut dup_bits);

        assert!(dup_bits.contains(7), "fragment must be flagged when pair present");
        assert!(!dup_bits.contains(42), "paired marker is not a fragment duplicate");
        assert_eq!(tracker.read_duplicates, 1);
        assert_eq!(tracker.reads_examined, 1); // marker does not count
    }

    #[test]
    fn paired_marker_alone_flags_nothing() {
        let mut tracker = SingleEndTracker::new();
        let mut dup_bits = BitVecDupSet::new();

        let key = SingleEndKey {
            library_idx: 0,
            barcode_hash: 0,
            ref_id: 0,
            unclipped_5prime: 1000,
            is_reverse: false,
        };
        tracker.add_read(
            key,
            ScoredSingle { score: 0, record_id: 1, is_paired_marker: true },
            &mut dup_bits,
        );
        tracker.flush(&mut dup_bits);
        assert_eq!(dup_bits.len(), 0);
        assert_eq!(tracker.reads_examined, 0);
    }

    // ---- Track A.1: barcode_hash dimension in keys ----

    fn paired_template() -> PairedEndKey {
        PairedEndKey {
            library_idx: 0,
            barcode_hash: 0,
            read1_barcode_hash: 0,
            read2_barcode_hash: 0,
            ref_id_lo: 0,
            pos_lo: 1000,
            is_reverse_lo: false,
            ref_id_hi: 0,
            pos_hi: 1500,
            is_reverse_hi: true,
        }
    }

    #[test]
    fn paired_key_distinct_barcode_hash_separate_groups() {
        let k1 = PairedEndKey { barcode_hash: 100, ..paired_template() };
        let k2 = PairedEndKey { barcode_hash: 200, ..paired_template() };
        assert_ne!(k1, k2);
    }

    #[test]
    fn paired_key_picard_missing_barcode_default_31_distinct_from_zero() {
        // BARCODE_TAG missing in Picard hashes to 31 (Objects.hash(null)).
        // The "no flag configured" baseline is 0. The two must not collide.
        let k_missing = PairedEndKey { barcode_hash: 31, ..paired_template() };
        let k_default = PairedEndKey { barcode_hash: 0, ..paired_template() };
        assert_ne!(k_missing, k_default);
    }

    #[test]
    fn paired_key_read1_vs_read2_barcode_distinguish_groups() {
        // (r1=A, r2=B) ≠ (r1=B, r2=A): A.2 must route firstOfPair correctly.
        let k_ab = PairedEndKey {
            read1_barcode_hash: 100,
            read2_barcode_hash: 200,
            ..paired_template()
        };
        let k_ba = PairedEndKey {
            read1_barcode_hash: 200,
            read2_barcode_hash: 100,
            ..paired_template()
        };
        assert_ne!(k_ab, k_ba);
    }

    #[test]
    fn single_key_distinct_barcode_hash_separate_groups() {
        let k1 = SingleEndKey {
            library_idx: 0,
            barcode_hash: 100,
            ref_id: 0,
            unclipped_5prime: 1000,
            is_reverse: false,
        };
        let k2 = SingleEndKey { barcode_hash: 200, ..k1.clone() };
        assert_ne!(k1, k2);
    }
}
