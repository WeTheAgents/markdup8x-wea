//! Pending mate tracking for paired-end reads in coordinate-sorted BAM.
//!
//! When we encounter read1 of a pair, its mate hasn't been seen yet.
//! Store minimal metadata keyed by QNAME hash until the mate arrives.

use rustc_hash::FxHashMap;

/// Minimal metadata for a pending mate. Packed to ~60 bytes (post Track A.1).
#[derive(Clone, Debug)]
pub struct PendingMate {
    pub name_hash: u64,
    pub check_hash: u32, // first 4 bytes of QNAME for collision detection
    pub ref_id: i32,
    pub unclipped_5prime: i64,
    pub is_reverse: bool,
    pub mate_ref_id: i32,
    pub mate_pos: i64, // aligned mate pos from MPOS field
    pub quality_sum: u32,
    pub record_id: u64,
    pub library_idx: u8,
    /// BARCODE_TAG hash carried from this record's RX (or equivalent). 0 when
    /// the feature is off; A.2 wires the extractor.
    pub barcode_hash: i32,
    /// READ_ONE_BARCODE_TAG or READ_TWO_BARCODE_TAG hash for this record —
    /// which one depends on the firstOfPair flag at pair-completion time.
    /// 0 when the feature is off.
    pub read_barcode_hash: i32,
}

/// Compute 64-bit FxHash of a mate-lookup key.
///
/// Picard (research §15) keys its mate lookup by `readGroupId + readName` so
/// reads with identical QNAMEs in different read groups do not accidentally
/// pair. This must stay scoped to the raw RG tag value, not the dedup
/// `library_idx`: different RGs may share the same LB and still must not be
/// eligible to mate with each other.
pub fn qname_hash(qname: &[u8], read_group: Option<&[u8]>) -> u64 {
    use std::hash::Hash;
    let mut hasher = rustc_hash::FxHasher::default();
    read_group.unwrap_or(b"").hash(&mut hasher);
    qname.hash(&mut hasher);
    std::hash::Hasher::finish(&hasher)
}

/// Compute a truncated check hash for collision detection inside one bucket.
pub fn check_hash(qname: &[u8], read_group: Option<&[u8]>) -> u32 {
    use std::hash::Hash;
    let mut hasher = rustc_hash::FxHasher::default();
    read_group.unwrap_or(b"").hash(&mut hasher);
    qname.hash(&mut hasher);
    std::hash::Hasher::finish(&hasher) as u32
}

/// Buffer for pending mates. Keyed by QNAME hash.
/// Handles collisions via check_hash + Vec for multiple entries per hash.
pub struct PendingMateBuffer {
    map: FxHashMap<u64, Vec<PendingMate>>,
    count: usize,
    peak: usize,
}

impl PendingMateBuffer {
    pub fn new() -> Self {
        Self {
            map: FxHashMap::default(),
            count: 0,
            peak: 0,
        }
    }

    /// Insert a pending mate. Returns None.
    pub fn insert(&mut self, mate: PendingMate) {
        self.map.entry(mate.name_hash).or_default().push(mate);
        self.count += 1;
        if self.count > self.peak {
            self.peak = self.count;
        }
    }

    /// Look up and remove a matching pending mate by QNAME hash + check_hash.
    /// Returns the mate if found.
    pub fn remove(&mut self, name_hash: u64, chk: u32) -> Option<PendingMate> {
        if let Some(entries) = self.map.get_mut(&name_hash) {
            if let Some(idx) = entries.iter().position(|m| m.check_hash == chk) {
                let mate = entries.swap_remove(idx);
                if entries.is_empty() {
                    self.map.remove(&name_hash);
                }
                self.count -= 1;
                return Some(mate);
            }
        }
        None
    }

    /// Number of pending mates currently stored.
    pub fn len(&self) -> usize {
        self.count
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Peak number of pending mates seen.
    pub fn peak(&self) -> usize {
        self.peak
    }

    /// Drain all remaining pending mates (orphans at EOF).
    pub fn drain_all(&mut self) -> Vec<PendingMate> {
        let mut orphans = Vec::new();
        for (_, entries) in self.map.drain() {
            orphans.extend(entries);
        }
        self.count = 0;
        orphans
    }
}

impl Default for PendingMateBuffer {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_pending(qname: &[u8], read_group: Option<&[u8]>, record_id: u64) -> PendingMate {
        PendingMate {
            name_hash: qname_hash(qname, read_group),
            check_hash: check_hash(qname, read_group),
            ref_id: 0,
            unclipped_5prime: 1000,
            is_reverse: false,
            mate_ref_id: 0,
            mate_pos: 1500,
            quality_sum: 500,
            record_id,
            library_idx: 0,
            barcode_hash: 0,
            read_barcode_hash: 0,
        }
    }

    #[test]
    fn insert_and_remove() {
        let mut buf = PendingMateBuffer::new();
        let qname = b"HISEQ:1:1:1000:2000";
        let mate = make_pending(qname, Some(b"rg1"), 0);
        let nh = mate.name_hash;
        let ch = mate.check_hash;

        buf.insert(mate);
        assert_eq!(buf.len(), 1);

        let found = buf.remove(nh, ch);
        assert!(found.is_some());
        assert_eq!(found.unwrap().record_id, 0);
        assert_eq!(buf.len(), 0);
    }

    #[test]
    fn remove_nonexistent() {
        let mut buf = PendingMateBuffer::new();
        assert!(buf.remove(12345, 0).is_none());
    }

    #[test]
    fn peak_tracking() {
        let mut buf = PendingMateBuffer::new();
        buf.insert(make_pending(b"READ_A", Some(b"rg1"), 0));
        buf.insert(make_pending(b"READ_B", Some(b"rg1"), 1));
        buf.insert(make_pending(b"READ_C", Some(b"rg1"), 2));
        assert_eq!(buf.peak(), 3);

        let nh = qname_hash(b"READ_A", Some(b"rg1"));
        let ch = check_hash(b"READ_A", Some(b"rg1"));
        buf.remove(nh, ch);
        assert_eq!(buf.len(), 2);
        assert_eq!(buf.peak(), 3); // Peak doesn't decrease
    }

    #[test]
    fn drain_orphans() {
        let mut buf = PendingMateBuffer::new();
        buf.insert(make_pending(b"ORPHAN_1", Some(b"rg1"), 0));
        buf.insert(make_pending(b"ORPHAN_2", Some(b"rg1"), 1));

        let orphans = buf.drain_all();
        assert_eq!(orphans.len(), 2);
        assert_eq!(buf.len(), 0);
    }

    #[test]
    fn collision_handling() {
        // Two different QNAMEs stored under same hash bucket
        let mut buf = PendingMateBuffer::new();
        let m1 = make_pending(b"READ_A", Some(b"rg1"), 0);
        let m2 = make_pending(b"READ_B", Some(b"rg1"), 1);

        // Even if they happen to share a hash, check_hash distinguishes them
        buf.insert(m1.clone());
        buf.insert(m2.clone());

        // Remove READ_A specifically
        let found = buf.remove(m1.name_hash, m1.check_hash);
        assert!(found.is_some());
        assert_eq!(found.unwrap().record_id, 0);

        // READ_B still there
        let found2 = buf.remove(m2.name_hash, m2.check_hash);
        assert!(found2.is_some());
        assert_eq!(found2.unwrap().record_id, 1);
    }

    #[test]
    fn same_qname_different_read_groups_use_different_lookup_keys() {
        let mut buf = PendingMateBuffer::new();
        let rg1 = make_pending(b"SHARED", Some(b"rg1"), 0);
        let rg2 = make_pending(b"SHARED", Some(b"rg2"), 1);

        assert_ne!(rg1.name_hash, rg2.name_hash);
        assert_ne!(rg1.check_hash, rg2.check_hash);

        buf.insert(rg1.clone());
        buf.insert(rg2.clone());

        let found = buf.remove(rg2.name_hash, rg2.check_hash).unwrap();
        assert_eq!(found.record_id, 1);

        let found = buf.remove(rg1.name_hash, rg1.check_hash).unwrap();
        assert_eq!(found.record_id, 0);
    }
}
