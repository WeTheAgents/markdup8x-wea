//! DupSet trait and implementations: BitVec vs FxHashSet.
//!
//! Both track which record IDs are duplicates. Benchmarked on yeast to pick winner.

use bitvec::prelude::*;
use rustc_hash::FxHashSet;

/// Trait for tracking duplicate record IDs.
pub trait DupSet {
    fn insert(&mut self, record_id: u64);
    fn contains(&self, record_id: u64) -> bool;
    fn len(&self) -> usize;
}

/// BitVec-backed dupset: 1 bit per record, O(1) lookup, ~25MB for 200M reads.
pub struct BitVecDupSet {
    bits: BitVec,
    count: usize,
}

impl BitVecDupSet {
    pub fn new() -> Self {
        Self {
            bits: BitVec::new(),
            count: 0,
        }
    }

    pub fn with_capacity(n: usize) -> Self {
        let mut bits = BitVec::with_capacity(n);
        bits.resize(n, false);
        Self { bits, count: 0 }
    }
}

impl DupSet for BitVecDupSet {
    fn insert(&mut self, record_id: u64) {
        let idx = record_id as usize;
        if idx >= self.bits.len() {
            self.bits.resize(idx + 1, false);
        }
        if !self.bits[idx] {
            self.count += 1;
        }
        self.bits.set(idx, true);
    }

    fn contains(&self, record_id: u64) -> bool {
        let idx = record_id as usize;
        idx < self.bits.len() && self.bits[idx]
    }

    fn len(&self) -> usize {
        self.count
    }
}

/// FxHashSet-backed dupset: stores only dup IDs, simpler but more memory.
pub struct HashDupSet {
    set: FxHashSet<u64>,
}

impl HashDupSet {
    pub fn new() -> Self {
        Self {
            set: FxHashSet::default(),
        }
    }
}

impl DupSet for HashDupSet {
    fn insert(&mut self, record_id: u64) {
        self.set.insert(record_id);
    }

    fn contains(&self, record_id: u64) -> bool {
        self.set.contains(&record_id)
    }

    fn len(&self) -> usize {
        self.set.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_dupset(ds: &mut dyn DupSet) {
        assert_eq!(ds.len(), 0);
        assert!(!ds.contains(0));
        assert!(!ds.contains(100));

        ds.insert(5);
        assert!(ds.contains(5));
        assert!(!ds.contains(4));
        assert_eq!(ds.len(), 1);

        ds.insert(100);
        ds.insert(200);
        assert_eq!(ds.len(), 3);
        assert!(ds.contains(100));
        assert!(ds.contains(200));

        // Double insert doesn't change count
        ds.insert(5);
        assert_eq!(ds.len(), 3);
    }

    #[test]
    fn bitvec_dupset() {
        let mut ds = BitVecDupSet::new();
        test_dupset(&mut ds);
    }

    #[test]
    fn hash_dupset() {
        let mut ds = HashDupSet::new();
        test_dupset(&mut ds);
    }

    #[test]
    fn bitvec_with_capacity() {
        let mut ds = BitVecDupSet::with_capacity(1000);
        ds.insert(999);
        assert!(ds.contains(999));
        // Beyond capacity — auto-grows
        ds.insert(2000);
        assert!(ds.contains(2000));
    }

    #[test]
    fn both_produce_same_results() {
        let ids: Vec<u64> = vec![0, 5, 10, 100, 999, 1_000_000];
        let mut bv = BitVecDupSet::new();
        let mut hs = HashDupSet::new();

        for &id in &ids {
            bv.insert(id);
            hs.insert(id);
        }

        for check in 0..1_000_010u64 {
            assert_eq!(
                bv.contains(check),
                hs.contains(check),
                "Mismatch at record_id {}",
                check
            );
        }
        assert_eq!(bv.len(), hs.len());
    }
}
