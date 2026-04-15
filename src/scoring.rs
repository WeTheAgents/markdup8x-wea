//! Quality-based scoring for duplicate group resolution.
//!
//! Picard's SUM_OF_BASE_QUALITIES strategy: sum all base qualities >= 15.

/// Minimum base quality to include in the score sum.
const MIN_BASE_QUALITY: u8 = 15;

/// Compute the quality score for a single read: sum of base qualities >= 15.
pub fn quality_sum(quals: &[u8]) -> u32 {
    quals
        .iter()
        .filter(|&&q| q >= MIN_BASE_QUALITY)
        .map(|&q| q as u32)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_high_quality() {
        assert_eq!(quality_sum(&[30, 25, 20, 35, 40]), 150);
    }

    #[test]
    fn mixed_quality() {
        assert_eq!(quality_sum(&[30, 10, 20, 5, 40]), 90);
    }

    #[test]
    fn all_low_quality() {
        assert_eq!(quality_sum(&[10, 8, 12, 14, 5]), 0);
    }

    #[test]
    fn empty_quals() {
        assert_eq!(quality_sum(&[]), 0);
    }

    #[test]
    fn edge_at_threshold() {
        assert_eq!(quality_sum(&[15, 14, 15, 14]), 30);
    }

    #[test]
    fn single_base() {
        assert_eq!(quality_sum(&[40]), 40);
    }

    #[test]
    fn max_quality() {
        assert_eq!(quality_sum(&[93, 93, 93]), 279);
    }
}
