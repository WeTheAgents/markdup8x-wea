//! Unclipped 5-prime position calculation for duplicate grouping.

use noodles::sam::alignment::record::cigar::op::Kind;

/// Compute the unclipped 5-prime position of a BAM record.
///
/// - Forward reads: pos - left_clips (soft + hard at 5' end)
/// - Reverse reads: alignment_end + right_clips (soft + hard at 3' genomic end = 5' read end)
///
/// `pos` is 0-based. Returns 0-based coordinate.
pub fn unclipped_5prime(pos: i64, cigar_ops: &[(Kind, usize)], is_reverse: bool) -> i64 {
    if cigar_ops.is_empty() {
        return pos;
    }

    if is_reverse {
        let mut end = pos;
        for &(kind, len) in cigar_ops {
            match kind {
                Kind::Match
                | Kind::Deletion
                | Kind::Skip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch => {
                    end += len as i64;
                }
                _ => {}
            }
        }
        let right_clips: i64 = cigar_ops
            .iter()
            .rev()
            .take_while(|(kind, _)| matches!(kind, Kind::SoftClip | Kind::HardClip))
            .map(|(_, len)| *len as i64)
            .sum();
        end + right_clips
    } else {
        let left_clips: i64 = cigar_ops
            .iter()
            .take_while(|(kind, _)| matches!(kind, Kind::SoftClip | Kind::HardClip))
            .map(|(_, len)| *len as i64)
            .sum();
        pos - left_clips
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forward_no_clip() {
        let ops = vec![(Kind::Match, 100)];
        assert_eq!(unclipped_5prime(1000, &ops, false), 1000);
    }

    #[test]
    fn forward_soft_clip() {
        let ops = vec![(Kind::SoftClip, 5), (Kind::Match, 95)];
        assert_eq!(unclipped_5prime(1000, &ops, false), 995);
    }

    #[test]
    fn forward_hard_and_soft_clip() {
        let ops = vec![(Kind::HardClip, 3), (Kind::SoftClip, 5), (Kind::Match, 92)];
        assert_eq!(unclipped_5prime(1000, &ops, false), 992);
    }

    #[test]
    fn reverse_soft_clip() {
        let ops = vec![(Kind::Match, 90), (Kind::SoftClip, 10)];
        assert_eq!(unclipped_5prime(1000, &ops, true), 1100);
    }

    #[test]
    fn reverse_hard_and_soft_clip() {
        let ops = vec![(Kind::Match, 85), (Kind::SoftClip, 10), (Kind::HardClip, 5)];
        assert_eq!(unclipped_5prime(1000, &ops, true), 1100);
    }

    #[test]
    fn reverse_spliced_with_clip() {
        let ops = vec![
            (Kind::Match, 50),
            (Kind::Skip, 1000),
            (Kind::Match, 50),
            (Kind::SoftClip, 5),
        ];
        assert_eq!(unclipped_5prime(1000, &ops, true), 2105);
    }

    #[test]
    fn empty_cigar() {
        assert_eq!(unclipped_5prime(500, &[], false), 500);
    }

    #[test]
    fn forward_hard_clip_only() {
        let ops = vec![(Kind::HardClip, 10), (Kind::Match, 90)];
        assert_eq!(unclipped_5prime(1000, &ops, false), 990);
    }

    #[test]
    fn reverse_no_right_clip() {
        let ops = vec![(Kind::Match, 100)];
        assert_eq!(unclipped_5prime(1000, &ops, true), 1100);
    }

    #[test]
    fn forward_with_insertion() {
        let ops = vec![(Kind::Match, 50), (Kind::Insertion, 2), (Kind::Match, 48)];
        assert_eq!(unclipped_5prime(1000, &ops, false), 1000);
    }

    #[test]
    fn reverse_with_deletion() {
        let ops = vec![(Kind::Match, 50), (Kind::Deletion, 5), (Kind::Match, 50)];
        assert_eq!(unclipped_5prime(1000, &ops, true), 1105);
    }
}
