//! Picard-format DuplicationMetrics output for MultiQC compatibility.

use anyhow::{Context, Result};
use std::io::Write;

/// Counters accumulated during Pass 1.
pub struct MetricsCounters {
    pub unpaired_reads_examined: u64,
    pub read_pairs_examined: u64,
    pub secondary_or_supplementary: u64,
    pub unmapped_reads: u64,
    pub unpaired_read_duplicates: u64,
    pub read_pair_duplicates: u64,
    pub group_sizes: Vec<u64>, // histogram: index = group_size, value = count
}

impl MetricsCounters {
    pub fn new() -> Self {
        Self {
            unpaired_reads_examined: 0,
            read_pairs_examined: 0,
            secondary_or_supplementary: 0,
            unmapped_reads: 0,
            unpaired_read_duplicates: 0,
            read_pair_duplicates: 0,
            group_sizes: Vec::new(),
        }
    }

    pub fn merge_histogram(&mut self, other: &[u64]) {
        if other.len() > self.group_sizes.len() {
            self.group_sizes.resize(other.len(), 0);
        }
        for (i, &count) in other.iter().enumerate() {
            self.group_sizes[i] += count;
        }
    }
}

/// Compute PERCENT_DUPLICATION per Picard formula.
pub fn percent_duplication(counters: &MetricsCounters) -> f64 {
    let numerator =
        counters.unpaired_read_duplicates as f64 + counters.read_pair_duplicates as f64 * 2.0;
    let denominator =
        counters.unpaired_reads_examined as f64 + counters.read_pairs_examined as f64 * 2.0;
    if denominator == 0.0 {
        0.0
    } else {
        numerator / denominator
    }
}

/// Estimate library size using Lander-Waterman equation.
/// Solve c/x = 1 - exp(-n/x) for x via Newton's method.
/// n = read_pairs_examined, c = unique pairs = n - duplicates.
pub fn estimate_library_size(pairs_examined: u64, pair_duplicates: u64) -> u64 {
    if pairs_examined == 0 {
        return 0;
    }
    let n = pairs_examined as f64;
    let c = (pairs_examined - pair_duplicates) as f64;

    if c == 0.0 {
        return 0;
    }
    if c == n {
        return n as u64; // No duplicates
    }

    // Newton's method: solve f(x) = c/x - 1 + exp(-n/x) = 0
    // f'(x) = -c/x^2 + n/x^2 * exp(-n/x)
    let mut x = c * c / (c - (n - c).max(1.0)); // Initial guess
    if x <= 0.0 {
        x = n; // Fallback
    }

    for _ in 0..100 {
        let exp_term = (-n / x).exp();
        let f = c / x - 1.0 + exp_term;
        let f_prime = -c / (x * x) + n / (x * x) * exp_term;

        if f_prime.abs() < 1e-20 {
            break;
        }

        let x_new = x - f / f_prime;
        if (x_new - x).abs() < 1.0 {
            x = x_new;
            break;
        }
        x = x_new;

        if x <= 0.0 {
            x = n;
            break;
        }
    }

    x.round().max(0.0) as u64
}

/// Write Picard-format metrics file.
pub fn write_metrics(
    path: &str,
    counters: &MetricsCounters,
    input_path: &str,
    output_path: &str,
) -> Result<()> {
    let mut f =
        std::fs::File::create(path).with_context(|| format!("Failed to create metrics: {}", path))?;

    let version = env!("CARGO_PKG_VERSION");
    let pct = percent_duplication(counters);
    let est = estimate_library_size(counters.read_pairs_examined, counters.read_pair_duplicates);

    // Header
    writeln!(f, "## htsjdk.samtools.metrics.StringHeader")?;
    writeln!(
        f,
        "# markdup-wea {} INPUT={} OUTPUT={}",
        version, input_path, output_path
    )?;

    // Metrics class line (TAB-separated — critical for MultiQC)
    writeln!(
        f,
        "## METRICS CLASS\tpicard.sam.markduplicates.DuplicationMetrics"
    )?;

    // Column headers
    writeln!(
        f,
        "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t\
         SECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\t\
         UNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t\
         READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\t\
         ESTIMATED_LIBRARY_SIZE"
    )?;

    // Data row (single library for now)
    writeln!(
        f,
        "default\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{:.6}\t{}",
        counters.unpaired_reads_examined,
        counters.read_pairs_examined,
        counters.secondary_or_supplementary,
        counters.unmapped_reads,
        counters.unpaired_read_duplicates,
        counters.read_pair_duplicates,
        pct,
        est
    )?;

    // Blank line before histogram
    writeln!(f)?;

    // Histogram
    writeln!(f, "## HISTOGRAM\tjava.lang.Double")?;
    writeln!(f, "BIN\tVALUE")?;
    for (i, &count) in counters.group_sizes.iter().enumerate() {
        if i >= 1 && count > 0 {
            writeln!(f, "{:.1}\t{}", i as f64, count)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percent_dup_normal() {
        let c = MetricsCounters {
            unpaired_reads_examined: 0,
            read_pairs_examined: 90481,
            secondary_or_supplementary: 0,
            unmapped_reads: 0,
            unpaired_read_duplicates: 0,
            read_pair_duplicates: 11688,
            group_sizes: Vec::new(),
        };
        let pct = percent_duplication(&c);
        assert!((pct - 0.129178).abs() < 0.001);
    }

    #[test]
    fn percent_dup_zero() {
        let c = MetricsCounters {
            unpaired_reads_examined: 0,
            read_pairs_examined: 50000,
            secondary_or_supplementary: 0,
            unmapped_reads: 0,
            unpaired_read_duplicates: 0,
            read_pair_duplicates: 0,
            group_sizes: Vec::new(),
        };
        assert_eq!(percent_duplication(&c), 0.0);
    }

    #[test]
    fn percent_dup_empty() {
        let c = MetricsCounters::new();
        assert_eq!(percent_duplication(&c), 0.0);
    }

    #[test]
    fn percent_dup_all_dups() {
        let c = MetricsCounters {
            unpaired_reads_examined: 0,
            read_pairs_examined: 100,
            secondary_or_supplementary: 0,
            unmapped_reads: 0,
            unpaired_read_duplicates: 0,
            read_pair_duplicates: 99,
            group_sizes: Vec::new(),
        };
        let pct = percent_duplication(&c);
        assert!((pct - 0.99).abs() < 0.001);
    }

    #[test]
    fn library_size_no_dups() {
        assert_eq!(estimate_library_size(50000, 0), 50000);
    }

    #[test]
    fn library_size_zero() {
        assert_eq!(estimate_library_size(0, 0), 0);
    }

    #[test]
    fn library_size_all_dups() {
        assert_eq!(estimate_library_size(100, 100), 0);
    }

    #[test]
    fn library_size_normal() {
        // With 90481 pairs and 11688 dups, library size should be a large positive number
        let est = estimate_library_size(90481, 11688);
        assert!(est > 0);
        assert!(est > 90481); // Library should be larger than observed pairs
    }

    #[test]
    fn metrics_file_format() {
        let c = MetricsCounters {
            unpaired_reads_examined: 100,
            read_pairs_examined: 5000,
            secondary_or_supplementary: 50,
            unmapped_reads: 20,
            unpaired_read_duplicates: 10,
            read_pair_duplicates: 500,
            group_sizes: vec![0, 4000, 400, 80, 15, 5],
        };

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("metrics.txt");
        let path_str = path.to_str().unwrap();

        write_metrics(path_str, &c, "input.bam", "output.bam").unwrap();

        let content = std::fs::read_to_string(&path).unwrap();

        // Check critical MultiQC header
        assert!(content.contains("## METRICS CLASS\tpicard.sam.markduplicates.DuplicationMetrics"));
        // Check column headers
        assert!(content.contains("LIBRARY\tUNPAIRED_READS_EXAMINED"));
        // Check histogram
        assert!(content.contains("## HISTOGRAM\tjava.lang.Double"));
        assert!(content.contains("BIN\tVALUE"));
        // Check data row starts with library name
        assert!(content.contains("default\t100\t5000"));
        // Optical dups = 0
        assert!(content.contains("\t0\t"));
    }
}
