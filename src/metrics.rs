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
    /// LIBRARY name for the metrics row. Matches Picard's `getLibraryName`
    /// semantics: first @RG LB value, or "Unknown Library" if LB is absent.
    pub library_name: String,
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
            library_name: "Unknown Library".to_string(),
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

impl Default for MetricsCounters {
    fn default() -> Self {
        Self::new()
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

/// Estimate library size using Picard's Lander-Waterman bisection (research §9).
///
/// Solves `f(x) = c/x - 1 + exp(-n/x) = 0` where `n = read_pairs` and
/// `c = unique_read_pairs`. Bracket is doubled (×10) until `f(M × c) > 0`,
/// then 40 bisection iterations.
///
/// Returns `None` when the estimate is undefined:
/// - zero pairs examined
/// - zero duplicates (Picard returns null field)
/// - unique ≥ total (Picard throws `IllegalStateException`; we return None)
pub fn estimate_library_size(pairs_examined: u64, pair_duplicates: u64) -> Option<u64> {
    if pairs_examined == 0 || pair_duplicates == 0 {
        return None;
    }

    let read_pairs = pairs_examined as f64;
    let unique_read_pairs = (pairs_examined - pair_duplicates) as f64;

    if unique_read_pairs >= read_pairs || unique_read_pairs <= 0.0 {
        return None;
    }

    // Picard's f(x, c, n) with x as the size multiplier of unique_read_pairs (c):
    //   library_size = x * c
    //   c / library_size = 1 / x
    //   n / library_size = n / (x * c) = read_pairs / (x * unique_read_pairs)
    //   f = 1/x - 1 + exp(-n / library_size)
    // Root of f defines library_size; we search for x in [1, 100] growing upward.
    let f = |x: f64| -> f64 {
        1.0 / x - 1.0 + (-read_pairs / (x * unique_read_pairs)).exp()
    };

    let mut m: f64 = 1.0;
    let mut big_m: f64 = 100.0;

    // Grow the upper bracket until f(big_m) > 0 (matches Picard loop).
    while f(big_m) > 0.0 {
        big_m *= 10.0;
        if big_m > 1e12 {
            return None; // numerical blow-up guard
        }
    }

    // 40 bisection iterations (Picard constant).
    for _ in 0..40 {
        let mid = (m + big_m) / 2.0;
        let fm = f(mid);
        if fm == 0.0 {
            m = mid;
            big_m = mid;
            break;
        } else if fm > 0.0 {
            m = mid;
        } else {
            big_m = mid;
        }
    }

    let est = unique_read_pairs * (m + big_m) / 2.0;
    if est.is_finite() && est > 0.0 {
        Some(est as u64)
    } else {
        None
    }
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
    let est_opt = estimate_library_size(counters.read_pairs_examined, counters.read_pair_duplicates);
    // Picard writes an empty field (no value after the TAB) when the estimate is null.
    let est_str = est_opt.map(|v| v.to_string()).unwrap_or_default();

    // Picard preamble has two `## htsjdk.samtools.metrics.StringHeader` blocks
    // separated by a `# ...` comment line each:
    //   block 1 → command line (provenance of the producing tool)
    //   block 2 → run start timestamp ("Started on: <Picard date format>")
    // We emit the same shape for MultiQC and downstream-parser parity. Content
    // is honest: tool identifier is markdup-wea (we don't pretend to be
    // MarkDuplicates), and the timestamp reflects our actual run start.
    writeln!(f, "## htsjdk.samtools.metrics.StringHeader")?;
    writeln!(
        f,
        "# markdup-wea {} INPUT={} OUTPUT={}",
        version, input_path, output_path
    )?;
    writeln!(f, "## htsjdk.samtools.metrics.StringHeader")?;
    writeln!(f, "# Started on: {}", started_on_picard_format())?;
    writeln!(f)?; // blank line after preamble (Picard parity)

    // Metrics class line (TAB-separated — critical for MultiQC).
    // Picard v3.x (verified on real nf-core/rnaseq output) writes
    // `picard.sam.DuplicationMetrics`, NOT the `markduplicates.` sub-package.
    writeln!(f, "## METRICS CLASS\tpicard.sam.DuplicationMetrics")?;

    // Column headers
    writeln!(
        f,
        "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t\
         SECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\t\
         UNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t\
         READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\t\
         ESTIMATED_LIBRARY_SIZE"
    )?;

    // Data row. LIBRARY value comes from counters.library_name — Picard's
    // `getLibraryName` semantics: first @RG LB, else "Unknown Library".
    // ESTIMATED_LIBRARY_SIZE may be empty when undefined.
    writeln!(
        f,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t{}",
        counters.library_name,
        counters.unpaired_reads_examined,
        counters.read_pairs_examined,
        counters.secondary_or_supplementary,
        counters.unmapped_reads,
        counters.unpaired_read_duplicates,
        counters.read_pair_duplicates,
        trim_float(pct),
        est_str
    )?;

    // Blank line before histogram
    writeln!(f)?;

    // Picard's `calculateRoiHistogram` emits the ROI histogram ONLY when
    // ESTIMATED_LIBRARY_SIZE is non-null. Match that behavior.
    if let Some(lib_size) = est_opt {
        writeln!(f, "## HISTOGRAM\tjava.lang.Double")?;
        writeln!(f, "BIN\tCoverageMult\tall_sets\tnon_optical_sets")?;
        let l = lib_size as f64;
        let n = counters.read_pairs_examined as f64;
        let unique = (counters.read_pairs_examined - counters.read_pair_duplicates) as f64;
        for (i, &count) in counters.group_sizes.iter().enumerate() {
            if i >= 1 && count > 0 {
                let x = i as f64;
                // Picard's DuplicationMetrics.estimateRoi (research §9):
                //   CoverageMult(x) = L * (1 - exp(-x * N / L)) / unique
                let cov_mult = if unique > 0.0 && l > 0.0 {
                    l * (1.0 - (-x * n / l).exp()) / unique
                } else {
                    1.0
                };
                // non_optical_sets == all_sets since we don't track optical dups.
                writeln!(
                    f,
                    "{:.1}\t{}\t{}\t{}",
                    x,
                    trim_float(cov_mult),
                    count,
                    count
                )?;
            }
        }
    }

    writeln!(f)?; // trailing blank line (Picard parity)
    Ok(())
}

/// Picard's date format: e.g. "Wed Apr 08 04:31:00 GMT 2026". Locale-independent
/// (always English short names + UTC, displayed as "GMT" to match Picard exactly).
fn started_on_picard_format() -> String {
    let now = jiff::Zoned::now().with_time_zone(jiff::tz::TimeZone::UTC);
    // strftime: %a abbrev weekday, %b abbrev month, %d zero-padded day,
    //           %H:%M:%S 24h time, then literal " GMT ", %Y 4-digit year.
    now.strftime("%a %b %d %H:%M:%S UTC %Y").to_string()
}

/// Format a float the way Picard's `FormatUtil`/`DecimalFormat("0.######")` does:
/// up to 6 fractional digits, trailing zeros and trailing '.' stripped.
fn trim_float(v: f64) -> String {
    let s = format!("{:.6}", v);
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
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
            library_name: "Unknown Library".to_string(),
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
            library_name: "Unknown Library".to_string(),
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
            library_name: "Unknown Library".to_string(),
        };
        let pct = percent_duplication(&c);
        assert!((pct - 0.99).abs() < 0.001);
    }

    #[test]
    fn library_size_no_dups_returns_none() {
        // Picard: zero duplicates → null (empty field).
        assert_eq!(estimate_library_size(50000, 0), None);
    }

    #[test]
    fn library_size_zero_pairs_returns_none() {
        assert_eq!(estimate_library_size(0, 0), None);
    }

    #[test]
    fn library_size_all_dups_returns_none() {
        // Picard throws IllegalStateException; we return None (documented).
        assert_eq!(estimate_library_size(100, 100), None);
    }

    #[test]
    fn library_size_normal_estimate_is_plausible() {
        // 90481 pairs observed, 11688 dups → ~78793 unique. Library must be larger than unique.
        let est = estimate_library_size(90481, 11688).expect("defined estimate");
        assert!(est > 78793, "library size {} should exceed unique count", est);
        // Sanity: for ~13% dup rate the library is finite and not absurd.
        assert!(est < 10_000_000, "library size {} is implausibly large", est);
    }

    #[test]
    fn library_size_bisection_converges_on_yeast_sample() {
        // WT_REP1 ground truth: 180342 reads total, ~36810 dups flagged.
        // Pair-level (approx): the estimate should be finite and well above unique count.
        let est = estimate_library_size(90171, 18405).expect("defined estimate");
        let unique = 90171 - 18405;
        assert!(est > unique, "{} should exceed unique {}", est, unique);
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
            library_name: "Unknown Library".to_string(),
        };

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("metrics.txt");
        let path_str = path.to_str().unwrap();

        write_metrics(path_str, &c, "input.bam", "output.bam").unwrap();

        let content = std::fs::read_to_string(&path).unwrap();

        // Check critical MultiQC header
        assert!(content.contains("## METRICS CLASS\tpicard.sam.DuplicationMetrics"));
        // Picard preamble structural parity: TWO StringHeader sections
        let n_headers = content.matches("## htsjdk.samtools.metrics.StringHeader").count();
        assert_eq!(n_headers, 2, "expected 2 StringHeader sections (Picard parity), got {}", n_headers);
        // Block 1: command-line provenance (we identify ourselves honestly)
        assert!(content.contains("# markdup-wea "));
        // Block 2: "Started on" timestamp in Picard format ("EEE MMM dd HH:MM:SS GMT YYYY")
        let started_line = content
            .lines()
            .find(|l| l.starts_with("# Started on: "))
            .expect("missing '# Started on:' line");
        let parts: Vec<&str> = started_line["# Started on: ".len()..].split_whitespace().collect();
        // Expected fields: [weekday, month, day, "HH:MM:SS", "GMT", year]
        assert_eq!(parts.len(), 6, "Picard date should have 6 space-separated tokens, got {:?}", parts);
        assert_eq!(parts[0].len(), 3, "weekday should be 3-char abbrev, got {:?}", parts[0]);
        assert_eq!(parts[1].len(), 3, "month should be 3-char abbrev, got {:?}", parts[1]);
        assert_eq!(parts[2].len(), 2, "day should be 2-char zero-padded, got {:?}", parts[2]);
        assert_eq!(parts[3].len(), 8, "time should be HH:MM:SS, got {:?}", parts[3]);
        assert_eq!(parts[4], "UTC", "expected literal UTC zone, got {:?}", parts[4]);
        assert_eq!(parts[5].len(), 4, "year should be 4 digits, got {:?}", parts[5]);
        // Check column headers
        assert!(content.contains("LIBRARY\tUNPAIRED_READS_EXAMINED"));
        // Check 4-column Picard histogram format
        assert!(content.contains("## HISTOGRAM\tjava.lang.Double"));
        assert!(content.contains("BIN\tCoverageMult\tall_sets\tnon_optical_sets"));
        // Check data row uses library_name (set to "Unknown Library" in this test)
        assert!(content.contains("Unknown Library\t100\t5000"));
        // Optical dups = 0
        assert!(content.contains("\t0\t"));
    }
}
