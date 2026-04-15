# Delta for picard-metrics

## ADDED Requirements

### Requirement: Picard-format DuplicationMetrics output
The system SHALL output a metrics file in the exact format expected by MultiQC's Picard module when the `-M` flag is provided.

#### Scenario: Metrics file header
- **WHEN** the user runs `markdup-wea input.bam -o out.bam -M metrics.txt`
- **THEN** `metrics.txt` starts with:
  ```
  ## htsjdk.samtools.metrics.StringHeader
  # markdup-wea {version} INPUT={input} OUTPUT={output}
  ## METRICS CLASS\tpicard.sam.markduplicates.DuplicationMetrics
  ```
- **AND** the third line uses a literal TAB character between "CLASS" and "picard..."

#### Scenario: Metrics data row
- **WHEN** the system writes the metrics data row
- **THEN** it contains TAB-separated columns in this exact order:
  LIBRARY, UNPAIRED_READS_EXAMINED, READ_PAIRS_EXAMINED, SECONDARY_OR_SUPPLEMENTARY_RDS, UNMAPPED_READS, UNPAIRED_READ_DUPLICATES, READ_PAIR_DUPLICATES, READ_PAIR_OPTICAL_DUPLICATES, PERCENT_DUPLICATION, ESTIMATED_LIBRARY_SIZE
- **AND** READ_PAIR_OPTICAL_DUPLICATES is 0 (optical detection not in MVP)

#### Scenario: MultiQC compatibility
- **GIVEN** a metrics file produced by markdup-wea
- **WHEN** the user runs `multiqc .` in the directory containing the file
- **THEN** MultiQC detects the file as Picard MarkDuplicates metrics
- **AND** generates a duplication rate plot with correct percentages

### Requirement: PERCENT_DUPLICATION calculation
The system SHALL compute PERCENT_DUPLICATION as `(unpaired_dups + read_pair_dups * 2) / (unpaired_examined + read_pairs_examined * 2)`.

#### Scenario: Typical RNA-seq duplication rate
- **GIVEN** 0 unpaired reads, 90,481 read pairs examined, 11,688 read pair duplicates
- **WHEN** the system computes PERCENT_DUPLICATION
- **THEN** the result is (0 + 11688 * 2) / (0 + 90481 * 2) = 0.129178 (approximately)

#### Scenario: Zero duplicates
- **GIVEN** 50,000 read pairs examined, 0 duplicates
- **WHEN** the system computes PERCENT_DUPLICATION
- **THEN** the result is 0.0

#### Scenario: All duplicates (degenerate)
- **GIVEN** 100 read pairs examined, 99 read pair duplicates (1 primary per group)
- **WHEN** the system computes PERCENT_DUPLICATION
- **THEN** the result is 198 / 200 = 0.99

### Requirement: ESTIMATED_LIBRARY_SIZE via Lander-Waterman
The system SHALL estimate library size using the Lander-Waterman equation: solve `c/x = 1 - exp(-n/x)` for x via Newton's method.

#### Scenario: Normal estimation
- **GIVEN** n = 90,481 read pairs examined, c = 78,793 unique pairs (n - duplicates)
- **WHEN** the system estimates library size
- **THEN** the result is a positive integer (the estimated number of unique molecules)
- **AND** the value is computed by Newton's method with convergence tolerance < 1.0

#### Scenario: No duplicates
- **GIVEN** n = 50,000 read pairs, c = 50,000 unique pairs
- **WHEN** the system estimates library size
- **THEN** the result equals n (trivial case, all unique)

#### Scenario: Zero examined
- **GIVEN** n = 0 read pairs examined
- **WHEN** the system estimates library size
- **THEN** the result is 0

### Requirement: Duplicate group size histogram
The system SHALL output a histogram of duplicate group sizes in the metrics file.

#### Scenario: Histogram format
- **WHEN** the metrics file is written
- **THEN** it contains a histogram section:
  ```
  ## HISTOGRAM\tjava.lang.Double
  BIN\tVALUE
  1.0\t{count of groups with 1 pair}
  2.0\t{count of groups with 2 pairs}
  ...
  ```
- **AND** BIN values are floating-point (1.0, 2.0, ...) for Picard compatibility

### Requirement: Per-library metrics
The system SHALL output one metrics row per library when multiple libraries are present.

#### Scenario: Multiple libraries
- **GIVEN** a BAM with reads from library "lib_A" and "lib_B"
- **WHEN** the metrics file is written
- **THEN** it contains two data rows, one for each library
- **AND** each row has independent counts and duplication rates

### Requirement: No metrics file when -M not specified
The system SHALL NOT create a metrics file unless the `-M` flag is provided.

#### Scenario: Metrics omitted
- **WHEN** the user runs `markdup-wea input.bam -o output.bam` (no -M flag)
- **THEN** no metrics file is created
- **AND** duplicate marking proceeds normally
