# Delta for sort-order-validation

## ADDED Requirements

### Requirement: Header-based sort order validation
The system SHALL validate that the input BAM header declares coordinate sort order via the @HD SO:coordinate field.

#### Scenario: Valid coordinate-sorted BAM
- **GIVEN** a BAM with @HD header containing SO:coordinate
- **WHEN** the system opens the file
- **THEN** processing proceeds normally

#### Scenario: Queryname-sorted BAM
- **GIVEN** a BAM with @HD SO:queryname
- **WHEN** the system opens the file
- **THEN** the system exits with a clear error: "Input BAM is sorted by queryname, not coordinate. markdup-wea requires coordinate-sorted input."
- **AND** exit code is non-zero

#### Scenario: Unsorted BAM
- **GIVEN** a BAM with @HD SO:unsorted or no SO field
- **WHEN** the system opens the file
- **THEN** the system exits with a clear error
- **AND** suggests: "Sort with: samtools sort input.bam -o sorted.bam"

#### Scenario: Override with --assume-sort-order
- **GIVEN** a BAM with @HD SO:unsorted
- **WHEN** the user runs with `--assume-sort-order coordinate`
- **THEN** the header check is bypassed
- **AND** runtime sort-order enforcement still applies

### Requirement: Runtime sort-order enforcement
The system SHALL verify coordinate ordering during Pass 1 by tracking the previous record's (ref_id, pos) and erroring if the current record violates ordering.

#### Scenario: Out-of-order records detected
- **GIVEN** a BAM where record at chr1:5000 is followed by chr1:3000
- **WHEN** the system encounters the out-of-order record
- **THEN** the system exits with error: "Input BAM is not coordinate-sorted: chr1:3000 follows chr1:5000"
- **AND** exit code is non-zero

#### Scenario: Chromosome boundary is valid
- **GIVEN** records: chr1:50000, chr2:100
- **WHEN** the system processes the chromosome boundary
- **THEN** processing continues normally (chr2:100 > chr1:anything by ref_id ordering)

### Requirement: Stdin input with temp file buffering
The system SHALL accept input from stdin (specified as `-`) by buffering to a temporary file for two-pass processing.

#### Scenario: Piped input
- **GIVEN** the user runs `samtools sort input.bam | markdup-wea - -o output.bam`
- **WHEN** the system detects stdin input
- **THEN** it creates a temporary file and copies stdin to it
- **AND** logs: "Reading from stdin; buffering to temp file"
- **AND** processes the temp file with the normal two-pass algorithm
- **AND** the temp file is automatically deleted on exit

#### Scenario: File input (normal path)
- **GIVEN** the user runs `markdup-wea input.bam -o output.bam`
- **WHEN** the system opens the input
- **THEN** no temp file is created
- **AND** the file is read directly in both passes
