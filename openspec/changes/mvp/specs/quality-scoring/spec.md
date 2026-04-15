# Delta for quality-scoring

## ADDED Requirements

### Requirement: Base quality sum scoring
The system SHALL score each read by summing base qualities that are >= 15 (Phred scale), matching Picard's SUM_OF_BASE_QUALITIES scoring strategy (the default).

#### Scenario: All high-quality bases
- **GIVEN** a read with base qualities [30, 25, 20, 35, 40]
- **WHEN** the system computes quality_sum
- **THEN** the result is 150 (all bases >= 15, sum all)

#### Scenario: Mixed quality bases
- **GIVEN** a read with base qualities [30, 10, 20, 5, 40]
- **WHEN** the system computes quality_sum
- **THEN** the result is 90 (only 30 + 20 + 40; bases with Q < 15 excluded)

#### Scenario: All low-quality bases
- **GIVEN** a read with base qualities [10, 8, 12, 14, 5]
- **WHEN** the system computes quality_sum
- **THEN** the result is 0 (no bases >= 15)

#### Scenario: Empty quality string
- **GIVEN** a read with no base quality data (qual = "*" in SAM)
- **WHEN** the system computes quality_sum
- **THEN** the result is 0

### Requirement: Paired-end combined scoring
The system SHALL score a read pair by summing the quality_sum of both mates.

#### Scenario: Paired scoring
- **GIVEN** read1 with quality_sum 500 and read2 with quality_sum 300
- **WHEN** the system computes the pair's combined_score
- **THEN** the result is 800

### Requirement: Duplicate group resolution by score
The system SHALL select the highest-scoring read (or pair) in each duplicate group as the primary, and flag all others as duplicates.

#### Scenario: Best pair wins
- **GIVEN** a group with pair A (score 800) and pair B (score 600)
- **WHEN** the group is resolved
- **THEN** pair A is primary (no FLAG 0x400 on either read)
- **AND** pair B is flagged (FLAG | 0x400 on both reads)

#### Scenario: Tie-breaking
- **GIVEN** a group with pair A (score 800) and pair B (score 800)
- **WHEN** the group is resolved
- **THEN** exactly one pair is primary and the other is flagged
- **AND** the pair with the lower record_id (first encountered in BAM) wins, matching Picard's index-based tie-breaking

#### Scenario: Group of one (no duplicates)
- **GIVEN** a group containing exactly one pair
- **WHEN** the group is resolved
- **THEN** the pair is primary (no FLAG 0x400)
- **AND** no duplicate is marked
