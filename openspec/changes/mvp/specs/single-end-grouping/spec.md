# Delta for single-end-grouping

## ADDED Requirements

### Requirement: Single-end grouping key
The system SHALL group single-end reads (and paired reads with unmapped mates) by `SingleEndKey = (library_idx, ref_id, unclipped_5prime, is_reverse)`.

#### Scenario: Two single-end duplicates
- **GIVEN** read A (unpaired) at chr1:1000(+) with quality_sum 500
- **AND** read B (unpaired) at chr1:1000(+) with quality_sum 400
- **WHEN** the system processes both reads
- **THEN** they are grouped by the same SingleEndKey
- **AND** read A (higher score) is kept as primary
- **AND** read B is flagged with FLAG | 0x400

#### Scenario: Same position, different strand
- **GIVEN** read A at chr1:1000(+) and read B at chr1:1000(-)
- **WHEN** the system computes grouping keys
- **THEN** they have different SingleEndKeys
- **AND** they are NOT grouped as duplicates

### Requirement: Mate-unmapped reads treated as single-end
The system SHALL route paired reads with unmapped mates (FLAG & 0x8) to single-end grouping.

#### Scenario: Paired read with unmapped mate
- **GIVEN** a paired read (FLAG & 0x1) at chr1:1000(+) with mate unmapped (FLAG & 0x8)
- **WHEN** the system classifies the read
- **THEN** it is routed to single-end grouping using SingleEndKey
- **AND** it is NOT stored in pending_mates

### Requirement: Consecutive grouping for single-end
The system SHALL exploit coordinate sort order to resolve single-end groups inline, without buffering.

#### Scenario: Inline resolution
- **GIVEN** reads at positions chr1:1000, chr1:1000, chr1:1000, chr1:1005
- **WHEN** the third read at chr1:1000 is processed and the next read is at chr1:1005
- **THEN** the group at chr1:1000 is resolved immediately (3 reads, best score wins)
- **AND** no buffering beyond the current group is required
