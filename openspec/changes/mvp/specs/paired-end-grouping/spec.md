# Delta for paired-end-grouping

## ADDED Requirements

### Requirement: Paired-end grouping key
The system SHALL group paired-end reads by `PairedEndKey = (library_idx, ref_id_lo, pos_lo, is_reverse_lo, ref_id_hi, pos_hi, is_reverse_hi)` where positions are unclipped 5-prime coordinates and "lo"/"hi" are ordered by (ref_id, position) ascending.

#### Scenario: Two duplicate pairs at same position
- **GIVEN** pair A with read1 at chr1:1000(+) and read2 at chr1:1500(-)
- **AND** pair B with read1 at chr1:1000(+) and read2 at chr1:1500(-)
- **WHEN** the system computes grouping keys
- **THEN** both pairs produce the identical PairedEndKey
- **AND** they are placed in the same duplicate group

#### Scenario: Same position, different orientation
- **GIVEN** pair A with read1 at chr1:1000(+) and read2 at chr1:1500(-)
- **AND** pair B with read1 at chr1:1000(-) and read2 at chr1:1500(+)
- **WHEN** the system computes grouping keys
- **THEN** the pairs produce different PairedEndKeys (orientation differs)
- **AND** they are NOT grouped as duplicates

#### Scenario: Encounter order independence
- **GIVEN** pair A where read2 appears before read1 in the BAM stream (read2 has lower coordinate)
- **WHEN** the system processes both reads
- **THEN** the PairedEndKey is identical regardless of encounter order
- **AND** "lo" always refers to the mate with the lower (ref_id, position)

### Requirement: Unclipped 5-prime position calculation
The system SHALL compute unclipped 5-prime positions accounting for soft and hard clips in the CIGAR string.

#### Scenario: Forward read with left soft clip
- **GIVEN** a forward-strand read at aligned position 1000 with CIGAR "5S95M"
- **WHEN** the system computes unclipped_5prime
- **THEN** the result is 995 (1000 - 5)

#### Scenario: Forward read with left hard clip
- **GIVEN** a forward-strand read at aligned position 1000 with CIGAR "3H5S92M"
- **WHEN** the system computes unclipped_5prime
- **THEN** the result is 992 (1000 - 5 - 3)

#### Scenario: Reverse read with right soft clip
- **GIVEN** a reverse-strand read at aligned position 1000 with CIGAR "90M10S"
- **WHEN** the system computes unclipped_5prime
- **THEN** the result is 1100 (1000 + 90 + 10, which is alignment_end + right_clips)

#### Scenario: Reverse read with right hard and soft clip
- **GIVEN** a reverse-strand read at aligned position 1000 with CIGAR "85M10S5H"
- **WHEN** the system computes unclipped_5prime
- **THEN** the result is 1100 (1000 + 85 + 10 + 5)

#### Scenario: Read with no clipping
- **GIVEN** a forward-strand read at position 1000 with CIGAR "100M"
- **WHEN** the system computes unclipped_5prime
- **THEN** the result is 1000

#### Scenario: Reverse read with spliced alignment
- **GIVEN** a reverse-strand read with CIGAR "50M1000N50M5S" at position 1000
- **WHEN** the system computes unclipped_5prime
- **THEN** the result is 2105 (1000 + 50 + 1000 + 50 + 5, end_pos + right_clips)

### Requirement: Mate tracking via QNAME hash
The system SHALL track unpaired mates using a FxHashMap keyed by 64-bit hash of the QNAME string.

#### Scenario: Normal mate resolution
- **GIVEN** read1 of pair "HISEQ:1:1:1000:2000" at chr1:500
- **AND** read2 of the same pair at chr1:800
- **WHEN** read1 is encountered first
- **THEN** read1 is stored in pending_mates with its metadata
- **AND** when read2 is encountered, read1 is found by QNAME hash
- **AND** the PairedEndKey is formed using both reads' data

#### Scenario: Hash collision detection
- **GIVEN** read A with QNAME "READ_A" and read B with QNAME "READ_B"
- **AND** FxHash("READ_A") == FxHash("READ_B") (collision)
- **WHEN** read B looks up its mate in pending_mates
- **THEN** the system detects the collision via check_hash mismatch (first 4 bytes of QNAME)
- **AND** read B is treated as a new pending mate, not paired with read A

#### Scenario: Orphan reads (mate never found)
- **GIVEN** a paired read whose mate is not present in the BAM
- **WHEN** the scan pass reaches EOF
- **THEN** the orphan read is logged as a warning to stderr
- **AND** the read passes through to output without FLAG 0x400
- **AND** processing does not abort

### Requirement: Cross-chromosome pair handling
The system SHALL detect duplicate pairs where mates map to different chromosomes.

#### Scenario: Inter-chromosomal duplicate pair
- **GIVEN** pair A with read1 at chr1:1000(+) and read2 at chr5:2000(-)
- **AND** pair B with read1 at chr1:1000(+) and read2 at chr5:2000(-)
- **WHEN** the system processes the entire BAM
- **THEN** pairs A and B are grouped together
- **AND** the lower-scoring pair is flagged as duplicate on both reads

#### Scenario: Cross-chromosome resolution timing
- **GIVEN** read1 of pair at chr1:1000 and read2 at chr22:500
- **WHEN** the stream reaches chr22
- **THEN** read2 resolves the pending mate from chr1
- **AND** the group containing this pair can be evaluated for resolution

### Requirement: Library-aware grouping
The system SHALL group reads by library identity derived from the @RG header.

#### Scenario: Multiple libraries
- **GIVEN** a BAM with @RG headers defining LB:lib1 and LB:lib2
- **AND** pair A from lib1 at chr1:1000(+)/chr1:1500(-)
- **AND** pair B from lib2 at chr1:1000(+)/chr1:1500(-)
- **WHEN** the system computes grouping keys
- **THEN** pairs A and B have different PairedEndKeys (library_idx differs)
- **AND** they are NOT grouped as duplicates

#### Scenario: No @RG header
- **GIVEN** a BAM with no @RG header lines
- **WHEN** the system processes reads
- **THEN** all reads are assigned to a single default library
- **AND** duplicate detection proceeds normally

### Requirement: Incremental group resolution
The system SHALL resolve duplicate groups incrementally as the coordinate stream advances, freeing memory for resolved groups.

#### Scenario: Group resolved when stream passes max position
- **GIVEN** a group with PairedEndKey where pos_hi = 5000
- **WHEN** the BAM stream advances to position 5001 on the same chromosome
- **THEN** the group is resolved (best pair selected, others flagged)
- **AND** the group's memory is freed

#### Scenario: Chromosome boundary resolution
- **GIVEN** pending groups on chr1
- **WHEN** the first record on chr2 is encountered
- **THEN** all groups for chr1 are resolved before processing chr2
