# Delta for duplicate-detection

## ADDED Requirements

### Requirement: Two-pass duplicate marking on coordinate-sorted BAM
The system SHALL read a coordinate-sorted BAM file in two passes: a scan pass to identify duplicates, and a write pass to apply duplicate flags and produce output.

#### Scenario: Standard paired-end RNA-seq BAM
- **GIVEN** a coordinate-sorted BAM with 180,342 paired-end reads (yeast WT_REP1)
- **WHEN** the user runs `markdup-wea input.bam -o output.bam`
- **THEN** the output BAM contains exactly 36,810 reads with FLAG bit 0x400 set
- **AND** the output BAM preserves all original records (same count, same order)
- **AND** all SAM tags (MD, NH, HI, AS, nM, etc.) are preserved unchanged

#### Scenario: BAM with mixed read types
- **GIVEN** a coordinate-sorted BAM containing paired-end, single-end, unmapped, secondary, and supplementary reads
- **WHEN** the system processes the BAM
- **THEN** only primary mapped reads participate in duplicate grouping
- **AND** secondary reads (FLAG & 0x100) pass through with original flags
- **AND** supplementary reads (FLAG & 0x800) pass through with original flags
- **AND** unmapped reads (FLAG & 0x4) pass through with original flags

### Requirement: Record classification
The system SHALL classify each BAM record into exactly one category for processing.

#### Scenario: Classification routing
- **WHEN** a record has FLAG & 0x4 (unmapped)
- **THEN** it is classified as "unmapped" and skipped for grouping
- **WHEN** a record has FLAG & 0x100 (secondary)
- **THEN** it is classified as "secondary" and skipped for grouping
- **WHEN** a record has FLAG & 0x800 (supplementary)
- **THEN** it is classified as "supplementary" and skipped for grouping
- **WHEN** a record has FLAG & 0x1 (paired) and FLAG & 0x8 (mate unmapped)
- **THEN** it is classified as "paired-mate-unmapped" and routed to single-end grouping
- **WHEN** a record has FLAG & 0x1 (paired) and !(FLAG & 0x4) and !(FLAG & 0x8)
- **THEN** it is classified as "paired-both-mapped" and routed to paired-end grouping
- **WHEN** a record has !(FLAG & 0x1) (unpaired) and !(FLAG & 0x4)
- **THEN** it is classified as "single-end" and routed to single-end grouping

### Requirement: Duplicate flag semantics
The system SHALL mark duplicates by setting SAM FLAG bit 0x400 (PCR or optical duplicate) on the BAM record.

#### Scenario: Flag bit is additive
- **GIVEN** a read with existing FLAG = 0x63 (paired, mapped, mate mapped, first in pair)
- **WHEN** the read is identified as a duplicate
- **THEN** the output FLAG is 0x463 (original | 0x400)
- **AND** no other FLAG bits are modified

#### Scenario: Already-flagged duplicates
- **GIVEN** a read with FLAG bit 0x400 already set in the input
- **WHEN** the system processes the read
- **THEN** the existing duplicate flag is cleared before grouping
- **AND** the flag is re-set only if the read is re-identified as a duplicate

### Requirement: Remove duplicates mode
The system SHALL support a `--remove-duplicates` flag that excludes duplicate reads from the output entirely.

#### Scenario: Remove duplicates
- **GIVEN** a BAM with 100 reads, 20 of which are duplicates
- **WHEN** the user runs with `--remove-duplicates`
- **THEN** the output BAM contains exactly 80 reads
- **AND** no read in the output has FLAG bit 0x400 set

### Requirement: Supplementary reads are NOT flagged in MVP
The system SHALL NOT set FLAG 0x400 on supplementary or secondary alignments, even if their primary alignment is marked as a duplicate.

#### Scenario: Supplementary of a duplicate primary
- **GIVEN** a primary read at chr1:1000 that is marked as duplicate
- **AND** a supplementary alignment of the same read at chr5:5000
- **WHEN** the system writes the output
- **THEN** the primary read has FLAG | 0x400
- **AND** the supplementary read retains its original flags (no 0x400 added)
