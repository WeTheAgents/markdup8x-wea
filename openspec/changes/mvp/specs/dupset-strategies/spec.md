# Delta for dupset-strategies

## ADDED Requirements

### Requirement: DupSet trait abstraction
The system SHALL implement a `DupSet` trait with two backing implementations (BitVec and FxHashSet) to enable benchmarking before committing to a strategy.

#### Scenario: BitVec mode
- **GIVEN** the system is configured to use BitVec dupset
- **WHEN** processing a BAM with 200M records and 20M duplicates
- **THEN** the dupset uses approximately 25MB of memory (1 bit per record)
- **AND** lookup by record_id is O(1)

#### Scenario: FxHashSet mode
- **GIVEN** the system is configured to use FxHashSet dupset
- **WHEN** processing a BAM with 200M records and 20M duplicates
- **THEN** the dupset stores only duplicate record entries (~230MB)
- **AND** insertion and lookup are O(1) amortized

### Requirement: Identical results regardless of strategy
The system SHALL produce byte-identical output BAM files regardless of which DupSet implementation is used.

#### Scenario: Cross-strategy validation
- **GIVEN** the same input BAM
- **WHEN** processed with BitVec mode and then with FxHashSet mode
- **THEN** both outputs have identical FLAG values on every record
- **AND** both outputs have identical metrics (if -M is used)

### Requirement: Strategy selection
The system SHALL select the DupSet strategy via compile-time feature flag or runtime flag.

#### Scenario: Default strategy
- **WHEN** no strategy is explicitly selected
- **THEN** the system uses the default strategy (determined after yeast benchmarking)
