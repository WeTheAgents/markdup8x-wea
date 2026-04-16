# Validation and Parity Gate

This repository treats **Picard MarkDuplicates 3.4.0** as the reference
implementation for the production pipeline we intend to replace.

The acceptance bar is deliberately strict:

- identical **decompressed BAM content**
- identical **metrics text**
- zero unexplained diffs on the agreed real-data corpus

Performance does not count as evidence. Speed only matters **after**
correctness has been demonstrated against the pinned Picard contract.

## What Is Pinned

The compatibility contract lives in a versioned TOML manifest:

- [scripts/validation/manifest.example.toml](/D:/GitHub/markdup-wea/scripts/validation/manifest.example.toml)

Before claiming parity, maintainers must replace placeholder values in a
checked-in manifest derived from the running Hetzner pipeline:

- exact Picard command line
- `picard.jar` path and SHA256
- Java version
- `samtools` version
- effective Picard defaults used by the pipeline
- corpus membership and input BAM SHA256 values
- expected Picard BAM/metrics outputs, either pre-generated or generated in-run

The initial corpus is the project’s existing proof target:

- 5 yeast samples: `WT_REP1`, `WT_REP2`, `RAP1_UNINDUCED_REP2`, `RAP1_UNINDUCED_REP1`, `RAP1_IAA_30M_REP1`
- 8 ENCODE samples: `K562_REP1`, `GM12878_REP1`, `GM12878_REP2`, `H1_REP1`, `H1_REP2`, `K562_REP2`, `MCF7_REP1`, `MCF7_REP2`

## Validation Commands

Local proof subset:

```powershell
python scripts/validation/run_local_gate.py
```

Local proof subset plus a real-data manifest:

```powershell
python scripts/validation/run_local_gate.py --manifest scripts/validation/hetzner.toml
```

Hetzner real-data parity batch:

```powershell
python scripts/validation/run_parity_batch.py --manifest scripts/validation/hetzner.toml
```

All commands emit machine-readable artifacts under `validation-results/`.

## Structured Comparator

The differential comparator is implemented as:

- [src/bin/validation_compare.rs](/D:/GitHub/markdup-wea/src/bin/validation_compare.rs)

It compares, in order:

- exact SAM header text
- record count and record order
- per-record core SAM fields
- duplicate flags
- auxiliary tags and values
- exact metrics text

The comparator stops at the **first mismatch** and writes JSON describing:

- sample ID
- mismatch class
- first mismatching record index, when applicable
- expected vs actual field values

## Release Gate

No rollout or “Picard replacement” claim is valid without:

1. `cargo test`
2. `cargo clippy --all-targets --all-features -- -D warnings`
3. green synthetic regression coverage for all in-scope parity blockers
4. green Hetzner parity batch against the pinned Picard 3.4.0 manifest

## Current Status

This repo now contains the **proof harness** and several blocker fixes, but it
is **not yet entitled to claim full Picard parity**.

Known high-importance gaps still requiring either implementation or explicit
proof of irrelevance to the pinned pipeline include:

- exact `@PG` / `PG:Z` parity
- multi-library metrics rows
- optical duplicate detection
- long-read score overflow parity if long reads are ever in scope
- any remaining metric-provenance text mismatches that prevent byte-identical output

Those gaps are acceptable only as open work, not as grounds for a parity claim.
