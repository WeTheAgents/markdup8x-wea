# markdup-wea

A Rust rewrite of Picard MarkDuplicates.

On 8 ENCODE RNA-seq BAMs (1.55 B records, 934 M duplicates) markdup-wea
produces **byte-identical** duplicate flags and metrics to
Picard MarkDuplicates 3.4.0 — **3.26× faster**, **54× less RAM**
(325 MB vs 17.6 GB peak on K562_REP1).

```
java -jar picard.jar MarkDuplicates I=in.bam O=out.bam M=metrics.txt
markdup-wea                        INPUT=in.bam OUTPUT=out.bam METRICS_FILE=metrics.txt
```

Same arguments (Picard-style `KEY=VALUE`), drop-in for nf-core/rnaseq's
`PICARD_MARKDUPLICATES` step. See `docs/deviations.md` for the short list
of intentional differences (optical duplicates, `@PG CL` text).

## UMI-aware duplicate marking

markdup-wea supports barcode-aware deduplication via `--barcode-tag <TAG>`
(equivalent to Picard's `BARCODE_TAG`). Two reads at the same 5'-unclipped
coordinate with the same barcode tag value are duplicates; different
barcode → independent fragments. Typical preparation:

```
umi_tools extract → STAR align → fgbio CopyUmiFromReadName → markdup-wea --barcode-tag RX
```

Validated against Picard 3.4.0 on 3 independent public UMI datasets
spanning 2 organisms and 3 chemistries: human UHRR 6bp SCRB-seq
(GSE75823), mouse testis 8bp (PRJNA416930), and mouse microglia
QuantSeq 6bp+TATA-spacer (GSE134031). **Zero flag divergence on a
combined 255.6M records, exact metrics match on all three.**
See [`docs/parity-evidence.md`](docs/parity-evidence.md) for full
tables, commands, and the reproduction harness.

Not yet supported: `DUPLEX_UMI` (mate-unclipped start via MC tag, planned
A.7), multi-library BAMs with distinct RG-level barcodes, optical-duplicate
pixel-distance detection. See `docs/deviations.md §8 / §10`.

## Testing & validation — where to read

Start at the top and drill in as needed:

| Doc | What's in it |
|-----|--------------|
| [`docs/parity-evidence.md`](docs/parity-evidence.md) | **Parity evidence pack** — single sheet summarizing every parity artifact: 8/8 ENCODE non-UMI byte equality, UMI feature arm (2 independent datasets, flag-set equality + exact metrics), A.4 regression gate, deviations summary, and the explicit "NOT tested" list. Start here. |
| [`docs/session-report.md`](docs/session-report.md) | **Parity gate result** — 8/8 ENCODE samples byte-identical vs Picard 3.4.0. Zero flag mismatches, zero metrics data divergence. Per-sample record counts, duplicate counts, the 7 fixes that got us to parity. |
| [`docs/perf.md`](docs/perf.md) | **Performance measurements** — K562_REP1 threading scaling (`-@ 1/4/8`), full 8-sample batch wall+RSS head-to-head, and Phase E post-markdup parallelism (Qualimap 4-wide, bedtools 2-wide). |
| [`docs/validation.md`](docs/validation.md) | **How the parity gate works** — acceptance bar (decompressed BAM + metrics must be byte-identical), the comparison harness, corpus convention. |
| [`docs/deviations.md`](docs/deviations.md) | **Known differences from Picard** — every place we deliberately diverge, with rationale and downstream-impact assessment. Read before adopting. |
| [`docs/picard-tuning-proposal.md`](docs/picard-tuning-proposal.md) | **nf-core/rnaseq memory-tuning proposal** — independent of markdup-wea. Measured sweep (heap × parallelism) showing nf-core's default `process_medium` (36 GB, `-Xmx28g`) is ~4× over-allocated for Picard on RNA-seq. One-line fix: demote to `process_low`. 3.05× end-to-end throughput on a 30 GB box, 8/8 byte-identical to the default config. |
| [`docs/SPEC.md`](docs/SPEC.md) | Algorithm + on-disk format spec. |

## Harness & scripts

- [`scripts/validation/run_parity_batch.py`](scripts/validation/run_parity_batch.py) — runs Picard + markdup-wea on the same inputs, diffs outputs
- [`scripts/validation/run_local_gate.py`](scripts/validation/run_local_gate.py) — local pre-commit gate
- [`tests/validation_compare.rs`](tests/validation_compare.rs) — in-tree BAM/metrics comparison tests
- [`src/bin/validation_compare.rs`](src/bin/validation_compare.rs) — standalone comparator binary

## Building

```
cargo build --release
./target/release/markdup-wea INPUT=in.bam OUTPUT=out.bam METRICS_FILE=m.txt
```

## License

MIT. See `LICENSE`.
