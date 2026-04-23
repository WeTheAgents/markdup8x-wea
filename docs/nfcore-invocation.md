# nf-core/rnaseq PICARD_MARKDUPLICATES — invocation and nf-test assertions

Source: `nf-core/modules` path `modules/nf-core/picard/markduplicates/` on `main`, inspected 2026-04-19 from a fresh clone.

## Module command

`modules/nf-core/picard/markduplicates/main.nf` runs Picard as:

```
picard \
    -Xmx${avail_mem}M \
    MarkDuplicates \
    ${args} \
    --INPUT ${reads} \
    --OUTPUT ${prefix}.${suffix} \
    ${reference} \
    --METRICS_FILE ${prefix}.metrics.txt
```

- CLI style is **GNU `--KEY VALUE`** (space-separated), not Picard-legacy `K=V`.
- `${args}` is `task.ext.args`, empty by default; configs can inject extra flags.
- `${reference}` is `--REFERENCE_SEQUENCE ${fasta}` when a fasta is supplied (CRAM path), otherwise empty.
- `-Xmx${avail_mem}M` is a JVM heap flag — our shim must silently swallow `-Xmx*`.

Version probe (from the same module):

```
picard MarkDuplicates --version 2>&1 | sed -n 's/.*Version://p'
```

Expected captured value: `3.4.0`. Our shim must emit a line containing `Version:3.4.0` on `--version`.

## nf-test assertions (real tests)

Three non-stub tests in `tests/main.nf.test`. Each asserts:

1. `process.success` — exit 0.
2. `snapshot(file(process.out.bam[0][1]).name, path(process.out.metrics.get(0).get(1)).readLines()[0..2], process.out.findAll { key, val -> key.startsWith("versions") }).match()`

That snapshot compares three things:

- **BAM/CRAM output filename** (string only, not md5). Must be `test.md.bam` (or `test.md.cram`).
- **First 3 lines of the metrics file** (exact string match against `tests/main.nf.test.snap`).
- **versions block** — `versions_picard = ["PICARD_MARKDUPLICATES", "picard", "3.4.0"]`.

**BAM content md5 is NOT checked.** Our per-record `PG:Z:MarkDuplicates` tag (which worried us) is irrelevant for nf-test. Big win.

### Metrics file first 3 lines (from snapshot)

Line 1 and line 3 are identical fixed strings:
```
## htsjdk.samtools.metrics.StringHeader
```

Line 2 is a **single long string** containing the Picard-style command line with **all default flags expanded**. Example (sarscov2 unsorted):

```
# MarkDuplicates --INPUT test.paired_end.bam --OUTPUT test.md.bam --METRICS_FILE test.md.metrics.txt --ASSUME_SORT_ORDER queryname --MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP 50000 --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 --SORTING_COLLECTION_SIZE_RATIO 0.25 --TAG_DUPLICATE_SET_MEMBERS false --REMOVE_SEQUENCING_DUPLICATES false --TAGGING_POLICY DontTag --CLEAR_DT true --DUPLEX_UMI false --FLOW_MODE false --FLOW_DUP_STRATEGY FLOW_QUALITY_SUM_STRATEGY --FLOW_USE_END_IN_UNPAIRED_READS false --FLOW_USE_UNPAIRED_CLIPPED_END false --FLOW_UNPAIRED_END_UNCERTAINTY 0 --FLOW_UNPAIRED_START_UNCERTAINTY 0 --FLOW_SKIP_FIRST_N_FLOWS 0 --FLOW_Q_IS_KNOWN_END false --FLOW_EFFECTIVE_QUALITY_THRESHOLD 15 --ADD_PG_TAG_TO_READS true --REMOVE_DUPLICATES false --ASSUME_SORTED false --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --PROGRAM_RECORD_ID MarkDuplicates --PROGRAM_GROUP_NAME MarkDuplicates --READ_NAME_REGEX <optimized capture of last three ':' separated fields as numeric values> --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 5 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false
```

The CRAM case additionally has `--REFERENCE_SEQUENCE genome.fasta` inserted after `--METRICS_FILE test.md.metrics.txt` and before `--MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP`.

**Consequence:** our metrics emitter cannot produce this line; Rust reimplementation doesn't know about `FLOW_UNPAIRED_END_UNCERTAINTY` etc. The cleanest fix is **post-processing in the `picard` wrapper**: after `markdup-wea` writes metrics, the wrapper rewrites the top 3 lines with a Picard-templated header derived from the observed CLI args, injecting the fixed default-flags suffix.

This is a shim concern, not a binary concern. The core `markdup-wea` stays clean.

## nf-test cases — success prediction

| Test | Input | Likely outcome | Reason |
|---|---|---|---|
| `sarscov2 [unsorted bam]` — real | queryname-sorted BAM | **FAIL** | binary rejects non-coordinate input (docs/deviations.md §2). Wrapper cannot fix without samtools or binary change. |
| `sarscov2 [sorted bam]` — real | coordinate-sorted BAM | **likely PASS** if wrapper produces correct metrics header | core parity path. |
| `homo_sapiens [cram]` — real | CRAM | **FAIL** | binary has no CRAM reader/writer. |
| three `- stub` tests | any | **PASS** trivially | touch empty files, md5s match empty-file md5 `d41d8cd98f00b204e9800998ecf8427e`. |

Realistic Track B outcome: **1/3 real tests pass, 3/3 stubs pass → 4/6 overall**, with two categorized failures (CRAM, queryname). Both failures map to known missing features already documented in `deviations.md` — not bugs. They are useful data points for the Phil/Jonathan response: "drop-in works for the coordinate-sorted BAM path that nf-core/rnaseq actually exercises (since STAR/HISAT output is coordinate-sorted via `samtools sort` upstream); CRAM and queryname-input paths would need additional engineering."

## Plan revisions driven by these findings

1. Wrapper parses `--KEY VALUE` (space-separated), not `K=V`. Simpler.
2. Wrapper post-processes metrics file: replace the top `## htsjdk.samtools.metrics.StringHeader` header + `# MarkDuplicates ...` line with a Picard-3.4.0-templated header. Template is the literal string above with `--INPUT`, `--OUTPUT`, `--METRICS_FILE`, optional `--REFERENCE_SEQUENCE` interpolated.
3. Wrapper `--version` emits `Version:3.4.0` (matching the sed capture).
4. Accept CRAM + unsorted-BAM failures as known; do not attempt to fix in this track.
5. BAM md5 irrelevant for nf-test → the `PG:Z per-record tag` concern from the original plan is **dropped**. No `--no-pg-per-record` flag needed.
