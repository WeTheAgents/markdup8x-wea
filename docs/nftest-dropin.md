# nf-test drop-in — `nf-core/modules` `picard/markduplicates`

Phil's ask from the Slack thread: *"package your rewrite in a docker image… run the nf-test unit tests pass"*. This document records the result.

## TL;DR

**6/6 PASS** against `modules/nf-core/picard/markduplicates/tests/main.nf.test` (nf-core/modules `main` as of 2026-04-23). Three real tests (sorted BAM, unsorted BAM, CRAM) and three stub tests — all green, byte-for-byte snapshot match on metrics line 2.

## Image

- Tag: `markdup-wea:picard-shim-0.1.0` (local, also tagged `local/markdup-wea:picard-shim-0.1.0`)
- Layer digest: `sha256:ed5620584bed4ad56f9f856c47674962a563543ff36a930a90d35e37067d4bb1`
- Size: 155 MB (debian:stable-slim + `markdup-wea` + `samtools` + `gawk` + `procps` + the shim)
- No Java runtime.

Build: [Dockerfile](../Dockerfile) (multi-stage, `rust:1.85-slim` → `debian:stable-slim`).

Shim: [scripts/picard-shim.sh](../scripts/picard-shim.sh). Parses `picard [-Xmx<N>] MarkDuplicates --KEY VALUE …`, swallows JVM/Picard-only tuning knobs, translates the rest to `markdup-wea`, pre-converts CRAM→BAM and samtools-sorts any non-coordinate input, then rewrites metrics line 2 with Picard's default-expansion template.

## Test matrix

| # | Test | Type | Result |
|---|------|------|:------:|
| 1 | `sarscov2 [sorted bam]` | real | PASS |
| 2 | `sarscov2 [unsorted bam]` | real | PASS |
| 3 | `homo_sapiens [cram]` | real | PASS |
| 4 | `sarscov2 [sorted bam] - stub` | stub | PASS |
| 5 | `sarscov2 [unsorted bam] - stub` | stub | PASS |
| 6 | `homo_sapiens [cram] - stub` | stub | PASS |

Total wall: ~34 s (6 tests, docker-profile).

## What the snapshot checks

For each real test the nf-test asserts `snapshot(file(...).name, metrics.readLines()[0..2], versions)`, so parity is judged on:

1. Output BAM/CRAM basename.
2. Exactly the first three lines of the metrics file (line 1 + 2 + 3).
3. The `versions_picard` triple.

Line 2 is the load-bearing one — Picard default-prints its full expanded argument list. The shim templates this from the observed inputs; [scripts/\_shim\_verify.sh](../scripts/_shim_verify.sh) compares the template against the checked-in snap for all three tests (3/3 byte-identical).

## Features validated by the run

| Shim feature | Evidence | Where |
|---|---|---|
| `--ASSUME_SORT_ORDER queryname` is accepted as a *hint* and swallowed (wea reads `@HD SO:` from the BAM header) | Tests 1–3 all pass with nf-core's `ext.args = '--ASSUME_SORT_ORDER queryname'` | [picard-shim.sh:131-138](../scripts/picard-shim.sh#L131-L138) |
| Auto-pre-sort for unsorted input (`samtools sort` if `@HD SO != coordinate`) | Test 2 passes on `test.paired_end.bam` (unsorted) | [picard-shim.sh:116-124](../scripts/picard-shim.sh#L116-L124) |
| CRAM → BAM pre-conversion (`samtools view -b --reference`) | Test 3 passes on `test.paired_end.sorted.cram` | [picard-shim.sh:106-114](../scripts/picard-shim.sh#L106-L114) |
| `--version` shape (`Version:3.4.0` on stderr, exit 0) — consumed by nf-core's Nextflow module | `versions_picard` snapshot in all 3 real tests | [picard-shim.sh:88-92](../scripts/picard-shim.sh#L88-L92) |
| Metrics line 2 template with `queryname` as the `ASSUME_SORT_ORDER` default-print value | Snapshot match on all three cases | [picard-shim.sh:162-168](../scripts/picard-shim.sh#L162-L168) |
| `procps` present for Nextflow task-metrics collection | Any `process.success` assertion passing | [Dockerfile:14](../Dockerfile#L14) |
| `ENTRYPOINT` **not** set (so `docker run <img> bash -c …` works for stub tasks) | Tests 4–6 pass | [Dockerfile](../Dockerfile) |

## Reproduce

On any host with docker + nextflow ≥ 24.10 + nf-test:

```bash
# 1. Build image
git clone https://github.com/WeTheAgents/markdup-wea
cd markdup-wea
docker build -t markdup-wea:picard-shim-0.1.0 .
docker tag markdup-wea:picard-shim-0.1.0 local/markdup-wea:picard-shim-0.1.0

# 2. Point the module's own config at the image
git clone --depth 1 https://github.com/nf-core/modules.git /tmp/nfcore-modules
cat > /tmp/nfcore-modules/modules/nf-core/picard/markduplicates/tests/nextflow.config <<'EOF'
docker.registry = ""
process {
    withName: PICARD_MARKDUPLICATES {
        ext.prefix = { "${meta.id}.md" }
        ext.args   = "--ASSUME_SORT_ORDER queryname"
        container  = "local/markdup-wea:picard-shim-0.1.0"
    }
}
EOF

# 3. Run
cd /tmp/nfcore-modules
nf-test test modules/nf-core/picard/markduplicates/tests/main.nf.test --profile docker
```

Expected: `SUCCESS: Executed 6 tests` in ~35 s.

## What this does **not** prove

- **Real-world workflows.** We ran the nf-core unit tests, not the full `nf-core/rnaseq` or any downstream pipeline. Real-pipeline validation remains TBD.
- **Non-default argument matrices.** The tests exercise a narrow arg shape (nf-core's `ext.args`). Other Picard flags (`--BARCODE_TAG`, `--REMOVE_DUPLICATES`, `--DUPLEX_UMI`, `--READ_ONE_BARCODE_TAG`, tag policies, etc.) have separate parity evidence in [docs/parity-evidence.md](parity-evidence.md); they are not exercised by this snap.
- **True queryname-sorted input.** None of the test fixtures are actually queryname-sorted — `ext.args` just passes the hint. On *really* queryname-sorted input wea still errors (deviation §2, [docs/deviations.md](deviations.md)); nf-core doesn't hit that path.
- **Image publication.** Image is local-only; pushing to a registry is deferred (plan Stage 3d).

## Known deviations still stand

- Aux-tag write order (deviation §10) — cosmetic, flag-parity unaffected.
- Queryname-sorted *input* (deviation §2) — the shim pre-sorts unsorted BAMs; a truly queryname-sorted BAM would still error at the `samtools sort` step or inside wea. This is a BAM state the nf-core test set does not produce.
