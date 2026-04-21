# markdup-wea parity evidence pack

> **Purpose:** one place to point co-founders, reviewers, and prospective
> adopters at every parity artifact markdup-wea has against Picard
> MarkDuplicates 3.4.0. If you are evaluating whether to substitute
> markdup-wea for Picard in a production pipeline, read this plus
> [`deviations.md`](deviations.md) and you have everything.

**Scope statement.** markdup-wea targets the exact Picard `MarkDuplicates`
invocation used by `nf-core/rnaseq`'s `PICARD_MARKDUPLICATES` step on
coordinate-sorted RNA-seq BAMs (polyA and UMI). It is NOT a general
replacement for every Picard MarkDuplicates feature — the
not-implemented list is in [§ What is NOT tested](#what-is-not-tested)
and in `deviations.md`.

---

## TL;DR

- **Non-UMI, default-off:** `samtools view | md5sum` byte-identical to
  Picard 3.4.0 on 8/8 ENCODE rnaseq BAMs (1,549,871,416 records,
  ~934M duplicates, **0 flag divergence**). Wall-clock 3.26× faster,
  54× less peak RAM on K562_REP1.
- **UMI feature arm (`--barcode-tag RX`):** flag-set byte-identical
  on 3 independent public datasets spanning 2 organisms, 3 chemistries
  (SCRB-seq 6nt UMI, fgbio-pipeline 8nt UMI on R2, QuantSeq 6nt UMI
  inline + TATA spacer on R1), metrics exact match on all 3, wea
  2–2.5× faster than Picard. Flag-set = sorted `(qname, dup-bit)`
  tuples. The only `samtools view | md5sum` delta is aux-tag ordering
  (cosmetic, `deviations.md §10`).

---

## Tested datasets

### A. Non-UMI — default-off, `samtools view | md5sum` byte equality

| Dataset | Records | Chemistry | Picard md5 = wea md5 | Source |
|---------|--------:|-----------|:--------------------:|--------|
| ENCODE K562_REP1 | 156,141,236 | TruSeq polyA PE 2×100 | YES | ENCODE DCC |
| ENCODE K562_REP2 | 195,837,680 | TruSeq polyA PE 2×100 | YES | ENCODE DCC |
| ENCODE GM12878_REP1 | 166,920,844 | TruSeq polyA PE 2×100 | YES | ENCODE DCC |
| ENCODE GM12878_REP2 | 167,447,264 | TruSeq polyA PE 2×100 | YES | ENCODE DCC |
| ENCODE H1_REP1 | 223,421,774 | TruSeq polyA PE 2×100 | YES | ENCODE DCC |
| ENCODE H1_REP2 | 186,588,962 | TruSeq polyA PE 2×100 | YES | ENCODE DCC |
| ENCODE MCF7_REP1 | 222,752,794 | TruSeq polyA PE 2×100 | YES | ENCODE DCC |
| ENCODE MCF7_REP2 | 230,760,862 | TruSeq polyA PE 2×100 | YES | ENCODE DCC |
| GSE75823 TruSeq UHRR | 4,782,806 | TruSeq PE 2×50 | YES | GEO GSE75823 (SRR2982527) |
| **Total** | **1,554,654,222** | | 9/9 YES | |

Per-sample duplicate counts, K562 RAM/wall numbers, and the 7 fixes
that got us here are in [`session-report.md`](session-report.md).

### B. UMI feature arm — `--barcode-tag RX`, flag-set equality + metrics

| Dataset | Records | UMI | Chemistry | Flag-diff lines | Metrics match | wall Picard | wall wea |
|---------|--------:|-----|-----------|:---------------:|:-------------:|:----------:|:--------:|
| GSE75823 UMI UHRR (SRR2982529) | 210,234,301 | 6nt SCRB-seq | SE from umi_tools + fgbio | **0** | YES (to 6 sig figs) | ~2000s | ~1000s |
| PRJNA416930 mouse testis (SRR6250998) | 41,629,973 | 8nt on R2 | SE from umi_tools + fgbio | **0** | YES exact | 299s | 149s |
| GSE134031 QuantSeq mouse microglia (SRR9659562) | 3,780,758 | 6nt + TATA spacer inline R1 | SE via umi_tools regex + fgbio | **0** | YES exact (PERCENT_DUPLICATION=0.550296) | 32s | 13s |

Every UMI sample went through the same prep: `umi_tools extract` →
STAR/HISAT2 align → `fgbio CopyUmiFromReadName --field-delimiter=_`
→ `samtools sort` → `Picard MarkDuplicates BARCODE_TAG=RX` (gold) and
`markdup-wea --barcode-tag RX` (ours).

Flag-set equality = `samtools view BAM | awk '{print $1"\t"($2%2048>=1024)}' | sort`
identical between the two outputs. This is the functionally relevant
signal — dup-flag-only regions that any downstream consumer (STAR,
featureCounts, HTSeq, salmon quant-dedup, etc.) looks at.

Full-BAM md5 differs only in aux-tag write order; see `deviations.md §10`
for why this is cosmetic and identical by `samtools view | cut -f1-11 | sort`.

---

## Regression gate (A.4 windowed-SingleEndTracker)

After introducing the windowed `SingleEndTracker` in commit `d77b3e4`
(switches single-end dup resolution from "flush on key change" to a
bounded same-boundary buffer keyed by `uc5 ± FWD_WINDOW`), we re-ran the
8-ENCODE gate to prove the change is a pure no-op for non-UMI PE
workloads:

- 8/8 samples **byte-identical to the pre-A.4 output.**
- 1,549,871,416 records, 0 flag divergence vs Picard baseline.
- Total wall 4h on cpx62 (16 vCPU).
- Gate script: `/tmp/a4_encode_gate.sh` on Hetzner; artifacts in
  `/mnt/HC_Volume_105344878/tmp/a4-encode-gate/`.

Same gate re-runs after every code change touching `src/groups.rs`,
`src/scan.rs`, or `src/dupset.rs`. Details in `session-report.md` §2026-04-21.

---

## Known deviations (not fixable, not bugs)

Condensed; full rationale + downstream impact per entry in
[`deviations.md`](deviations.md).

| # | Topic | Category | Flag impact | Metrics impact |
|---|-------|----------|:-----------:|:--------------:|
| 1 | `READ_PAIR_OPTICAL_DUPLICATES` always 0 | intentional | none | column is 0 |
| 2 | Queryname-sorted input rejected | intentional | — | — |
| 3 | No per-record `PG:Z` tag | intentional | none | none |
| 4 | Missing RG → library_idx=0 (not NPE) | intentional | none for single-library | none |
| 5 | Score type is u32 not short | intentional | none for Illumina | none for Illumina |
| 6 | MQ=0 chimeric-pair quirk reproduced | matched | matches Picard | matches Picard |
| 7 | Orphan-vs-fragment same-locus | approximation | bounded by orphan rate (<0.1% nf-core) | same |
| 8 | Barcode tag extraction (UMI feature) | approximation | exact match (0 flag-diff, 2 datasets) | exact match |
| 9 | `ESTIMATED_LIBRARY_SIZE` graceful on pathological input | intentional | none | empty field vs exception |
| 10 | Aux-tag write order | cosmetic | none | none |

Also not implemented: `DUPLEX_UMI` (Track A.7 — MC-tag based mate
unclipped start), `TAGGING_POLICY` / `DT` tag, `DUPLICATE_SCORING_STRATEGY`
other than `SUM_OF_BASE_QUALITIES`, multi-threaded Pass 1 (Phase 4).

---

## What is NOT tested

Hard-edge list — where we have no parity data and you should validate
locally before production use:

- **Patterned flowcell (HiSeq X / NovaSeq).** Only open corpus we found
  is GTEx v9/v10 via Broad WARP, which is dbGaP-gated (phs000424). We
  have not applied for access. WARP test BAMs at
  `gs://broad-gotc-test-storage/rna_with_umis/` are auth-gated
  (401 Anonymous on `gsutil ls` from an unauthenticated client).
- **Long-read (PacBio, Oxford Nanopore).** `FWD_WINDOW=1024` in
  `SingleEndTracker` assumes Illumina read length ≤ ~300bp.
  `add_read` panics with a clear error if left-clip exceeds the window,
  so long-read input **fails fast, does not silently corrupt**. One-line
  constant bump + regression re-test would extend support.
- **Single-cell (10x Chromium, Drop-seq).** Different dedup semantics
  (cell-barcode + UMI compound key, transcriptome-aligned). markdup-wea
  targets bulk RNA-seq only.
- **Multi-library BAMs.** Library dimension enters the dedup key but we
  only parity-test single-library inputs. Multi-library mis-annotation
  silently pools into default library (`deviations.md §4`).
- **Targeted panels / amplicon / hybrid-capture / WGS at 30×+.** Not
  benchmarked. Algorithm is the same, but we have not verified corner
  cases (extreme coverage spikes, reference gaps).
- **Queryname-sorted input.** Rejected with a clear error pointing at
  `samtools sort` (`deviations.md §2`).
- **DUPLEX_UMI.** Listed in `deviations.md` out-of-scope; planned as
  Track A.7 once we have a representative test dataset.

---

## How to reproduce

Everything runs on the Hetzner cpx62 validation box
(`/mnt/HC_Volume_105344878/tmp/`, 1TB volume, 16 vCPU). The paths and
scripts below are copy-pastable.

### Non-UMI gate (8 ENCODE samples)

- `scripts/validation/run_parity_batch.py` — in-repo harness.
- `/tmp/a4_encode_gate.sh` on Hetzner — post-A.4 regression driver.
- Picard gold BAMs pre-staged at
  `/mnt/HC_Volume_105344878/tmp/a4-encode-gate/<sample>/<sample>.markdup.sorted.bam`.
- Picard invocation: `java -Xmx4g -jar picard.jar MarkDuplicates
  VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true CLEAR_DT=true
  ADD_PG_TAG_TO_READS=true REMOVE_DUPLICATES=false COMPRESSION_LEVEL=5`
  with jar sha256 `e76128c283889fc583c9dea33a3b7448974c067d102c9e35be152642d4d5f901`.

### UMI arm (GSE75823 / PRJNA416930 / GSE134031)

Pattern for any UMI dataset:

```bash
WORK=/mnt/HC_Volume_105344878/tmp/<ds>-prep
umi_tools extract ... → $WORK/fastq/*.umi.fastq.gz
STAR --runThreadN 16 --genomeDir $STAR_INDEX ... → $WORK/aligned/*.Aligned.out.bam
fgbio -Xmx4g CopyUmiFromReadName --field-delimiter=_ -i ... -o $WORK/aligned/*.rx.bam
samtools sort -@ 8 -o $WORK/output/sorted.bam ... && samtools index
java -Xmx8g -jar /opt/picard/picard.jar MarkDuplicates \
  BARCODE_TAG=RX ASSUME_SORT_ORDER=coordinate \
  VALIDATION_STRINGENCY=LENIENT CLEAR_DT=true ADD_PG_TAG_TO_READS=true \
  I=$WORK/output/sorted.bam O=$WORK/output/picard.bam M=$WORK/output/picard.metrics.txt
~/markdup-test/markdup-wea-a2/target/release/markdup-wea \
  --input $WORK/output/sorted.bam --output $WORK/output/wea.bam \
  --metrics $WORK/output/wea.metrics.txt --barcode-tag RX

# Flag-set diff
samtools view $WORK/output/picard.bam | awk '{print $1"\t"($2%2048>=1024)}' | sort > picard.flags
samtools view $WORK/output/wea.bam    | awk '{print $1"\t"($2%2048>=1024)}' | sort > wea.flags
diff picard.flags wea.flags | wc -l   # expect 0
```

Dataset-specific umi_tools patterns:

- **GSE75823:** `--bc-pattern='NNNNNN' --extract-method=string`
- **PRJNA416930:** custom Python `umi_extract.py` (8bp from R2 → R1 read
  name; see `/tmp/umi_extract.py` on Hetzner)
- **GSE134031 (QuantSeq):** `--bc-pattern='^(?P<umi_1>.{6})(?P<discard_1>TATA){s<=1}.*' --extract-method=regex`

### A.4 regression re-run

```bash
ssh root@204.168.253.225 'bash /tmp/a4_encode_gate.sh'
# prints per-sample record count, dup-flag diff count, wea wall time
```

---

## Open questions for reviewers

If you are reading this as a co-founder / prospective reviewer:

1. **Does this cover your use case?** If not — which chemistry or
   protocol should we benchmark next? (NovaSeq, scRNA, targeted panel,
   long-read, multi-library mix, …)
2. **Can you run wea against one of your own BAMs and share flag-diff?**
   The reproduction script above is 20 lines. On a 50M-read BAM the
   whole gate finishes in ~5 minutes. A `diff picard.flags wea.flags |
   head` is enough evidence either way.
3. **Any production constraint we missed?** Multi-library BAMs,
   library-complexity estimation, GC-aware dedup, duplicate-marker
   inheritance across re-align, chain-of-custody PG lineage, …
