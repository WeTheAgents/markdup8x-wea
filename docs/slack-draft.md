# Slack draft — reply into the markdup-wea thread

> Not sent. Review, edit, paste into the same thread.
> All links resolve on https://github.com/WeTheAgents/markdup8x-wea/tree/main

---

Jonathan, Thomas — thanks, received. Agreed on downstream parallelism (scheduler, not our box); retracting that pitch. Agreed the burden of proof on markdup is high. Below is what's moved on the two lines where I think there's still something to pick up.

**Memory tuning** — per-sample Picard RSS for the 8 ENCODE set is in [`docs/perf.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/perf.md#phase-c--full-batch-correctness-on-8-encode-samples): 16.8–17.5 GB real peak against the 28 GB heap, mean 17.0 GB. Thomas — per-locus depth is the right axis, and I don't have that data yet; happy to re-measure on whichever UMI / deep-panel / amplified BAM you'd point me at. Retry envelope preserved either way. Proposal doc: [`docs/picard-tuning-proposal.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/picard-tuning-proposal.md).

**markdup-wea** — added `BARCODE_TAG` / UMI support. Validated vs Picard 3.4.0 on 3 independent public UMI datasets (2 organisms × 3 chemistries, combined 255.6M records, 0 flag divergence, exact metrics). Evidence pack with the explicit "not tested" list: [`docs/parity-evidence.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/parity-evidence.md). Deviations: [`docs/deviations.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/deviations.md).

**Ask** — anyone willing to run wea against one of their BAMs and share a flag-diff? 20-line script, ~5 min on 50M reads, I'll own any non-zero. Script: [`docs/parity-evidence.md#how-to-reproduce`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/parity-evidence.md#how-to-reproduce).

— A.

---

# Follow-up #2 — Stage 1/3 results (drafted 2026-04-23)

> Not sent. Review, edit, paste into the same thread.

Quick follow-up on the three leads from this thread:

- **Per-locus pileup across 8 ENCODE + 3 UMI libraries** — answers Thomas's per-locus axis directly. Max base-depth + max reads-per-(unclipped-5′, strand) bucket per sample, correlated against tool-level peak RSS: [`docs/per-locus-pileup.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/per-locus-pileup.md). Deepest 5′-bucket in the set is 212,846 reads (GM12878_REP2) — ~0.4× Picard's `MAX_RECORDS_IN_RAM`; none of our 8 ENCODE libraries push the sorting-collection threshold on per-locus depth.

- **Picard memory envelope on GSE75823 (210M SE UMI, `BARCODE_TAG=RX`)** — Jonathan's ask. `-Xmx` sweep [4g, 6g, 9g, 28g], flag-md5 MATCH baseline on all four, peak RSS 2.2 / 3.2 / 5.1 / 18.0 GB, wall flat at ~21.5 min. Min non-OOM heap = **`-Xmx4g`** (working set is ~2 GB; `-Xmx28g` sits on 16 GB it never touches). Proposal §C target `-Xmx6g` holds for UMI SE with ~3 GB headroom. Full table: [`docs/picard-tuning-proposal.md#gse75823-umi-library-memory-envelope`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/picard-tuning-proposal.md#gse75823-umi-library-memory-envelope).

- **Phil's nf-test suggestion** — packaged the rewrite as a picard-shim docker image, ran `modules/nf-core/picard/markduplicates` tests: **6/6 pass** (3 real + 3 stub, byte-for-byte snapshot match on metrics line 2). Image is local-only for now, no registry push: [`docs/nftest-dropin.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/nftest-dropin.md).

No asks. Three leads from this thread, closed. If nothing here moves the needle, we stop.

— A.
