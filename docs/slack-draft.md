# Slack draft — reply into the markdup-wea thread

> Not sent. Review, edit, paste into the same thread.
> All links resolve on https://github.com/WeTheAgents/markdup8x-wea/tree/main

---

Jonathan, Thomas — thanks, received. Agreed on downstream parallelism (scheduler, not our box); retracting that pitch. Agreed the burden of proof on markdup is high. Below is what's moved on the two lines where I think there's still something to pick up.

**Memory tuning** — per-sample Picard RSS for the 8 ENCODE set is in [`docs/perf.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/perf.md#phase-c--full-batch-correctness-on-8-encode-samples): 16.8–17.5 GB real peak against the 28 GB heap, mean 17.0 GB. Thomas — per-locus depth is the right axis, and I don't have that data yet; happy to re-measure on whichever UMI / deep-panel / amplified BAM you'd point me at. Retry envelope preserved either way. Proposal doc: [`docs/picard-tuning-proposal.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/picard-tuning-proposal.md).

**markdup-wea** — added `BARCODE_TAG` / UMI support. Validated vs Picard 3.4.0 on 3 independent public UMI datasets (2 organisms × 3 chemistries, combined 255.6M records, 0 flag divergence, exact metrics). Evidence pack with the explicit "not tested" list: [`docs/parity-evidence.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/parity-evidence.md). Deviations: [`docs/deviations.md`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/deviations.md).

**Ask** — anyone willing to run wea against one of their BAMs and share a flag-diff? 20-line script, ~5 min on 50M reads, I'll own any non-zero. Script: [`docs/parity-evidence.md#how-to-reproduce`](https://github.com/WeTheAgents/markdup8x-wea/blob/main/docs/parity-evidence.md#how-to-reproduce).

— A.
