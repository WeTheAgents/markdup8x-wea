# Per-locus pileup — 8 ENCODE + 3 UMI libraries

Thomas's correction in the Slack thread was precise: *"the memory challenge does not come from the total number of reads (alone), but from the number of reads mapping to any given locus"*. We conceded the gap in the reply draft ([docs/slack-draft.md](slack-draft.md)). This document closes it: for each of our 11 bench libraries we report two locus-depth metrics alongside the existing tool-level peak-RSS numbers.

## Metrics

- **M1 — Max base coverage**: `samtools depth -a` over all mapped bases, chromosome/position/count of the single deepest base.
- **M2 — Max reads per (unclipped-5′ coordinate, strand) bucket**: the quantity MarkDuplicates actually keys on (`ReadEnds` is hashed by unclipped 5′ + strand + library), derived from a strand-aware CIGAR parse. This is the direct memory proxy for Picard's coordinate path.

Primary data on Hetzner: `/mnt/HC_Volume_105344878/tmp/pileup/<sample>/{records,m1,m2}.tsv`.

## Data

| Sample | Records | M1 max base depth (loc) | M2 max 5′-bucket (loc:strand → count) | Picard peak RSS | wea peak RSS |
|---|---:|---|---|---:|---:|
| ENCODE GM12878_REP1 | 166,920,844 | 577,366 @ 21:9827270 | 21:9827249:F → 86,874 | 17.0 GB¹ | 0.31 GB¹ |
| ENCODE GM12878_REP2 | 167,447,264 | 1,065,704 @ 21:9827271 | 21:9827249:F → 212,846 | 17.0 GB¹ | 0.31 GB¹ |
| ENCODE H1_REP1 | 223,421,774 | 1,302,187 @ MT:6979 | MT:5917:F → 115,724 | 17.0 GB¹ | 0.31 GB¹ |
| ENCODE H1_REP2 | 186,588,962 | 586,280 @ 21:9827090 | 21:9827361:R → 138,807 | 17.0 GB¹ | 0.31 GB¹ |
| ENCODE K562_REP1 | 156,141,236 | 932,849 @ 21:9827269 | 21:9827249:F → 202,731 | 17.0 GB¹ | 0.31 GB¹ |
| ENCODE K562_REP2 | 195,837,680 | 932,783 @ 21:9827269 | 14:50329559:R → 100,056 | 17.0 GB¹ | 0.31 GB¹ |
| ENCODE MCF7_REP1 | 222,752,794 | 1,097,230 @ MT:1939 | MT:1683:F → 193,345 | 17.0 GB¹ | 0.31 GB¹ |
| ENCODE MCF7_REP2 | 230,760,862 | 1,299,723 @ MT:1932 | MT:1683:F → 209,791 | 17.0 GB¹ | 0.31 GB¹ |
| UMI GSE75823 (SCRB-seq) | 210,234,301 | 323,838 @ 21:9827335 | 15:69746012:F → 63,902 | 2.23 GB² | 0.15 GB² |
| UMI PRJNA416930 | 41,629,973 | 663,954 @ chr12:69408103 | chr9:78082620:F → 79,235 | — | — |
| UMI GSE134031 | 3,780,758 | 39,347 @ chrM:1870 | chrM:1850:F → 17,307 | — | — |

¹ ENCODE tool-level RSS numbers from [docs/perf.md](perf.md) and [docs/picard-tuning-proposal.md](picard-tuning-proposal.md); mean across the 8-sample sweep.
² GSE75823 RSS from the `-Xmx4g` run of the Stage 2 sweep (picard at `-Xmx4g` peaks at 2.23 GB; wea at 147 MB on the same input). PRJNA416930 / GSE134031 RSS not measured in this pass.

### Note on MCF7_REP1 M2

The in-memory awk counter OOM'd on 222 M reads (too many distinct 5′-bucket keys to fit in `awk` associative memory on a 30 GB / 0-swap host). Refilled via a disk-backed `sort | uniq -c` pipeline: MT:1683:F → 193,345 reads. M1 was unaffected by the OOM.

## Observations

**Ranking by M2 (the MarkDuplicates-relevant axis).** On the ENCODE 8-pack the top three are GM12878_REP2 (212,846), MCF7_REP2 (209,791), K562_REP1 (202,731) — all coordinate-sorted RNA-seq. Four of the eight top buckets sit at `21:9827249:F` (a chr21 ribosomal hotspot in that test region); the rest are chrMT/MT or scattered.

**Ranking by M1 (raw coverage).** H1_REP1 tops at 1.3 M bases at a single MT position; four of the eight peaks are on the mitochondrial contig. That is a coverage measurement, not a memory measurement — mtDNA depth is high but reads collapse into few 5′ buckets, which is why M2 ≪ M1 everywhere.

**Correlation with tool RSS.** The 8 ENCODE libraries occupy Picard peaks 16.8–17.5 GB against a 28 GB heap (mean 17.0 GB), independent of the M2 spread (86k → 213k, ~2.5×). Picard's memory floor on these libraries is dominated by the sorting-collection structures, not the per-locus bucket size — which matches Picard's own architecture (pre-allocated on-disk collector with `MAX_RECORDS_IN_RAM=500000`). M2 would only drive RSS up if it exceeded the collection threshold, which it does not on ENCODE.

For UMI libraries the correlation is the open question — Stage 2 sweep on GSE75823 will land the data.

**Thomas's axis, answered directly.** The deepest 5′-bucket we've seen in this set is 212,846 reads (GM12878_REP2). That is ~0.4× `MAX_RECORDS_IN_RAM`. No ENCODE sample pushes Picard past the sorting-collection threshold on the basis of per-locus depth. If there is a pathological library in the fleet, it is not in this 11-sample set — and we would want to see it before recommending a heap reduction below what these libraries already establish as safe.

## Method

Driver: `/mnt/HC_Volume_105344878/tmp/pileup/run_pileup.sh` (on Hetzner). Four-wide via `xargs -P 4`, each worker running `samtools` with `--threads 4` (16 vCPU total). Strand-aware uc5 calculation:

- Forward (`flag & 16 == 0`): `uc5 = pos − leading-soft-clip`.
- Reverse: `uc5 = pos + ref-consuming-ops + trailing-soft-clip − 1`.

`-F 0x104` filters unmapped + secondary; supplementary are kept (they contribute to the duplicate-marking pool). CIGAR `*` entries skipped.

## Reproduce

On any host with `samtools` ≥ 1.18 and `gawk`:

```bash
BAM=/path/to/sample.sorted.bam

# M1
samtools depth -a --threads 8 "$BAM" \
  | awk 'BEGIN{m=0} $3>m{m=$3;c=$1;p=$2} END{print c"\t"p"\t"m}'

# M2 (awk in-memory; see MCF7 note for disk-backed fallback)
samtools view -F 0x104 --threads 4 "$BAM" \
  | gawk '
    $6=="*"{next}
    {
      flag=$2; ref=$3; pos=$4; cigar=$6
      rev=(and(flag,16)!=0)
      if(!rev){ lead=0; if(match(cigar,/^[0-9]+S/)) lead=substr(cigar,1,RLENGTH-1)+0; uc5=pos-lead; strand="F" }
      else{ rest=cigar; refspan=0
        while(match(rest,/[0-9]+[MIDNSHP=X]/)){ tok=substr(rest,RSTART,RLENGTH); op=substr(tok,length(tok)); n=substr(tok,1,length(tok)-1)+0; if(op~/[MDN=X]/) refspan+=n; rest=substr(rest,RSTART+RLENGTH) }
        trail=0; if(match(cigar,/[0-9]+S$/)) trail=substr(cigar,RSTART,RLENGTH-1)+0; uc5=pos+refspan+trail-1; strand="R" }
      key=ref":"uc5":"strand
      if(++c[key]>m){ m=c[key]; mk=key }
    }
    END{ print mk"\t"m }'
```
