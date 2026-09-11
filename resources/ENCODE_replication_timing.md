# ENCODE Replication-Timing Candidates

This file records ENCODE Repli-seq datasets evaluated for the BTEP structural-context sensitivity analysis. The current BTEP model requires an hg38 four-column interval file:

```text
chr  start  end  replication_timing_score
```

## Recommended Path

Use a paired early-S/late-S ENCODE replication-timing series on GRCh38, calculate a binned timing score such as $\log_2((E + \epsilon)/(L + \epsilon))$, and export it as an hg38 bedGraph. The best verified paired GRCh38 candidate is the LNCAP series:

| Phase | Experiment | Alignment file |
|---|---|---|
| Early S | [ENCSR385QAX](https://www.encodeproject.org/experiments/ENCSR385QAX/) | [ENCFF312QYV](https://www.encodeproject.org/files/ENCFF312QYV/@@download/ENCFF312QYV.bam) |
| Late S | [ENCSR089VDE](https://www.encodeproject.org/experiments/ENCSR089VDE/) | [ENCFF403WZE](https://www.encodeproject.org/files/ENCFF403WZE/@@download/ENCFF403WZE.bam) |

These experiments belong to [Replication Timing Series ENCSR888VII](https://www.encodeproject.org/replication-timing-series/ENCSR888VII/). This is an exploratory shared genomic-context proxy: LNCAP is a prostate cancer cell line, not liver, skeletal muscle, or adipose.

## Legacy Processed Signals

ENCODE has convenient released bigWig tracks for HepG2, IMR-90, BJ, and primary keratinocyte, but the verified tracks are on **hg19**. Do not set them directly as `replication_timing_url` in BTEP, because BTEP probe coordinates are hg38. They require a documented bigWig-to-bedGraph conversion and hg19-to-hg38 liftover first.

HepG2 is the closest available liver-related proxy, but it is a hepatocellular carcinoma cell line. Its processed track is [ENCFF001GPC](https://www.encodeproject.org/files/ENCFF001GPC/@@download/ENCFF001GPC.bigWig), from [ENCSR000CXG](https://www.encodeproject.org/experiments/ENCSR000CXG/).

## Use in BTEP

After creating an hg38 bedGraph, set a local file rather than a generic URL:

```yaml
replication_timing_bed: "/data/RBL_NCI/ccrrbl20/annotations/encode_lncaps_early_late_hg38.bedGraph.gz"
```

Record the ENCODE experiment and file accessions, the genome build, bin size, normalization, score definition, and any liftover steps in the analysis provenance.