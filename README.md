# Quadron Methylation and Aging

[main.nf](main.nf) handles all computation for the paper. It uses one [Nextflow module](modules/btep_g4_methylation.nf) and the commented R calculation scripts in `bin/`. [G4_methylation_paper.Rmd](G4_methylation_paper.Rmd) renders the paper directly; there is no separate replication framework or rendering wrapper.

## Run and Render

From BTEP, with the existing prepared data and packages available:

```bash
nextflow run main.nf -profile standard
Rscript -e 'rmarkdown::render("G4_methylation_paper.Rmd", output_format = "all")'
```

This fits CpG age slopes in the three tissues, runs M0/M1, and verifies all six result tables. Outputs go to `results/paper`. Rendering uses the frozen paper tables and original figures, never silently replaces them with a new run, and produces [G4_methylation_paper.md](G4_methylation_paper.md) and [G4_methylation_paper.html](G4_methylation_paper.html). HTML display resources are embedded; internet access is needed at render time for KaTeX, but the finished HTML can be read offline.

Resume a managed run or execute on SLURM with the same entry point:

```bash
nextflow run main.nf -profile standard -resume
nextflow run main.nf -profile slurm --slurm_partition norm -resume
```

Use a new `--outdir` for a separate run. Populated, unmanaged output directories and paths overlapping protected inputs are rejected. Nextflow owns task scheduling, caching, command logs and resume; scientific scripts are staged inputs so code changes invalidate the corresponding cached tasks. Local jobs are serialized by the standard profile to limit memory use.

## Models

| GEO series | Tissue | Prepared specimens |
|---|---|---:|
| GSE61257 | Adipose | 32 |
| GSE61258 | Liver | 79 |
| GSE61259 | Skeletal muscle | 26 |

These are specimen counts, not independent donor counts. Each CpG first uses `beta ~ age + sex + bmi`. Its signed age coefficient is the response for:

```r
# M0: mean-methylation adjusted, without G-richness correction.
slope ~ has_stable_g4 + has_unstable_g4 + mean_beta

# M1: only the two strand-oriented G fractions are added.
slope ~ has_stable_g4 + has_unstable_g4 + mean_beta +
				probe_strand_g_fraction + opposite_strand_g_fraction
```

No earlier regression family, chromatin correction, density stratification, or elastic-net analysis is invoked. The stable and unstable indicators can overlap; neither model has their interaction. `mean_beta` is the same-tissue specimen mean, not an independent baseline. BH adjustments are across three tissues separately for each term and for the direct stable-minus-unstable contrasts. Per-decade effects are yearly coefficients multiplied by ten; the approximate 95% intervals are not multiplicity-adjusted.

## Inputs and Entry Points

All defaults are in [nextflow.config](nextflow.config). The default `--data_dir` is the existing `results/g4_methylation_age_quadron_100bp` directory. Ignored input data and historical outputs are not deleted by this consolidation.

| Entry point | Required inputs under `data_dir` |
|---|---|
| `--start_from prepared` (default) | Three `prepared_inputs/GSE*_tissue_beta.rds` and matching metadata CSVs; merged and unmerged Quadron checkpoints; external `--probe_manifest` |
| `--start_from slopes` | Three `per_tissue/GSE*_tissue_probe_slopes.csv.gz` files with stored mean beta; unmerged Quadron checkpoint |
| `--start_from geo` | Merged/unmerged checkpoints and manifest; processed GEO matrices imported by the pipeline |

Exactly the three named cohort files are selected, not arbitrary extra files in the directory. Beta matrices use probe IDs as rows and sample IDs as columns. Metadata require `sample_id`, `age`, `sex`, and `bmi`. The manifest must provide mapped hg38 coordinates. Both windows are named `g4_motifs_100bp_merged.bed` and `g4_motifs_100bp_unmerged.bed` under `g4_windows`.

```bash
# Faster replay from saved annotated slopes.
nextflow run main.nf -profile standard --start_from slopes --outdir results/paper_slopes

# Explicit processed-GEO import. Without series_matrix_dir this downloads data.
nextflow run main.nf -profile standard --start_from geo --outdir results/paper_geo

# Optional offline GEO imports: directory contains GSE*_series_matrix.txt.gz.
nextflow run main.nf -profile standard --start_from geo \
	--series_matrix_dir /path/to/frozen/matrices --outdir results/paper_geo_cached
```

To rebuild G4 neighborhoods, supply both `--g4_stable` and `--g4_unstable` with the original class-separated Quadron BEDs. Otherwise the existing checkpoints are used. G4 padding is fixed at **100 bp**; the separate probe-centered sequence flank is **1,000 bp**, controlled by `--sequence_flank_bp`. They are not interchangeable.

The original Quadron cutoff/version and raw-IDAT preprocessing history are not reconstructed here. The processed-GEO/original-BED route requires those source inputs and was not tested on the unavailable original BEDs. Existing coordinate and strand conventions are preserved for numerical reproduction; they are not certified by a matching replay.

## Dependencies

Validated with Nextflow 26.04.1, R 4.5.2, data.table 1.17.8, ggplot2 4.0.0, rmarkdown 2.30, knitr 1.51 and Pandoc 3.9. Use a Java version supported by the installed Nextflow. Packages can be installed explicitly in a suitable R library:

```r
install.packages(c("data.table", "optparse", "ggplot2", "rmarkdown", "knitr", "xml2", "BiocManager"))
BiocManager::install(c("IlluminaHumanMethylation450kanno.ilmn12.hg19",
											 "BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "IRanges",
											 "GenomeInfoDb", "Biostrings", "BSgenome"), update = FALSE, ask = FALSE)
```

The genome package is large. Both models retain the original script's genome/annotation helper dependencies, even where a computed helper is not part of M0's formula. No installer or lockfile is generated automatically.

## Verification and Layout

The final Nextflow task checks every column of the six regenerated model tables against the frozen paper bundle: row identities, formulas, counts, estimates, SEs, test statistics, P/q values and direct contrasts. Counts and labels are exact; estimates permit `1e-12 + 1e-7 * abs(reference)` floating-point differences. P/q use relative tolerance only, so materially different tiny probabilities cannot match through an absolute floor. A mismatch fails the run; expected results are never updated to force agreement. `--verify false` explicitly skips this reference comparison for a different analysis, not for a claimed reproduction.

Both real-data routes (`prepared` and `slopes`) were run successfully through Nextflow after consolidation on 2026-09-16, including the six-table verification. Fresh plots reproduce the estimates but may differ cosmetically from the original paper PNGs. The preserved literature review remains in [AI_Project.Rmd](AI_Project.Rmd).

```text
main.nf                       # Only computational entry point
nextflow.config               # Input/settings and local/SLURM profiles
modules/btep_g4_methylation.nf # Import, windows, slopes, M0/M1, verification
bin/                          # Four commented calculation kernels + checker
G4_methylation_paper.Rmd      # Direct render to Markdown and offline HTML
G4_methylation_paper.bib      # References
G4_methylation_paper_files/   # Frozen tables and the two original figures
```

New outputs live under `results/paper/{per_tissue,minimal_model,minimal_model_strand_g_richness,verification}`. Imported matrices and rebuilt windows are published when requested. Nextflow commands and task logs remain in the work directory; the execution trace is `.nextflow-paper-trace.txt`. All run artifacts are ignored by Git rather than checked into another replication-report tree.
