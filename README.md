# BTEP G4 Methylation-Aging Pipeline

## Reproduce the Paper

For the two final Quadron models in [G4_methylation_paper.md](G4_methylation_paper.md), use the [curated replication guide](replication/README.md). It provides commented scripts, explicit M0/M1 settings, checkpoint and prepared-input routes, package checks, and numerical verification against the manuscript tables. The older exploratory workflow below is not the paper's replication entry point.

This Nextflow DSL2 workflow tests whether human CpG methylation age-slopes are larger within annotated G4Hunter motifs plus 100 bp flanks than outside those coordinates.

The default experiment uses the public GSE61256 methylation subseries:

| Series | Tissue | Samples | Assay |
|---|---:|---:|---|
| GSE61257 | adipose | 32 | Illumina 450K methylation |
| GSE61258 | liver | 79 | Illumina 450K methylation |
| GSE61259 | muscle | 26 | Illumina 450K methylation |

## Run

From this directory:

```bash
nextflow run main.nf -params-file params_GSE61256.yaml -profile standard
```

On the cluster, submit the provided SLURM wrapper:

```bash
sbatch run_gse61256_cluster.sh
```

To run the independent Quadron annotation sensitivity, which writes to a separate result directory, use:

```bash
sbatch run_gse61256_cluster.sh params_GSE61256_quadron.yaml
```

The wrapper accepts an optional params file, profile, and extra Nextflow arguments:

```bash
sbatch run_gse61256_cluster.sh params_GSE61256.yaml slurm --g4_long_min_bp 50000
```

The workflow can also run from local beta/metadata files:

```bash
nextflow run main.nf --samplesheet samplesheet.template.csv -profile standard
```

## Inputs

- `../G4Hunter/G4Hunter.hg38.stable.bed`
- `../G4Hunter/G4Hunter.hg38.unstable.bed`
- `../results/submission_v3/reference_manifest/probe_gene_manifest.csv.gz`

The local-input samplesheet must contain:

```text
cohort,tissue,beta_matrix,metadata,age_col,sample_id_col,covariates
```

Beta matrices should have CpG probe IDs as row names and sample IDs as columns. Metadata must include the sample ID column and chronological age column.

## Coordinate Overlap Policy

G4Hunter has dense and overlapping motifs, even after each motif is expanded by 100 bp. The primary analysis therefore uses a non-overlapping coordinate-of-interest set:

1. Read stable and unstable G4Hunter BED files.
2. Expand every motif by 100 bp on both sides.
3. Merge overlapping or adjacent expanded windows across both classes.
4. Flag unusually long merged clusters using `max(99th percentile width, median + 3 * MAD, optional minimum bp)`.
5. Retain annotation columns with motif counts and class labels: `stable_only`, `unstable_only`, `stable_unstable_overlap`, or `long_g4_cluster`.

This makes every CpG eligible to be counted once in the primary test. For secondary motif-level analyses, use `g4_motifs_100bp_unmerged.bed`, but assign overlapping CpGs by nearest motif, priority class, or fractional weights to avoid double counting.

Long clusters are isolated because they likely represent runs of consecutive or highly dense G4 motifs rather than ordinary single-motif neighborhoods. Keeping them as `long_g4_cluster` prevents the bulk G4-window estimate from being dominated by unusually broad regions. Tune this behavior with `g4_long_quantile`, `g4_long_mad_multiplier`, and `g4_long_min_bp`.

## Main Outputs

- `g4_windows/g4_motifs_100bp_merged.bed`: primary non-overlapping G4 plus 100 bp windows.
- `g4_windows/g4_motifs_100bp_long_clusters.bed`: unusually long merged G4 clusters isolated from the bulk classes.
- `g4_windows/g4_motifs_100bp_unmerged.bed`: motif-level expanded windows for secondary analyses.
- `per_tissue/*_probe_slopes.csv.gz`: probe-level age slopes and G4-window labels.
- `summary/btep_g4_cross_tissue_summary.csv`: per-tissue summaries by G4 class.
- `summary/btep_g4_cross_tissue_tests.csv`: G4-window versus outside-CpG tests.
- `per_tissue/*_g4_regression_effects.csv`: regression estimates using overlapping predictors: `has_stable_g4`, `has_unstable_g4`, their interaction, `is_long_cluster`, and motif density.

## Structural-Context Sensitivity Analysis

The pipeline automatically annotates probes with hg38 RepeatMasker class and segmental-duplication overlap using UCSC tables. On Biowulf it uses the centrally maintained UCSC goldenPath mirror at `/fdb/genomebrowser/goldenPath/hg38/database/` first. It falls back to UCSC internet downloads only when those files are unavailable. The cluster location is configurable through `ucsc_goldenpath_dir`. It then runs an additional G4 sensitivity model adjusted for these annotations. Outputs are written to:

- `sensitivity_annotations/probe_structure_annotations.csv.gz`
- `sensitivity_annotations/probe_structure_annotation_provenance.csv`
- `sensitivity_annotations/probe_structure_annotation_qc.csv`
- `structure_sensitivity/structure_sensitivity_g4_terms.csv`
- `structure_sensitivity/structure_sensitivity_stable_minus_unstable.csv`
- `structure_sensitivity/structure_sensitivity_g4_terms.png`

Optional paths can add further annotations to this sensitivity model:

```text
rloop_bed: "/path/to/hg38_rloop_regions.bed"
rloop_url: "https://example.org/hg38_rloop_regions.bed.gz"
replication_timing_bed: "/path/to/hg38_replication_timing.bed"
replication_timing_url: "https://hgdownload.soeucsc.edu/goldenPath/hg38/database/your_replication_track.bed.gz"
cross_reactive_probes: "/path/to/cross_reactive_450k_probe_ids.txt"
cross_reactive_probes_url: "https://example.org/cross_reactive_450k_probe_ids.txt"
```

URL parameters are downloaded once into the annotation task directory and recorded in the provenance table. UCSC supplies the automatic RepeatMasker and segmental-duplication tables. Select a specific UCSC replication-timing track appropriate for the tissue/model before setting `replication_timing_url`; UCSC does not provide one universal timing profile. R-loop maps and array cross-reactive probe lists should use a cited external source or a locally curated track. R-loop and segmental-duplication tracks are modeled as overlap indicators. Replication timing is modeled as a numeric value. When a cross-reactive probe list is supplied, those probes are excluded before fitting the sensitivity model.

## Regression Interpretation

For the main biological question, use the `age_slope_delta_beta_per_year` response in `*_g4_regression_effects.csv`. This is equivalent to asking whether G4 annotations modify the age effect, without forcing probes into mutually exclusive stable/unstable bins. In full-model notation, the hypothesis is represented by interaction terms such as:

```text
methylation_beta ~ age + sex + bmi + age:has_stable_g4 + age:has_unstable_g4 + age:has_stable_g4:has_unstable_g4
```

Because G4 status is a probe-level feature, the efficient implementation first estimates each probe's age slope adjusted for sex and BMI, then regresses those slopes on overlapping G4 predictors. The `has_stable_g4:has_unstable_g4` term captures regions where stable and unstable annotations overlap rather than discarding them or forcing priority labels.

For the primary stable-versus-unstable contrast, `include_g4_architecture` is `false`: motif density, long-cluster status, and strand-specific motif counts are derived from the same G4 annotation and are therefore treated as exploratory sensitivities rather than confounders. Set `include_g4_architecture: true` only for that secondary sensitivity.
