# Reproduce the Quadron methylation paper

This is the curated entry point for the results in [../G4_methylation_paper.md](../G4_methylation_paper.md). It runs **only the two final Quadron minimal models**. The older Nextflow defaults, other annotation methods, descriptive tests, chromatin sensitivities, and elastic-net analyses are not part of this replication route.

**Verified locally:** both a saved-slope replay and a reconstruction from prepared beta matrices reproduced all six manuscript result tables. The prepared route also matched each of the three complete slope/annotation checkpoints. The scientific source files and manuscript figures were not overwritten. See [validation_summary.md](validation_summary.md) for the recorded checks and remaining limits.

## 1. What is being reproduced?

The first stage estimates a signed methylation-versus-age slope for each CpG separately within each tissue:

```r
methylation_beta ~ age + sex + bmi
```

The final models analyze those estimated slopes:

```r
# M0: no G-richness correction, but mean methylation is included.
slope ~ has_stable_g4 + has_unstable_g4 + mean_beta

# M1: exactly the two additional base-composition covariates.
slope ~ has_stable_g4 + has_unstable_g4 + mean_beta +
        probe_strand_g_fraction + opposite_strand_g_fraction
```

The two figure/result prefixes are:

- `quadron_minimal_stable_unstable_mean_beta`
- `quadron_minimal_strand_g_richness`

Each produces coefficients, G4 terms, and a direct stable-minus-unstable contrast. These six CSVs are compared with the unchanged paper bundle. Fresh PNG/PDF plots are generated separately; the original two manuscript PNGs remain the reference display assets.

## 2. Curated script inventory

There is one maintained implementation of each scientific step. The files below were commented in place rather than copied into a divergent second pipeline.

| Script | Purpose | Human review points |
|---|---|---|
| [run.R](run.R) | Explicit orchestration and safe output handling | Exact three cohorts, model flags, staged inputs, nonempty-output refusal, command logs, and checksums |
| [config.yaml](config.yaml) | Paths and declared settings | Paths resolve relative to BTEP; G4 padding and sequence flank are distinct |
| [../bin/download_geo_series_matrix.R](../bin/download_geo_series_matrix.R) | Import a frozen or explicitly downloaded processed GEO series matrix | Not IDAT normalization; sample-ID alignment; numeric metadata conversion; `--series-matrix` enables offline replay |
| [../bin/prepare_g4_windows.R](../bin/prepare_g4_windows.R) | Pad and merge class-separated Quadron intervals | Input-file labels are not calculated thresholds; overlapping/adjacent windows merge; historical coordinates are preserved |
| [../bin/analyze_g4_methylation_age.R](../bin/analyze_g4_methylation_age.R) | Fit CpG age slopes and attach merged-window annotations | Sample filtering, rank/precision, mean beta, inclusive overlap behavior; `--slopes-only` skips older results |
| [../bin/compare_g4_effects_adjusted.R](../bin/compare_g4_effects_adjusted.R) | Fit M0/M1 and generate statistics/plots | Minimal branch, overlapping indicators, strand-oriented fractions, equal CpG weights, direct contrasts, and BH families |
| [compare_results.R](compare_results.R) | Numerical reproduction gate | Semantic row keys, all columns, exact probe counts/formulas, small-P comparison, protected report writes |
| [collect_validation.R](collect_validation.R) | Package diagnostic summaries from the two documented runs | Checks successful comparisons before exporting reviewable summaries; no scientific model fitting |
| [install_packages.R](install_packages.R) | Check or explicitly install dependencies | Check-only default; no automatic upgrades; separates current tested environment from historical provenance |
| [../render_G4_methylation_paper.R](../render_G4_methylation_paper.R) | Render Rmd to Markdown and offline HTML | Does not fit scientific models; embeds the original figures and KaTeX resources |
| [../validate_G4_methylation_paper.R](../validate_G4_methylation_paper.R) | Validate the manuscript artifacts | Source hashes, exact figure paths, scope, tables, citations, links, and embedded resources |

Tests: [tests/test_inputs.R](tests/test_inputs.R), [tests/test_comparison.R](tests/test_comparison.R), and [tests/test_paths.R](tests/test_paths.R). Each uses small temporary fixtures and makes no external data request.

## 3. Environment

Run from the workspace root in the examples below; from BTEP, omit the initial `BTEP/` in commands.

```bash
# Check only: no installation or data download.
Rscript BTEP/replication/install_packages.R

# Optional: explicitly install missing CRAN/Bioconductor packages.
Rscript BTEP/replication/install_packages.R --install
```

[environment_validated.csv](environment_validated.csv) records the successful replay environment, including R 4.5.2, data.table 1.17.8, ggplot2 4.0.0, and hg38 BSgenome 1.4.5. It is a tested version inventory, **not a dependency lockfile or a recovered historical environment**. The installer selects a Bioconductor release compatible with the active R and does not upgrade already installed packages. Re-run the numerical gate after any environment change. A dedicated R library set through `R_LIBS_USER` avoids modifying a shared library.

The final script still computes some shared helper annotations that are absent from the minimal formulas. Consequently, **even M0** requires the 450K annotation package, hg38 genome package, Bioconductor range/sequence packages, and unmerged G4 checkpoint. The runner retains these dependencies rather than refactoring the code that generated the results. Neither Nextflow nor SLURM is required for this curated Rscript route.

## 4. Preferred route: saved-slope replay

By default [config.yaml](config.yaml) points to the existing `results/g4_methylation_age_quadron_100bp` tree under BTEP. It requires:

```text
per_tissue/GSE61257_adipose_probe_slopes.csv.gz
per_tissue/GSE61258_liver_probe_slopes.csv.gz
per_tissue/GSE61259_muscle_probe_slopes.csv.gz
g4_windows/g4_motifs_100bp_unmerged.bed
```

The slope tables must include stored `mean_beta`, coordinates, motif counts, and cohort/tissue labels. Requiring stored means avoids the older script's fallback to averaging potentially different prepared specimens. Exactly the three files above are staged, so an extra file in the source directory cannot silently enter the regression.

```bash
# Read-only input/package/reference-table validation.
Rscript BTEP/replication/run.R --mode check

# Refit only M0/M1 and verify all six resulting tables.
# This example uses a new directory; reruns must choose another empty path.
Rscript BTEP/replication/run.R --mode models \
  --outdir results/quadron_paper_reproduction_new

# Verify an existing completed run without refitting.
Rscript BTEP/replication/run.R --mode verify \
  --outdir results/quadron_paper_reproduction_new
```

The default output directory already contains the verified replay from this curation session. Use `--mode verify` on that run, or select a new `--outdir` for another numerical run. Inputs, scripts, the paper bundle, and populated run directories are protected. Traversal through nonexistent parents and symlinked report destinations are rejected. Failed runs are deliberately retained with logs; do not reuse their output directories.

The sequence flank is explicitly set to **1,000 bp**. That setting reproduced the reference statistics in the tested environment, but this is empirical replay evidence, not a recovered historical command. Changing the flank creates a new candidate run that must pass the same comparison; do not change the expected tables or relax tolerances to hide differences.

## 5. Reconstruct from prepared beta matrices

This route was also verified. It additionally requires the three prepared beta RDS/metadata CSV pairs, the merged G4 checkpoint, and the shared probe manifest configured in YAML.

```bash
Rscript BTEP/replication/run.R --mode prepared \
  --outdir results/quadron_paper_reproduction_prepared_new
```

The runner checks unique probe IDs and mapped hg38 coordinates against the saved slope checkpoints. Unmapped/nonstandard entries are retained in the original manifest but excluded by the existing analysis script; no coordinates are invented or shifted. In the locally checked manifest, 381,179 of 657,261 entries map to the standard chromosomes. Every saved checkpoint coordinate matched. This establishes compatibility with the retained artifacts, not the origin of the original manifest.

Age/sex/BMI columns, sample-ID uniqueness, and specimen-level design rank are checked. The original first-stage script is called with `--slopes-only`, so no older Wilcoxon or architecture-heavy regression side products are run. Each regenerated **full** slope/annotation table is compared with its saved checkpoint before M0/M1 are fitted. A mismatch stops the run and writes a per-tissue report. This makes it possible to distinguish an upstream reconstruction failure from a second-stage model failure.

## 6. Optional processed-source route

This route is implemented and fixture-tested, but was **not validated with the original full source inputs** because the original class-separated Quadron BEDs and their generation provenance are not available locally.

Supply the original stable/unstable BED paths in the configuration, together with frozen processed GEO matrix files where available. The Quadron inputs must have the historical schema expected by the window builder; do not substitute a current prediction export without checking it.

```bash
# With all three processed GEO files specified locally:
Rscript BTEP/replication/run.R --mode source \
  --config BTEP/replication/config.yaml \
  --outdir results/quadron_paper_reproduction_source_new

# Only if fresh GEO retrieval is intended, add --allow-downloads.
```

`--mode source` rebuilds the windows, imports processed GEO matrices, estimates slopes, and fits M0/M1. It still uses the saved slope checkpoints and six paper tables as comparison targets. GEO downloads require explicit authorization through `--allow-downloads`; input and code hashes, commands and logs are recorded. New downloads may differ from historical files.

**Not provided or claimed:** raw-IDAT normalization, original probe-QC decisions, training/running Quadron, or reconstruction of its original stability cutoff. The supplied stable/unstable file labels are inherited, not recomputed. Source-level replay cannot establish those missing provenance facts by itself.

## 7. Outputs and success criteria

Each new run has the following organization:

```text
provenance/
  config_requested.yaml
  run_resolved.yaml
  input_manifest.csv
  package_versions.csv
  session_info.txt
  commands.txt
  source_used/                 # Copies of the code/config actually executed
logs/
staged_slopes/                 # Exact three-cohort input selection
minimal_model/
minimal_model_strand_g_richness/
verification/
  numerical_comparison.csv
```

Prepared/source runs add their reconstructed inputs and per-tissue slope-parity reports. Commands record every option explicitly, including the disabling of chromatin, architecture, density stratification and cross-reactive filtering in the final models. No structure-annotation file is supplied. Matching the saved statistics confirms the effective replay, not whether a historical redundant filter was once used.

The gate checks all columns of all six CSVs, including exact formulas, row identities, probe counts, estimates, SEs, test statistics, P values, BH values, and contrasts. Ordinary numerical fields use absolute tolerance `1e-12` plus relative tolerance `1e-7`; P/q values use **relative tolerance only**, preventing a tiny-P mismatch from being hidden by an absolute floor. NaN/NA patterns and nonfinite values must agree. A failure returns nonzero status and preserves diagnostics.

New figures reproduce the plotted estimates and tests, not necessarily byte-identical images. Device/font/package differences affect graphics, and the current plotter includes direct-comparison brackets that were absent from one retained figure. The paper's original PNGs remain unchanged. Re-rendering the manuscript uses those original images, not the new replay images.

## 8. Rebuild and check the documents

```bash
# Network access may be required to embed KaTeX; the finished HTML is offline.
Rscript BTEP/render_G4_methylation_paper.R
Rscript BTEP/validate_G4_methylation_paper.R

# Fast offline tests of the computational entry points and safety gates.
Rscript BTEP/replication/tests/test_inputs.R
Rscript BTEP/replication/tests/test_comparison.R
Rscript BTEP/replication/tests/test_paths.R
```

Document validation is **not** model re-estimation. A rendered paper can be correct while a new scientific run fails numerical parity; these are intentionally separate checks.

## 9. Human double-check checklist

- Are all three specimens/metadata files from the named GEO subseries, aligned by sample ID rather than position? Shared subject IDs do not imply 137 independent people.
- Is `mean_beta` the mean over the same matched specimens as the slopes, rather than a young reference, a new normalization, or an imputed mean?
- Are G4 windows padded by 100 bp while sequence composition uses the separate 1,000-bp probe-centered flank?
- Were the original supplied coordinate numbers retained? Standard BED is zero-based/half-open; the historical `foverlaps` implementation is closed and sequence extraction uses one-based `GRanges`. This replay does not silently fix the unresolved coordinate contract.
- Do stable and unstable indicators both remain true for mixed merged neighborhoods? Motif classes must not be forced into exclusive priority bins.
- Do the two formulas contain exactly the documented terms? G fractions are base composition, not strand-specific motif counts or a chromatin correction.
- Do the models retain equal CpG weights and conventional OLS uncertainty? Correlated probes and first-stage slope precision are not accounted for by this reproduction.
- Is the stable-minus-unstable difference tested using the coefficient covariance? A significant term versus a non-significant term is not evidence that the classes differ.
- Are BH adjustments made across three tissues separately for each term and separately for contrasts within each model?
- Are plotted effects multiplied by 10 for beta/decade, and by another 100 only for percentage points? Whiskers are ordinary approximate 95% intervals, not BH-adjusted intervals.
- Do both numerical gates pass before a rerun is described as a reproduction? Preserve all mismatch reports and do not overwrite the paper's expected tables.

The extensive comments in the curated scripts point to these exact calculation and control-flow locations. A corrected coordinate system, dependence-aware standard errors, different normalization, additional covariates, or a new class cutoff would be a separately labelled scientific extension, not an invisible replication repair.