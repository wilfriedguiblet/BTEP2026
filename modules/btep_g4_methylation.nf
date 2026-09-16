// Quote filesystem paths individually so spaces or apostrophes never become
// shell syntax. Fixed cohort/model names are defined in main.nf, not read from
// arbitrary file globs. R scripts are staged inputs and participate in caching.
def shellQuote(value) {
    "'" + value.toString().replace("'", "'\\''") + "'"
}

process IMPORT_GEO {
    tag "${meta.cohort}_${meta.tissue}"
    publishDir "${params.outdir}/prepared_inputs", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(series_matrix)
    path calculation

    output:
    tuple val(meta),
          path("${meta.cohort}_${meta.tissue}_beta.rds"),
          path("${meta.cohort}_${meta.tissue}_metadata.csv"),
          emit: prepared

    script:
    def cached = series_matrix ? "--series-matrix ${shellQuote(series_matrix)}" : ''
    """
    Rscript ${shellQuote(calculation)} \
        --accession ${meta.cohort} --tissue ${meta.tissue} \
        --prefix ${meta.cohort}_${meta.tissue} \
        ${cached} --outdir .
    """
}

process PREPARE_G4_WINDOWS {
    tag 'Quadron +/-100 bp'
    publishDir "${params.outdir}/g4_windows", mode: 'copy', overwrite: true

    input:
    path stable_bed
    path unstable_bed
    path calculation

    output:
    path 'g4_motifs_100bp_unmerged.bed', emit: unmerged
    path 'g4_motifs_100bp_merged.bed', emit: merged
    path 'g4_window_overlap_summary.csv', emit: summary

    script:
    """
    Rscript ${shellQuote(calculation)} \
        --stable ${shellQuote(stable_bed)} --unstable ${shellQuote(unstable_bed)} \
        --flank 100 --long-quantile 0.99 --long-mad-multiplier 3 --long-min-bp 0 --outdir .
    """
}

process FIT_AGE_SLOPES {
    tag "${meta.cohort}_${meta.tissue}"
    publishDir "${params.outdir}/per_tissue", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(beta_matrix), path(metadata)
    path probe_manifest
    path g4_windows
    path calculation

    output:
    tuple val(meta), path("${meta.cohort}_${meta.tissue}_probe_slopes.csv.gz"), emit: slopes

    script:
    """
    Rscript ${shellQuote(calculation)} \
        --beta ${shellQuote(beta_matrix)} --metadata ${shellQuote(metadata)} \
        --manifest ${shellQuote(probe_manifest)} --g4-windows ${shellQuote(g4_windows)} \
        --age-col age --sample-col sample_id --covariates 'sex|bmi' --min-samples 2 \
        --cohort ${meta.cohort} --tissue ${meta.tissue} \
        --prefix ${meta.cohort}_${meta.tissue} --slopes-only --outdir .
    """
}

process FIT_PAPER_MODEL {
    tag "${model.id}"
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    val model
    path slope_files
    path g4_unmerged
    path calculation

    output:
    tuple val(model), path("${model.directory}"), emit: results

    script:
    """
    Rscript ${shellQuote(calculation)} \
        --slopes-dir . --g4-unmerged ${shellQuote(g4_unmerged)} \
        --minimal-model true --include-strand-g-richness ${model.richness} \
        --include-chromatin false --include-g4-architecture false \
        --stratify-g4-density false --exclude-cross-reactive false \
        --gc-flank-bp ${params.sequence_flank_bp} \
        --outdir ${model.directory} --prefix ${model.prefix}
    """
}

process VERIFY_PAPER_RESULTS {
    tag 'six paper tables'
    label 'process_low'
    publishDir "${params.outdir}/verification", mode: 'copy', overwrite: true

    input:
    path result_directories
    path expected_tables, stageAs: 'reference'
    path comparison

    output:
    path 'numerical_comparison.csv', emit: report

    script:
    """
    Rscript -e 'library(data.table); source("${comparison}"); verify_paper_results(".", "reference", "numerical_comparison.csv")'
    """
}