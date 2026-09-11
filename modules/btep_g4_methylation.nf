process DOWNLOAD_GEO_SERIES_MATRIX {
    tag "${meta.cohort}_${meta.tissue}"
    label 'process_high'

    publishDir "${params.outdir}/prepared_inputs", mode: 'copy'

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta),
          path("${meta.cohort}_${meta.tissue}_beta.rds"),
          path("${meta.cohort}_${meta.tissue}_metadata.csv"),
          emit: prepared

    script:
    def g4_architecture_arg = params.include_g4_architecture ?: false
    """
    Rscript ${projectDir}/bin/download_geo_series_matrix.R \
        --accession ${accession} \
        --tissue '${meta.tissue}' \
        --prefix ${meta.cohort}_${meta.tissue} \
        --outdir .
    """
}

process PREPARE_G4_WINDOWS {
    tag 'g4_100bp_windows'
    label 'process_medium'

    publishDir "${params.outdir}/g4_windows", mode: 'copy'

    input:
    path stable_bed
    path unstable_bed

    output:
    path 'g4_motifs_core.bed', emit: core_motifs
    path 'g4_motifs_100bp_unmerged.bed', emit: unmerged_windows
    path 'g4_motifs_100bp_merged.bed', emit: windows
    path 'g4_motifs_100bp_long_clusters.bed', emit: long_clusters
    path 'g4_window_overlap_summary.csv', emit: overlap_summary

    script:
    """
    Rscript ${projectDir}/bin/prepare_g4_windows.R \
        --stable ${stable_bed} \
        --unstable ${unstable_bed} \
        --flank ${params.g4_flank_bp} \
        --long-quantile ${params.g4_long_quantile} \
        --long-mad-multiplier ${params.g4_long_mad_multiplier} \
        --long-min-bp ${params.g4_long_min_bp} \
        --outdir .
    """
}

process BUILD_PROBE_STRUCTURE_ANNOTATIONS {
    tag 'probe_structure_annotations'
    label 'process_medium'

    publishDir "${params.outdir}/sensitivity_annotations", mode: 'copy'

    input:
    path probe_manifest

    output:
    path 'probe_structure_annotations.csv.gz', emit: annotations
    path 'probe_structure_annotation_provenance.csv', emit: provenance
    path 'probe_structure_annotation_qc.csv', emit: qc

    script:
    def rloop_arg = params.rloop_bed ?: ''
    def goldenpath_arg = params.ucsc_goldenpath_dir ?: '/fdb/genomebrowser/goldenPath/hg38'
    def rloop_url_arg = params.rloop_url ?: ''
    def timing_arg = params.replication_timing_bed ?: ''
    def timing_url_arg = params.replication_timing_url ?: ''
    def cross_arg = params.cross_reactive_probes ?: ''
    def cross_url_arg = params.cross_reactive_probes_url ?: ''
    """
    Rscript ${projectDir}/bin/build_probe_structure_annotations.R \
        --manifest ${probe_manifest} \
        --ucsc-goldenpath-dir '${goldenpath_arg}' \
        --rloop-bed '${rloop_arg}' \
        --rloop-url '${rloop_url_arg}' \
        --replication-timing-bed '${timing_arg}' \
        --replication-timing-url '${timing_url_arg}' \
        --cross-reactive-probes '${cross_arg}' \
        --cross-reactive-probes-url '${cross_url_arg}' \
        --outdir .
    """
}

process ANALYZE_G4_METHYLATION_AGE {
    tag "${meta.cohort}_${meta.tissue}"
    label 'process_high'

    publishDir "${params.outdir}/per_tissue", mode: 'copy'

    input:
    tuple val(meta), path(beta_matrix), path(metadata), path(probe_manifest), path(g4_windows)

    output:
    tuple val(meta), path("${meta.cohort}_${meta.tissue}_probe_slopes.csv.gz"), emit: slopes
    tuple val(meta), path("${meta.cohort}_${meta.tissue}_g4_age_summary.csv"), emit: summary
    tuple val(meta), path("${meta.cohort}_${meta.tissue}_g4_age_tests.csv"), emit: tests
    tuple val(meta), path("${meta.cohort}_${meta.tissue}_g4_regression_effects.csv"), emit: regression_effects

    script:
    """
    Rscript ${projectDir}/bin/analyze_g4_methylation_age.R \
        --beta ${beta_matrix} \
        --metadata ${metadata} \
        --manifest ${probe_manifest} \
        --g4-windows ${g4_windows} \
        --age-col '${meta.age_col}' \
        --sample-col '${meta.sample_id_col}' \
        --covariates '${meta.covariates}' \
        --min-samples ${params.min_samples} \
        --cohort '${meta.cohort}' \
        --tissue '${meta.tissue}' \
        --prefix ${meta.cohort}_${meta.tissue} \
        --outdir .
    """
}

process RUN_STRUCTURE_SENSITIVITY {
    tag 'structure_sensitivity'
    label 'process_high'

    publishDir "${params.outdir}/structure_sensitivity", mode: 'copy'

    input:
    path slope_files
    path g4_unmerged
    path structure_annotations

    output:
    path 'structure_sensitivity_g4_terms.csv'
    path 'structure_sensitivity_stable_minus_unstable.csv'
    path 'structure_sensitivity_g4_terms.png'
    path 'structure_sensitivity_g4_terms.pdf'

    script:
    def g4_architecture_arg = params.include_g4_architecture ?: false
    """
    Rscript ${projectDir}/bin/compare_g4_effects_adjusted.R \
        --slopes-dir . \
        --g4-unmerged ${g4_unmerged} \
        --structure-annotations ${structure_annotations} \
        --exclude-cross-reactive true \
        --include-g4-architecture ${g4_architecture_arg} \
        --outdir . \
        --prefix structure_sensitivity
    """
}

process SUMMARIZE_G4_TISSUES {
    tag 'cross_tissue_summary'
    label 'process_low'

    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path summary_files
    path test_files
    path regression_files
    path overlap_summary

    output:
    path 'btep_g4_cross_tissue_summary.csv', emit: summary
    path 'btep_g4_cross_tissue_tests.csv', emit: tests
    path 'btep_g4_cross_tissue_regression_effects.csv', emit: regression_effects
    path 'btep_g4_experiment_readme.md', emit: readme

    script:
    def summaries_arg = summary_files.collect { it.toString() }.join(',')
    def tests_arg = test_files.collect { it.toString() }.join(',')
    def regressions_arg = regression_files.collect { it.toString() }.join(',')
    """
    Rscript ${projectDir}/bin/summarize_g4_tissues.R \
        --summaries '${summaries_arg}' \
        --tests '${tests_arg}' \
        --regressions '${regressions_arg}' \
        --overlap-summary ${overlap_summary} \
        --outdir .
    """
}