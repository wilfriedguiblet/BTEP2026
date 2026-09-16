#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * The single computational entry point for the Quadron methylation paper.
 * prepared: beta/metadata -> CpG age slopes -> M0/M1 -> numerical verification.
 * slopes:   replay the two models from existing annotated CpG slopes.
 * geo:      explicitly import processed GEO matrices before fitting slopes.
 * Rendering is separate: knit G4_methylation_paper.Rmd after the run succeeds.
 */

include {
    IMPORT_GEO
    PREPARE_G4_WINDOWS
    FIT_AGE_SLOPES
    FIT_PAPER_MODEL
    VERIFY_PAPER_RESULTS
} from './modules/btep_g4_methylation'

def requiredFile(value, label) {
    if (!value) error "Missing --${label}"
    def candidate = file(value)
    if (!candidate.exists()) error "${label} not found: ${value}"
    return candidate
}

// Resolve existing symlinks before checking destinations. A missing ancestor
// followed by '..' is rejected: creating it later could redirect output to data.
def physicalPath(value) {
    def candidate = file(value).toFile()
    if (candidate.exists()) return candidate.toPath().toRealPath()
    if (java.nio.file.Files.isSymbolicLink(candidate.toPath())) error "Dangling symlink: ${value}"
    if (candidate.name in ['.', '..']) error "Unresolved parent traversal: ${value}"
    return physicalPath(candidate.parent).resolve(candidate.name)
}

workflow {
    if (params.help) {
        log.info '''
Quadron methylation paper

  nextflow run main.nf --start_from prepared
  nextflow run main.nf --start_from slopes --outdir results/paper_slopes
  nextflow run main.nf --start_from geo --outdir results/paper_geo
  nextflow run main.nf -profile slurm -resume

Inputs (defaults in nextflow.config):
  --data_dir PATH          Prepared inputs, annotated slopes and G4 checkpoints
  --probe_manifest PATH    The shared hg38 probe manifest (prepared/geo routes)
  --g4_stable PATH         Optional original Quadron stable BED
  --g4_unstable PATH       Optional original Quadron unstable BED; supply both
  --series_matrix_dir PATH Optional frozen GEO matrix files for the geo route

Settings:
  --sequence_flank_bp 1000 Probe-centered sequence flank; G4 padding is fixed at 100 bp
  --outdir PATH           New output directory; use -resume to continue a managed run
  --verify false         Skip comparison with the six frozen paper tables

Only the two additive minimal models are run. No chromatin, density,
RepeatMasker, or older regression families are included. The geo route may
download processed matrices; it does not reproduce raw-IDAT normalization.
'''.stripIndent()
    } else {
        if (!(params.start_from in ['prepared', 'slopes', 'geo'])) error 'start_from must be prepared, slopes, or geo'
        if (!(params.sequence_flank_bp.toString() ==~ /\d+/)) error 'sequence_flank_bp must be a nonnegative integer'
        if ((params.g4_stable as Boolean) != (params.g4_unstable as Boolean)) error 'Supply both Quadron class BEDs or neither'

        def destination = physicalPath(params.outdir)
        def protectedPaths = [params.data_dir, params.expected_tables, "${projectDir}/bin",
                              "${projectDir}/modules", "${projectDir}/G4_methylation_paper_files"]
        protectedPaths.each { value ->
            def protectedPath = physicalPath(value)
            if (destination.startsWith(protectedPath) || protectedPath.startsWith(destination)) {
                error "Output overlaps an input or source directory: ${destination}"
            }
        }
        def marker = destination.resolve('.quadron-paper-run').toFile()
        def directory = destination.toFile()
        if (directory.exists() && directory.list()?.size() && !(workflow.resume && marker.exists())) {
            error "Output is not empty: ${destination}. Choose a new --outdir or -resume a managed run."
        }
        directory.mkdirs()
        if (!marker.exists()) marker.text = 'Managed by BTEP/main.nf\n'
        params.outdir = destination.toString()

        def cohorts = [[cohort: 'GSE61257', tissue: 'adipose'],
                       [cohort: 'GSE61258', tissue: 'liver'],
                       [cohort: 'GSE61259', tissue: 'muscle']]
        def merged
        def unmerged
        if (params.g4_stable) {
            PREPARE_G4_WINDOWS(channel.value(requiredFile(params.g4_stable, 'g4_stable')),
                               channel.value(requiredFile(params.g4_unstable, 'g4_unstable')),
                               channel.value(requiredFile("${projectDir}/bin/prepare_g4_windows.R", 'window script')))
            merged = PREPARE_G4_WINDOWS.out.merged
            unmerged = PREPARE_G4_WINDOWS.out.unmerged
        } else {
            unmerged = channel.value(requiredFile("${params.data_dir}/g4_windows/g4_motifs_100bp_unmerged.bed", 'unmerged Quadron checkpoint'))
            if (params.start_from != 'slopes') {
                merged = channel.value(requiredFile("${params.data_dir}/g4_windows/g4_motifs_100bp_merged.bed", 'merged Quadron checkpoint'))
            }
        }

        def slopes
        if (params.start_from == 'slopes') {
            slopes = channel.fromList(cohorts).map { meta ->
                requiredFile("${params.data_dir}/per_tissue/${meta.cohort}_${meta.tissue}_probe_slopes.csv.gz", 'annotated slope checkpoint')
            }
        } else {
            def prepared
            if (params.start_from == 'geo') {
                def imports = channel.fromList(cohorts).map { meta ->
                    def matrix = params.series_matrix_dir ? requiredFile("${params.series_matrix_dir}/${meta.cohort}_series_matrix.txt.gz", 'series matrix') : []
                    tuple(meta, matrix)
                }
                IMPORT_GEO(imports, channel.value(requiredFile("${projectDir}/bin/download_geo_series_matrix.R", 'GEO script')))
                prepared = IMPORT_GEO.out.prepared
            } else {
                prepared = channel.fromList(cohorts).map { meta ->
                    tuple(meta,
                          requiredFile("${params.data_dir}/prepared_inputs/${meta.cohort}_${meta.tissue}_beta.rds", 'beta matrix'),
                          requiredFile("${params.data_dir}/prepared_inputs/${meta.cohort}_${meta.tissue}_metadata.csv", 'metadata'))
                }
            }
            FIT_AGE_SLOPES(prepared, channel.value(requiredFile(params.probe_manifest, 'probe_manifest')),
                           merged, channel.value(requiredFile("${projectDir}/bin/analyze_g4_methylation_age.R", 'age-slope script')))
            slopes = FIT_AGE_SLOPES.out.slopes.map { meta, path -> path }
        }

        // Stable/unstable are overlapping ADDITIVE predictors. Only the two
        // strand-oriented G fractions distinguish M1 from M0; no other model
        // branches are invoked. Sort the collected inputs for repeatable replay.
        def models = channel.of(
            [id: 'M0', directory: 'minimal_model', prefix: 'quadron_minimal_stable_unstable_mean_beta', richness: false],
            [id: 'M1', directory: 'minimal_model_strand_g_richness', prefix: 'quadron_minimal_strand_g_richness', richness: true]
        )
        FIT_PAPER_MODEL(models, slopes.collect().map { paths -> paths.sort { it.name } },
                         unmerged, channel.value(requiredFile("${projectDir}/bin/compare_g4_effects_adjusted.R", 'model script')))

        if (params.verify) {
            VERIFY_PAPER_RESULTS(FIT_PAPER_MODEL.out.results.map { model, path -> path }.collect(),
                                  channel.value(requiredFile(params.expected_tables, 'expected_tables')),
                                  channel.value(requiredFile("${projectDir}/bin/verify_paper_results.R", 'comparison script')))
        }
    }
}