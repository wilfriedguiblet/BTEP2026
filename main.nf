#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * BTEP G4 methylation-aging workflow
 *
 * Tests whether human CpG methylation age-slopes are larger within annotated
 * G4 motifs plus a 100 bp surrounding window. The default input is the public
 * GSE61256 methylation subseries: adipose, liver and muscle.
 */

include {
    DOWNLOAD_GEO_SERIES_MATRIX
    PREPARE_G4_WINDOWS
    BUILD_PROBE_STRUCTURE_ANNOTATIONS
    ANALYZE_G4_METHYLATION_AGE
    RUN_STRUCTURE_SENSITIVITY
    SUMMARIZE_G4_TISSUES
} from './modules/btep_g4_methylation'

def boolParam(value) {
    return value instanceof Boolean ? value : value.toString().toBoolean()
}

def requiredFile(value, label) {
    if (!value) {
        error "ERROR: --${label} is required"
    }
    def candidate = file(value)
    if (!candidate.exists()) {
        error "ERROR: ${label} not found: ${value}"
    }
    return candidate
}

def helpText() {
    return """
BTEP G4 methylation-aging workflow

Default experiment:
  Downloads processed GEO Series Matrix files for GSE61257, GSE61258 and
  GSE61259, representing adipose, liver and muscle Illumina 450K methylation.

Required reference inputs:
  --probe_manifest PATH     Probe coordinates with probe_id, chr_hg38/start_hg38/end_hg38
  --g4_stable PATH          G4Hunter stable BED, hg38
  --g4_unstable PATH        G4Hunter unstable BED, hg38

Optional local-input mode:
  --samplesheet PATH        CSV with columns cohort,tissue,beta_matrix,metadata,age_col,sample_id_col,covariates

Example:
  nextflow run main.nf -profile standard
  nextflow run main.nf --samplesheet samplesheet.template.csv
""".stripIndent()
}

params.help = false
params.analysis_id = 'gse61256_g4_100bp'
params.outdir = 'results/g4_methylation_age_100bp'
params.samplesheet = null
params.geo_series = 'GSE61257:adipose,GSE61258:liver,GSE61259:muscle'
params.g4_stable = '../G4Hunter/G4Hunter.hg38.stable.bed'
params.g4_unstable = '../G4Hunter/G4Hunter.hg38.unstable.bed'
params.probe_manifest = '../results/submission_v3/reference_manifest/probe_gene_manifest.csv.gz'
params.g4_flank_bp = 100
params.g4_long_quantile = 0.99
params.g4_long_mad_multiplier = 3
params.g4_long_min_bp = 0
params.ucsc_goldenpath_dir = '/fdb/genomebrowser/goldenPath/hg38'
params.rloop_bed = ''
params.rloop_url = ''
params.replication_timing_bed = ''
params.replication_timing_url = ''
params.cross_reactive_probes = ''
params.cross_reactive_probes_url = ''
params.include_g4_architecture = false
params.min_samples = 2
params.default_covariates = 'sex|bmi'

workflow {
    if (boolParam(params.help)) {
        log.info helpText()
    } else {
        def stable_g4 = requiredFile(params.g4_stable, 'g4_stable')
        def unstable_g4 = requiredFile(params.g4_unstable, 'g4_unstable')
        def probe_manifest = requiredFile(params.probe_manifest, 'probe_manifest')

        log.info """
        ╔══════════════════════════════════════════════════════════════╗
        ║             BTEP G4 Methylation Aging Pipeline              ║
        ╠══════════════════════════════════════════════════════════════╣
        ║ Output:          ${params.outdir}
        ║ G4 flank:        ${params.g4_flank_bp} bp
        ║ Long clusters:   >= max(q${params.g4_long_quantile}, median + ${params.g4_long_mad_multiplier} MAD, ${params.g4_long_min_bp} bp)
        ║ Stable BED:      ${stable_g4}
        ║ Unstable BED:    ${unstable_g4}
        ║ Probe manifest:  ${probe_manifest}
        ║ Samplesheet:     ${params.samplesheet ?: 'GEO defaults: ' + params.geo_series}
        ╚══════════════════════════════════════════════════════════════╝
        """.stripIndent()

        PREPARE_G4_WINDOWS(
            channel.value(stable_g4),
            channel.value(unstable_g4)
        )
        BUILD_PROBE_STRUCTURE_ANNOTATIONS(channel.value(probe_manifest))

        def ch_samples
        if (params.samplesheet) {
            def samplesheet_file = requiredFile(params.samplesheet, 'samplesheet')
            def sampleLines = new File(samplesheet_file.toString())
                .readLines('UTF-8')
                .findAll { line -> line.trim() }
            if (sampleLines.size() < 2) {
                error 'ERROR: samplesheet must contain a header and at least one sample row'
            }
            def headers = sampleLines[0].split(',', -1) as List
            def requiredColumns = ['cohort', 'tissue', 'beta_matrix', 'metadata', 'age_col', 'sample_id_col']
            def missingHeaders = requiredColumns - headers
            if (missingHeaders) {
                error "ERROR: samplesheet missing columns: ${missingHeaders.join(', ')}"
            }
            def sampleTuples = []
            sampleLines.drop(1).eachWithIndex { line, rowIndex ->
                def values = line.split(',', -1) as List
                if (values.size() != headers.size()) {
                    error "ERROR: samplesheet row ${rowIndex + 2} has ${values.size()} fields; expected ${headers.size()}"
                }
                def row = [headers, values].transpose().collectEntries()
                requiredColumns.each { column ->
                    if (!row[column]?.toString()?.trim()) {
                        error "ERROR: samplesheet row ${rowIndex + 2} is missing '${column}'"
                    }
                }
                def beta = file(row.beta_matrix)
                def metadata = file(row.metadata)
                if (!beta.exists()) error "ERROR: beta_matrix not found: ${row.beta_matrix}"
                if (!metadata.exists()) error "ERROR: metadata not found: ${row.metadata}"
                def meta = [
                    cohort: row.cohort.toString(),
                    tissue: row.tissue.toString(),
                    age_col: row.age_col.toString(),
                    sample_id_col: row.sample_id_col.toString(),
                    covariates: row.covariates?.toString()?.trim() ?: params.default_covariates.toString()
                ]
                sampleTuples << tuple(meta, beta, metadata)
            }
            ch_samples = channel.fromList(sampleTuples)
        } else {
            def geoTuples = params.geo_series.toString().split(',').collect { item ->
                def fields = item.split(':', -1)
                if (fields.size() != 2 || !fields[0] || !fields[1]) {
                    error "ERROR: malformed --geo_series item '${item}'. Use ACCESSION:tissue"
                }
                tuple([
                    cohort: fields[0].toString(),
                    tissue: fields[1].toString(),
                    age_col: 'age',
                    sample_id_col: 'sample_id',
                    covariates: params.default_covariates.toString()
                ], fields[0].toString())
            }
            DOWNLOAD_GEO_SERIES_MATRIX(channel.fromList(geoTuples))
            ch_samples = DOWNLOAD_GEO_SERIES_MATRIX.out.prepared
        }

        ch_analysis_inputs = ch_samples
            .combine(channel.value(probe_manifest))
            .combine(PREPARE_G4_WINDOWS.out.windows)
            .map { meta, beta, metadata, manifest, windows ->
                tuple(meta, beta, metadata, manifest, windows)
            }

        ANALYZE_G4_METHYLATION_AGE(ch_analysis_inputs)

        RUN_STRUCTURE_SENSITIVITY(
            ANALYZE_G4_METHYLATION_AGE.out.slopes.map { meta, slopes -> slopes }.collect(),
            PREPARE_G4_WINDOWS.out.unmerged_windows,
            BUILD_PROBE_STRUCTURE_ANNOTATIONS.out.annotations
        )

        SUMMARIZE_G4_TISSUES(
            ANALYZE_G4_METHYLATION_AGE.out.summary.map { meta, summary -> summary }.collect(),
            ANALYZE_G4_METHYLATION_AGE.out.tests.map { meta, tests -> tests }.collect(),
            ANALYZE_G4_METHYLATION_AGE.out.regression_effects.map { meta, regression -> regression }.collect(),
            PREPARE_G4_WINDOWS.out.overlap_summary
        )
    }
}