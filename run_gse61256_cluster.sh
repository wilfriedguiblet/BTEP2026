#!/bin/bash
#SBATCH --job-name=btep_g4
#SBATCH --output=btep_g4_%j.out
#SBATCH --error=btep_g4_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=2

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARAMS_INPUT="${1:-${SCRIPT_DIR}/params_GSE61256.yaml}"
PARAMS="$(cd "$(dirname "${PARAMS_INPUT}")" && pwd)/$(basename "${PARAMS_INPUT}")"
PROFILE="${2:-slurm}"
EXTRA_ARGS=("${@:3}")

module load nextflow/24.04
module load R/4.3

Rscript -e '
required <- c("data.table", "optparse")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
    stop("Missing required cluster R packages: ", paste(missing, collapse = ", "))
}
' >/dev/null

mkdir -p "${SCRIPT_DIR}/reports"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

echo "============================================"
echo "BTEP G4 Methylation Aging Pipeline"
echo "  Nextflow: $(nextflow -version 2>&1 | head -n 1)"
echo "  Params:   ${PARAMS}"
echo "  Profile:  ${PROFILE}"
echo "  Work dir: ${SCRIPT_DIR}/work"
echo "============================================"

cd "${SCRIPT_DIR}"

nextflow -C "${SCRIPT_DIR}/nextflow.config" \
    run "${SCRIPT_DIR}/main.nf" \
    -params-file "${PARAMS}" \
    -profile "${PROFILE}" \
    -with-report "${SCRIPT_DIR}/reports/btep_g4_report_${TIMESTAMP}.html" \
    -with-timeline "${SCRIPT_DIR}/reports/btep_g4_timeline_${TIMESTAMP}.html" \
    -with-trace "${SCRIPT_DIR}/reports/btep_g4_trace_${TIMESTAMP}.txt" \
    -with-dag "${SCRIPT_DIR}/reports/btep_g4_dag_${TIMESTAMP}.html" \
    "${EXTRA_ARGS[@]}" \
    -resume
