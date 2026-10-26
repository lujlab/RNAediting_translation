#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

SAMPLE_LIST="${SAMPLE_LIST:-${script_dir}/pausepred_sample_list_Dsec_merged.tsv}"
ALLED_IT_CDS="${ALLED_IT_CDS:-/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed/sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv}"
PAUSEPRED_DIR="${PAUSEPRED_DIR:-/data/dongxn/RNAediting/codex/02_stalling/2pausepred/pausepred_results_new}"
OUT_DIR="${OUT_DIR:-${script_dir}/pausepred_enrichment_Dsec_merged_alledit_cds}"
CDS_INTERVALS="${CDS_INTERVALS:-/path/to/RNA_editing_DATAprocessing/05.psite/dsec_all_r1.3_with_UTR_CDSstop_UTR_position_CDS}"
MIN_READS="${MIN_READS:-2}"
MIN_COVERAGE="${MIN_COVERAGE:-20}"
PAUSE_CUTOFFS="${PAUSE_CUTOFFS:-5,10}"
WINDOWS="${WINDOWS:-0,2}"
MAX_JOBS="${MAX_JOBS:-2}"
PYTHON_SCRIPT="${PYTHON_SCRIPT:-${script_dir}/pausepred_alledit_cds_to_full_transcript.py}"

for required_file in "$SAMPLE_LIST" "$ALLED_IT_CDS" "$CDS_INTERVALS" "$PYTHON_SCRIPT"; do
    if [[ ! -f "$required_file" ]]; then
        echo "ERROR: required file not found: $required_file" >&2
        exit 1
    fi
done
if [[ ! -d "$PAUSEPRED_DIR" ]]; then
    echo "ERROR: PausePred directory not found: $PAUSEPRED_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "========== Dsec merged PausePred enrichment =========="
echo "SAMPLE_LIST:   $SAMPLE_LIST"
echo "ALLED_IT_CDS:  $ALLED_IT_CDS"
echo "PAUSEPRED_DIR: $PAUSEPRED_DIR"
echo "OUT_DIR:       $OUT_DIR"
echo "CDS_INTERVALS: $CDS_INTERVALS"
echo "MIN_READS:     $MIN_READS"
echo "MIN_COVERAGE:  $MIN_COVERAGE"
echo "PAUSE_CUTOFFS: $PAUSE_CUTOFFS"
echo "WINDOWS:       $WINDOWS nt"
echo "MAX_JOBS:      $MAX_JOBS"
echo "======================================================"

mapfile -t samples < <(tail -n +2 "$SAMPLE_LIST" | awk -F '\t' 'NF && $1 !~ /^#/ {print $1}')
if [[ "${#samples[@]}" -eq 0 ]]; then
    echo "ERROR: no sample names found in $SAMPLE_LIST" >&2
    exit 1
fi

for sample in "${samples[@]}"; do
    pausepred_file="${PAUSEPRED_DIR}/${sample}.pausepred.csv"
    if [[ ! -f "$pausepred_file" ]]; then
        echo "ERROR: expected PausePred file not found: $pausepred_file" >&2
        exit 1
    fi
done

process_sample() {
    local sample="$1"
    local pausepred_file="${PAUSEPRED_DIR}/${sample}.pausepred.csv"
    local output_file="${OUT_DIR}/${sample}.pausepred_enrichment.tsv"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing $sample"
    python3 "$PYTHON_SCRIPT" \
        --pausepred "$pausepred_file" \
        --alledit-cds "$ALLED_IT_CDS" \
        --cds-intervals "$CDS_INTERVALS" \
        --out "$output_file" \
        --min-reads "$MIN_READS" \
        --min-coverage "$MIN_COVERAGE" \
        --pause-cutoffs "$PAUSE_CUTOFFS" \
        --windows "$WINDOWS" \
        --write-intermediate
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed $sample"
}
export -f process_sample
export PAUSEPRED_DIR ALLED_IT_CDS OUT_DIR CDS_INTERVALS MIN_READS MIN_COVERAGE \
    PAUSE_CUTOFFS WINDOWS PYTHON_SCRIPT

printf '%s\n' "${samples[@]}" | xargs -r -n 1 -P "$MAX_JOBS" bash -c 'process_sample "$1"' _
echo "All merged Dsec samples processed. Results: $OUT_DIR"
