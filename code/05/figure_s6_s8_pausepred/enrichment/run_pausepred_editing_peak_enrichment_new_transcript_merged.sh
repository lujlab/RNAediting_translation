#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

SAMPLE_LIST="${SAMPLE_LIST:?Set SAMPLE_LIST to a one-column sample list}"
EDIT_SITES_FILE="${EDIT_SITES_FILE:?Set EDIT_SITES_FILE to the mapped editing-site TSV}"
PAUSEPRED_DIR="${PAUSEPRED_DIR:-/data/dongxn/RNAediting/codex/02_stalling/2pausepred/pausepred_results_new}"
OUT_DIR="${OUT_DIR:-${PAUSEPRED_DIR}/enrichment_alledit_cds}"
MIN_READS="${MIN_READS:-2}"
MIN_COVERAGE="${MIN_COVERAGE:-20}"
PAUSE_CUTOFFS="${PAUSE_CUTOFFS:-5,10}"
WINDOWS="${WINDOWS:-0,2}"
MAX_JOBS="${MAX_JOBS:-2}"
TARGET_REGION="${TARGET_REGION:-all}"
WRITE_INTERMEDIATE="${WRITE_INTERMEDIATE:-true}"
SCRIPT="${SCRIPT:-${script_dir}/pausepred_editing_peak_enrichment_new_transcript_merged.py}"

for required_file in "$SAMPLE_LIST" "$EDIT_SITES_FILE" "$SCRIPT"; do
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

mapfile -t samples < <(tail -n +2 "$SAMPLE_LIST" | awk -F '\t' 'NF && $1 !~ /^#/ {gsub(/\r/,"",$1); print $1}')
if [[ "${#samples[@]}" -eq 0 ]]; then
    echo "ERROR: no samples found in $SAMPLE_LIST" >&2
    exit 1
fi

for sample in "${samples[@]}"; do
    f="${PAUSEPRED_DIR}/${sample}.pausepred.csv"
    if [[ ! -f "$f" ]]; then
        echo "ERROR: expected PausePred file not found: $f" >&2
        exit 1
    fi
done

process_sample() {
    local sample="$1"
    local pausepred_file="${PAUSEPRED_DIR}/${sample}.pausepred.csv"
    local output_file="${OUT_DIR}/${sample}.pausepred_enrichment.tsv"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing $sample"
    local cmd=(
        python3 "$SCRIPT"
        --pausepred "$pausepred_file"
        --edit-sites "$EDIT_SITES_FILE"
        --target-region "$TARGET_REGION"
        --out "$output_file"
        --min-reads "$MIN_READS"
        --min-coverage "$MIN_COVERAGE"
        --pause-cutoffs "$PAUSE_CUTOFFS"
        --windows "$WINDOWS"
    )
    case "${WRITE_INTERMEDIATE,,}" in
        true | 1 | yes) cmd+=(--write-intermediate) ;;
        false | 0 | no) ;;
        *)
            echo "ERROR: WRITE_INTERMEDIATE must be true/false" >&2
            return 2
            ;;
    esac
    "${cmd[@]}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed $sample"
}
export -f process_sample
export PAUSEPRED_DIR EDIT_SITES_FILE OUT_DIR MIN_READS MIN_COVERAGE \
    PAUSE_CUTOFFS WINDOWS TARGET_REGION WRITE_INTERMEDIATE SCRIPT

printf '%s\n' "${samples[@]}" | xargs -r -n 1 -P "$MAX_JOBS" bash -c 'process_sample "$1"' _
echo "All samples processed. Results: $OUT_DIR"
