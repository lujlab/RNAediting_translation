#!/usr/bin/env bash
# Run the rigorous flank30 edgeR analysis with sample-specific RNA scaling.
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
analysis_root="$script_dir"
source_root=$(cd -- "$analysis_root/.." && pwd)
publication_root=$(cd -- "$script_dir/../../../.." && pwd)

rscript_bin=${RSCRIPT_BIN:-Rscript}
master_manifest=${MASTER_MANIFEST:-$source_root/stalling_master_tables/master_table_manifest.tsv}
result_dir=${EDGER_V2_DIR:-$analysis_root/results}
summary_dir="$result_dir/run_summaries"
jobs=${EDGER_JOBS:-1}
r_threads=${R_THREADS:-4}
sample_fraction=${SAMPLE_FRACTION:-1}
sample_seed=${SAMPLE_SEED:-20260805}
run_label=${RUN_LABEL:-peak_only_CDS150_m6A_20261002}
normalization_method=${NORMALIZATION_METHOD:-TMM}
min_stalling_log2fc=${MIN_STALLING_LOG2FC:-0}
min_transcript_median_rna=${MIN_TRANSCRIPT_MEDIAN_RNA:-0}
min_transcript_median_rpf=${MIN_TRANSCRIPT_MEDIAN_RPF:-0}
transcript_median_rule=${TRANSCRIPT_MEDIAN_RULE:-any_both}
project_r_library=${PROJECT_R_LIBRARY:-$source_root/R_library}
groups_csv=${GROUP_IDS:-${GROUPS_CSV:-}}
exclude_head_replicates=${EXCLUDE_HEAD_REPLICATES:-0}
rna_scale_manifest=${RNA_SCALE_MANIFEST:-$publication_root/config/rna_count_scaling.tsv}

if [[ -d "$project_r_library" ]]; then
    export R_LIBS_USER="$project_r_library${R_LIBS_USER:+:$R_LIBS_USER}"
fi
if [[ "$jobs" -lt 1 || "$r_threads" -lt 1 ]]; then
    echo "EDGER_JOBS and R_THREADS must be positive" >&2
    exit 2
fi
if [[ ! -s "$master_manifest" ]]; then
    echo "Missing master manifest: $master_manifest" >&2
    exit 2
fi
if [[ ! -s "$rna_scale_manifest" ]]; then
    echo "Missing RNA scale manifest: $rna_scale_manifest" >&2
    exit 2
fi

mkdir -p "$summary_dir"
task_file=$(mktemp)
trap 'rm -f "$task_file"' EXIT

if [[ -n "$groups_csv" ]]; then
    IFS=',' read -r -a groups_array <<< "$groups_csv"
    printf '%s\n' "${groups_array[@]}" > "$task_file"
else
    # ISO1 has one replicate and is excluded. Head samples remain replicates.
    awk -F '\t' 'NR > 1 && $6 >= 2 && !seen[$1]++ { print $1 }' \
        "$master_manifest" > "$task_file"
fi

if [[ ! -s "$task_file" ]]; then
    echo "No groups selected" >&2
    exit 2
fi

run_group() {
    local group_id=$1
    local summary_file="$summary_dir/${group_id}.tsv"
    local args=(
        "$script_dir/edgeR_peak_only.R"
        --master-manifest "$master_manifest"
        --output-dir "$result_dir"
        --rna-scale-manifest "$rna_scale_manifest"
        --focal-reference-manifest "${FOCAL_REFERENCE_MANIFEST:?FOCAL_REFERENCE_MANIFEST is required}"
        --groups "$group_id"
        --flank-sizes "30"
        --min-focal-rna-per-rep "5"
        --min-flank-rna-per-rep "100"
        --min-flank-n-per-rep "10"
        --min-focal-rpf-sum "0"
        --normalization-method "$normalization_method"
        --fdr "0.05"
        --min-stalling-log2fc "$min_stalling_log2fc"
        --min-transcript-median-rna "$min_transcript_median_rna"
        --min-transcript-median-rpf "$min_transcript_median_rpf"
        --transcript-median-rule "$transcript_median_rule"
        --sample-fraction "$sample_fraction"
        --sample-seed "$sample_seed"
        --threads "$r_threads"
        --run-label "$run_label"
        --summary-file "$summary_file"
    )
    if [[ ${OVERWRITE_EDGER:-0} == 1 ]]; then
        args+=(--overwrite)
    fi
    if [[ "$exclude_head_replicates" == 1 ]]; then
        args+=(--exclude-head-replicates)
    fi
    echo "[start] $group_id" >&2
    "$rscript_bin" "${args[@]}"
    echo "[done] $group_id" >&2
}

export -f run_group
export rscript_bin script_dir master_manifest result_dir summary_dir
export rna_scale_manifest FOCAL_REFERENCE_MANIFEST
export min_stalling_log2fc OVERWRITE_EDGER
export r_threads sample_fraction sample_seed run_label normalization_method
export min_transcript_median_rna min_transcript_median_rpf transcript_median_rule
export exclude_head_replicates

xargs -a "$task_file" -n 1 -P "$jobs" bash -c 'run_group "$1"' _

if compgen -G "$summary_dir/*.tsv" > /dev/null; then
    "$rscript_bin" "$script_dir/aggregate_peaks.R" \
        "$result_dir" "$master_manifest"
fi
