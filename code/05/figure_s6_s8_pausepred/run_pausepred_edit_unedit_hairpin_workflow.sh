#!/usr/bin/env bash
set -euo pipefail

# Reproducible post-PausePred workflow.
#
# It performs, in order:
#   1. all-edit CDS enrichment for every sample in the 4-species sample sheets;
#   2. unedited-A-base enrichment, retaining editing sites but restricting only
#      the unedited background to reference A bases;
#   3. merged tables and the standard odds-ratio/fraction PDF for each scope;
#   4. the corresponding CDS-coordinate hairpin analysis and two PDFs.
#
# The script waits for the safe PausePred runner to finish successfully unless
# --no-wait is supplied. Results are written to a run-tagged directory, so an
# updated sequence run does not overwrite the previous enrichment results.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$SCRIPT_DIR/../../.." && pwd)}"
PAUSEPRED_DIR="${PAUSEPRED_DIR:-${PIPELINE_DIR}/results/05_stalling/pausepred_results_new}"
SAMPLE_SHEETS_DIR="${SAMPLE_SHEETS_DIR:-${PIPELINE_DIR}/results/05_stalling/psite_2times_merged_bam/pausepred_sample_sheets}"
RESCUE_ROOT="${RESCUE_ROOT:-/path/to/Proj1_editing/codex/rescue_3R_18907675_1v1_pipeline_clean}"
IPOS_DIR="${IPOS_DIR:-${SCRIPT_DIR}/enrichment}"
HAIRPIN_DIR="${HAIRPIN_DIR:-/path/to/Proj1_editing/codex/02_stalling/3flank/hairpin}"

# These are the plotting parameters for the analysis.
WINDOWS="${WINDOWS:-2}"
OFFSET="${OFFSET:-5}"

# The enrichment/hairpin tables retain both 0-nt and the requested window by
# default. The plot selects --windows and --offset below.
PAUSE_CUTOFFS="${PAUSE_CUTOFFS:-5,8,10}"
ENRICH_WINDOWS="${ENRICH_WINDOWS:-}"
ENRICH_JOBS="${ENRICH_JOBS:-6}"
HAIRPIN_JOBS="${HAIRPIN_JOBS:-4}"
WAIT_FOR_PAUSEPRED="${WAIT_FOR_PAUSEPRED:-1}"
WAIT_SECONDS="${WAIT_SECONDS:-600}"
ALLOW_OUTPUT_ONLY="${ALLOW_OUTPUT_ONLY:-0}"
RUN_TAG="${RUN_TAG:-seq90_up90_down90}"

usage() {
    cat << 'EOF'
Usage:
  bash run_pausepred_edit_unedit_hairpin_workflow.sh [options]

Main options:
  --windows N       Plot window in nt (default: 2)
  --offset N        Pause_score cutoff used for the plot (default: 5)
  --run-tag TAG     Output tag (default: seq90_up90_down90)
  --enrich-jobs N   Parallel enrichment samples (default: 6)
  --hairpin-jobs N  Parallel hairpin samples (default: 4)
  --wait-seconds N  Status-check interval while PausePred is running (default: 600)
  --no-wait         Do not wait; require successful PausePred outputs immediately
  --allow-output-only
                    Accept non-empty PausePred CSVs without checking runner logs
  -h, --help        Show this help

Examples:
  bash run_pausepred_edit_unedit_hairpin_workflow.sh \
      --windows 2 --offset 5

  ENRICH_JOBS=3 HAIRPIN_JOBS=2 RUN_TAG=seq90_retry \
      bash run_pausepred_edit_unedit_hairpin_workflow.sh \
      --windows 2 --offset 5
EOF
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

require_file() {
    local path="$1"
    [[ -s "$path" ]] || die "Required file is missing or empty: $path"
}

is_number() {
    [[ "$1" =~ ^[0-9]+([.][0-9]+)?$ ]]
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --windows)
            [[ $# -ge 2 ]] || die "--windows requires a value"
            WINDOWS="$2"
            shift 2
            ;;
        --offset)
            [[ $# -ge 2 ]] || die "--offset requires a value"
            OFFSET="$2"
            shift 2
            ;;
        --run-tag)
            [[ $# -ge 2 ]] || die "--run-tag requires a value"
            RUN_TAG="$2"
            shift 2
            ;;
        --enrich-jobs)
            [[ $# -ge 2 ]] || die "--enrich-jobs requires a value"
            ENRICH_JOBS="$2"
            shift 2
            ;;
        --hairpin-jobs)
            [[ $# -ge 2 ]] || die "--hairpin-jobs requires a value"
            HAIRPIN_JOBS="$2"
            shift 2
            ;;
        --wait-seconds)
            [[ $# -ge 2 ]] || die "--wait-seconds requires a value"
            WAIT_SECONDS="$2"
            shift 2
            ;;
        --no-wait)
            WAIT_FOR_PAUSEPRED=0
            shift
            ;;
        --allow-output-only)
            ALLOW_OUTPUT_ONLY=1
            shift
            ;;
        -h | --help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1; use --help"
            ;;
    esac
done

is_number "$WINDOWS" || die "--windows must be numeric: $WINDOWS"
is_number "$OFFSET" || die "--offset must be numeric: $OFFSET"
[[ "$ENRICH_JOBS" =~ ^[1-9][0-9]*$ ]] || die "ENRICH_JOBS must be a positive integer"
[[ "$HAIRPIN_JOBS" =~ ^[1-9][0-9]*$ ]] || die "HAIRPIN_JOBS must be a positive integer"
[[ "$WAIT_SECONDS" =~ ^[1-9][0-9]*$ ]] || die "WAIT_SECONDS must be a positive integer"
[[ "$RUN_TAG" =~ ^[A-Za-z0-9_.-]+$ ]] || die "RUN_TAG contains unsafe characters: $RUN_TAG"
command -v python3 > /dev/null || die "python3 not found"
command -v awk > /dev/null || die "awk not found"
command -v pgrep > /dev/null || die "pgrep not found"
python3 -c 'import numpy, pandas, matplotlib, scipy' || die "Python needs numpy, pandas, matplotlib and scipy"

COMPLETION_PY="${SCRIPT_DIR}/check_pausepred_completion.py"
PAUSEPRED_QC_DIR="${PAUSEPRED_QC_DIR:-$(dirname "$PAUSEPRED_DIR")/qc}"

NORMAL_RUNNER="${IPOS_DIR}/run_pausepred_editing_peak_enrichment_new_transcript_merged.sh"
NORMAL_PY="${IPOS_DIR}/pausepred_editing_peak_enrichment_new_transcript_merged.py"
DSEC_RUNNER="${IPOS_DIR}/run_pausepred_enrichment_Dsec_merged_alledit_cds.sh"
DSEC_PY="${IPOS_DIR}/pausepred_alledit_cds_to_full_transcript.py"
A_BASE_MAP_PY="${IPOS_DIR}/build_pausepred_a_base_map.py"
A_BASE_STATS_PY="${IPOS_DIR}/pausepred_editing_peak_enrichment_unedit_A_base.py"
MERGE_PY="${IPOS_DIR}/merge_pausepred_enrichment_no_bh.py"
PLOT_PY="${IPOS_DIR}/02_03_04_plot_pausepred_combined_odds_ratio_fraction.py"
HAIRPIN_PY="${SCRIPT_DIR}/enrichment/pausepred_hairpin_edit_unedit_enrichment.py"
HAIRPIN_PLOT_PY="${SCRIPT_DIR}/enrichment/plot_pausepred_hairpin_combined.py"

for path in \
    "$NORMAL_RUNNER" "$NORMAL_PY" "$DSEC_RUNNER" "$DSEC_PY" \
    "$A_BASE_MAP_PY" "$A_BASE_STATS_PY" "$MERGE_PY" "$PLOT_PY" \
    "$HAIRPIN_PY" "$HAIRPIN_PLOT_PY" "$COMPLETION_PY"; do
    require_file "$path"
done
[[ -d "$PAUSEPRED_DIR" ]] || die "PausePred directory not found: $PAUSEPRED_DIR"
[[ -d "$SAMPLE_SHEETS_DIR" ]] || die "Sample-sheet directory not found: $SAMPLE_SHEETS_DIR"
[[ -d "$HAIRPIN_DIR" ]] || die "Hairpin directory not found: $HAIRPIN_DIR"

ALLEDIT_OUT="${PAUSEPRED_DIR}/enrichment_${RUN_TAG}_alledit_cds"
ABASE_OUT="${PAUSEPRED_DIR}/enrichment_${RUN_TAG}_unedit_A_base"
HAIRPIN_OUT="${PAUSEPRED_DIR}/hairpin_edit_unedit_enrichment_${RUN_TAG}"
WORK_DIR="${PAUSEPRED_DIR}/workflow_${RUN_TAG}"
mkdir -p "$WORK_DIR" "$ALLEDIT_OUT" "$ABASE_OUT" "$HAIRPIN_OUT"

if [[ -z "$ENRICH_WINDOWS" ]]; then
    if [[ "$WINDOWS" == "0" ]]; then
        ENRICH_WINDOWS="0"
    else
        ENRICH_WINDOWS="0,${WINDOWS}"
    fi
elif [[ "$WINDOWS" == "0" ]]; then
    ENRICH_WINDOWS="0"
elif [[ "$ENRICH_WINDOWS" == "0,0" ]]; then
    ENRICH_WINDOWS="0"
fi

echo "[CONFIG] PAUSEPRED_DIR=$PAUSEPRED_DIR"
echo "[CONFIG] RUN_TAG=$RUN_TAG"
echo "[CONFIG] PAUSE_CUTOFFS=$PAUSE_CUTOFFS"
echo "[CONFIG] ENRICH_WINDOWS=$ENRICH_WINDOWS"
echo "[CONFIG] plot windows=$WINDOWS; plot offset=$OFFSET"
echo "[CONFIG] ENRICH_JOBS=$ENRICH_JOBS; HAIRPIN_JOBS=$HAIRPIN_JOBS"
echo "[CONFIG] ALLEDIT_OUT=$ALLEDIT_OUT"
echo "[CONFIG] ABASE_OUT=$ABASE_OUT"
echo "[CONFIG] HAIRPIN_OUT=$HAIRPIN_OUT"

# The safe runner writes one generated sheet per species. Reading these sheets
# keeps this post-processing workflow synchronized with the exact sample set
# used by the PausePred run, excluding the configured analysis samples.
PAUSEPRED_EXCLUDE_SAMPLES="${PAUSEPRED_EXCLUDE_SAMPLES-AdarG_F,AdarG_M,AdarS_F,AdarS_M}"
echo "[CONFIG] excluded samples=$PAUSEPRED_EXCLUDE_SAMPLES"
declare -a ALL_SAMPLES=()
declare -A SAMPLE_SPECIES=()
declare -A SEEN_SAMPLE=()

read_samples_from_sheet() {
    local species="$1"
    local sheet="$2"
    require_file "$sheet"
    while IFS= read -r sample; do
        [[ -n "$sample" ]] || continue
        case ",$PAUSEPRED_EXCLUDE_SAMPLES," in
            *",$sample,"*) continue ;;
        esac
        if [[ -n "${SEEN_SAMPLE[$sample]+yes}" ]]; then
            continue
        fi
        SEEN_SAMPLE["$sample"]=1
        SAMPLE_SPECIES["$sample"]="$species"
        ALL_SAMPLES+=("$sample")
    done < <(awk -F '\t' 'NR > 1 && NF && $1 !~ /^#/ {gsub(/\r/, "", $1); print $1}' "$sheet")
}

read_samples_from_sheet Dmel "${SAMPLE_SHEETS_DIR}/pausepred_sample_list.dmel.tsv"
read_samples_from_sheet Dere "${SAMPLE_SHEETS_DIR}/pausepred_sample_list.dere.tsv"
read_samples_from_sheet Dsec "${SAMPLE_SHEETS_DIR}/pausepred_sample_list.dsec.tsv"
read_samples_from_sheet Dsim "${SAMPLE_SHEETS_DIR}/pausepred_sample_list.dsim.tsv"
((${#ALL_SAMPLES[@]} > 0)) || die "No samples were found in the generated PausePred sample sheets"

active_pausepred() {
    pgrep -f '[r]un_pausepred_4species_safe\.sh|[o]ffline_pausepred_no_flag_filter(_mergefix)?\.pl' > /dev/null 2>&1
}

sample_success() {
    local sample="$1"
    if [[ "$ALLOW_OUTPUT_ONLY" == "1" ]]; then
        [[ -s "${PAUSEPRED_DIR}/${sample}.pausepred.csv" ]]
        return
    fi
    python3 "$COMPLETION_PY" --csv-dir "$PAUSEPRED_DIR" \
        --sample "$sample" --qc-dir "$PAUSEPRED_QC_DIR" > /dev/null
}

pausepred_status() {
    local pending=0
    local sample
    local missing=()
    for sample in "${ALL_SAMPLES[@]}"; do
        if ! sample_success "$sample"; then
            pending=$((pending + 1))
            missing+=("$sample")
        fi
    done
    echo "$pending:${missing[*]-}"
}

wait_for_pausepred() {
    local status pending names
    while true; do
        status="$(pausepred_status)"
        pending="${status%%:*}"
        names="${status#*:}"
        if [[ "$pending" == "0" ]]; then
            if active_pausepred; then
                echo "[WAIT] all expected CSV/log pairs exist, but PausePred processes are still active; waiting"
            else
                echo "[OK] all ${#ALL_SAMPLES[@]} PausePred samples passed the completion check"
                return 0
            fi
        elif active_pausepred; then
            echo "[WAIT] pending=$pending samples: $names; next check in ${WAIT_SECONDS}s"
        else
            echo "[ERROR] PausePred is no longer running, but these samples did not finish successfully: $names" >&2
            echo "[ERROR] Inspect ${PAUSEPRED_DIR}/logs/<sample>.pausepred.log and the species_runner_logs directory." >&2
            return 1
        fi
        [[ "$WAIT_FOR_PAUSEPRED" == "1" ]] || return 1
        sleep "$WAIT_SECONDS"
    done
}

wait_for_pausepred

write_header_list() {
    local path="$1"
    printf 'sample\n' > "$path"
}

DMEL_LIST="${WORK_DIR}/samples_Dmel_full.tsv"
WBG_LIST="${WORK_DIR}/samples_WBg.tsv"
DERE_LIST="${WORK_DIR}/samples_Dere.tsv"
DSEC_LIST="${WORK_DIR}/samples_Dsec.tsv"
DSIM_LIST="${WORK_DIR}/samples_Dsim.tsv"
for list in "$DMEL_LIST" "$WBG_LIST" "$DERE_LIST" "$DSEC_LIST" "$DSIM_LIST"; do
    write_header_list "$list"
done

for sample in "${ALL_SAMPLES[@]}"; do
    species="${SAMPLE_SPECIES[$sample]}"
    case "$species" in
        Dmel)
            if [[ "$sample" == WBg_* ]]; then
                printf '%s\n' "$sample" >> "$WBG_LIST"
            else
                printf '%s\n' "$sample" >> "$DMEL_LIST"
            fi
            ;;
        Dere) printf '%s\n' "$sample" >> "$DERE_LIST" ;;
        Dsec) printf '%s\n' "$sample" >> "$DSEC_LIST" ;;
        Dsim) printf '%s\n' "$sample" >> "$DSIM_LIST" ;;
        *) die "Unsupported species for sample $sample: $species" ;;
    esac
done

EDIT_ROOT="${EDIT_ROOT:-${RESCUE_ROOT}/03_cdspos/result/cdspos_for_alledit_rescue_3R_18907675_sourceCDS_keepTargetUTR3_fixed}"
ISO1_FULL="${EDIT_ROOT}/ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv"
ISO1_CDSPOS="${EDIT_ROOT}/ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"
ERE_FULL="${EDIT_ROOT}/ere.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv"
SIM_FULL="${EDIT_ROOT}/sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.full_with_longest_projection.tsv"
SEC_CDSPOS="${EDIT_ROOT}/sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"
DSEC_INTERVALS="${DSEC_INTERVALS:-/path/to/RNA_editing_DATAprocessing/05.psite/dsec_all_r1.3_with_UTR_CDSstop_UTR_position_CDS}"
for path in "$ISO1_FULL" "$ISO1_CDSPOS" "$ERE_FULL" "$SIM_FULL" "$SEC_CDSPOS" "$DSEC_INTERVALS"; do
    require_file "$path"
done

run_normal_group() {
    local label="$1"
    local sample_list="$2"
    local edit_sites="$3"
    local n
    n=$(($(wc -l < "$sample_list") - 1))
    ((n > 0)) || {
        echo "[SKIP] no samples in $label"
        return 0
    }
    echo "========== all-edit $label: $n sample(s) =========="
    env \
        SAMPLE_LIST="$sample_list" \
        EDIT_SITES_FILE="$edit_sites" \
        PAUSEPRED_DIR="$PAUSEPRED_DIR" \
        OUT_DIR="$ALLEDIT_OUT" \
        MIN_READS=2 \
        MIN_COVERAGE=20 \
        PAUSE_CUTOFFS="$PAUSE_CUTOFFS" \
        WINDOWS="$ENRICH_WINDOWS" \
        MAX_JOBS="$ENRICH_JOBS" \
        TARGET_REGION=all \
        WRITE_INTERMEDIATE=true \
        SCRIPT="$NORMAL_PY" \
        bash "$NORMAL_RUNNER"
}

# WBg uses its historical cdspos-for-alledit table; the other Dmel samples
# use the full longest-transcript projection table.
run_normal_group Dmel "$DMEL_LIST" "$ISO1_FULL"
run_normal_group WBg "$WBG_LIST" "$ISO1_CDSPOS"
run_normal_group Dere "$DERE_LIST" "$ERE_FULL"
run_normal_group Dsim "$DSIM_LIST" "$SIM_FULL"

if (($(wc -l < "$DSEC_LIST") > 1)); then
    echo "========== all-edit Dsec =========="
    env \
        SAMPLE_LIST="$DSEC_LIST" \
        ALLED_IT_CDS="$SEC_CDSPOS" \
        PAUSEPRED_DIR="$PAUSEPRED_DIR" \
        OUT_DIR="$ALLEDIT_OUT" \
        CDS_INTERVALS="$DSEC_INTERVALS" \
        MIN_READS=2 \
        MIN_COVERAGE=20 \
        PAUSE_CUTOFFS="$PAUSE_CUTOFFS" \
        WINDOWS="$ENRICH_WINDOWS" \
        MAX_JOBS="$ENRICH_JOBS" \
        PYTHON_SCRIPT="$DSEC_PY" \
        bash "$DSEC_RUNNER"
fi

merge_and_plot() {
    local enrichment_dir="$1"
    local plot_dir="${enrichment_dir}/plots"
    mkdir -p "$plot_dir"
    python3 "$MERGE_PY" \
        --indir "$enrichment_dir" \
        --outdir "$plot_dir" \
        --offset "$OFFSET" \
        --samples "$(
            IFS=,
            printf '%s' "${ALL_SAMPLES[*]}"
        )"
    python3 "$PLOT_PY" \
        --indir "$plot_dir" \
        --outdir "$plot_dir" \
        --windows "$WINDOWS" \
        --offset "$OFFSET" \
        --significance-column p_value
}

echo "========== merge/plot all-edit CDS =========="
merge_and_plot "$ALLEDIT_OUT"

# Build reference-A maps once per species. These are nucleotide-level CDS
# maps, not A-containing-codon maps: only unedited positions on reference A
# bases are retained in the second analysis.
DMEL_CDS="/path/to/RNA_editing_DATAprocessing_new/ref_genome/dmel-r6.12_GTF.CDS_position.tsv"
DMEL_FASTA="/path/to/RNA_editing_DATAprocessing_new/ref_genome/dmel_r6.12_longestCDS_transcriptome.fa"
DSIM_CDS="/path/to/RNA_editing_DATAprocessing_new/ref_genome/dsim-all-r2.02_GTF.CDS_position.tsv"
DSIM_FASTA="/path/to/RNA_editing_DATAprocessing_new/ref_genome/dsim_r2.02_longestCDS_transcriptome.fa"
DSEC_CDS="/path/to/RNA_editing_DATAprocessing/05.psite/dsec_all_r1.3_with_UTR_CDSstop_UTR_position_CDS"
DSEC_FASTA="/path/to/RNA_editing_DATAprocessing_new/ref_genome/dsec-all-transcript-r1.3.with_UTR.simple_longestCDS.fa"
DERE_CDS="/path/to/RNA_editing_DATAprocessing_new/ref_genome/dere-all-r1.05_GTF.CDS_position.tsv"
DERE_FASTA="/path/to/RNA_editing_DATAprocessing_new/ref_genome/dere_r1.05_longestCDS_transcriptome.fa"

MAP_DIR="${ABASE_OUT}/a_base_maps"
mkdir -p "$MAP_DIR"

build_a_map() {
    local species="$1"
    local cds="$2"
    local fasta="$3"
    local output="${MAP_DIR}/${species}.A_base_map.tsv.gz"
    require_file "$cds"
    require_file "$fasta"
    if [[ -s "$output" ]]; then
        echo "[REUSE] A-base map: $output"
    else
        echo "[BUILD] A-base map: $species"
        python3 "$A_BASE_MAP_PY" \
            --cds-intervals "$cds" \
            --reference-fasta "$fasta" \
            --out "$output"
    fi
}

build_a_map Dmel "$DMEL_CDS" "$DMEL_FASTA"
build_a_map Dere "$DERE_CDS" "$DERE_FASTA"
build_a_map Dsec "$DSEC_CDS" "$DSEC_FASTA"
build_a_map Dsim "$DSIM_CDS" "$DSIM_FASTA"

A_BASE_JOBS="${WORK_DIR}/unedit_A_base_jobs.tsv"
: > "$A_BASE_JOBS"
for sample in "${ALL_SAMPLES[@]}"; do
    species="${SAMPLE_SPECIES[$sample]}"
    printf '%s\t%s\n' "$sample" "${MAP_DIR}/${species}.A_base_map.tsv.gz" >> "$A_BASE_JOBS"
done

process_a_base() {
    local sample="$1"
    local a_map="$2"
    local pause_file="${PAUSEPRED_DIR}/${sample}.pausepred.csv"
    local edit_file="${ALLEDIT_OUT}/${sample}.pausepred_enrichment.editing_sites_used_in_universe.tsv"
    local output_file="${ABASE_OUT}/${sample}.pausepred_enrichment.tsv"
    require_file "$pause_file"
    require_file "$edit_file"
    require_file "$a_map"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] unedited-A $sample"
    python3 "$A_BASE_STATS_PY" \
        --pausepred "$pause_file" \
        --editing-sites-used-in-universe "$edit_file" \
        --a-base-map "$a_map" \
        --out "$output_file" \
        --min-reads 2 \
        --min-coverage 20 \
        --pause-cutoffs "$PAUSE_CUTOFFS" \
        --windows "$ENRICH_WINDOWS" \
        --background-mode unedited_A_base \
        --write-intermediate
}
export -f process_a_base require_file die
export PAUSEPRED_DIR ALLEDIT_OUT ABASE_OUT A_BASE_STATS_PY PAUSE_CUTOFFS ENRICH_WINDOWS

echo "========== unedited-A-base enrichment: ${#ALL_SAMPLES[@]} sample(s) =========="
xargs -a "$A_BASE_JOBS" -r -n 2 -P "$ENRICH_JOBS" \
    bash -c 'process_a_base "$1" "$2"' _

echo "========== merge/plot unedited-A-base =========="
merge_and_plot "$ABASE_OUT"

SAMPLE_CSV="$(
    IFS=,
    printf '%s' "${ALL_SAMPLES[*]}"
)"
HAIRPIN_SUMMARY="${HAIRPIN_OUT}/hairpin_pausepred_edit_unedit_enrichment_summary.tsv"
echo "========== hairpin enrichment: ${#ALL_SAMPLES[@]} sample(s) =========="
python3 "$HAIRPIN_PY" \
    --all-edit-dir "$ALLEDIT_OUT" \
    --a-base-dir "$ABASE_OUT" \
    --hairpin-dir "$HAIRPIN_DIR" \
    --out-dir "$HAIRPIN_OUT" \
    --samples "$SAMPLE_CSV" \
    --pause-cutoffs "$PAUSE_CUTOFFS" \
    --windows "$ENRICH_WINDOWS" \
    --jobs "$HAIRPIN_JOBS" \
    --executor process

[[ -s "$HAIRPIN_SUMMARY" ]] || die "Hairpin summary was not created: $HAIRPIN_SUMMARY"
HAIRPIN_PLOT_DIR="${HAIRPIN_OUT}/plots"
mkdir -p "$HAIRPIN_PLOT_DIR"
python3 "$HAIRPIN_PLOT_PY" \
    --summary "$HAIRPIN_SUMMARY" \
    --out-dir "$HAIRPIN_PLOT_DIR" \
    --windows "$WINDOWS" \
    --cutoff "$OFFSET"

CONFIG_TSV="${WORK_DIR}/workflow_config.tsv"
{
    printf 'parameter\tvalue\n'
    printf 'run_tag\t%s\n' "$RUN_TAG"
    printf 'pausepred_dir\t%s\n' "$PAUSEPRED_DIR"
    printf 'pause_score_cutoffs_calculated\t%s\n' "$PAUSE_CUTOFFS"
    printf 'windows_calculated\t%s\n' "$ENRICH_WINDOWS"
    printf 'plot_windows\t%s\n' "$WINDOWS"
    printf 'plot_offset\t%s\n' "$OFFSET"
    printf 'min_reads\t2\n'
    printf 'min_coverage\t20\n'
    printf 'all_edit_output\t%s\n' "$ALLEDIT_OUT"
    printf 'unedited_A_base_output\t%s\n' "$ABASE_OUT"
    printf 'hairpin_output\t%s\n' "$HAIRPIN_OUT"
    printf 'excluded_samples\t%s\n' "$PAUSEPRED_EXCLUDE_SAMPLES"
    printf 'sample_count\t%s\n' "${#ALL_SAMPLES[@]}"
    printf 'samples\t%s\n' "$SAMPLE_CSV"
} > "$CONFIG_TSV"

echo
echo "[DONE] PausePred post-processing completed."
echo "[DONE] all-edit table/plot:       $ALLEDIT_OUT/plots"
echo "[DONE] unedited-A table/plot:     $ABASE_OUT/plots"
echo "[DONE] hairpin table/plots:       $HAIRPIN_OUT"
echo "[DONE] run configuration:         $CONFIG_TSV"
echo "[DONE] selected PDFs:"
echo "       $ALLEDIT_OUT/plots/odds_ratio_fraction_combined.window_${WINDOWS}.cutoff_${OFFSET}.pdf"
echo "       $ABASE_OUT/plots/odds_ratio_fraction_combined.window_${WINDOWS}.cutoff_${OFFSET}.pdf"
echo "       $HAIRPIN_OUT/plots/odds_ratio_fraction_combined.all_edit_cds.window_${WINDOWS}.cutoff_${OFFSET}.pdf"
echo "       $HAIRPIN_OUT/plots/odds_ratio_fraction_combined.unedit_A_base.window_${WINDOWS}.cutoff_${OFFSET}.pdf"
