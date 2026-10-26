#!/usr/bin/env bash

# Run PausePred for Dmel, Dere, Dsec and Dsim concurrently.
#
# The old batch runner starts several large BAM scans inside each species.
# This wrapper supports bounded sample-level parallelism within each species.
# Set SPECIES_JOBS and MAX_JOBS according to available memory/CPU.

set -uo pipefail

MODE="run"
case "${1:-}" in
    "" | --run) MODE="run" ;;
    --review) MODE="review" ;;
    *)
        echo "Usage: bash $0 [--review|--run]" >&2
        exit 2
        ;;
esac

THIS_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$THIS_DIR/../../.." && pwd)}"
PAUSE_RUNTIME="${PAUSE_RUNTIME:-$PIPELINE_DIR/.runtime/06}"
MERGED_BAM_DIR="${MERGED_BAM_DIR:-$PIPELINE_DIR/results/05_stalling/psite_2times_merged_bam}"
MASTER_SHEET="${MASTER_SHEET:-${MERGED_BAM_DIR}/pausepred_sample_sheet.offset_plus3.tsv}"
SHEET_DIR="${SHEET_DIR:-${MERGED_BAM_DIR}/pausepred_sample_sheets}"
SCRIPT="${SCRIPT:-${THIS_DIR}/offline_pausepred_no_flag_filter_mergefix.pl}"
REF_DIR="${REF_DIR:-/path/to/RNA_editing_DATAprocessing_new/ref_genome}"
OUTDIR="${OUTDIR:-${PIPELINE_DIR}/results/05_stalling/pausepred_results_new}"

SPECIES_SELECTION="${SPECIES_SELECTION:-dmel dere dsec dsim}"
SPECIES_JOBS="${SPECIES_JOBS:-4}"
MAX_JOBS="${MAX_JOBS:-1}"

WINDOW_SIZE="${WINDOW_SIZE:-500}"
# A zero cutoff retains almost every covered position and can make the Perl
# script consume hundreds of GB while sorting its in-memory candidate list.
# Override explicitly with PAUSE_CUTOFF=0 only if that is intended.
PAUSE_CUTOFF="${PAUSE_CUTOFF:-0}"
WINDOW_COVERAGE="${WINDOW_COVERAGE:-10}"
UPSTREAM_SEQ="${UPSTREAM_SEQ:-180}"
DOWNSTREAM_SEQ="${DOWNSTREAM_SEQ:-180}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"

# PausePred calls samtools directly and requires Bio::DB::Fasta.
export PATH="${PAUSE_RUNTIME}/bin:${PATH}"
export PERL5LIB="${THIS_DIR}/lib:${PAUSE_RUNTIME}/lib${PERL5LIB:+:${PERL5LIB}}"

declare -A SAMPLE_REGEX=(
    [dmel]='^(AdarG_|AdarS_|WBg_|Dmel_)'
    [dere]='^Dere_'
    [dsec]='^Dsec_'
    [dsim]='^Dsim_'
)
declare -A FASTA=(
    [dmel]="${REF_DIR}/dmel_r6.12_longestCDS_transcriptome.fa"
    [dere]="${REF_DIR}/dere_r1.05_longestCDS_transcriptome.fa"
    [dsec]="${REF_DIR}/dsec-all-transcript-r1.3.with_UTR.simple_longestCDS.fa"
    [dsim]="${REF_DIR}/dsim_r2.02_longestCDS_transcriptome.fa"
)

read -r -a SPECIES_GROUPS <<< "${SPECIES_SELECTION}"

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || die "Missing or empty file: ${path}"
}

validate_integer() {
    local name="$1"
    local value="$2"
    [[ "${value}" =~ ^[0-9]+$ && "${value}" -ge 1 ]] || die "${name} must be a positive integer: ${value}"
}

prepare_sheets() {
    require_file "${MASTER_SHEET}"
    require_file "${SCRIPT}"
    command -v perl > /dev/null || die "perl not found"
    command -v samtools > /dev/null || die "samtools not found"
    perl -MBio::DB::Fasta -e 1 2> /dev/null || die "Perl module Bio::DB::Fasta is unavailable"

    awk -F '\t' '
        NR == 1 {
            if ($1 != "sample" || $2 != "bam" || $3 != "read_lengths" || $4 != "offsets") {
                print "[ERROR] Unexpected sample-sheet header" > "/dev/stderr"; exit 1
            }
            next
        }
        NF != 4 { print "[ERROR] Expected 4 tab-separated columns at line " NR > "/dev/stderr"; bad=1 }
        $4 ~ /(^|,)NA(,|$)/ { print "[ERROR] NA offset for sample " $1 > "/dev/stderr"; bad=1 }
        n_lengths=split($3, lengths, ","); n_offsets=split($4, offsets, ",")
        n_lengths != n_offsets {
            print "[ERROR] read_lengths/offsets count mismatch for " $1 > "/dev/stderr"; bad=1
        }
        END { exit bad }
    ' "${MASTER_SHEET}" || die "Master sample sheet validation failed"

    mkdir -p "${SHEET_DIR}" "${OUTDIR}/logs" "${OUTDIR}/species_runner_logs"

    for group in "${SPECIES_GROUPS[@]}"; do
        [[ -n "${SAMPLE_REGEX[${group}]+yes}" ]] || die "Unknown species group: ${group}"
        require_file "${FASTA[${group}]}"

        local sheet="${SHEET_DIR}/pausepred_sample_list.${group}.tsv"
        local tmp="${sheet}.tmp.$$"
        awk -F '\t' -v OFS='\t' -v regex="${SAMPLE_REGEX[${group}]}" \
            'NR == 1 || $1 ~ regex' "${MASTER_SHEET}" > "${tmp}" || die "Could not create ${sheet}"
        mv -f "${tmp}" "${sheet}"

        local n_samples
        n_samples="$(($(wc -l < "${sheet}") - 1))"
        ((n_samples > 0)) || die "No samples selected for ${group}"
        echo "[PLAN] group=${group} samples=${n_samples} fasta=${FASTA[${group}]}"

        while IFS=$'\t' read -r sample bam read_lengths offsets; do
            [[ "${sample}" == "sample" ]] && continue
            require_file "${bam}"
            if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
                die "BAM index missing: ${bam}"
            fi
            samtools quickcheck -v "${bam}" < /dev/null || die "BAM quickcheck failed: ${bam}"
        done < "${sheet}"
    done
}

run_one_sample() {
    local group="$1"
    local fasta="$2"
    local sample="$3"
    local bam="$4"
    local read_lengths="$5"
    local offsets="$6"
    local out_csv="${OUTDIR}/${sample}.pausepred.csv"
    local log="${OUTDIR}/logs/${sample}.pausepred.log"
    local tmp_out=""
    local rc

    # The log is intentionally replaced for a rerun, while an existing
    # non-empty result is preserved when SKIP_EXISTING=1.
    {
        echo "[INFO] Start PausePred: $(date)"
        echo "[INFO] group=${group}"
        echo "[INFO] sample=${sample}"
        echo "[INFO] bam=${bam}"
        echo "[INFO] output=${out_csv}"
        echo "[INFO] window=${WINDOW_SIZE} cutoff=${PAUSE_CUTOFF} coverage=${WINDOW_COVERAGE} upstream=${UPSTREAM_SEQ} downstream=${DOWNSTREAM_SEQ}"

        if [[ "${SKIP_EXISTING}" == "1" && -s "${out_csv}" ]]; then
            echo "[INFO] Skip existing output: ${out_csv}"
            return 0
        fi

        local n_len n_off
        n_len="$(awk -F',' '{print NF}' <<< "${read_lengths}")"
        n_off="$(awk -F',' '{print NF}' <<< "${offsets}")"
        if [[ "${n_len}" -ne "${n_off}" ]]; then
            echo "[ERROR] read_lengths count != offsets count"
            echo "[ERROR] read_lengths=${read_lengths}"
            echo "[ERROR] offsets=${offsets}"
            return 1
        fi

        if ! samtools quickcheck -v "${bam}"; then
            echo "[ERROR] samtools quickcheck failed: ${bam}"
            return 1
        fi

        if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
            echo "[INFO] BAM index not found; creating it"
            if ! samtools index -@ 2 "${bam}"; then
                echo "[ERROR] samtools index failed: ${bam}"
                return 1
            fi
        fi

        # PausePred obtains the output path from STDIN.  Use a temporary path
        # so a killed Perl process cannot leave a false successful CSV.
        tmp_out="$(mktemp "${out_csv}.tmp.XXXXXX")" || {
            echo "[ERROR] Could not create temporary output path"
            return 1
        }
        echo "[INFO] temporary_output=${tmp_out}"

        if printf '%s\n' "${tmp_out}" | perl "${SCRIPT}" \
            "${bam}" \
            "${WINDOW_SIZE}" \
            "${PAUSE_CUTOFF}" \
            "${fasta}" \
            "${read_lengths}" \
            "${WINDOW_COVERAGE}" \
            "${UPSTREAM_SEQ}" \
            "${DOWNSTREAM_SEQ}" \
            "${offsets}"; then
            rc=0
        else
            rc=$?
        fi

        if [[ "${rc}" -ne 0 ]]; then
            echo "[ERROR] PausePred exited with status ${rc}"
            rm -f "${tmp_out}"
            tmp_out=""
            return "${rc}"
        fi
        if [[ ! -s "${tmp_out}" ]]; then
            echo "[ERROR] PausePred returned success but temporary CSV is empty"
            rm -f "${tmp_out}"
            tmp_out=""
            return 1
        fi
        mv -f "${tmp_out}" "${out_csv}" || {
            echo "[ERROR] Could not install output: ${out_csv}"
            rm -f "${tmp_out}"
            tmp_out=""
            return 1
        }
        tmp_out=""
        echo "[INFO] End PausePred: $(date)"
        echo "[INFO] Done: ${sample}"
    } > "${log}" 2>&1

    rc=$?
    if [[ -n "${tmp_out}" ]]; then
        rm -f "${tmp_out}"
    fi
    return "${rc}"
}

run_species() {
    local group="$1"
    local sheet="${SHEET_DIR}/pausepred_sample_list.${group}.tsv"
    local species_log="${OUTDIR}/species_runner_logs/${group}.runner.log"
    local failed=0
    local total=0
    local done_count=0
    local sample bam read_lengths offsets rc

    {
        echo "[INFO] group=${group}"
        echo "[INFO] sheet=${sheet}"
        echo "[INFO] fasta=${FASTA[${group}]}"
        echo "[INFO] max sample jobs=${MAX_JOBS}"
        echo "[INFO] started=$(date)"

        declare -a sample_pids=()
        declare -a sample_names=()

        reap_one() {
            local pid="$1"
            local finished_sample="$2"
            local rc
            if wait "${pid}"; then
                done_count=$((done_count + 1))
                echo "[DONE] ${finished_sample}"
            else
                rc=$?
                failed=$((failed + 1))
                echo "[FAILED] ${finished_sample} exit=${rc} log=${OUTDIR}/logs/${finished_sample}.pausepred.log"
                tail -n 20 "${OUTDIR}/logs/${finished_sample}.pausepred.log" || true
            fi
        }

        while IFS=$'\t' read -r sample bam read_lengths offsets; do
            [[ -z "${sample}" || "${sample}" == "sample" || "${sample}" =~ ^# ]] && continue
            total=$((total + 1))
            run_one_sample "${group}" "${FASTA[${group}]}" "${sample}" "${bam}" "${read_lengths}" "${offsets}" &
            sample_pids+=("$!")
            sample_names+=("${sample}")

            if ((${#sample_pids[@]} >= MAX_JOBS)); then
                reap_one "${sample_pids[0]}" "${sample_names[0]}"
                sample_pids=("${sample_pids[@]:1}")
                sample_names=("${sample_names[@]:1}")
            fi
        done < "${sheet}"

        for i in "${!sample_pids[@]}"; do
            reap_one "${sample_pids[$i]}" "${sample_names[$i]}"
        done

        echo "[SUMMARY] group=${group} total=${total} done=${done_count} failed=${failed} finished=$(date)"
        if [[ "${failed}" -ne 0 ]]; then
            return 1
        fi
    } > "${species_log}" 2>&1
}

validate_integer "SPECIES_JOBS" "${SPECIES_JOBS}"
validate_integer "MAX_JOBS" "${MAX_JOBS}"

if [[ "${PAUSE_CUTOFF}" == "0" ]]; then
    echo "[WARN] PAUSE_CUTOFF=0 can create a very large in-memory candidate list; PAUSE_CUTOFF=2 or 5 is safer" >&2
fi

prepare_sheets
samtools --version | head -1
echo "[INFO] SPECIES_SELECTION=${SPECIES_SELECTION}"
echo "[INFO] SPECIES_JOBS=${SPECIES_JOBS}"
echo "[INFO] MAX_JOBS=${MAX_JOBS}"
echo "[INFO] OUTDIR=${OUTDIR}"
echo "[INFO] PAUSE_CUTOFF=${PAUSE_CUTOFF}"

if [[ "${MODE}" == "review" ]]; then
    echo "[REVIEW ONLY] No PausePred jobs were started."
    exit 0
fi

failed_groups=0
declare -a pids=()
declare -a pid_groups=()

for group in "${SPECIES_GROUPS[@]}"; do
    echo "[SUBMIT] group=${group} log=${OUTDIR}/species_runner_logs/${group}.runner.log"
    run_species "${group}" &
    pids+=("$!")
    pid_groups+=("${group}")

    if ((${#pids[@]} >= SPECIES_JOBS)); then
        for i in "${!pids[@]}"; do
            if wait "${pids[$i]}"; then
                echo "[SPECIES DONE] ${pid_groups[$i]}"
            else
                echo "[SPECIES FAILED] ${pid_groups[$i]} inspect ${OUTDIR}/species_runner_logs/${pid_groups[$i]}.runner.log" >&2
                failed_groups=$((failed_groups + 1))
            fi
        done
        pids=()
        pid_groups=()
    fi
done

for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
        echo "[SPECIES DONE] ${pid_groups[$i]}"
    else
        echo "[SPECIES FAILED] ${pid_groups[$i]} inspect ${OUTDIR}/species_runner_logs/${pid_groups[$i]}.runner.log" >&2
        failed_groups=$((failed_groups + 1))
    fi
done

if [[ "${failed_groups}" -ne 0 ]]; then
    echo "[ERROR] ${failed_groups} species group(s) failed" >&2
    exit 1
fi

echo "[DONE] All species groups completed: $(date)"
