#!/usr/bin/env bash
# Adapted from run_Dmel_T25Z3.sh; prediction and interval parsing/merging are unchanged.
set -euo pipefail
export LC_ALL=C

species="${RNALFOLD_SPECIES:-Dmel}"
case "$species" in Dmel | Dsim | Dsec | Dere) ;; *)
    echo "Unknown species: $species" >&2
    exit 2
    ;;
esac
input_fasta="${RNALFOLD_INPUT_FASTA:?Set RNALFOLD_INPUT_FASTA to an unedited longest-CDS FASTA}"
rnalfold_command="${RNALFOLD_BIN:-RNALfold}"
rnalfold_bin="$(command -v "$rnalfold_command")" || {
    echo "RNALfold executable not found: $rnalfold_command; set RNALFOLD_BIN." >&2
    exit 2
}
[[ -x "$rnalfold_bin" ]] || {
    echo "Not executable: $rnalfold_bin" >&2
    exit 2
}
[[ -s "$input_fasta" ]] || {
    echo "Missing FASTA: $input_fasta" >&2
    exit 2
}
input_fasta="$(realpath "$input_fasta")"
rnalfold_bin="$(realpath "$rnalfold_bin")"
records_per_chunk="${RECORDS_PER_CHUNK:-100}"
parallel_jobs="${RNALFOLD_JOBS:-2}"
[[ "$records_per_chunk" =~ ^[1-9][0-9]*$ && "$parallel_jobs" =~ ^[1-9][0-9]*$ ]] || {
    echo 'RECORDS_PER_CHUNK and RNALFOLD_JOBS must be positive integers.' >&2
    exit 2
}
output_dir="${RNALFOLD_OUT_DIR:?Set RNALFOLD_OUT_DIR}"
mkdir -p "$output_dir"
output_dir="$(cd -- "$output_dir" && pwd)"
run_root="${output_dir}/${species}_T25Z3_run"
chunk_root="${run_root}/chunks"
raw_root="${run_root}/raw"
raw_output="${output_dir}/${species}.fa.T25Z3"
one_tsv="${output_dir}/${species}.fa.T25Z3.1tsv"
two_tsv="${output_dir}/${species}.fa.T25Z3.2tsv"
sorted_tsv="${run_root}/${species}.fa.T25Z3.2.sorted.tsv"
merged_tsv="${output_dir}/${species}.fa.T25Z3.2.merged.tsv"

# Resume completed chunks only for the same FASTA, executable, script and chunk size.
# Use a new output directory if any of these inputs change.
signature="$(sha256sum "$input_fasta" "$rnalfold_bin" "${BASH_SOURCE[0]}")
records_per_chunk=$records_per_chunk
parameters=--noLP --temp=25 --zscore=-3"
if [[ -f "$run_root/input_signature.txt" ]]; then
    [[ "$(cat "$run_root/input_signature.txt")" == "$signature" ]] || {
        echo "RNALfold inputs changed; use a new RNALFOLD_OUT_DIR: $run_root" >&2
        exit 2
    }
elif [[ -d "$run_root" && -n "$(find "$run_root" -mindepth 1 -print -quit)" ]]; then
    echo "Cannot validate existing RNALfold cache; use a new output directory: $run_root" >&2
    exit 2
fi
mkdir -p "$chunk_root" "$raw_root"
printf '%s\n' "$signature" > "$run_root/input_signature.txt"

if [[ -z "$(find "$chunk_root" -maxdepth 1 -type f -name '*.fa' -print -quit)" ]]; then
    awk -v chunk_prefix="${chunk_root}/chunk" -v per_chunk="$records_per_chunk" '
        /^>/ {
            if (record_count % per_chunk == 0) {
                if (out_file != "")
                    close(out_file)
                chunk_number++
                out_file = sprintf("%s.%05d.fa", chunk_prefix, chunk_number)
            }
            record_count++
        }
        {
            if (out_file == "")
                next
            print > out_file
        }
        END {
            if (out_file != "")
                close(out_file)
        }
    ' "$input_fasta"
fi

chunk_count="$(find "$chunk_root" -maxdepth 1 -type f -name '*.fa' | wc -l)"
if [[ "$chunk_count" -eq 0 ]]; then
    printf 'No FASTA chunks were generated.\n' >&2
    exit 1
fi

find "$chunk_root" -maxdepth 1 -type f -name '*.fa' -print0 |
    sort -z |
    xargs -0 -r -n 1 -P "$parallel_jobs" bash -c '
        raw_dir="$1"
        rnalfold="$2"
        chunk_path="$3"
        chunk_name="${chunk_path##*/}"
        chunk_stem="${chunk_name%.fa}"
        raw_path="${raw_dir}/${chunk_stem}.fold"
        log_path="${raw_dir}/${chunk_stem}.log"
        done_path="${raw_dir}/${chunk_stem}.done"

        if [[ -f "$done_path" && -s "$raw_path" ]]; then
            exit 0
        fi

        temporary_path="${raw_path}.tmp.$$"
        if "$rnalfold" --noLP --temp=25 --zscore=-3 -i "$chunk_path" > "$temporary_path" 2> "$log_path"; then
            mv "$temporary_path" "$raw_path"
            : > "$done_path"
        else
            status=$?
            rm -f "$temporary_path"
            exit "$status"
        fi
    ' _ "$raw_root" "$rnalfold_bin"

raw_count="$(find "$raw_root" -maxdepth 1 -type f -name '*.fold' | wc -l)"
if [[ "$raw_count" -ne "$chunk_count" ]]; then
    printf 'Expected %s raw chunk outputs but found %s.\n' "$chunk_count" "$raw_count" >&2
    exit 1
fi

find "$raw_root" -maxdepth 1 -type f -name '*.fold' -print0 |
    sort -z |
    xargs -0 -r cat > "$raw_output"

awk '
BEGIN {
    OFS = "\t"
    print "transcript", "start", "end", "length", "MFE", "zscore", "structure"
}

/^>/ {
    header = substr($0, 2)
    split(header, a, /[[:space:]]+/)
    transcript = a[1]
    next
}

/^[.()]+[[:space:]]+\(-?[0-9.]+\)[[:space:]]+[0-9]+[[:space:]]+z=/ {
    structure = $1

    mfe = $2
    gsub(/[()]/, "", mfe)

    start = $3
    zscore = $5

    len = length(structure)
    end = start + len - 1

    print transcript, start, end, len, mfe, zscore, structure
}
' "$raw_output" > "$one_tsv"

awk '
BEGIN {
    OFS = "\t"
}

/^transcript/ {
    print
    next
}

{
    structure = $7
    depth = 0
    hairpin_loops = 0
    valid = 1

    delete contains_pair

    for (i = 1; i <= length(structure); i++) {
        c = substr(structure, i, 1)

        if (c == "(") {
            if (depth > 0)
                contains_pair[depth] = 1

            depth++
            contains_pair[depth] = 0
        }
        else if (c == ")") {
            if (depth <= 0) {
                valid = 0
                break
            }

            if (contains_pair[depth] == 0)
                hairpin_loops++

            delete contains_pair[depth]
            depth--
        }
    }

    if (depth != 0)
        valid = 0

    if (valid && hairpin_loops == 1)
        print
}
' "$one_tsv" > "$two_tsv"

awk 'NR > 1' "$two_tsv" |
    LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n > "$sorted_tsv"

awk -F '\t' -v OFS='\t' '
function begin_group(transcript_value, start_value, end_value, mfe_value, zscore_value, structure_value) {
    current_transcript = transcript_value
    current_start = start_value + 0
    current_end = end_value + 0
    hit_count = 1
    best_mfe = mfe_value + 0
    best_zscore = zscore_value + 0
    best_structure = structure_value
    have_group = 1
}

function emit_group() {
    print current_transcript, current_start, current_end,
          current_end - current_start + 1, hit_count,
          best_mfe, best_zscore, best_structure
}

BEGIN {
    print "transcript", "start", "end", "length", "n_hits", "best_MFE", "best_zscore", "best_structure"
}

{
    transcript = $1
    start = $2 + 0
    end = $3 + 0
    mfe = $5 + 0
    zscore = $6 + 0
    structure = $7

    if (!have_group) {
        begin_group(transcript, start, end, mfe, zscore, structure)
        next
    }

    if (transcript == current_transcript && start <= current_end) {
        if (end > current_end)
            current_end = end
        hit_count++
        if (mfe < best_mfe) {
            best_mfe = mfe
            best_structure = structure
        }
        if (zscore < best_zscore)
            best_zscore = zscore
    }
    else {
        emit_group()
        begin_group(transcript, start, end, mfe, zscore, structure)
    }
}

END {
    if (have_group)
        emit_group()
}
' "$sorted_tsv" > "$merged_tsv"

input_transcripts="$(awk '/^>/{n++} END{print n+0}' "$input_fasta")"
output_transcripts="$(awk '/^>/{n++} END{print n+0}' "$raw_output")"
one_rows="$(awk 'NR>1{n++} END{print n+0}' "$one_tsv")"
two_rows="$(awk 'NR>1{n++} END{print n+0}' "$two_tsv")"
merged_rows="$(awk 'NR>1{n++} END{print n+0}' "$merged_tsv")"

printf 'input_transcripts=%s\n' "$input_transcripts"
printf 'raw_output_transcripts=%s\n' "$output_transcripts"
printf 'predictions_1tsv=%s\n' "$one_rows"
printf 'hairpin_rows_2tsv=%s\n' "$two_rows"
printf 'merged_regions=%s\n' "$merged_rows"

if [[ "$input_transcripts" -ne "$output_transcripts" ]]; then
    printf 'Transcript count mismatch between input and RNALfold output.\n' >&2
    exit 1
fi
