#!/bin/bash
set -euo pipefail

# ============================================================
# DIA-NN 2.6.1 search for selected Orbitrap Astral RAW files
# Generic version for multiple species.
#
# Usage:
#   bash run_diann_astral_multi_species.sh \
#       SPECIES RUN_TAG LIB FASTA PREFIX1 [PREFIX2 ...]
#
# Arguments:
#   SPECIES : dmel, dsim, dsec, dere, or another species tag
#   RUN_TAG : name of this DIA-NN run
#   LIB     : predicted/empirical spectral library path
#   FASTA   : matching protein FASTA path; use NONE to omit FASTA
#   PREFIX  : one or more RAW file prefixes
#
# Example:
#   bash run_diann_astral_multi_species.sh \
#       dmel WBg_ISO1 \
#       /data/dongxn/MS/DIA_NN_RUN/Dmel/dmel.predicted.speclib \
#       /data/dongxn/MS/dmel-all-CDS-r6.12.longest_CDS_per_gene.protein.fasta \
#       XQH_20260320_WBg_ XQH_20260320_ISO1_
#
# Optional environment variables:
#   RAW_DIR       RAW file directory
#   OUT_ROOT_BASE Base output directory
#   THREADS       DIA-NN threads
#   MS1_ACC       MS1 mass accuracy
#   MS2_ACC       MS2 mass accuracy
#   SCAN_WINDOW   scan window radius
# ============================================================

if [ "$#" -lt 5 ]; then
    echo "Usage:" >&2
    echo "  bash $0 SPECIES RUN_TAG LIB FASTA PREFIX1 [PREFIX2 ...]" >&2
    echo >&2
    echo "Example:" >&2
    echo "  bash $0 dmel WBg_ISO1 /path/dmel.predicted.speclib /path/dmel.protein.fasta XQH_20260320_WBg_ XQH_20260320_ISO1_" >&2
    echo >&2
    echo "Use FASTA=NONE only when intentionally searching without a FASTA database." >&2
    exit 1
fi

SPECIES_RAW="$1"
RUN_TAG="$2"
LIB="$3"
FASTA="$4"
shift 4
PREFIXES=("$@")

SPECIES=$(printf '%s' "$SPECIES_RAW" | tr '[:upper:]' '[:lower:]')
case "$SPECIES" in
    dmel) SPECIES_DIR="Dmel" ;;
    dsim) SPECIES_DIR="Dsim" ;;
    dsec) SPECIES_DIR="Dsec" ;;
    dere) SPECIES_DIR="Dere" ;;
    *)
        # Preserve a readable directory name for additional species.
        SPECIES_DIR="${SPECIES_RAW}"
        ;;
esac

# ---------- DIA-NN ----------
DIANN_DIR="${DIANN_DIR:-$HOME/software/diann-2.6.1/}"
DIANN="${DIANN_DIR%/}/diann-linux"
DIANN_DIR="${DIANN_DIR%/}/"

# ---------- .NET ----------
DOTNET_DIR="${DOTNET_DIR:-$HOME/software/dotnet8}"
export DOTNET_ROOT="$DOTNET_DIR"
export DOTNET_ROOT_X64="$DOTNET_DIR"
export DOTNET_HOST_PATH="$DOTNET_DIR/dotnet"
export DOTNET_MULTILEVEL_LOOKUP=0
export PATH="$DOTNET_DIR:$PATH"
export LD_LIBRARY_PATH="${DIANN_DIR}:${LD_LIBRARY_PATH:-}"

# ---------- Input/output ----------
RAW_DIR="${RAW_DIR:-/data/dongxn/MS/Rawdata_Astral}"
OUT_ROOT_BASE="${OUT_ROOT_BASE:-/data/dongxn/RNAediting/codex/03_MS/DIA_NN_run_trypsin}"
OUT_ROOT="${OUT_ROOT_BASE}/${SPECIES_DIR}"
OUT_DIR="${OUT_ROOT}/${RUN_TAG}"
TEMP_DIR="${OUT_DIR}/temp"
RAW_LIST="${OUT_DIR}/selected_raw_files.txt"
COMMAND_FILE="${OUT_DIR}/diann_command.txt"

mkdir -p "$OUT_DIR"
rm -rf "$TEMP_DIR"
mkdir -p "$TEMP_DIR"

# ---------- Orbitrap Astral parameters ----------
MS1_ACC="${MS1_ACC:-4}"
MS2_ACC="${MS2_ACC:-10}"
SCAN_WINDOW="${SCAN_WINDOW:-5}"
THREADS="${THREADS:-20}"

# ---------- Basic checks ----------
if [ ! -x "$DIANN" ]; then
    echo "ERROR: DIA-NN is missing or not executable:" >&2
    echo "  $DIANN" >&2
    exit 1
fi

if [ ! -r "$LIB" ]; then
    echo "ERROR: spectral library is missing or unreadable:" >&2
    echo "  $LIB" >&2
    exit 1
fi

USE_FASTA=1
case "$FASTA" in
    NONE | none | -)
        USE_FASTA=0
        ;;
    *)
        if [ ! -r "$FASTA" ]; then
            echo "ERROR: FASTA is missing or unreadable:" >&2
            echo "  $FASTA" >&2
            exit 1
        fi
        ;;
esac

if [ ! -d "$RAW_DIR" ]; then
    echo "ERROR: RAW directory does not exist:" >&2
    echo "  $RAW_DIR" >&2
    exit 1
fi

# ---------- Find RAW files by prefix ----------
: > "$RAW_LIST"

echo "Searching RAW files in:"
echo "  $RAW_DIR"
echo

for prefix in "${PREFIXES[@]}"; do
    prefix_count=0

    while IFS= read -r -d '' raw_file; do
        printf '%s\n' "$raw_file" >> "$RAW_LIST"
        prefix_count=$((prefix_count + 1))
    done < <(
        find "$RAW_DIR" \
            -maxdepth 1 \
            -type f \
            -name "${prefix}*.raw" \
            -print0
    )

    echo "Prefix: $prefix"
    echo "Matched RAW files: $prefix_count"
    echo
done

sort -u "$RAW_LIST" -o "$RAW_LIST"
N_RAW=$(wc -l < "$RAW_LIST" | tr -d '[:space:]')

if [ "$N_RAW" -eq 0 ]; then
    echo "ERROR: no RAW files matched the supplied prefixes." >&2
    echo "Prefixes:" >&2
    printf '  %s\n' "${PREFIXES[@]}" >&2
    exit 1
fi

echo "============================================"
echo "Species: $SPECIES"
echo "Run tag: $RUN_TAG"
echo "Spectral library: $LIB"
if [ "$USE_FASTA" -eq 1 ]; then
    echo "FASTA: $FASTA"
else
    echo "FASTA: omitted"
fi
echo "Unique RAW files selected: $N_RAW"
echo "RAW list: $RAW_LIST"
echo "Output directory: $OUT_DIR"
echo "MS1 accuracy: ${MS1_ACC} ppm"
echo "MS2 accuracy: ${MS2_ACC} ppm"
echo "Scan window radius: $SCAN_WINDOW"
echo "Threads: $THREADS"
echo "============================================"
echo
cat "$RAW_LIST"
echo

# ---------- Build DIA-NN arguments ----------
RAW_ARGS=()
while IFS= read -r raw_file; do
    if [ ! -r "$raw_file" ]; then
        echo "ERROR: RAW file is unreadable:" >&2
        echo "  $raw_file" >&2
        exit 1
    fi
    RAW_ARGS+=(--f "$raw_file")
done < "$RAW_LIST"

FASTA_ARGS=()
if [ "$USE_FASTA" -eq 1 ]; then
    FASTA_ARGS+=(--fasta "$FASTA")
fi

DIANN_ARGS=(
    --parent "$DIANN_DIR"
    "${RAW_ARGS[@]}"
    --lib "$LIB"
    "${FASTA_ARGS[@]}"
    --temp "$TEMP_DIR"
    --mass-acc-ms1 "$MS1_ACC"
    --mass-acc "$MS2_ACC"
    --window "$SCAN_WINDOW"
    --threads "$THREADS"
    --verbose 2
    --out "$OUT_DIR/report.parquet"
    --qvalue 0.05
    --matrix-qvalue 0.01
    --matrices
    --proteoforms
    --reanalyse
    --pg-level 0
    --ids-to-names
    --unimod4
    --var-mods 1
    --var-mod "UniMod:35,15.994915,M"
    --no-skyline
)

# Save a reproducible shell-escaped command before execution.
{
    printf '%q ' "$DIANN"
    printf '%q ' "${DIANN_ARGS[@]}"
    printf '\n'
} > "$COMMAND_FILE"

echo "Installed .NET runtimes:"
dotnet --list-runtimes
echo

echo "DIA-NN command saved to:"
echo "  $COMMAND_FILE"
echo

cd "$DIANN_DIR"
"$DIANN" "${DIANN_ARGS[@]}"

echo
echo "============================================"
echo "DIA-NN finished"
echo "Species: $SPECIES"
echo "Selected RAW files: $RAW_LIST"
echo "Main report: $OUT_DIR/report.parquet"
echo "Precursor matrix: $OUT_DIR/report.pr_matrix.tsv"
echo "Protein-group matrix: $OUT_DIR/report.pg_matrix.tsv"
echo "Stats report: $OUT_DIR/report.stats.tsv"
echo "Command: $COMMAND_FILE"
echo "============================================"
