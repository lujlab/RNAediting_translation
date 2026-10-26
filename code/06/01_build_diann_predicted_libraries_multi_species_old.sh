#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Build DIA-NN predicted spectral libraries for multiple species
#
# Default: process all supported species found in FASTA_DIR.
# Optional: specify species on command line, e.g.
#   bash build_diann_predicted_libraries_multi_species.sh dmel dsim
#
# Expected FASTA names, for example:
#   dmel-all-CDS-r6.12.longest_CDS_per_gene.protein.fasta
#   dsim-all-CDS-r2.02.longest_CDS_per_gene.protein.fasta
#   dsec-all-CDS-r1.3.longest_CDS_per_gene.protein.fasta
#   dere-all-CDS-r1.05.longest_CDS_per_gene.protein.fasta
# ============================================================

# ---------- Paths ----------
DIANN_DIR="${DIANN_DIR:-$HOME/software/diann-2.6.1/}"
DIANN_DIR="${DIANN_DIR%/}/"
DIANN="${DIANN:-${DIANN_DIR}diann-linux}"

FASTA_DIR="${FASTA_DIR:-/data/dongxn/MS}"
OUT_ROOT="${OUT_ROOT:-/data/dongxn/MS/DIA_NN_RUN}"

# ---------- DIA-NN settings ----------
THREADS=24
MISSED_CLEAVAGES=1
MIN_PEP_LEN=7
MAX_PEP_LEN=40
MIN_PR_CHARGE=2
MAX_PR_CHARGE=4

# Fixed Cys carbamidomethylation and variable Met oxidation.
MAX_VAR_MODS=1

# ---------- Supported species ----------
SUPPORTED_SPECIES=(dmel dsim dsec dere)

species_dir_name() {
    case "$1" in
        dmel) printf '%s\n' "Dmel" ;;
        dsim) printf '%s\n' "Dsim" ;;
        dsec) printf '%s\n' "Dsec" ;;
        dere) printf '%s\n' "Dere" ;;
        *)
            echo "ERROR: unsupported species: $1" >&2
            return 1
            ;;
    esac
}

is_supported_species() {
    local query="$1"
    local species
    for species in "${SUPPORTED_SPECIES[@]}"; do
        if [[ "$query" == "$species" ]]; then
            return 0
        fi
    done
    return 1
}

# ---------- Basic checks ----------
if [[ ! -x "$DIANN" ]]; then
    echo "ERROR: DIA-NN is missing or not executable:" >&2
    echo "  $DIANN" >&2
    exit 1
fi

if [[ ! -d "$FASTA_DIR" ]]; then
    echo "ERROR: FASTA directory does not exist:" >&2
    echo "  $FASTA_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_ROOT"

# ---------- Decide which species to process ----------
if [[ "$#" -gt 0 ]]; then
    SPECIES_TO_RUN=("$@")
else
    SPECIES_TO_RUN=("${SUPPORTED_SPECIES[@]}")
fi

for species in "${SPECIES_TO_RUN[@]}"; do
    if ! is_supported_species "$species"; then
        echo "ERROR: unsupported species '$species'." >&2
        echo "Supported species: ${SUPPORTED_SPECIES[*]}" >&2
        exit 1
    fi
done

# ---------- Build one library per species ----------
for species in "${SPECIES_TO_RUN[@]}"; do
    out_species=$(species_dir_name "$species")
    out_dir="${OUT_ROOT}/${out_species}"
    mkdir -p "$out_dir"

    # There should normally be exactly one matching longest-CDS protein FASTA.
    mapfile -t fasta_matches < <(
        find "$FASTA_DIR" \
            -maxdepth 1 \
            -type f \
            -name "${species}-all-CDS-*.longest_CDS_per_gene.protein.fasta" \
            -print | sort
    )

    if [[ "${#fasta_matches[@]}" -eq 0 ]]; then
        echo "WARNING: no FASTA found for $species; skipping." >&2
        echo "Expected pattern:" >&2
        echo "  ${FASTA_DIR}/${species}-all-CDS-*.longest_CDS_per_gene.protein.fasta" >&2
        continue
    fi

    if [[ "${#fasta_matches[@]}" -gt 1 ]]; then
        echo "ERROR: multiple FASTA files found for $species:" >&2
        printf '  %s\n' "${fasta_matches[@]}" >&2
        echo "Please remove old versions or specify one version in the script." >&2
        exit 1
    fi

    fasta="${fasta_matches[0]}"
    fasta_base=$(basename "$fasta" .fasta)
    output_lib="${out_dir}/${fasta_base}.predicted.speclib"
    log_file="${out_dir}/${fasta_base}.prediction.log"

    echo "============================================================"
    echo "Species:        $species"
    echo "FASTA:          $fasta"
    echo "Output library: $output_lib"
    echo "Log:            $log_file"
    echo "============================================================"

    if [[ ! -r "$fasta" ]]; then
        echo "ERROR: FASTA is not readable: $fasta" >&2
        exit 1
    fi

    # Remove a previous library with exactly the same output name.
    predicted_lib="${output_lib%.speclib}.predicted.speclib"
    rm -f "$output_lib" "$predicted_lib"

    "$DIANN" \
        --parent "$DIANN_DIR" \
        --fasta "$fasta" \
        --fasta-search \
        --predictor \
        --gen-spec-lib \
        --out-lib "$output_lib" \
        --cut 'K*,R*' \
        --missed-cleavages "$MISSED_CLEAVAGES" \
        --min-pep-len "$MIN_PEP_LEN" \
        --max-pep-len "$MAX_PEP_LEN" \
        --min-pr-charge "$MIN_PR_CHARGE" \
        --max-pr-charge "$MAX_PR_CHARGE" \
        --unimod4 \
        --var-mods "$MAX_VAR_MODS" \
        --var-mod 'UniMod:35,15.994915,M' \
        --ids-to-names \
        --pg-level 0 \
        --threads "$THREADS" \
        --verbose 2 \
        2>&1 | tee "$log_file"

    # DIA-NN 2.6.1 writes .predicted.predicted.speclib for this --out-lib value.
    if [[ -s "$predicted_lib" ]]; then
        ln -sfn "$(basename "$predicted_lib")" "$output_lib"
    fi
    if [[ ! -s "$output_lib" ]]; then
        echo "ERROR: predicted library was not created or is empty:" >&2
        echo "  $output_lib" >&2
        exit 1
    fi

    echo "Finished: $output_lib"
    echo

done

echo "All requested species finished."
