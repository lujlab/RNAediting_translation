#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# RNA/Ribo codon counts, editing-level comparisons, and Tables 1, S4 and S5.
# Sections 01-09 run sequentially after configuration in section 00.
# Section 10 provides optional MS comparison commands.
# Existing unscaled alledit tables can be supplied to section 06.

# =============================================================================
# 00. Paths and environment
# =============================================================================
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
export SCRIPT_DIR="$PIPELINE_DIR/code/04"
SITE_CODE_DIR="$PIPELINE_DIR/code/04"
EDGE_CODE="$SITE_CODE_DIR/Table1"
SCALE_CODE="$SITE_CODE_DIR/Table1"
export REF_ROOT="${REF_ROOT:-$PIPELINE_DIR/..}"
export SOURCE_PROJECT="${SOURCE_PROJECT:-/path/to/Proj1_editing/codex/rescue_3R_18907675_1v1_pipeline_clean}"
export CANDIDATE_ROOT="${CANDIDATE_ROOT:-$PIPELINE_DIR/results/03_editing_sites_test2_union}"
export LONGEST_MODE="${LONGEST_MODE:-counter}"
export MS_ROOT="${MS_ROOT:-$SOURCE_PROJECT/../03_MS}"
OUT_ROOT="${OUT_ROOT:-$PIPELINE_DIR/results/04_RNA_Ribo_alledit_Table1}"
export CDSPOS_DIR="${CDSPOS_DIR:-$OUT_ROOT/01_cdspos}"
EDITED_CDS_DIR="${EDITED_CDS_DIR:-$OUT_ROOT/02_edited_CDS}"
COUNT_DIR="${COUNT_DIR:-$OUT_ROOT/02_counts}"
export ALLEDIT_DIR="${ALLEDIT_DIR:-$OUT_ROOT/03_alledit}"
PAIR_DIR="${PAIR_DIR:-$REF_ROOT/10.stalling}"
RNA_COUNT_DIR="${RNA_COUNT_DIR:-$REF_ROOT/RNA_processing_JWZ/06_codon_count}"
RIBO_COUNT_DIR="${RIBO_COUNT_DIR:-$REF_ROOT/05.psite/psite_2times}"
MAX_JOBS="${MAX_JOBS:-12}"
export PYTHON="${PYTHON:-python3}"
export RSCRIPT="${RSCRIPT:-Rscript}"
export R_LIBS_USER="${R05_LIBRARY:-$PIPELINE_DIR/.runtime/05/R_library}:${R_LIBS_USER:-/path/to/Proj1_editing/codex/02_stalling/3flank/peak_find/R_library}"
PLOT_PYTHON="${PLOT_PYTHON:-/usr/bin/python3}"
TABLE1_WORK_ROOT="${TABLE1_WORK_ROOT:-$PIPELINE_DIR/results/04_Table1_TableS4S5_RNAdiv2}"
TABLE1_NAME=Table1_edgeR_QL_Fisher_minrep2_pairedcoverage
export TABLE1_OUT="${TABLE1_OUT:-$OUT_ROOT/05_Table1_edgeR_QL_Fisher_minrep2_pairedcoverage}"
RNA_SCALE_MANIFEST="${RNA_SCALE_MANIFEST:-$PIPELINE_DIR/config/rna_count_scaling.tsv}"
COMPARE_THREADS="${COMPARE_THREADS:-8}"
MIN_REPLICATES="${MIN_REPLICATES:-2}"

# =============================================================================
# 01. Project candidates onto representative transcripts
# =============================================================================
# Use workflow 03 species-union candidates without additional DP/MQ thresholds.
# Rank transcripts by CDS length, transcript length, then transcript ID.
# Exclude failed coordinate projections; count-level source-transcript fallback is separate.
CDSOUT="$CDSPOS_DIR"
AUDIT_DIR="$OUT_ROOT/01_cdspos_audit"
TMP="$OUT_ROOT/01_cdspos_tmp"
LONGEST_LIST_DIR="${LONGEST_LIST_DIR:-${AUDIT_DIR}/lists}"
PROJECT_GROUPS="${PROJECT_GROUPS:-ISO1 sim sec ere WBG AG AS}"
case "$LONGEST_MODE" in counter | legacy) ;; *)
    echo 'LONGEST_MODE=counter or legacy' >&2
    exit 2
    ;;
esac
mkdir -p "$CDSOUT" "$AUDIT_DIR" "$LONGEST_LIST_DIR"
for group in $PROJECT_GROUPS; do
    case "$group" in
        ISO1 | WBG | AG | AS)
            sp=dmel
            tag=dmel_r6.12
            txinfo=dmel_r6.12
            gtf=dmel-all-r6.12
            ;;
        sim)
            sp=dsim
            tag=dsim_r2.02
            txinfo=dsim_all_r2.02
            gtf=dsim-all-r2.02
            ;;
        sec)
            sp=dsec
            tag=dsec_r1.3
            txinfo=dsec_all_r1.3
            gtf=dsec-all-r1.3
            ;;
        ere)
            sp=dere
            tag=dere_r1.05
            txinfo=dere_all_r1.05
            gtf=dere-all-r1.05
            ;;
        *)
            echo "Unknown group $group" >&2
            exit 2
            ;;
    esac
    input="$CANDIDATE_ROOT/raw_rows_by_group/$group.candidate.AtoG.annotations.long.tsv"
    [[ -s "$input" ]] || {
        echo "Missing current candidate rows: $input" >&2
        exit 1
    }
    longest="$LONGEST_LIST_DIR/$tag.ids.txt"
    if [[ "$LONGEST_MODE" == counter ]]; then
        # Same ties as the existing 01/02 counters and the previous 05 default.
        awk 'BEGIN{FS=OFS="\t"} NR>1 {
       gene=$1;tx=$2;txlen=$5+0;cds=$6+0;if(cds<=0)next
       if(!(gene in seen)||cds>best[gene]||(cds==best[gene]&&txlen>len[gene])||(cds==best[gene]&&txlen==len[gene]&&tx<sel[gene]))
       {seen[gene]=1;best[gene]=cds;len[gene]=txlen;sel[gene]=tx}}
       END{for(g in sel)print sel[g]}' "$REF_ROOT/05.psite/${txinfo}_txinfo.biotype.tsv" | sort > "$longest"
    else
        # Exact longest list referenced in run_sites.sh.
        cp "$REF_ROOT/05.psite/${tag}_longestCDS_transcript_ids.txt" "$longest"
    fi
    prefix="$CDSOUT/$group.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3"
    "$PYTHON" "$SITE_CODE_DIR/rebuild_candidate_cds_outputs_strict.py" \
        --input "$input" --group "$group" --gtf "$REF_ROOT/ref_genome/$gtf.gtf" \
        --longest-list "$longest" --fallback exclude \
        --map-script "$SITE_CODE_DIR/map_transcript_coordinates_v6.py" \
        --longest-script "$SITE_CODE_DIR/add_longest_and_deduplicate.py" \
        --output-cdspos "$prefix.targetCDS_only.tmp.tsv" \
        --output-full "$prefix.full_with_longest_projection.tsv" \
        --summary "$prefix.strict_projection.summary.tsv" --tmp-dir "${TMP:-$AUDIT_DIR/tmp}" \
        > "$AUDIT_DIR/$group.projection.log" 2>&1
    "$PYTHON" "$SITE_CODE_DIR/make_cdspos_keep_target_region_from_full_projection.py" \
        --input-full "$prefix.full_with_longest_projection.tsv" \
        --output "$prefix.cdspos_for_alledit.tsv" --summary "$prefix.keepTargetRegion.summary.tsv"
    echo "Completed projection: $group"
done

# =============================================================================
# 02. Generate edited CDS sequences and codon annotations
# =============================================================================
# Invalid or skipped edits stop the workflow before downstream counting.
EDIT_GROUPS="${EDIT_GROUPS:-ISO1 sim sec ere}"
# Check before creating any FASTA; run the original Biopython editor without --force.
"$PYTHON" -c 'from Bio import SeqIO; from Bio.Seq import Seq'
command -v perl > /dev/null
mkdir -p "$EDITED_CDS_DIR"
for group in $EDIT_GROUPS; do
    case "$group" in
        ISO1) name=dmel-all-CDS-r6.12.fasta ;;
        sim) name=dsim-all-CDS-r2.02.fasta ;;
        sec) name=dsec-all-CDS-r1.3.fasta ;;
        ere) name=dere-all-CDS-r1.05.fasta ;;
        *)
            echo "Unknown edit group: $group" >&2
            exit 2
            ;;
    esac
    fasta="$REF_ROOT/ref_genome/$name"
    sites="$CDSPOS_DIR/$group.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"
    output="$EDITED_CDS_DIR/$name.edited"
    [[ -s "$fasta" && -s "$sites" ]] || {
        echo "Missing FASTA or sites for $group" >&2
        exit 1
    }
    "$PYTHON" "$SITE_CODE_DIR/apply_editing_to_fasta_new_format.py" "$fasta" "$sites" "$output" \
        --coord-column new_cds_pos_1based > "$EDITED_CDS_DIR/$group.edit.log" 2>&1
    # Keep the original Python -> sed -> Perl codon extraction order.
    "$PYTHON" "$SITE_CODE_DIR/split_CDS.py" "$output" | sed 's/,/\t/1' | sed 's/;//1' > "$output.tab"
    perl "$SITE_CODE_DIR/split_protein_codons.pl" "$output.tab" > "$output.tab.split"
    "$PYTHON" "$SITE_CODE_DIR/ensure_edited_codons.py" --reference "$fasta" --edited "$output" \
        --sites "$sites" --codons "$output.tab.split"
    echo "Completed edited FASTA and codon extraction: $group"
done

# =============================================================================
# 03. Export site paths and prepare four-column RNA/Ribo counts
# =============================================================================
"$PYTHON" "$SITE_CODE_DIR/write_sites_handoff.py" \
    --cdspos-dir "$CDSPOS_DIR" \
    --edited-dir "$EDITED_CDS_DIR" \
    --candidate-root "$CANDIDATE_ROOT" \
    --alledit-dir "$ALLEDIT_DIR" \
    --out "$OUT_ROOT/sites_for_06.env"

# Read existing counts from workflows 01 and 02 without rounding.
# Exclude UTR coordinates and write four-column count files and a sample manifest.
"$PYTHON" "$SCRIPT_DIR/prepare_alledit_count_inputs.py" \
    --out-dir "$COUNT_DIR" \
    --groups as wbg ag oldDmel sim sec ere \
    --site-types Asite Psite Esite \
    --pair-dir "$PAIR_DIR" \
    --rna-dir "$RNA_COUNT_DIR" \
    --ribo-dir "$RIBO_COUNT_DIR"

# Alternative input: an explicit sample-pair manifest.
# "$PYTHON" "$SCRIPT_DIR/prepare_alledit_count_inputs.py" \
#   --out-dir "$COUNT_DIR" --groups as wbg ag oldDmel sim sec ere \
#   --site-types Asite Psite Esite --manifest /absolute/path/pairs.tsv

# =============================================================================
# 04. Quantify editing states for each sample and A/P/E position
# =============================================================================
# Use target-transcript RNA coverage, with source-transcript fallback where required.
# Write unfiltered alledit tables.
[[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]]
mkdir -p "$ALLEDIT_DIR" "$OUT_ROOT/logs/alledit"
pids=()
while IFS=$'\t' read -r group rna ribo site rf bf rn bn; do
    [[ "$group" == group ]] && continue
    case "$group" in
        as | wbg | ag | oldDmel) sp=ISO1 ;;
        sim | sec | ere) sp="$group" ;;
        *)
            echo "Unsupported group: $group" >&2
            exit 1
            ;;
    esac
    cdspos="$CDSPOS_DIR/$sp.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv"
    [[ -f "$rn" && -f "$bn" && -f "$cdspos" ]]
    "$RSCRIPT" "$SCRIPT_DIR/recompute_alledit_with_original_tx_fallback.R" \
        "$rn" "$bn" "$cdspos" "$ALLEDIT_DIR/${ribo}.${site}" \
        > "$OUT_ROOT/logs/alledit/${rna}__${ribo}.${site}.log" 2>&1 &
    pids+=("$!")
    if ((${#pids[@]} >= MAX_JOBS)); then
        failed=0
        for pid in ${pids[@]+"${pids[@]}"}; do wait "$pid" || failed=1; done
        [[ "$failed" == 0 ]]
        pids=()
    fi
done < "$COUNT_DIR/count_manifest.tsv"
# Bash 4.2 + nounset: the last batch may leave an empty PID array.
failed=0
for pid in ${pids[@]+"${pids[@]}"}; do wait "$pid" || failed=1; done
[[ "$failed" == 0 ]]

# =============================================================================
# 05. Generate I-position summaries and violins from unscaled alledit counts
# =============================================================================
# I-position panels retain RNA>=3, Ribo>=3 and edited RNA>=1.
# Scatter panels use the Table1 filter workflow below.
"$RSCRIPT" "$SCRIPT_DIR/replot_multispecies_alledit_counter_rule.R" \
    ALL "$OUT_ROOT/04_alledit_scatter"

# =============================================================================
# 06. Prepare RNA-count-corrected inputs for Tables 1, S4 and S5
# =============================================================================
# Sequential execution uses the alledit tables generated in section 04.
# Optional input override for existing unscaled alledit tables:
# INPUT_ALLEDIT_DIR="$SOURCE_PROJECT/04_alledit_new/results"
# Inputs must contain unscaled counts to avoid applying the RNA correction twice.
# Sections 07-09 share the corrected input tables.
INPUT_ALLEDIT_DIR="${INPUT_ALLEDIT_DIR:-$ALLEDIT_DIR}"
"$PYTHON" "$SCALE_CODE/prepare_rna_scaled_alledit.py" \
    --input "$INPUT_ALLEDIT_DIR" \
    --manifest "$RNA_SCALE_MANIFEST" \
    --out "$TABLE1_WORK_ROOT"

# =============================================================================
# 07. Table 1: edgeR inference and Fisher sensitivity comparisons
# =============================================================================
# Divide paired-end RNA allele counts by two; retain single-end RNA and Ribo counts.
# Require matched RNA/Ribo coverage>=5 in >=2 replicates and pooled RNA editing>=0.05.
# The editing-level threshold applies to pooled counts, not individual replicates.
# Omit uninformative rounded Fisher tables; combine each direction across informative replicates.
# Record significance in both directions as direction_conflict.
# Test A and P separately; AP is their union. Dmel combines ISO1 and WBg.
# The R script exports method-overlap plots.
mkdir -p "$TABLE1_WORK_ROOT/logs"
# Paired-end RNA correction.
table_out="$TABLE1_WORK_ROOT/RNAdiv2/$TABLE1_NAME"
"$RSCRIPT" --vanilla "$EDGE_CODE/Table1_edgeR_QL_Fisher.R" \
    --input-dir "$TABLE1_WORK_ROOT/RNAdiv2/alledit" \
    --output-dir "$table_out" \
    --threads "$COMPARE_THREADS" \
    --min-replicates "$MIN_REPLICATES" \
    --min-rna-per-rep-equiv 5 \
    --min-ribo-per-rep-equiv 5 \
    --min-lm 0.05 \
    --fdr 0.05 \
    > "$TABLE1_WORK_ROOT/logs/RNAdiv2.Table1_edgeR_QL_Fisher.log" 2>&1

"$RSCRIPT" --vanilla "$EDGE_CODE/export_edgeR_manuscript_Table1.R" "$table_out"
mkdir -p "$(dirname "$TABLE1_OUT")"
ln -sfnT "$TABLE1_WORK_ROOT/RNAdiv2/$TABLE1_NAME" "$TABLE1_OUT"

# =============================================================================
# 07b. Table1 filter: scatter panels using the Table 1 minrep=2 eligibility rules
# =============================================================================
RSCRIPT="$RSCRIPT" TABLE1_SCALED_ROOT="$TABLE1_WORK_ROOT/RNAdiv2" \
    TABLE1_RESULT_DIR="$table_out" TABLE1_FILTER_OUT="$OUT_ROOT/04_Table1_filter" \
    bash "$SCRIPT_DIR/table1_filter/run.sh"

# =============================================================================
# 08. Compare N/S composition of edgeR A-site I>A and remaining sites
# Plot P values are from two-sided Fisher tests of N/S composition.
# =============================================================================
# Paired-end RNA correction.
table_out="$TABLE1_WORK_ROOT/RNAdiv2/$TABLE1_NAME"
"$RSCRIPT" --vanilla "$EDGE_CODE/prepare_NS_plot.R" \
    "$table_out/Table1_edgeR_final.tsv" \
    "$table_out/31_ns_ratio_plot_data.edgeR.tsv"

"$PLOT_PYTHON" "$EDGE_CODE/plot_NS.py" \
    --data "$table_out/31_ns_ratio_plot_data.edgeR.tsv" \
    --y-axis-max 16 \
    --output "$table_out/31_ns_ratio_grouped_bar_chart_reference_lines.edgeR.pdf"

# =============================================================================
# 09. Tables S4/S5: coverage summaries from RNA-corrected counts
# =============================================================================
# Count covered sites independently of edgeR or Fisher significance.
# Use the species-specific CDS-position tables generated in section 01.
# Paired-end RNA correction.
s4_out="$TABLE1_WORK_ROOT/RNAdiv2/TableS4S5"
mkdir -p "$s4_out" "$TABLE1_WORK_ROOT/logs"
"$PYTHON" "$SCALE_CODE/count_riboseq_editing_sites_new_annotation.py" \
    --site-dir "$TABLE1_WORK_ROOT/RNAdiv2/alledit" \
    --annotation-dir "$CDSPOS_DIR" \
    --dmel-annotation "$CDSPOS_DIR/ISO1.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv" \
    --dsim-annotation "$CDSPOS_DIR/sim.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv" \
    --dsec-annotation "$CDSPOS_DIR/sec.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv" \
    --dere-annotation "$CDSPOS_DIR/ere.candidate.rescue_3R_18907675.sourceCDS.keepTargetUTR3.cdspos_for_alledit.tsv" \
    --output-prefix "$s4_out/riboseq_editing_site_counts" \
    > "$TABLE1_WORK_ROOT/logs/RNAdiv2.TableS4S5.log" 2>&1

"$PYTHON" "$SCALE_CODE/cal_ribo_rna_coverage_codon_num_percentage.py" \
    --input "$s4_out/riboseq_editing_site_counts.formatted.tsv" \
    --output "$s4_out/riboseq_editing_site_counts.CDS_percentage.tsv" \
    >> "$TABLE1_WORK_ROOT/logs/RNAdiv2.TableS4S5.log" 2>&1

mkdir -p "$OUT_ROOT"
ln -sfnT "$TABLE1_WORK_ROOT/RNAdiv2/TableS4S5" "$OUT_ROOT/06_TableS4S5"

# =============================================================================
# 10. Optional MS Lp versus Lm/Lr plots using workflow 06 outputs
# =============================================================================
# CODE_DIR="$PIPELINE_DIR/code/06" OUTPUT_ROOT="$OUT_ROOT/08_Lp_Lm_Lr_figures" \
#   bash "$PIPELINE_DIR/code/06/run_Lp_Lm_Lr_counter_rule_figures.sh"
