#!/usr/bin/env bash
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"

# RNA-seq processing workflow for one sample.
# Replace all values enclosed in <...> before running the script.
#
# Two library types were used in this study:
#   1) first-generation single-end RNA-seq;
#   2) second-generation strand-specific paired-end RNA-seq.
#
# Step 01 provides cutadapt examples for both library types. Steps 02-08 use
# the second-generation paired-end library as the runnable example.

# ============================================================================
# 00. Input files and references
# ============================================================================

# Use the sample ID listed in code/03/snpeff_samples.tsv (annotation path).
SAMPLE="${SAMPLE:-<sample>}"

# First-generation single-end RNA-seq input (used in the Step 01A example).
SE_READ="${SE_READ:-<sample>.fq.gz}"

# Second-generation strand-specific paired-end RNA-seq inputs.
R1="${R1:-<sample>_1.fq.gz}"
R2="${R2:-<sample>_2.fq.gz}"

OUT_DIR="${OUT_DIR:-${PIPELINE_DIR}/results/${SAMPLE}}"
TRNA_INDEX="${TRNA_INDEX:-<path/to/tRNA_bowtie2_index>}"
RRNA_INDEX="${RRNA_INDEX:-<path/to/rRNA_bowtie2_index>}"
STAR_INDEX="${STAR_INDEX:-<path/to/STAR_genome_index>}"
GENOME_FASTA="${GENOME_FASTA:-<path/to/genome.fa>}"
GTF="${GTF:-<path/to/genome_annotation.gtf>}"
TXINFO="${TXINFO:-<path/to/transcript_info.tsv>}"

# snpEff settings used for variant annotation.
SNPEFF_JAR="${SNPEFF_JAR:-<path/to/snpEff.jar>}"
SNPEFF_CONFIG="${SNPEFF_CONFIG:-<path/to/snpEff.config>}"
SNPEFF_DATABASE="${SNPEFF_DATABASE:-<snpEff_database_name>}"

SCRIPT_DIR="${PIPELINE_DIR}/code/02"

THREADS=16
WORKERS=4
MIN_MAPQ=0
MIN_BASEQ=13

mkdir -p \
    "${OUT_DIR}/00_cutadapt" \
    "${OUT_DIR}/01_fastqc" \
    "${OUT_DIR}/02_non_rRNA_tRNA" \
    "${OUT_DIR}/03_alignment" \
    "${OUT_DIR}/05_variants" \
    "${OUT_DIR}/06_codon_counts" \
    "${OUT_DIR}/reference" \
    "${OUT_DIR}/logs"

# ============================================================================
# 01. Adapter and quality trimming for the two RNA-seq library types
# ============================================================================

# ----------------------------------------------------------------------------
# 01A. First-generation single-end RNA-seq
#
# The first-generation libraries were single-end. The following command is an
# alternative to Step 01B. Uncomment it when processing a single-end library.
# ----------------------------------------------------------------------------

# cutadapt \
#     -j "${THREADS}" \
#     -a TGGAATTCTCGGGTGCCAAGG \
#     --minimum-length 18 \
#     -o "${OUT_DIR}/00_cutadapt/${SAMPLE}.SE.clean.fq.gz" \
#     "${SE_READ}" \
#     > "${OUT_DIR}/logs/${SAMPLE}.SE.cutadapt.log" 2>&1

# ----------------------------------------------------------------------------
# 01B. Second-generation strand-specific paired-end RNA-seq
#
# This is the library type used in the runnable example below. R1 and R2 were
# trimmed simultaneously. A complete read pair was discarded if either mate
# was shorter than 35 nt after trimming (--pair-filter=any).
# ----------------------------------------------------------------------------

cutadapt \
    -j "${THREADS}" \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 20,20 \
    --trim-n \
    --minimum-length 35 \
    --pair-filter=any \
    -o "${OUT_DIR}/00_cutadapt/${SAMPLE}.R1.clean.fq.gz" \
    -p "${OUT_DIR}/00_cutadapt/${SAMPLE}.R2.clean.fq.gz" \
    "${R1}" "${R2}" \
    > "${OUT_DIR}/logs/${SAMPLE}.cutadapt.log" 2>&1

# ============================================================================
# 02. Assess read quality (second-generation paired-end example)
# ============================================================================

fastqc \
    --threads "${THREADS}" \
    --outdir "${OUT_DIR}/01_fastqc" \
    "${OUT_DIR}/00_cutadapt/${SAMPLE}.R1.clean.fq.gz" \
    "${OUT_DIR}/00_cutadapt/${SAMPLE}.R2.clean.fq.gz"

# ============================================================================
# 03. Remove reads mapped to tRNA and rRNA (paired-end example)
# ============================================================================

bowtie2 \
    -p "${THREADS}" \
    -x "${TRNA_INDEX}" \
    -1 "${OUT_DIR}/00_cutadapt/${SAMPLE}.R1.clean.fq.gz" \
    -2 "${OUT_DIR}/00_cutadapt/${SAMPLE}.R2.clean.fq.gz" \
    --un-conc-gz "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.no_tRNA.fq.gz" \
    -S /dev/null \
    2> "${OUT_DIR}/logs/${SAMPLE}.tRNA.log"

bowtie2 \
    -p "${THREADS}" \
    -x "${RRNA_INDEX}" \
    -1 "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.no_tRNA.fq.1.gz" \
    -2 "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.no_tRNA.fq.2.gz" \
    --un-conc-gz "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.clean.fq.gz" \
    --al-conc-gz "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.rRNA.fq.gz" \
    -S /dev/null \
    2> "${OUT_DIR}/logs/${SAMPLE}.rRNA.log"

mv "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.clean.fq.1.gz" \
    "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.R1.fq.gz"
mv "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.clean.fq.2.gz" \
    "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.R2.fq.gz"

# ============================================================================
# 04. Align RNA-seq reads with STAR (paired-end example)
# ============================================================================

STAR \
    --runMode alignReads \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn \
    "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.R1.fq.gz" \
    "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.R2.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${OUT_DIR}/03_alignment/${SAMPLE}." \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --quantMode TranscriptomeSAM \
    --outFilterMultimapNmax 1

samtools index -@ "${THREADS}" \
    "${OUT_DIR}/03_alignment/${SAMPLE}.Aligned.sortedByCoord.out.bam"

# Gene-level fragment counting is performed at the start of workflow 08.

# ============================================================================
# 05A. Call variants from genome-aligned RNA-seq reads (paired-end example)
# ============================================================================

bcftools mpileup \
    -Ou \
    -f "${GENOME_FASTA}" \
    --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
    --max-depth 100000 \
    --threads "${THREADS}" \
    "${OUT_DIR}/03_alignment/${SAMPLE}.Aligned.sortedByCoord.out.bam" |
    bcftools call \
        --keep-alts \
        --multiallelic-caller \
        --ploidy 1 \
        -Oz \
        -o "${OUT_DIR}/05_variants/${SAMPLE}.vcf.gz"

bcftools index --tbi "${OUT_DIR}/05_variants/${SAMPLE}.vcf.gz"

# ============================================================================
# 05B. Remove records with missing REF or ALT alleles
# ============================================================================

# Preserve VCF headers and records with REF != "." and ALT != ".".
# No additional QUAL, DP or MQ threshold is applied at this step.
bcftools view -Ov "${OUT_DIR}/05_variants/${SAMPLE}.vcf.gz" |
    perl "${SCRIPT_DIR}/FilterVCF.pl" \
        > "${OUT_DIR}/05_variants/${SAMPLE}.filtered.vcf"

# ============================================================================
# 06. Annotate variants with snpEff (paired-end example)
# ============================================================================

java -Xmx8g -jar "${SNPEFF_JAR}" ann \
    -c "${SNPEFF_CONFIG}" \
    "${SNPEFF_DATABASE}" \
    "${OUT_DIR}/05_variants/${SAMPLE}.filtered.vcf" \
    > "${OUT_DIR}/05_variants/${SAMPLE}.annotated.vcf" \
    2> "${OUT_DIR}/logs/${SAMPLE}.snpEff.log"

# ANN annotations are stored in the VCF INFO column; no separate .ann file is needed.
# Workflow 03 reads results/<sample>/05_variants/<sample>.annotated.vcf through
# code/03/snpeff_samples.tsv; its default ANNOTATION_ROOT is code_upload/results.
# Continue with 04: parse -> filter/A>G -> liftover to dm6 -> selection -> plotting.
# Inputs must contain the required INFO/AD and INFO/MQ fields and use the intended
# upstream variant filters.

# ============================================================================
# 07. Generate transcript-coordinate CDS positions (paired-end example)
# ============================================================================

# This reference-preparation step only needs to be run once per annotation.
python3 "${SCRIPT_DIR}/gtf_to_transcript_CDS_position.py" \
    --gtf "${GTF}" \
    --out-prefix "${OUT_DIR}/reference/annotation"

# ============================================================================
# 08. Generate longest-CDS RNA codon counts (paired-end example)
# ============================================================================

python3 "${SCRIPT_DIR}/transcriptome_bam_codon_counter_longestCDS.py" \
    --bam "${OUT_DIR}/03_alignment/${SAMPLE}.Aligned.toTranscriptome.out.bam" \
    --cds-pos "${OUT_DIR}/reference/annotation.CDS_position.tsv" \
    --txinfo "${TXINFO}" \
    --select-mode longest_cds \
    --cds-start-base 1 \
    --out-prefix "${OUT_DIR}/06_codon_counts/${SAMPLE}.longestCDS" \
    --threads "${THREADS}" \
    --workers "${WORKERS}" \
    --min-mapq "${MIN_MAPQ}" \
    --min-baseq "${MIN_BASEQ}" \
    > "${OUT_DIR}/logs/${SAMPLE}.codon_count.stdout.log" \
    2> "${OUT_DIR}/logs/${SAMPLE}.codon_count.stderr.log"

echo "RNA-seq processing completed: ${SAMPLE}"
