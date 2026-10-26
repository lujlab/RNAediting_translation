#!/usr/bin/env bash
set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"

# Ribo-seq processing workflow for one sample.
# Forward-only transcriptome alignments for P/A/E extraction.
# Replace all values enclosed in <...> before running the script.

# ============================================================================
# 00. Input files and references
# ============================================================================

SAMPLE="${SAMPLE:-<sample>}"
R1="${R1:-<sample>_1.fq.gz}"
R2="${R2:-<sample>_2.fq.gz}"

OUT_DIR="${OUT_DIR:-${PIPELINE_DIR}/results/${SAMPLE}}"
TRNA_INDEX="${TRNA_INDEX:-<path/to/tRNA_bowtie2_index>}"
RRNA_INDEX="${RRNA_INDEX:-<path/to/rRNA_bowtie2_index>}"
STAR_INDEX="${STAR_INDEX:-<path/to/STAR_genome_index>}"
GTF="${GTF:-<path/to/genome_annotation.gtf>}"
TRANSCRIPT_FASTA="${TRANSCRIPT_FASTA:-<path/to/transcripts.fa>}"
TXINFO="${TXINFO:-<path/to/transcript_info.tsv>}"
SCRIPT_DIR="${PIPELINE_DIR}/code/01"

THREADS=16
MIN_LENGTH=25
RPF_MIN_LENGTH=26
RPF_MAX_LENGTH=34

mkdir -p \
    "${OUT_DIR}/00_umi" \
    "${OUT_DIR}/01_fastp" \
    "${OUT_DIR}/02_non_rRNA_tRNA" \
    "${OUT_DIR}/03_alignment" \
    "${OUT_DIR}/04_deduplicated" \
    "${OUT_DIR}/05_psite" \
    "${OUT_DIR}/06_counts" \
    "${OUT_DIR}/logs"

# ============================================================================
# 01. Extract UMIs from paired-end reads
# ============================================================================

RA3="TGGAATTCTCGGGTGCCAAGG"
RA5_RC="GATCGTCGGACTGTAGAACTCTGAAC"
RA5_UMI_BLOCK="(?:A.{3}TAC.{3}|C.{3}AGT.{3})"

# Define regexes separately: their quantifier braces must not terminate a
# Bash ${variable:-default} expansion. Explicit caller overrides are retained.
umi_default_r1="^(?:(?P<discard_ra5_umi_part>.{0,40}?${RA5_UMI_BLOCK}))?(?P<insert_seq>.{0,140}?)(?P<umi_1>.{3})(?P<discard_s7>.{0,3})(?P<umi_2>.{3})(?P<discard_s8>.{0,3})(?P<umi_3>.{3})(?P<discard_3p>${RA3}).*"
umi_default_r2="^(?:(?P<discard_ra3_umi_part>.{0,40}?A.{3}GAC.{3}))?(?P<insert_seq>.{0,140}?)(?P<umi_4>.{3})(?P<discard_s3>.{0,3})(?P<umi_5>.{3})(?P<discard_s4>.{0,3})(?P<umi_6>.{3})(?P<discard_5p_rc>${RA5_RC}).*"
BC_R1="${BC_R1:-$umi_default_r1}"
BC_R2="${BC_R2:-$umi_default_r2}"

umi_tools extract \
    --extract-method=regex \
    --ignore-read-pair-suffixes \
    --bc-pattern="${BC_R1}" \
    --bc-pattern2="${BC_R2}" \
    --stdin="${R1}" \
    --read2-in="${R2}" \
    --stdout >(pigz -p "${THREADS}" > "${OUT_DIR}/00_umi/${SAMPLE}.R1.extracted.fq.gz") \
    --read2-out >(pigz -p "${THREADS}" > "${OUT_DIR}/00_umi/${SAMPLE}.R2.extracted.fq.gz") \
    --filtered-out >(pigz -p "${THREADS}" > "${OUT_DIR}/00_umi/${SAMPLE}.R1.unmatched.fq.gz") \
    --filtered-out2 >(pigz -p "${THREADS}" > "${OUT_DIR}/00_umi/${SAMPLE}.R2.unmatched.fq.gz") \
    --log="${OUT_DIR}/logs/${SAMPLE}.umi_extract.log"

# ============================================================================
# 02. Recover unmatched R1 reads as single-end reads
# ============================================================================

umi_tools extract \
    --extract-method=regex \
    --bc-pattern="${BC_R1}" \
    --stdin="${OUT_DIR}/00_umi/${SAMPLE}.R1.unmatched.fq.gz" \
    --stdout >(awk 'NR % 4 == 1 {$0=$0 "/SE"} 1' |
        pigz -p "${THREADS}" > "${OUT_DIR}/00_umi/${SAMPLE}.R1.SE.fq.gz") \
    --log="${OUT_DIR}/logs/${SAMPLE}.umi_recover_R1.log"

# ============================================================================
# 03. Adapter trimming and read-quality filtering
# ============================================================================

# Strictly extracted paired-end reads
fastp \
    -i "${OUT_DIR}/00_umi/${SAMPLE}.R1.extracted.fq.gz" \
    -I "${OUT_DIR}/00_umi/${SAMPLE}.R2.extracted.fq.gz" \
    -o "${OUT_DIR}/01_fastp/${SAMPLE}.R1.clean.fq.gz" \
    -O "${OUT_DIR}/01_fastp/${SAMPLE}.R2.clean.fq.gz" \
    --adapter_sequence AACTCCAGTCAC \
    --adapter_sequence_r2 GTGTAGATCTCGGTGGTCGCAGTATCATT \
    --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --length_required "${MIN_LENGTH}" --trim_poly_g --trim_poly_x \
    --thread "${THREADS}" \
    --html "${OUT_DIR}/01_fastp/${SAMPLE}.PE.fastp.html" \
    --json "${OUT_DIR}/01_fastp/${SAMPLE}.PE.fastp.json"

# Recovered single-end R1 reads
fastp \
    -i "${OUT_DIR}/00_umi/${SAMPLE}.R1.SE.fq.gz" \
    -o "${OUT_DIR}/01_fastp/${SAMPLE}.R1.SE.clean.fq.gz" \
    --adapter_sequence AACTCCAGTCAC \
    --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --length_required "${MIN_LENGTH}" --trim_poly_g --trim_poly_x \
    --thread "${THREADS}" \
    --html "${OUT_DIR}/01_fastp/${SAMPLE}.SE.fastp.html" \
    --json "${OUT_DIR}/01_fastp/${SAMPLE}.SE.fastp.json"

# ============================================================================
# 04. Remove reads mapped to tRNA and rRNA
# ============================================================================

# Paired-end library: only R1 is used for downstream RPF analysis.
bowtie2 -p "${THREADS}" --local --no-unal \
    -x "${TRNA_INDEX}" \
    -U "${OUT_DIR}/01_fastp/${SAMPLE}.R1.clean.fq.gz" \
    --un-gz "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.PE.no_tRNA.fq.gz" \
    -S /dev/null \
    2> "${OUT_DIR}/logs/${SAMPLE}.PE.tRNA.log"

bowtie2 -p "${THREADS}" --local --no-unal \
    -x "${RRNA_INDEX}" \
    -U "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.PE.no_tRNA.fq.gz" \
    --un-gz "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.PE.no_tRNA_rRNA.fq.gz" \
    -S /dev/null \
    2> "${OUT_DIR}/logs/${SAMPLE}.PE.rRNA.log"

# Recovered single-end library
bowtie2 -p "${THREADS}" --local --no-unal \
    -x "${TRNA_INDEX}" \
    -U "${OUT_DIR}/01_fastp/${SAMPLE}.R1.SE.clean.fq.gz" \
    --un-gz "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.SE.no_tRNA.fq.gz" \
    -S /dev/null \
    2> "${OUT_DIR}/logs/${SAMPLE}.SE.tRNA.log"

bowtie2 -p "${THREADS}" --local --no-unal \
    -x "${RRNA_INDEX}" \
    -U "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.SE.no_tRNA.fq.gz" \
    --un-gz "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.SE.no_tRNA_rRNA.fq.gz" \
    -S /dev/null \
    2> "${OUT_DIR}/logs/${SAMPLE}.SE.rRNA.log"

# ============================================================================
# 05. Align RPF reads to the genome and transcriptome with STAR
# ============================================================================

# Paired-end-derived R1 reads
STAR \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.PE.no_tRNA_rRNA.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${OUT_DIR}/03_alignment/${SAMPLE}.PE." \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --outSAMprimaryFlag AllBestScore \
    --outMultimapperOrder Random \
    --outFilterMultimapNmax 10 \
    --twopassMode Basic \
    --quantMode TranscriptomeSAM \
    --outSAMunmapped Within \
    --runRNGseed 777

# Recovered single-end R1 reads
STAR \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${OUT_DIR}/02_non_rRNA_tRNA/${SAMPLE}.SE.no_tRNA_rRNA.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${OUT_DIR}/03_alignment/${SAMPLE}.SE." \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --outSAMprimaryFlag AllBestScore \
    --outMultimapperOrder Random \
    --outFilterMultimapNmax 10 \
    --twopassMode Basic \
    --quantMode TranscriptomeSAM \
    --outSAMunmapped Within \
    --runRNGseed 777

# ============================================================================
# 06. Retain uniquely mapped reads and remove UMI duplicates
# ============================================================================

# Paired-end-derived R1 reads
samtools view -h -@ "${THREADS}" \
    "${OUT_DIR}/03_alignment/${SAMPLE}.PE.Aligned.sortedByCoord.out.bam" |
    awk '$0 ~ /^@/ || $0 ~ /\tNH:i:1(\t|$)/' |
    samtools sort -@ "${THREADS}" \
        -o "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.unique.bam" -

samtools index -@ "${THREADS}" \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.unique.bam"

umi_tools dedup \
    -I "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.unique.bam" \
    -S "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.genome.dedup.bam" \
    --method=directional \
    --spliced-is-unique \
    --umi-separator=_ \
    --output-stats="${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.dedup_stats"

# Recovered single-end R1 reads
samtools view -h -@ "${THREADS}" \
    "${OUT_DIR}/03_alignment/${SAMPLE}.SE.Aligned.sortedByCoord.out.bam" |
    awk '$0 ~ /^@/ || $0 ~ /\tNH:i:1(\t|$)/' |
    samtools sort -@ "${THREADS}" \
        -o "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.unique.bam" -

samtools index -@ "${THREADS}" \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.unique.bam"

umi_tools dedup \
    -I "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.unique.bam" \
    -S "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.genome.dedup.bam" \
    --method=directional \
    --spliced-is-unique \
    --umi-separator=_ \
    --output-stats="${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.dedup_stats"

# ============================================================================
# 07. Apply the genome deduplication result to transcriptome BAM files
# ============================================================================

# Paired-end-derived R1 reads
samtools view "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.genome.dedup.bam" |
    cut -f1 | sort -u \
    > "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.kept_read_names.txt"

samtools view -@ "${THREADS}" \
    -N "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.kept_read_names.txt" \
    -u "${OUT_DIR}/03_alignment/${SAMPLE}.PE.Aligned.toTranscriptome.out.bam" |
    samtools sort -@ "${THREADS}" \
        -o "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.transcript.dedup.bam" -

# Recovered single-end R1 reads
samtools view "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.genome.dedup.bam" |
    cut -f1 | sort -u \
    > "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.kept_read_names.txt"

samtools view -@ "${THREADS}" \
    -N "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.kept_read_names.txt" \
    -u "${OUT_DIR}/03_alignment/${SAMPLE}.SE.Aligned.toTranscriptome.out.bam" |
    samtools sort -@ "${THREADS}" \
        -o "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.transcript.dedup.bam" -

# ============================================================================
# 08. Merge paired-end-derived and recovered single-end alignments
# ============================================================================

samtools merge -@ "${THREADS}" -u - \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.genome.dedup.bam" \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.genome.dedup.bam" |
    samtools sort -@ "${THREADS}" \
        -o "${OUT_DIR}/04_deduplicated/${SAMPLE}.genome.merged.sorted.bam" -

samtools index -@ "${THREADS}" \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.genome.merged.sorted.bam"

samtools merge -@ "${THREADS}" -u - \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.PE.transcript.dedup.bam" \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.SE.transcript.dedup.bam" |
    samtools sort -@ "${THREADS}" \
        -o "${OUT_DIR}/04_deduplicated/${SAMPLE}.transcript.merged.sorted.bam" -

samtools index -@ "${THREADS}" \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.transcript.merged.sorted.bam"

# ============================================================================
# 09. Predict P-sites and count P/A/E-site codons
# ============================================================================

python3 "${SCRIPT_DIR}/train_mod.py" \
    -i -t cds_longest \
    "${TRANSCRIPT_FASTA}" \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.transcript.merged.sorted.bam" \
    "${OUT_DIR}/05_psite/${SAMPLE}" \
    "${TXINFO}"

python3 "${SCRIPT_DIR}/predict_mod.py" \
    --rlen_min "${RPF_MIN_LENGTH}" \
    --rlen_max "${RPF_MAX_LENGTH}" \
    -i \
    "${TRANSCRIPT_FASTA}" \
    "${OUT_DIR}/04_deduplicated/${SAMPLE}.transcript.merged.sorted.bam" \
    "${OUT_DIR}/05_psite/${SAMPLE}.gbt.pickle" \
    "${OUT_DIR}/05_psite/${SAMPLE}.tx_startcodon.tsv" \
    "${OUT_DIR}/05_psite/${SAMPLE}.RPF.transcriptome.tag.bam"

python3 "${SCRIPT_DIR}/extract_PAE_codon_from_tagged_bam_forward.py" \
    "${OUT_DIR}/05_psite/${SAMPLE}.RPF.transcriptome.tag.bam" \
    "${TXINFO}" \
    "${OUT_DIR}/05_psite/${SAMPLE}.PAE_codon.tsv" \
    8 50000

perl "${SCRIPT_DIR}/count_PAE_codon.pl" \
    "${OUT_DIR}/05_psite/${SAMPLE}" \
    < "${OUT_DIR}/05_psite/${SAMPLE}.PAE_codon.tsv"

# ============================================================================
# 10A. Figure S3: prepare transcript annotations
# ============================================================================
# Configure the genome GTF directory; this script prepares all four species.
export REF_DIR="${REF_DIR:-<path/to/reference_gtf_directory>}"
export FS3_ANNOTATION_DIR="${FS3_ANNOTATION_DIR:-${PIPELINE_DIR}/results/01/FS3_annotation}"
"${RSCRIPT:-Rscript}" "${PIPELINE_DIR}/code/01/ribowaltz_ann.R"

# ============================================================================
# 10B. Figure S3: QC plots from existing PS/XF-tagged transcriptome BAMs
# ============================================================================
# Run separately for each species after preparing its annotation and sample table.
# Sample table columns: sample, group, bam; one biological replicate per row.
FS3_SAMPLES="${FS3_SAMPLES:-<path/to/species_samples.tsv>}"
FS3_ANNOTATION="${FS3_ANNOTATION:-<path/to/species_annotation.tsv>}"
FS3_OUTPUT_DIR="${FS3_OUTPUT_DIR:-${PIPELINE_DIR}/results/01/FS3}"

"${RSCRIPT:-Rscript}" "${PIPELINE_DIR}/code/01/FS3/run_psite_tagged.R" \
    "${FS3_SAMPLES}" "${FS3_ANNOTATION}" "${FS3_OUTPUT_DIR}"
