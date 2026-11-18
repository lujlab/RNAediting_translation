#!/bin/bash
#
# Exon and intron read coverage from bedtools genomecov, bedtools intersect,
# and coverage.py.
# Required variables:
#   GTF_NO_EXON: exon-coordinate GTF
#   GTF_INTRON: intron-coordinate GTF
#   GENOME_SIZES: chromosome sizes table
#   BAM_PATH: directory containing BAM files

# Temporary directory for intermediate files.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TEMP=$(mktemp -d "${TMPDIR:-/tmp}/druid.XXXXXX")

# create temp strand specific gtf files
GTF_NO_EXON_NEG=${TEMP}/exon.neg.gtf
GTF_NO_EXON_POS=${TEMP}/exon.pos.gtf
cat $GTF_NO_EXON | grep -P '\s-\s' | sort -k1,1 -k4,4n > $GTF_NO_EXON_NEG
cat $GTF_NO_EXON | grep -P '\s\+\s' | sort -k1,1 -k4,4n > $GTF_NO_EXON_POS

GTF_INTRON_NEG=${TEMP}/intron.neg.gtf
GTF_INTRON_POS=${TEMP}/intron.pos.gtf
cat $GTF_INTRON | grep -P '\s-\s' | sort -k1,1 -k4,4n > $GTF_INTRON_NEG
cat $GTF_INTRON | grep -P '\s\+\s' | sort -k1,1 -k4,4n > $GTF_INTRON_POS

mkdir intersect

# Coverage by BAM file.
for BAM in ${BAM_PATH}/*bam; do
    BASE=$(basename $BAM)

    # create temporary bam file containing only unique mappers
    samtools view -h $BAM |
        grep -P "^@|NH:i:1" |
        samtools view -bS - > ${TEMP}/${BASE}

    # create a bedgraph file for both positive and negative strands
    bedtools genomecov -ibam ${TEMP}/${BASE} -bg -split -strand - |
        sort -k1,1 -k2,2n > ${TEMP}/${BASE%bam}neg.bedgraph
    bedtools genomecov -ibam ${TEMP}/${BASE} -bg -split -strand + |
        sort -k1,1 -k2,2n > ${TEMP}/${BASE%bam}pos.bedgraph
    rm ${TEMP}/${BASE}

    # determine intersection between bedgraph and gtfs
    bedtools intersect -a $GTF_NO_EXON_NEG \
        -b ${TEMP}/${BASE%bam}pos.bedgraph -wo -sorted \
        -g $GENOME_SIZES > ${TEMP}/${BASE}_exon_pos.intersect
    "${PYTHON_BIN:-python3}" "$SCRIPT_DIR/coverage.py" ${TEMP}/${BASE}_exon_pos.intersect > \
        intersect/${BASE%.bam}_exon_pos.csv
    rm ${TEMP}/${BASE}_exon_pos.intersect
    bedtools intersect -a $GTF_NO_EXON_POS \
        -b ${TEMP}/${BASE%bam}neg.bedgraph -wo -sorted \
        -g $GENOME_SIZES > ${TEMP}/${BASE}_exon_neg.intersect
    "${PYTHON_BIN:-python3}" "$SCRIPT_DIR/coverage.py" ${TEMP}/${BASE}_exon_neg.intersect > \
        intersect/${BASE%.bam}_exon_neg.csv
    rm ${TEMP}/${BASE}_exon_neg.intersect

    bedtools intersect -a $GTF_INTRON_NEG \
        -b ${TEMP}/${BASE%bam}pos.bedgraph -wo -sorted \
        -g $GENOME_SIZES > ${TEMP}/${BASE}_intron_pos.intersect
    "${PYTHON_BIN:-python3}" "$SCRIPT_DIR/coverage.py" ${TEMP}/${BASE}_intron_pos.intersect > \
        intersect/${BASE%.bam}_intron_pos.csv
    rm ${TEMP}/${BASE}_intron_pos.intersect
    bedtools intersect -a $GTF_INTRON_POS \
        -b ${TEMP}/${BASE%bam}neg.bedgraph -wo -sorted \
        -g $GENOME_SIZES > ${TEMP}/${BASE}_intron_neg.intersect
    "${PYTHON_BIN:-python3}" "$SCRIPT_DIR/coverage.py" ${TEMP}/${BASE}_intron_neg.intersect > \
        intersect/${BASE%.bam}_intron_neg.csv
    rm ${TEMP}/${BASE}_intron_neg.intersect

    rm ${TEMP}/${BASE%bam}pos.bedgraph ${TEMP}/${BASE%bam}neg.bedgraph
done

# remove strand specific gtf files
rm $GTF_NO_EXON_NEG $GTF_NO_EXON_POS $GTF_INTRON_NEG $GTF_INTRON_POS
