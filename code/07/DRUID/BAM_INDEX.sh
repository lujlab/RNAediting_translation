#!/bin/bash

BAM_DIR="/path/to/RNA_editing_DATAprocessing/RNA_DRUID_analysis/Analysis/HEK293_R1/BAM"
for bam in ${BAM_DIR}/SRR*.bam; do
    if [[ -f "$bam" ]]; then
        echo "Indexing: $(basename $bam)"
        samtools index "$bam"
    fi
done
