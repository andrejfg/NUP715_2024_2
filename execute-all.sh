#!/bin/bash



bash scripts/get-reads.sh
bash scripts/get-qualities.sh
bash scripts/filter-reads.sh data/reads/SRR14457781.fastq
bash scripts/get-qualities.sh data/filtered_reads/SRR14457781/good/SRR14457781_filtered_Q30_30.fastq
bash scripts/index-reference.sh
bash scripts/align_reads.sh data/reads/SRR1570792.fastq data/filtered_reads/SRR14457781/good/SRR14457781_filtered_Q30_30.fastq 
bash scripts/normalize-bams.sh
bash scripts/bam-analysis.sh
bash scripts/identify-variants.sh
bash scripts/isec-variants.sh
