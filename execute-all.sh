#!/bin/bash

bash scripts/install-dependencies.sh
bash scripts/get-reads.sh
bash scripts/get-qualities.sh
bash scripts/filter-reads.sh
bash scripts/index-reference.sh
bash scripts/align_reads.sh
bash scripts/normalize-bams.sh
bash scripts/bam-analysis.sh
bash scripts/identify-variants.sh
bash scripts/isec-variants.sh
