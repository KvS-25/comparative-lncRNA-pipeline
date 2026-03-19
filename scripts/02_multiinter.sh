#!/bin/bash
# ============================================================
# Script 02: Multi-sample interval comparison using bedtools
#            multiinter and filtering into region categories
#
# Usage: bash scripts/02_multiinter.sh
# Requires: config/config.yaml, results from 01_align.sh
# Environment: alignment (bedtools)
# ============================================================

set -euo pipefail

# ---- Load config ----
CONFIG="config/config.yaml"

OUT_BED=$(grep "paf:" $CONFIG | awk '{print $2}' | sed 's/paf/bed/')
OUT_BASE=$(grep "base:" $CONFIG | awk '{print $2}')
SAMPLES=($(grep -A 10 "names:" $CONFIG | grep "^    - " | awk '{print $2}'))

# ---- Activate environment ----
export PATH="/path/to/micromamba/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate alignment

echo "[$(date)] Running bedtools multiinter"

# ---- Multi-sample interval comparison ----
bedtools multiinter \
    -i ${OUT_BED}/${SAMPLES[0]}.bed \
       ${OUT_BED}/${SAMPLES[1]}.bed \
       ${OUT_BED}/${SAMPLES[2]}.bed \
       ${OUT_BED}/${SAMPLES[3]}.bed \
    -names ${SAMPLES[@]} \
    > ${OUT_BASE}/multiinter_output.bed

echo "[$(date)] Filtering regions by category"

# ---- Filter by category ----
# Conserved in all 4 samples
awk '$4 == 4' ${OUT_BASE}/multiinter_output.bed \
    > ${OUT_BASE}/conserved.bed

# Sample specific
for i in 0 1 2 3; do
    SAMPLE=${SAMPLES[$i]}
    awk -v s="$SAMPLE" '$4 == 1 && $5 == s' \
        ${OUT_BASE}/multiinter_output.bed \
        > ${OUT_BASE}/${SAMPLE}_specific.bed
    echo "  $SAMPLE specific: $(wc -l < ${OUT_BASE}/${SAMPLE}_specific.bed) regions"
done

# Tissue specific (columns 6-9 are binary presence/absence)
awk '$6==1 && $7==0 && $8==1 && $9==0' \
    ${OUT_BASE}/multiinter_output.bed > ${OUT_BASE}/needle_specific.bed

awk '$6==0 && $7==1 && $8==0 && $9==1' \
    ${OUT_BASE}/multiinter_output.bed > ${OUT_BASE}/root_specific.bed

# Condition specific
awk '$6==1 && $7==1 && $8==0 && $9==0' \
    ${OUT_BASE}/multiinter_output.bed > ${OUT_BASE}/cold_specific.bed

awk '$6==0 && $7==0 && $8==1 && $9==1' \
    ${OUT_BASE}/multiinter_output.bed > ${OUT_BASE}/drought_specific.bed

echo "[$(date)] Region counts:"
wc -l ${OUT_BASE}/needle_specific.bed \
       ${OUT_BASE}/root_specific.bed \
       ${OUT_BASE}/cold_specific.bed \
       ${OUT_BASE}/drought_specific.bed

echo "[$(date)] Done"