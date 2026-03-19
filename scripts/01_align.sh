#!/bin/bash
#SBATCH --job-name=lncRNA_align
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=errors/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL

# ============================================================
# Script 01: Align candidate lncRNA transcripts to reference
#            transcriptome using minimap2 with splice preset
#
# Usage: sbatch scripts/01_align.sh
# Requires: config/config.yaml to be filled in
# Environment: alignment (minimap2, bedtools)
# ============================================================

set -euo pipefail

# ---- Load config ----
CONFIG="config/config.yaml"

GENOME=$(grep "reference:" $CONFIG | awk '{print $2}')
SAMPLES_DIR=$(grep "directory:" $CONFIG | awk '{print $2}')
OUT_PAF=$(grep "paf:" $CONFIG | awk '{print $2}')
OUT_BED=$(grep "bed:" $CONFIG | awk '{print $2}')
THREADS=$(grep "threads:" $CONFIG | awk '{print $2}')
INDEX_SIZE=$(grep "index_size:" $CONFIG | awk '{print $2}' | tr -d '"')
PRESET=$(grep "preset:" $CONFIG | awk '{print $2}' | tr -d '"')
ACCOUNT=$(grep "account:" $CONFIG | awk '{print $2}' | tr -d '"')
EMAIL=$(grep "email:" $CONFIG | awk '{print $2}' | tr -d '"')

# ---- SLURM settings from config ----
#SBATCH -A $ACCOUNT
#SBATCH -p main
#SBATCH -n $THREADS
#SBATCH -t 2-00:00:00
#SBATCH --array=1-4%1
#SBATCH --mail-user=$EMAIL

# ---- Activate environment ----
export PATH="/path/to/micromamba/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate alignment

# ---- Get sample name ----
SAMPLES=($(grep -A 10 "names:" $CONFIG | grep "^    - " | awk '{print $2}'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

mkdir -p "$OUT_PAF" "$OUT_BED" logs errors

echo "[$(date)] Starting alignment for $SAMPLE"

# ---- Align with minimap2 ----
minimap2 \
    -t "$THREADS" \
    -I "$INDEX_SIZE" \
    -x "$PRESET" \
    "$GENOME" \
    "${SAMPLES_DIR}/candidate_transcript_${SAMPLE}.fasta" \
    > "${OUT_PAF}/${SAMPLE}.paf"

if [ ! -s "${OUT_PAF}/${SAMPLE}.paf" ]; then
    echo "ERROR: PAF file is empty for $SAMPLE" >&2
    exit 1
fi

echo "[$(date)] Converting PAF to BED for $SAMPLE"

# ---- Convert PAF to BED ----
awk 'BEGIN{OFS="\t"} {print $6, $8, $9, $1, $12}' \
    "${OUT_PAF}/${SAMPLE}.paf" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -c 4,5 -o distinct,mean \
    > "${OUT_BED}/${SAMPLE}.bed"

echo "[$(date)] Finished $SAMPLE"