#!/bin/bash
# ============================================================
# Script 01: minimap2 alignment — SLURM array job
#
# Reads all settings from config/config.yaml.
# One array task per sample — submits automatically via SLURM.
#
# Usage:
#   sbatch scripts/01_align.sh
# ============================================================

#SBATCH --job-name=lncrna_align
#SBATCH --account=u2015037
#SBATCH --partition=main
#SBATCH --ntasks=20
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

CONFIG="config/config.yaml"

# ---- Parse config -------------------------------------------
SAMPLES_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['samples']['directory'])")
GENOME=$(python3      -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['genome']['reference'])")
OUT_PAF=$(python3     -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output']['paf'])")
OUT_BED=$(python3     -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output']['bed'])")
LOGS=$(python3        -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['slurm']['logs'])")
ERRORS=$(python3      -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['slurm']['errors'])")
THREADS=$(python3     -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['slurm']['threads'])")
INDEX_SIZE=$(python3  -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['minimap2']['index_size'])")
PRESET=$(python3      -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['minimap2']['preset'])")

mapfile -t SAMPLES < <(python3 -c "
import yaml
c = yaml.safe_load(open('$CONFIG'))
for s in c['samples']['names']:
    print(s)
")

TOTAL=${#SAMPLES[@]}

# If not running as array task, re-submit as one
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Submitting as SLURM array with $TOTAL tasks..."
    sbatch --array=0-$((TOTAL - 1))%1 "$0"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# ---- Set up logging -----------------------------------------
mkdir -p "$LOGS" "$ERRORS" "$OUT_PAF" "$OUT_BED"
exec > >(tee -a "$LOGS/align_${SAMPLE}.log") 2> >(tee -a "$ERRORS/align_${SAMPLE}.err" >&2)

echo "[$(date)] Starting alignment: $SAMPLE"

INPUT_FASTA="${SAMPLES_DIR}/candidate_transcript_${SAMPLE}.fasta"

[ -f "$INPUT_FASTA" ] || { echo "ERROR: Input FASTA not found: $INPUT_FASTA" >&2; exit 1; }
[ -f "$GENOME" ]      || { echo "ERROR: Reference not found: $GENOME" >&2; exit 1; }

# ---- Activate environment -----------------------------------
eval "$(micromamba shell hook --shell bash)"
micromamba activate alignment

# ---- minimap2 -----------------------------------------------
echo "[$(date)] Running minimap2"
minimap2 \
    -x "$PRESET" \
    -I "$INDEX_SIZE" \
    -t "$THREADS" \
    "$GENOME" \
    "$INPUT_FASTA" \
    > "${OUT_PAF}/${SAMPLE}.paf"

[ -s "${OUT_PAF}/${SAMPLE}.paf" ] || { echo "ERROR: Empty PAF for $SAMPLE" >&2; exit 1; }
echo "[$(date)] PAF: $(wc -l < "${OUT_PAF}/${SAMPLE}.paf") alignments"

# ---- PAF to BED ---------------------------------------------
echo "[$(date)] Converting to BED"
awk 'BEGIN{OFS="\t"} $5 != "*" {
    score = int($12 * 1000 / 60);
    if (score > 1000) score = 1000;
    print $6, $8, $9, $1, score, $5
}' "${OUT_PAF}/${SAMPLE}.paf" \
    | sort -k1,1 -k2,2n \
    > "${OUT_BED}/${SAMPLE}.bed"

echo "[$(date)] BED: $(wc -l < "${OUT_BED}/${SAMPLE}.bed") regions"
echo "[$(date)] Done: $SAMPLE"