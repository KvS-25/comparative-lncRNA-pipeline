#!/bin/bash
#SBATCH -A u2015037
#SBATCH -p main -n 20
#SBATCH -t 2-00:00:00
#SBATCH --job-name=minimap2_spruce
#SBATCH --array=1-6%1
#SBATCH --output=ThesisData/logs/spruce/%x_%A_%a.out
#SBATCH --error=ThesisData/errors/spruce/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kasi0057@student.umu.se

set -euo pipefail

export PATH="/mnt/picea/home/ksingh/.local/bin:$PATH"
eval "$(/mnt/picea/home/ksingh/.local/bin/micromamba shell hook --shell bash)"
micromamba activate alignment

SAMPLES=(SCN SCR SDN SDR SSE SZE)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

GENOME="/mnt/picea/home/ksingh/ThesisData/spruce_ref/merge_transcripts_spruce.fasta"
FASTA="/mnt/picea/home/ksingh/ThesisData/spruce_samples/candidate_transcript_${SAMPLE}.fasta"
OUT="ThesisData/output/spruce"

mkdir -p "$OUT"/{paf,bed}
mkdir -p ThesisData/logs/spruce ThesisData/errors/spruce

echo "[$(date)] Starting alignment for $SAMPLE"

minimap2 -t 20 -I 4G -x splice \
    "$GENOME" "$FASTA" \
    > "$OUT/paf/${SAMPLE}.paf"

if [ ! -s "$OUT/paf/${SAMPLE}.paf" ]; then
    echo "ERROR: PAF file is empty for $SAMPLE" >&2
    exit 1
fi

echo "[$(date)] Converting PAF to BED for $SAMPLE"

awk 'BEGIN{OFS="\t"} {print $6, $8, $9, $1, $12}' \
    "$OUT/paf/${SAMPLE}.paf" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -c 4,5 -o distinct,mean \
    > "$OUT/bed/${SAMPLE}.bed"

echo "[$(date)] Finished $SAMPLE"
