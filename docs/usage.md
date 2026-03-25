# Usage

## Prerequisites

Before running the pipeline, make sure you have:
- Micromamba or Conda installed
- Access to a SLURM cluster (for alignment step)
- Candidate lncRNA FASTA files (one per sample)
- A reference transcriptome FASTA file
- eggNOG-mapper annotation file for your species

## Input File Format

### Candidate lncRNA FASTAs
```
>MSTRG.1.1
ATCGATCGATCGATCGATCG...
>MSTRG.2.1
GCTAGCTAGCTAGCTAGCTA...
```

Headers must follow StringTie MSTRG naming convention: `MSTRG.gene.transcript`

### Reference Transcriptome
Standard FASTA format. Used as the alignment target for minimap2.

### eggNOG-mapper Annotation
Tab-separated file from eggNOG-mapper v2 with GO terms in column 10 and
KEGG pathways in column 13.

## Configuration

All parameters are set in `config/config.yaml`. Copy the template first:
```bash
cp config/config.yaml.template config/config.yaml
```

Key parameters to set:

| Parameter | Description |
|-----------|-------------|
| `genome.reference` | Path to reference transcriptome FASTA |
| `genome.annotation` | Path to eggNOG-mapper annotation file |
| `samples.directory` | Directory containing candidate lncRNA FASTAs |
| `samples.names` | List of sample names |
| `slurm.account` | Your SLURM account |
| `slurm.email` | Email for job notifications |

## Running the Pipeline

### Manual step-by-step
```bash
# Step 1: Align (SLURM)
sbatch scripts/01_align.sh

# Step 2: Multi-sample comparison
bash scripts/02_multiinter.sh

# Step 3: GO enrichment
micromamba activate goanalysis
Rscript scripts/03_go_analysis.R

# Step 4: KEGG enrichment
Rscript scripts/04_kegg_analysis.R

# Step 5: Generate plots
Rscript scripts/05_plots.R
```

### Automated with Snakemake
```bash
micromamba activate snakemake

# Dry run first
snakemake --dry-run --cores 4

# Run on SLURM
snakemake --profile profiles/slurm --use-conda
```

## Test Run

A small synthetic test dataset is provided in `test_data/`. To test:
```bash
# Update config to point to test data
cp config/config.yaml.template config/config.yaml
# Edit config.yaml:
#   genome.reference: "test_data/reference/test_reference.fasta"
#   samples.directory: "test_data/samples"
#   samples.names: [TEST1, TEST2]

snakemake --dry-run --cores 2
```
## Adapting for Other Experimental Designs

The pipeline was developed for a 2×2 factorial design (2 stress conditions × 2 tissues).
For other designs such as developmental time series or embryogenesis comparisons,
the main thing to update is the awk filter logic in `scripts/02_multiinter.sh`.

### Example: zygotic vs somatic embryo (2 samples)
```bash
bedtools multiinter \
    -i SZE.bed SSE.bed \
    -names SZE SSE \
    > multiinter_output.bed

# Zygotic-specific
awk '$4==1 && $5=="SZE"' multiinter_output.bed > zygotic_specific.bed

# Somatic-specific
awk '$4==1 && $5=="SSE"' multiinter_output.bed > somatic_specific.bed

# Conserved between both
awk '$4==2' multiinter_output.bed > conserved.bed
```

All downstream GO and KEGG steps work without changes regardless of experimental design.

### Example: multi-species combined run (8 samples)
```bash
bedtools multiinter \
    -i results_pine/bed/*.bed results_spruce/bed/*.bed \
    -names PCN PCR PDN PDR SCN SCR SDN SDR \
    > multiinter_combined.bed

# Conserved across all 8 samples (both species, both conditions)
awk '$4==8' multiinter_combined.bed > conserved_all.bed

# Pine-specific (present in pine samples but not spruce)
awk '$6==1 && $7==1 && $8==1 && $9==1 && $10==0 && $11==0 && $12==0 && $13==0' \
    multiinter_combined.bed > pine_specific.bed

# Spruce-specific
awk '$6==0 && $7==0 && $8==0 && $9==0 && $10==1 && $11==1 && $12==1 && $13==1' \
    multiinter_combined.bed > spruce_specific.bed
```