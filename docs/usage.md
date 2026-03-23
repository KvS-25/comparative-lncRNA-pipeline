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
