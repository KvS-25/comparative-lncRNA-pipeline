# Usage

## Prerequisites

Before running the pipeline, make sure you have:
- Micromamba or Conda installed
- Access to a SLURM cluster (for alignment step)
- Candidate lncRNA FASTA files (one per sample)
- A reference transcriptome FASTA file
- eggNOG-mapper annotation file for your species

## Installation

**1. Clone the repository:**
```bash
git clone https://github.com/KvS-25/comparative-lncRNA-pipeline.git
cd comparative-lncRNA-pipeline
```

**2. Install Snakemake and the SLURM executor plugin:**
```bash
micromamba create -n snakemake -c conda-forge -c bioconda snakemake
micromamba activate snakemake
pip install snakemake-executor-plugin-slurm
```

**3. Ensure conda ≥ 24.7.1 is available in the snakemake environment:**
```bash
micromamba install "conda>=24.7.1" -c conda-forge
```

**4. Create conda environments:**
```bash
micromamba env create -f envs/alignment.yaml
micromamba env create -f envs/goanalysis.yaml
```

**5. Set up config:**
```bash
cp config/config.yaml.template config/config.yaml
nano config/config.yaml  # fill in your paths
```

## Configuration

All parameters are set in `config/config.yaml`. Copy the template first:
```bash
cp config/config.yaml.template config/config.yaml
```

Key parameters to set for each species block:

| Parameter | Description |
|-----------|-------------|
| `genome.reference` | Path to reference transcriptome FASTA |
| `genome.annotation` | Path to eggNOG-mapper annotation file |
| `samples.directory` | Directory containing candidate lncRNA FASTAs |
| `samples.names` | List of sample names |
| `slurm.account` | Your SLURM account |
| `slurm.partition` | Your SLURM partition |

The default species to run is set via `default_species` at the top of the config.
You can override this at runtime with `--config species=spruce`.

## Input File Format

### Candidate lncRNA FASTAs
```
>MSTRG.1.1
ATCGATCGATCGATCGATCG...
>MSTRG.2.1
GCTAGCTAGCTAGCTAGCTA...
```

Headers must follow StringTie MSTRG naming convention (`MSTRG.gene.transcript`)
or reference gene ID format (`PA_chrXX_GXXXXXX.mRNA.X`). Both are supported.

### Reference Transcriptome
Standard FASTA format. Used as the alignment target for minimap2.

### eggNOG-mapper Annotation
Tab-separated file from eggNOG-mapper v2 with GO terms in column 10 and
KEGG pathways in column 13.

## Running the Pipeline

### On a SLURM cluster (recommended)
```bash
micromamba activate snakemake

# Dry run first — always do this
snakemake --config species=pine --dry-run --cores 4

# Run for pine
snakemake --config species=pine --profile profiles/slurm --use-conda

# Run for spruce
snakemake --config species=spruce --profile profiles/slurm --use-conda
```

### Locally (no SLURM)
```bash
snakemake --config species=pine --cores 4 --use-conda
```

## Monitoring Jobs

Check running SLURM jobs:
```bash
squeue -u $USER
```

Watch live:
```bash
watch -n 10 squeue -u $USER
```

Follow a specific log:
```bash
tail -f ThesisData/logs/align_SCN.log
```

## Troubleshooting

### conda solver errors
If you see `LibMambaUnsatisfiableError`, check your channel priority:
```bash
conda config --show channel_priority
```

If set to `strict`, change it to `flexible`:
```bash
conda config --set channel_priority flexible
```

Then clear the Snakemake conda cache and retry:
```bash
rm -rf .snakemake/conda/
snakemake --config species=pine --profile profiles/slurm --use-conda
```

### Snakemake lock errors
If Snakemake exits unexpectedly and leaves a lock:
```bash
snakemake --config species=pine --profile profiles/slurm --use-conda --unlock
```

### Rerunning incomplete jobs
```bash
snakemake --config species=pine --profile profiles/slurm --use-conda --rerun-incomplete
```

## Test Run

A small synthetic test dataset is provided in `test_data/`. To test:
```bash
cp config/config.yaml.template config/config.yaml
# Edit config.yaml:
#   genome.reference: "test_data/reference/test_reference.fasta"
#   samples.directory: "test_data/samples"
#   samples.names: [TEST1, TEST2]

snakemake --dry-run --cores 2
```

## Adapting for Other Experimental Designs

The pipeline was developed for a 2×2 factorial design (2 stress conditions × 2 tissues).
The `scripts/02_multiinter.sh` awk filters are automatically generated based on
sample name conventions in the config — samples ending in `N` are treated as needle,
`R` as root, second character `C` as cold, `D` as drought, and `SSE`/`SZE` as embryo.

For other designs, update the filter logic in `scripts/02_multiinter.sh` to match
your sample columns. See below for examples.

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

### Example: multi-species combined run (8 samples)
```bash
bedtools multiinter \
    -i results_pine/bed/*.bed results_spruce/bed/*.bed \
    -names PCN PCR PDN PDR SCN SCR SDN SDR \
    > multiinter_combined.bed

# Conserved across all 8 samples
awk '$4==8' multiinter_combined.bed > conserved_all.bed

# Pine-specific
awk '$6==1 && $7==1 && $8==1 && $9==1 && $10==0 && $11==0 && $12==0 && $13==0' \
    multiinter_combined.bed > pine_specific.bed

# Spruce-specific
awk '$6==0 && $7==0 && $8==0 && $9==0 && $10==1 && $11==1 && $12==1 && $13==1' \
    multiinter_combined.bed > spruce_specific.bed
```

All downstream GO and KEGG steps work without changes regardless of experimental design.