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

The config has a separate block for each species (`pine`, `spruce`). You specify
which species to run at the command line — you do not need to edit the config
each time you switch species.

Key parameters to set in each species block:

| Parameter | Description |
|-----------|-------------|
| `genome.reference` | Path to reference transcriptome FASTA |
| `genome.annotation` | Path to eggNOG-mapper annotation file |
| `samples.directory` | Directory containing candidate lncRNA FASTAs |
| `samples.names` | List of sample names |
| `slurm.account` | Your SLURM account |

## Running the Pipeline

### Recommended: Snakemake (automated)

Snakemake handles all steps in the correct order, including running GO and KEGG
enrichment once per category automatically.
```bash
micromamba activate snakemake

# Dry run first — always do this to check everything looks right
snakemake --config species=pine --dry-run --cores 4

# Run on SLURM (recommended for full datasets)
snakemake --config species=pine --profile profiles/slurm --use-conda

# For spruce
snakemake --config species=spruce --profile profiles/slurm --use-conda

# Run locally (small datasets / test data only)
snakemake --config species=pine --cores 4 --use-conda
```

Snakemake will automatically:
1. Align each sample with minimap2 (submitted as SLURM array jobs)
2. Run bedtools multiinter and classify regions
3. Build GO and KEGG maps once from the annotation file
4. Run GO and KEGG enrichment for each region category
5. Generate all plots

### Manual step-by-step

Only needed for debugging individual steps. In normal use, run Snakemake instead.
```bash
# Step 1: Align (submits SLURM array job, one task per sample)
sbatch scripts/01_align.sh

# Step 2: Multi-sample region comparison (run on login node after Step 1 finishes)
bash scripts/02_multiinter.sh

# Step 3: GO enrichment (run on login node)
# First, build the gene-to-GO map once from the annotation file:
micromamba activate goanalysis
Rscript scripts/03_go_analysis.R build_maps \
    /path/to/annotation.tsv.gz \
    /path/to/output/GO_analysis/gene_to_GO.txt \
    /path/to/output/GO_analysis/mstrg_to_refgene.txt

# Then run enrichment for each region category:
Rscript scripts/03_go_analysis.R enrich \
    /path/to/output/needle_specific.bed \
    /path/to/output/GO_analysis/gene_to_GO.txt \
    /path/to/output/GO_analysis/mstrg_to_refgene.txt \
    /path/to/output/GO_analysis/needle_specific_refgenes.txt \
    /path/to/output/GO_analysis/needle_specific_BP_GO_enrichment.txt
# Repeat the enrich command for each category (root_specific, cold_specific, etc.)

# Step 4: KEGG enrichment (run on login node)
# First, build the KEGG map once:
Rscript scripts/04_kegg_analysis.R build_maps \
    /path/to/annotation.tsv.gz \
    /path/to/output/KEGG/gene_to_KEGG.txt \
    /path/to/output/KEGG/kegg_pathway_names.txt

# Then run enrichment for each category:
Rscript scripts/04_kegg_analysis.R enrich \
    /path/to/output/GO_analysis/needle_specific_refgenes.txt \
    /path/to/output/KEGG/gene_to_KEGG.txt \
    /path/to/output/KEGG/kegg_pathway_names.txt \
    /path/to/output/KEGG/needle_specific_KEGG_enrichment.txt
# Repeat the enrich command for each category

# Step 5: Generate plots (run on login node after Steps 3 and 4 finish)
Rscript scripts/05_plots.R
```

## Test Run

A small synthetic test dataset is provided in `test_data/`. To test:
```bash
cp config/config.yaml.template config/config.yaml
# Edit config.yaml:
#   Under the pine block:
#     genome.reference: "test_data/reference/test_reference.fasta"
#     samples.directory: "test_data/samples"
#     samples.names: [TEST1, TEST2]

snakemake --config species=pine --dry-run --cores 2
```

## Adapting for Other Experimental Designs

The pipeline was developed for a 2×2 factorial design (2 stress conditions × 2 tissues).
`scripts/02_multiinter.sh` auto-detects which categories to produce based on your sample
names — needle vs. root from names ending in N or R, cold vs. drought from the second
character being C or D, and embryo from SSE/SZE being present. For most new designs,
updating `samples.names` in the config is sufficient.

For completely custom designs not covered by the naming convention, update the awk
filter logic in `scripts/02_multiinter.sh`. All downstream GO and KEGG steps work
without changes regardless of experimental design.

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