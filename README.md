# comparative-lncRNA-pipeline

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Version](https://img.shields.io/badge/version-1.5.0-blue)
![Snakemake](https://img.shields.io/badge/workflow-Snakemake-1a9c3e?logo=snakemake)
![Language: R](https://img.shields.io/badge/language-R-276DC3?logo=r)
![Language: Bash](https://img.shields.io/badge/language-Bash-4EAA25?logo=gnubash)
![Platform: HPC](https://img.shields.io/badge/platform-HPC%20SLURM-orange)
![Species: Conifer](https://img.shields.io/badge/species-pine%20%7C%20spruce-2e8b57)

A comparative genomics pipeline for identifying and characterising tissue-specific and stress-responsive long non-coding RNAs (lncRNAs) in conifers using minimap2, bedtools, and GO/KEGG enrichment analysis.

Developed as part of an MSc thesis at Umeå University, initially applied to *Pinus sylvestris* under cold and drought stress conditions across needle and root tissues. Being extended to *Picea abies* (Norway spruce) including stress response and embryogenesis samples.

**Contact Information:**
- Email: kvs.ms.2512@gmail.com
- GitHub: [KvS-25](https://github.com/KvS-25)

---

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Automated Workflow (Snakemake)](#automated-workflow-snakemake)
- [Output Structure](#output-structure)
- [Multi-species Analysis](#multi-species-analysis)
- [Citation](#citation)
- [License](#license)

---

## Pipeline Overview

![Pipeline Flowchart](images/pipeline_flowchart.svg)

---

## Requirements

- [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) or Conda
- SLURM workload manager (for alignment step)
- Internet access (for KEGG pathway name download)

---

## Installation

**1. Clone the repository:**
```bash
git clone https://github.com/KvS-25/comparative-lncRNA-pipeline.git
cd comparative-lncRNA-pipeline
```

**2. Create environments:**
```bash
micromamba env create -f envs/alignment.yaml
micromamba env create -f envs/goanalysis.yaml

# Install remaining R packages
micromamba activate goanalysis
Rscript -e "install.packages(c('yaml', 'UpSetR'), repos='https://cloud.r-project.org')"
```

**3. Set up config:**
```bash
cp config/config.yaml.template config/config.yaml
nano config/config.yaml  # fill in your paths
```

---

## Usage

Run scripts in order:
```bash
# Step 1: Align (SLURM)
sbatch scripts/01_align.sh

# Step 2: Multi-sample comparison (login node)
bash scripts/02_multiinter.sh

# Step 3: GO enrichment (login node)
micromamba activate goanalysis
Rscript scripts/03_go_analysis.R

# Step 4: KEGG enrichment (login node)
Rscript scripts/04_kegg_analysis.R

# Step 5: Generate plots (login node)
Rscript scripts/05_plots.R
```

---

## Automated Workflow (Snakemake)
```bash
# Install Snakemake
micromamba create -n snakemake -c conda-forge -c bioconda snakemake
micromamba activate snakemake

# Dry run to check workflow
snakemake --dry-run --cores 4

# Run locally
snakemake --cores 4 --use-conda

# Run on SLURM cluster
snakemake --profile profiles/slurm --use-conda
```

---

## Output Structure
```
results/
├── paf/                    # minimap2 alignment output
├── bed/                    # converted BED files
├── fasta/                  # extracted FASTA sequences
├── GO_analysis/            # GO enrichment results
│   ├── gene_to_GO.txt
│   ├── mstrg_to_refgene.txt
│   ├── *_refgenes.txt
│   └── *_GO_enrichment.txt
├── KEGG/                   # KEGG pathway results
│   ├── gene_to_KEGG.txt
│   ├── kegg_pathway_names.txt
│   └── *_KEGG_enrichment.txt
├── plots/                  # all figures
│   ├── upset_plot.png
│   ├── region_counts_bar.png
│   ├── GO_bar_*.png
│   └── KEGG_bubble_plot.png
├── multiinter_output.bed
├── conserved.bed
├── needle_specific.bed
├── root_specific.bed
├── cold_specific.bed
└── drought_specific.bed
```

---

## Multi-species Analysis

The pipeline is designed to be species-agnostic. It has been applied to *Pinus sylvestris* and is being extended to *Picea abies* (Norway spruce) for both stress response and embryogenesis comparisons.

### Running for a new species

1. Obtain candidate lncRNA FASTAs using [Plant LncRNA Pipeline v2](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2)
2. Obtain a reference transcriptome and eggNOG-mapper annotation for your species
3. Copy and update the config:
```bash
cp config/config.yaml.template config/config.yaml
# Update species, genome paths, sample names and output directory
```

4. Run as normal — the pipeline requires no other changes

### Suggested sample naming convention

| Code | Meaning |
|------|---------|
| PCN | Pine Cold Needle |
| PCR | Pine Cold Root |
| PDN | Pine Drought Needle |
| PDR | Pine Drought Root |
| SCN | Spruce Cold Needle |
| SCR | Spruce Cold Root |
| SDN | Spruce Drought Needle |
| SDR | Spruce Drought Root |
| SZE | Spruce Zygotic Embryo |
| SSE | Spruce Somatic Embryo |

> **Note for embryogenesis or other experimental designs**: The `scripts/02_multiinter.sh`
> awk filters are currently designed for a 2×2 stress/tissue design. For other designs
> update the filter logic to match your sample columns. See `docs/usage.md` for details.

### Cross-species comparison

To compare pine and spruce results:
- Run the pipeline separately for each species with separate output directories
- GO and KEGG enrichment results can be compared directly between species
- For a combined multi-sample analysis, use all 8 BED files in `scripts/02_multiinter.sh`
```bash
bedtools multiinter \
    -i results_pine/bed/*.bed results_spruce/bed/*.bed \
    -names PCN PCR PDN PDR SCN SCR SDN SDR \
    > results_combined/multiinter_output.bed
```

---

## Citation

Please see [CITATIONS.md](CITATIONS.md) for full citation information.

---

## License

MIT License — free to use and modify with attribution.