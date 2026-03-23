# comparative-lncRNA-pipeline

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Version](https://img.shields.io/badge/version-1.0.0-blue)
![Snakemake](https://img.shields.io/badge/workflow-Snakemake-1a9c3e?logo=snakemake)
![Language: R](https://img.shields.io/badge/language-R-276DC3?logo=r)
![Language: Bash](https://img.shields.io/badge/language-Bash-4EAA25?logo=gnubash)
![Platform: HPC](https://img.shields.io/badge/platform-HPC%20SLURM-orange)
![Species: Conifer](https://img.shields.io/badge/species-conifer-2e8b57)

A comparative genomics pipeline for identifying and characterising tissue-specific and stress-responsive long non-coding RNAs (lncRNAs) in conifers using minimap2, bedtools, and GO/KEGG enrichment analysis.

Developed as part of an MSc thesis at Umeå University, initially applied to *Pinus sylvestris* under cold and drought stress conditions across needle and root tissues.

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

## Citation

Please see [CITATIONS.md](CITATIONS.md) for full citation information.

---

## License

MIT License — free to use and modify with attribution.