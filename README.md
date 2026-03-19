# comparative-lncRNA-pipeline

A comparative genomics pipeline for identifying and characterising tissue-specific and stress-responsive long non-coding RNAs (lncRNAs) in conifers using minimap2, bedtools, and GO/KEGG enrichment analysis.

Developed as part of an MSc thesis at Umeå University, initially applied to *Pinus sylvestris* under cold and drought stress conditions across needle and root tissues.

---

## Overview

This pipeline takes candidate lncRNA transcript sequences (e.g. from the [Plant LncRNA Pipeline v2](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2)) and performs:

1. Alignment to a reference transcriptome using minimap2
2. PAF to BED format conversion
3. Multi-sample interval comparison to identify conserved and condition/tissue-specific regions
4. Gene Ontology (GO) enrichment analysis using topGO
5. KEGG pathway enrichment analysis
6. Publication-ready visualisation plots

---

## Pipeline Overview
```
Candidate lncRNA FASTAs
        │
        ▼
01_align.sh          ── minimap2 splice alignment ──► PAF files
        │
        ▼
        │            ── PAF to BED conversion    ──► BED files
        │
        ▼
02_multiinter.sh     ── bedtools multiinter      ──► Region categories
        │                                            (conserved, tissue-specific,
        │                                             condition-specific)
        ▼
03_go_analysis.R     ── topGO enrichment         ──► GO results (BP/MF/CC)
        │
        ▼
04_kegg_analysis.R   ── Fisher's exact test      ──► KEGG pathway results
        │
        ▼
05_plots.R           ── ggplot2 / UpSetR         ──► Publication figures
```

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

To run the entire pipeline automatically:
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

If you use this pipeline please cite:

- **minimap2**: Li, H. (2018). Bioinformatics, 34(18), 3094–3100.
- **bedtools**: Quinlan, A.R. and Hall, I.M. (2010). Bioinformatics, 26(6), 841–842.
- **topGO**: Alexa, A., Rahnenführer, J. and Lengauer, T. (2006). Bioinformatics, 22(13), 1600–1607.
- **StringTie**: Pertea, M. et al. (2015). Nature Biotechnology, 33(3), 290–295.
- **eggNOG-mapper**: Cantalapiedra, C.P. et al. (2021). Molecular Biology and Evolution, 38(12), 5825–5829.
- **ggplot2**: Wickham, H. (2016). Springer-Verlag, New York.

---

## Author

Karanveer Singh — MSc Molecular Biology, Umeå University  
GitHub: [KvS-25](https://github.com/KvS-25)

---

## License

MIT License — free to use and modify with attribution.
