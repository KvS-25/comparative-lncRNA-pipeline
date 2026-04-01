# Changelog

## v1.0.0 — 2026-03-19

### Initial release

- minimap2 alignment with splice preset
- PAF to BED conversion
- bedtools multiinter for multi-sample comparison
- Region categorisation (conserved, tissue-specific, condition-specific)
- GO enrichment analysis using topGO (BP, MF, CC)
- KEGG pathway enrichment using Fisher's exact test
- Visualisation with ggplot2 (bubble plots, bar plots, heatmap, UpSet plot)
- Snakemake workflow with SLURM integration
- Synthetic test dataset
- Config template for easy setup

## v1.5.0 - 2026-04-01

### Addition of spruce

-Addition of spruce dataset rules and outputs to config and template
-Updated Snakefile to allow choosing of spieces to run
-Updated scripts with revisied directory structure and logic