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

## v1.5.0 — 2026-04-01

### Multi-species support

- Config restructured into per-species blocks (`pine`, `spruce`) with `default_species` key
- Runtime species selection: `snakemake --config species=spruce`
- Spruce sample set added: SCN, SCR, SDN, SDR (stress) and SSE, SZE (embryogenesis)
- `scripts/02_multiinter.sh` now auto-detects tissue/condition/embryo categories from sample names — no manual edits needed for new designs
- `scripts/03_go_analysis.R` split into `build_maps` (runs once) and `enrich` (per-category) modes
- `scripts/04_kegg_analysis.R` split into `build_maps` and `enrich` modes; KEGG API fetched once, saved as RDS
- `scripts/05_plots.R` extended with embryo-specific colour palette entries
- `Snakefile` rewritten: `build_go_maps` and `build_kegg_maps` rules run once; `go_analysis` and `kegg_analysis` wildcard rules run per category; embryo categories supported