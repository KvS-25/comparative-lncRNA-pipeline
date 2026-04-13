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

## v1.5.1 — 2026-04-02

### Bug fixes and compatibility

- Fixed `envs/goanalysis.yaml`: removed `r-reshape2` (caused unsatisfiable dependency conflicts), added `r-yaml` and `r-upsetr`, pinned `r-base`, `libdeflate`, and `libtiff` to resolve `libtiff`/`libdeflate` solver conflicts
- Fixed `Snakefile`: added `conda: "envs/alignment.yaml"` directive to `align`, `extract_fasta`, and `multiinter` rules — minimap2 and bedtools were not found in SLURM job environments without this
- Fixed `multiinter` rule: converted from `run:` block to `shell:` block so the conda environment is correctly activated
- Fixed `scripts/03_go_analysis.R`: added automatic BED format detection to handle both chromosome-based (pine) and transcript-based (spruce) reference alignments; added `strip_isoform()` helper to correctly match gene IDs against eggNOG annotation
- Fixed `scripts/04_kegg_analysis.R`: added `strip_isoform()` helper to match refgene IDs against KEGG annotation
- Fixed `scripts/05_plots.R`: updated config parsing to use multi-species block structure (`config[[SPECIES]]`) instead of top-level keys; added `tryCatch` to handle empty GO/KEGG result files gracefully
- Fixed `profiles/slurm/config.yaml`: replaced template placeholders with actual SLURM account
## v1.5.2 — 2026-04-13

### GO enrichment improvements

- Added non-plant GO term filter to `03_go_analysis.R` — excludes animal/mammalian-specific terms (immune, neurological, hormonal) that bleed through from eggNOG cross-kingdom annotations
- Extended GO enrichment to all three ontology domains (BP, MF, CC) — previously only BP was run
- Updated Snakefile `go_analysis` rule to output BP, MF, and CC files per category
- Updated `rule all` and `rule plots` input blocks to include MF and CC outputs
- Updated `05_plots.R` GO bar plots to three-panel faceted style (BP/MF/CC), top 5 terms per domain, coloured by domain
- Fixed `05_plots.R` to respect species passed via command line argument rather than always reading `default_species` from config
