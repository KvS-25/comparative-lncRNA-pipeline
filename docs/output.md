# Output

All results are saved in the `results/` directory as specified in `config/config.yaml`.

## Directory Structure
```
results/
├── paf/                          # minimap2 PAF alignment files
│   └── {SAMPLE}.paf
├── bed/                          # BED format alignment files
│   └── {SAMPLE}.bed
├── fasta/                        # Extracted FASTA sequences
│   ├── conserved.fasta
│   ├── {CATEGORY}_specific.fasta
│   └── tissue_condition/
├── GO_analysis/                  # GO enrichment results
│   ├── gene_to_GO.txt            # Gene to GO term mapping
│   ├── mstrg_to_refgene.txt      # MSTRG to reference gene mapping
│   ├── {CATEGORY}_refgenes.txt   # Reference genes per category
│   ├── {CATEGORY}_GO.txt         # GO annotated genes
│   └── {CATEGORY}_{DOMAIN}_GO_enrichment.txt
├── KEGG/                         # KEGG pathway results
│   ├── gene_to_KEGG.txt          # Gene to KEGG mapping
│   ├── kegg_pathway_names.txt    # KEGG pathway name lookup
│   └── {CATEGORY}_KEGG_enrichment.txt
├── plots/                        # Publication figures
│   ├── upset_plot.png
│   ├── region_counts_bar.png
│   ├── GO_bar_{CATEGORY}.png
│   └── KEGG_bubble_plot.png
├── multiinter_output.bed         # Full multi-sample comparison
├── conserved.bed                 # Regions in all samples
├── needle_specific.bed           # Needle tissue specific
├── root_specific.bed             # Root tissue specific
├── cold_specific.bed             # Cold stress specific
├── drought_specific.bed          # Drought stress specific
├── somatic_embryo_specific.bed   # Somatic embryo specific (SSE only, if present)
└── zygotic_embryo_specific.bed   # Zygotic embryo specific (SZE only, if present)
```

## File Descriptions

### PAF files
Pairwise Alignment Format output from minimap2. Contains alignment
coordinates, mapping quality, and CIGAR strings.

### BED files
Converted and merged alignment intervals. Columns:
1. Sequence name
2. Start position
3. End position
4. Query transcript IDs (comma-separated if merged)
5. Mean mapping quality

### GO enrichment files
Tab-separated results from topGO. Columns:
- GO.ID — Gene Ontology term ID
- Term — GO term description
- Annotated — total genes annotated to this term
- Significant — significant genes in your set
- Expected — expected by chance
- weight01 — p-value from Fisher's exact test

### KEGG enrichment files
Tab-separated results. Columns:
- Pathway_ID — KEGG pathway ID (map##### format)
- Pathway_Name — pathway description
- Annotated — total genes in pathway
- Significant — significant genes in your set
- Expected — expected by chance
- p_value — Fisher's exact test p-value
