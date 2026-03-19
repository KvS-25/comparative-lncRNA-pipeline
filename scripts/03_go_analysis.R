# ============================================================
# Script 03: GO enrichment analysis using topGO
#
# Usage: Rscript scripts/03_go_analysis.R
# Requires: config/config.yaml, results from 02_multiinter.sh
# Environment: goanalysis (R, topGO)
# ============================================================

library(topGO)
library(yaml)

# ---- Load config ----
config <- yaml.load_file("config/config.yaml")

ANNOTATION  <- config$genome$annotation
OUT_BED     <- config$output$bed
OUT_GO      <- config$output$go
PVAL_CUTOFF <- config$go$pvalue_cutoff
TOP_NODES   <- config$go$top_nodes
ALGORITHM   <- config$go$algorithm
STATISTIC   <- config$go$statistic

dir.create(OUT_GO, recursive=TRUE, showWarnings=FALSE)

# ---- Build gene to GO mapping ----
cat("[", format(Sys.time()), "] Building gene to GO mapping\n")

gene_to_go_file <- file.path(OUT_GO, "gene_to_GO.txt")

if (!file.exists(gene_to_go_file)) {
  cmd <- paste0(
    "zcat ", ANNOTATION, " | grep -v '^#' | ",
    "awk -F'\\t' '$10 != \"-\" {split($1, a, \".mRNA\"); print a[1] \"\\t\" $10}' > ",
    gene_to_go_file
  )
  system(cmd)
}

geneID2GO <- readMappings(gene_to_go_file)
universe  <- names(geneID2GO)
cat("Universe size:", length(universe), "genes\n")

# ---- Build MSTRG to reference gene mapping ----
mstrg_file <- file.path(OUT_GO, "mstrg_to_refgene.txt")
mstrg_map  <- read.table(mstrg_file, header=FALSE, sep="\t",
                          col.names=c("mstrg", "refgene"))

# ---- Run enrichment for each category ----
categories <- c("needle_specific", "root_specific",
                 "cold_specific", "drought_specific")
domains    <- c("BP", "MF", "CC")

for (category in categories) {
  cat("\n=== Processing:", category, "===\n")

  # Load and map gene IDs
  bed_file  <- file.path(dirname(OUT_GO), paste0(category, ".bed"))
  mstrg_ids <- sub("\\.[0-9]+$", "",
                    unique(read.table(bed_file, header=FALSE)$V1))
  ref_genes <- unique(mstrg_map$refgene[mstrg_map$mstrg %in% mstrg_ids])
  ref_genes <- ref_genes[ref_genes %in% universe]

  cat("Genes with GO annotation:", length(ref_genes), "\n")

  # Save reference genes
  write.table(ref_genes,
              file.path(OUT_GO, paste0(category, "_refgenes.txt")),
              row.names=FALSE, col.names=FALSE, quote=FALSE)

  gene_list        <- factor(as.integer(universe %in% ref_genes))
  names(gene_list) <- universe

  for (domain in domains) {
    cat("  Running", domain, "enrichment...\n")

    GOdata <- new("topGOdata",
                  ontology  = domain,
                  allGenes  = gene_list,
                  annot     = annFUN.gene2GO,
                  gene2GO   = geneID2GO)

    result    <- runTest(GOdata, algorithm=ALGORITHM, statistic=STATISTIC)
    res_table <- GenTable(GOdata,
                          weight01  = result,
                          orderBy   = "weight01",
                          topNodes  = TOP_NODES)

    res_table$weight01 <- as.numeric(res_table$weight01)
    sig_results        <- res_table[res_table$weight01 < PVAL_CUTOFF, ]

    cat("  Significant GO terms:", nrow(sig_results), "\n")

    out_file <- file.path(OUT_GO, paste0(category, "_", domain, "_GO_enrichment.txt"))
    write.table(sig_results, out_file, sep="\t", row.names=FALSE, quote=FALSE)
  }
}

cat("\n[", format(Sys.time()), "] GO analysis complete. Results saved to", OUT_GO, "\n")