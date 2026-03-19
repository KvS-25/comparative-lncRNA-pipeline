# ============================================================
# Script 04: KEGG pathway enrichment analysis
#
# Usage: Rscript scripts/04_kegg_analysis.R
# Requires: config/config.yaml, results from 03_go_analysis.R
# Environment: goanalysis (R)
# ============================================================

library(yaml)

# ---- Load config ----
config <- yaml.load_file("config/config.yaml")

ANNOTATION  <- config$genome$annotation
OUT_GO      <- config$output$go
OUT_KEGG    <- config$output$kegg
PVAL_CUTOFF <- config$kegg$pvalue_cutoff

dir.create(OUT_KEGG, recursive=TRUE, showWarnings=FALSE)

# ---- Build gene to KEGG mapping ----
cat("[", format(Sys.time()), "] Building gene to KEGG mapping\n")

kegg_file <- file.path(OUT_KEGG, "gene_to_KEGG.txt")

if (!file.exists(kegg_file)) {
  cmd <- paste0(
    "zcat ", ANNOTATION, " | grep -v '^#' | ",
    "awk -F'\\t' '$13 != \"-\" {split($1, a, \".mRNA\"); print a[1] \"\\t\" $13}' > ",
    kegg_file
  )
  system(cmd)
}

kegg_data <- read.table(kegg_file, header=FALSE, sep="\t",
                         col.names=c("gene", "pathways"), quote="")

# ---- Build lookup tables ----
gene_to_path <- list()
for (i in 1:nrow(kegg_data)) {
  gene  <- kegg_data$gene[i]
  paths <- unlist(strsplit(as.character(kegg_data$pathways[i]), ","))
  paths <- paths[grepl("^map", paths)]
  gene_to_path[[gene]] <- paths
}

path_to_genes <- list()
for (gene in names(gene_to_path)) {
  for (path in gene_to_path[[gene]]) {
    path_to_genes[[path]] <- c(path_to_genes[[path]], gene)
  }
}

universe <- names(gene_to_path)
cat("Universe size:", length(universe), "genes\n")
cat("Total pathways:", length(path_to_genes), "\n\n")

# ---- Download pathway names ----
pathway_names_file <- file.path(OUT_KEGG, "kegg_pathway_names.txt")

if (!file.exists(pathway_names_file)) {
  cat("Downloading KEGG pathway names\n")
  pathway_url  <- "https://rest.kegg.jp/list/pathway"
  pathway_names <- tryCatch(
    read.table(pathway_url, sep="\t", header=FALSE, quote=""),
    error=function(e) NULL
  )
  if (!is.null(pathway_names)) {
    colnames(pathway_names) <- c("Pathway_ID", "Name")
    pathway_names$Pathway_ID <- gsub("path:", "", pathway_names$Pathway_ID)
    write.table(pathway_names, pathway_names_file,
                sep="\t", row.names=FALSE, quote=FALSE)
  }
} else {
  pathway_names <- read.table(pathway_names_file,
                               header=TRUE, sep="\t", quote="")
}

# ---- Run enrichment for each category ----
categories <- c("needle_specific", "root_specific",
                 "cold_specific", "drought_specific")

for (category in categories) {
  cat("=== Processing:", category, "===\n")

  refgene_file    <- file.path(OUT_GO, paste0(category, "_refgenes.txt"))
  genes_of_interest <- readLines(refgene_file)
  genes_of_interest <- genes_of_interest[genes_of_interest %in% universe]

  cat("Genes with KEGG annotation:", length(genes_of_interest), "\n")

  results <- data.frame()

  for (path in names(path_to_genes)) {
    path_genes <- path_to_genes[[path]]
    a <- sum(genes_of_interest %in% path_genes)
    b <- length(genes_of_interest) - a
    c <- sum(universe %in% path_genes) - a
    d <- length(universe) - a - b - c
    if (a == 0) next
    mat  <- matrix(c(a, b, c, d), nrow=2)
    test <- fisher.test(mat, alternative="greater")
    results <- rbind(results, data.frame(
      Pathway_ID  = path,
      Annotated   = length(path_genes),
      Significant = a,
      Expected    = round(length(genes_of_interest) * length(path_genes) / length(universe), 2),
      p_value     = test$p.value
    ))
  }

  results     <- results[order(results$p_value), ]
  sig_results <- results[results$p_value < PVAL_CUTOFF, ]

  # Add pathway names
  sig_results$Pathway_Name <- pathway_names$Name[
    match(sig_results$Pathway_ID, pathway_names$Pathway_ID)]
  sig_results <- sig_results[, c("Pathway_ID", "Pathway_Name",
                                  "Annotated", "Significant",
                                  "Expected", "p_value")]

  cat("Significant pathways:", nrow(sig_results), "\n\n")

  out_file <- file.path(OUT_KEGG, paste0(category, "_KEGG_enrichment.txt"))
  write.table(sig_results, out_file, sep="\t", row.names=FALSE, quote=FALSE)
}

cat("[", format(Sys.time()), "] KEGG analysis complete. Results saved to", OUT_KEGG, "\n")