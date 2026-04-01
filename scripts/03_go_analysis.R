# ============================================================
# Script 03: GO enrichment analysis using topGO
#
# Reads config/config.yaml to find output paths and annotation.
# Iterates over all category BED files that are non-empty.
# Produces BP enrichment results per category.
#
# Usage:
#   micromamba activate goanalysis
#   Rscript scripts/03_go_analysis.R
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(topGO)
})

# ---- Load config --------------------------------------------
config    <- yaml.load_file("config/config.yaml")
OUT_BASE  <- config$output$base
OUT_GO    <- config$output$go
ANNOT     <- config$genome$annotation
P_CUTOFF  <- config$go$pvalue_cutoff
TOP_NODES <- config$go$top_nodes
ALGORITHM <- config$go$algorithm
STATISTIC <- config$go$statistic

dir.create(OUT_GO, recursive = TRUE, showWarnings = FALSE)
dir.create(config$slurm$logs, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(config$slurm$logs, "go_analysis.log")
con <- file(log_file, open = "at")
msg <- function(...) {
  m <- paste0("[", format(Sys.time()), "] ", ...)
  message(m); cat(m, "\n", file = con)
}

# ---- All possible categories --------------------------------
all_cats <- c(
  "needle_specific", "root_specific",
  "cold_specific",   "drought_specific",
  "somatic_embryo_specific", "zygotic_embryo_specific"
)

cats <- Filter(function(cat) {
  f <- file.path(OUT_BASE, paste0(cat, ".bed"))
  file.exists(f) && file.info(f)$size > 0
}, all_cats)

if (length(cats) == 0) {
  msg("No non-empty category BED files found — nothing to do.")
  quit(status = 0)
}
msg("Categories: ", paste(cats, collapse = ", "))

# ---- Load annotation ----------------------------------------
if (is.null(ANNOT) || nchar(ANNOT) == 0 || !file.exists(ANNOT)) {
  stop("Annotation file not set or not found. Set genome.annotation in config.yaml.")
}
msg("Loading annotation: ", ANNOT)
annot <- read.table(gzfile(ANNOT), sep = "\t", header = FALSE,
                    comment.char = "#", quote = "", fill = TRUE)

gene2GO <- list()
for (i in seq_len(nrow(annot))) {
  gene <- as.character(annot[i, 1])
  gos  <- as.character(annot[i, 10])   # col 10 = GO terms in eggNOG output
  if (nchar(gos) > 0 && gos != "-") {
    gene2GO[[gene]] <- unlist(strsplit(gos, ","))
  }
}

write.table(
  data.frame(gene = names(gene2GO),
             GO   = sapply(gene2GO, paste, collapse = ",")),
  file = file.path(OUT_GO, "gene_to_GO.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
msg("gene_to_GO.txt: ", length(gene2GO), " genes")

# ---- Process each category ----------------------------------
for (cat in cats) {
  msg("Processing: ", cat)

  bed      <- read.table(file.path(OUT_BASE, paste0(cat, ".bed")),
                         header = FALSE, sep = "\t")
  raw_ids  <- as.character(bed[, 4])
  cat_genes <- unique(unlist(strsplit(raw_ids, ",")))

  # Separate MSTRG (novel) from reference IDs
  mstrg_ids <- cat_genes[grepl("^MSTRG", cat_genes)]
  ref_genes  <- cat_genes[!grepl("^MSTRG", cat_genes)]

  write.table(
    data.frame(MSTRG = mstrg_ids,
               ref_gene = NA_character_),   # filled if orthologue mapping available
    file = file.path(OUT_GO, "mstrg_to_refgene.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  writeLines(ref_genes, file.path(OUT_GO, paste0(cat, "_refgenes.txt")))

  all_genes <- names(gene2GO)
  gene_list <- factor(as.integer(all_genes %in% ref_genes))
  names(gene_list) <- all_genes

  n_annotated <- sum(gene_list == 1)
  if (n_annotated == 0) {
    msg("  No annotated reference genes for ", cat, " — skipping")
    file.create(file.path(OUT_GO, paste0(cat, "_BP_GO_enrichment.txt")))
    next
  }
  msg("  Annotated genes: ", n_annotated)

  go_data <- new("topGOdata",
    ontology = "BP",
    allGenes = gene_list,
    gene2GO  = gene2GO,
    annot    = annFUN.gene2GO,
    nodeSize = 10
  )

  result    <- runTest(go_data, algorithm = ALGORITHM, statistic = STATISTIC)
  go_table  <- GenTable(go_data, pvalue = result,
                        topNodes = TOP_NODES, numChar = 100)
  go_table$pvalue <- as.numeric(go_table$pvalue)
  sig_table <- go_table[!is.na(go_table$pvalue) & go_table$pvalue < P_CUTOFF, ]

  out_file <- file.path(OUT_GO, paste0(cat, "_BP_GO_enrichment.txt"))
  write.table(sig_table, file = out_file,
              sep = "\t", row.names = FALSE, quote = FALSE)
  msg("  Significant GO terms: ", nrow(sig_table), " -> ", out_file)
}

close(con)
message("GO analysis complete.")