# ============================================================
# Script 03: GO enrichment analysis (topGO)
#
# Two modes, called by Snakemake:
#
#   Mode 1 — build shared maps (runs once):
#     Rscript scripts/03_go_analysis.R \
#       build_maps <annotation> <out_gene_go> <out_mstrg_map>
#
#   Mode 2 — per-category enrichment:
#     Rscript scripts/03_go_analysis.R \
#       enrich <bed_file> <gene_go> <mstrg_map> <out_refgenes> <out_result>
# ============================================================

suppressPackageStartupMessages({
  library(topGO)
  library(yaml)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
MODE <- args[1]

# ---- GO parameters (from config if present, else defaults) ----
GO_PVAL <- 0.05; TOP_N <- 50; ALGO <- "weight01"; STAT <- "fisher"
if (file.exists("config/config.yaml")) {
  cfg     <- yaml.load_file("config/config.yaml")
  GO_PVAL <- cfg$go$pvalue_cutoff %||% 0.05
  TOP_N   <- cfg$go$top_nodes     %||% 50
  ALGO    <- cfg$go$algorithm     %||% "weight01"
  STAT    <- cfg$go$statistic     %||% "fisher"
}

# ============================================================
# Shared helper: load eggNOG annotation
# ============================================================
load_annotation <- function(annot_file) {
  cat("[", format(Sys.time()), "] Loading annotation:", annot_file, "\n")
  annot <- read.table(
    gzfile(annot_file), sep = "\t", header = FALSE,
    comment.char = "#", stringsAsFactors = FALSE, fill = TRUE
  )
  col_names <- c("query","seed_ortholog","evalue","score",
                 "eggNOG_OGs","max_annot_lvl","COG_category",
                 "Description","Preferred_name","GOs",
                 "EC","KEGG_ko","KEGG_Pathway","KEGG_Module",
                 "KEGG_Reaction","KEGG_rclass","BRITE",
                 "KEGG_TC","CAZy","BiGG_Reaction","PFAMs")
  n <- min(ncol(annot), length(col_names))
  colnames(annot)[1:n] <- col_names[1:n]
  annot$gene_id <- sub("\\.(mRNA|t)\\.[0-9]+$", "", annot$query)
  annot$gene_id <- sub("\\.mRNA[0-9]+$",         "", annot$gene_id)
  annot
}

# ============================================================
# MODE 1: build_maps
# Args: build_maps <annotation> <out_gene_go> <out_mstrg_map>
# ============================================================
if (MODE == "build_maps") {
  if (length(args) < 4)
    stop("Usage: build_maps <annotation> <out_gene_go> <out_mstrg_map>")

  ANNOT_FILE  <- args[2]
  OUT_GENE_GO <- args[3]
  OUT_MSTRG   <- args[4]

  dir.create(dirname(OUT_GENE_GO), recursive = TRUE, showWarnings = FALSE)

  annot <- load_annotation(ANNOT_FILE)

  # Build gene → GO map
  annot_go <- annot[!is.na(annot$GOs) & annot$GOs != "" & annot$GOs != "-", ]
  gene_go_list <- list()
  for (i in seq_len(nrow(annot_go))) {
    gid <- annot_go$gene_id[i]
    gos <- trimws(strsplit(annot_go$GOs[i], ",")[[1]])
    gos <- gos[grepl("^GO:", gos)]
    if (length(gos) > 0)
      gene_go_list[[gid]] <- unique(c(gene_go_list[[gid]], gos))
  }
  gene_go_df <- do.call(rbind, lapply(names(gene_go_list), function(g) {
    data.frame(gene = g, GO = paste(gene_go_list[[g]], collapse = ","),
               stringsAsFactors = FALSE)
  }))
  write.table(gene_go_df, OUT_GENE_GO, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("  gene→GO map:", nrow(gene_go_df), "entries written\n")

  # Write empty MSTRG map placeholder (populated per-category in enrich mode)
  write.table(data.frame(mstrg = character(), ref = character()),
              OUT_MSTRG, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  MSTRG map placeholder written\n")
  cat("[", format(Sys.time()), "] build_maps done\n")

# ============================================================
# MODE 2: enrich
# Args: enrich <bed> <gene_go> <mstrg_map> <out_refgenes> <out_result>
# ============================================================
} else if (MODE == "enrich") {
  if (length(args) < 6)
    stop("Usage: enrich <bed> <gene_go> <mstrg_map> <out_refgenes> <out_result>")

  BED_FILE     <- args[2]
  GENE_GO_FILE <- args[3]
  MSTRG_FILE   <- args[4]
  OUT_REFGENES <- args[5]
  OUT_RESULT   <- args[6]

  dir.create(dirname(OUT_REFGENES), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(OUT_RESULT),   recursive = TRUE, showWarnings = FALSE)

  cat("[", format(Sys.time()), "] GO enrichment for:", BED_FILE, "\n")

  # Load BED
  bed     <- read.table(BED_FILE, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE)
  all_ids <- unique(unlist(strsplit(bed$V4, ",")))
  mstrg_ids <- all_ids[grepl("^MSTRG", all_ids)]
  known_ids <- all_ids[!grepl("^MSTRG", all_ids)]
  cat("  Regions:", nrow(bed), "| MSTRG IDs:", length(mstrg_ids),
      "| Known IDs:", length(known_ids), "\n")

  # Build MSTRG → ref gene map from BED col 4 co-occurrences
  mstrg_map <- data.frame(mstrg = character(), ref = character(),
                          stringsAsFactors = FALSE)
  for (col4 in bed$V4) {
    ids        <- strsplit(col4, ",")[[1]]
    m_here     <- ids[grepl("^MSTRG", ids)]
    k_here     <- ids[!grepl("^MSTRG", ids)]
    if (length(m_here) > 0 && length(k_here) > 0) {
      for (m in m_here) for (k in k_here)
        mstrg_map <- rbind(mstrg_map,
                           data.frame(mstrg = m, ref = k,
                                      stringsAsFactors = FALSE))
    }
  }
  mstrg_map <- unique(mstrg_map)

  # Resolve ref genes
  resolved   <- unique(c(known_ids, mstrg_map$ref[mstrg_map$mstrg %in% mstrg_ids]))
  resolved   <- resolved[nchar(resolved) > 0]
  write(resolved, OUT_REFGENES)
  cat("  Resolved ref genes:", length(resolved), "\n")

  # Load gene→GO map
  gene_go_df   <- read.table(GENE_GO_FILE, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE)
  gene_go_list <- setNames(
    lapply(gene_go_df$GO, function(x) strsplit(x, ",")[[1]]),
    gene_go_df$gene
  )

  # topGO enrichment
  all_genes <- names(gene_go_list)
  gene_list <- factor(as.integer(all_genes %in% resolved))
  names(gene_list) <- all_genes

  if (sum(gene_list == 1) == 0) {
    cat("  WARNING: No genes of interest in GO annotation — writing empty result\n")
    write.table(data.frame(), OUT_RESULT, sep = "\t",
                row.names = FALSE, quote = FALSE)
    quit(status = 0)
  }

  go_data <- new("topGOdata",
                 ontology = "BP",
                 allGenes = gene_list,
                 annot    = annFUN.gene2GO,
                 gene2GO  = gene_go_list,
                 nodeSize = 5)

  result   <- runTest(go_data, algorithm = ALGO, statistic = STAT)
  go_table <- GenTable(go_data, pvalue = result,
                       topNodes = min(TOP_N, length(score(result))),
                       numChar  = 100)
  go_table$FDR <- p.adjust(as.numeric(go_table$pvalue), method = "BH")
  go_sig       <- go_table[go_table$FDR < GO_PVAL, ]

  cat("  Significant GO terms (FDR <", GO_PVAL, "):", nrow(go_sig), "\n")
  write.table(go_sig, OUT_RESULT, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("[", format(Sys.time()), "] enrich done\n")

} else {
  stop("Unknown mode '", MODE, "'. Use 'build_maps' or 'enrich'.")
}