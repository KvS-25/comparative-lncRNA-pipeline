# ============================================================
# Script 04: KEGG pathway enrichment analysis
#
# Two modes, called by Snakemake:
#
#   Mode 1 — build shared maps (runs once):
#     Rscript scripts/04_kegg_analysis.R \
#       build_maps <annotation> <out_gene_kegg> <out_pathway_names>
#
#   Mode 2 — per-category enrichment:
#     Rscript scripts/04_kegg_analysis.R \
#       enrich <refgenes_file> <gene_kegg> <pathway_names> <out_result>
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

args      <- commandArgs(trailingOnly = TRUE)
MODE      <- args[1]

KEGG_PVAL <- 0.05
if (file.exists("config/config.yaml")) {
  cfg       <- yaml.load_file("config/config.yaml")
  KEGG_PVAL <- cfg$kegg$pvalue_cutoff %||% 0.05
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
# Shared helper: strip isoform suffix from gene IDs
# ============================================================
strip_isoform <- function(ids) {
  ids <- sub("\\.(mRNA|t)\\.[0-9]+$", "", ids)
  ids <- sub("\\.mRNA[0-9]+$",         "", ids)
  unique(ids)
}

# ============================================================
# MODE 1: build_maps
# Args: build_maps <annotation> <out_gene_kegg> <out_pathway_names>
# ============================================================
if (MODE == "build_maps") {
  if (length(args) < 4)
    stop("Usage: build_maps <annotation> <out_gene_kegg> <out_pathway_names>")

  ANNOT_FILE        <- args[2]
  OUT_GENE_KEGG     <- args[3]
  OUT_PATHWAY_NAMES <- args[4]

  dir.create(dirname(OUT_GENE_KEGG),     recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(OUT_PATHWAY_NAMES), recursive = TRUE, showWarnings = FALSE)

  annot <- load_annotation(ANNOT_FILE)

  # gene → KEGG KO
  annot_kegg <- annot[!is.na(annot$KEGG_ko) &
                        annot$KEGG_ko != "" & annot$KEGG_ko != "-", ]
  gene_ko <- unique(data.frame(
    gene = annot_kegg$gene_id,
    KO   = annot_kegg$KEGG_ko,
    stringsAsFactors = FALSE
  ))
  gene_ko <- gene_ko[nchar(gene_ko$gene) > 0, ]
  write.table(gene_ko, OUT_GENE_KEGG, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("  gene→KEGG map:", nrow(gene_ko), "entries written\n")

  # pathway → gene map
  annot_path <- annot[!is.na(annot$KEGG_Pathway) &
                        annot$KEGG_Pathway != "" & annot$KEGG_Pathway != "-", ]
  pathway_genes <- list()
  for (i in seq_len(nrow(annot_path))) {
    gid      <- annot_path$gene_id[i]
    pathways <- trimws(strsplit(annot_path$KEGG_Pathway[i], ",")[[1]])
    pathways <- pathways[grepl("^ko", pathways)]
    for (p in pathways)
      pathway_genes[[p]] <- unique(c(pathway_genes[[p]], gid))
  }

  # Fetch pathway names from KEGG REST API
  cat("[", format(Sys.time()), "] Fetching pathway names from KEGG API\n")
  pw_ids   <- names(pathway_genes)
  pw_names <- data.frame(Pathway = pw_ids, Name = NA_character_,
                         stringsAsFactors = FALSE)
  for (i in seq_along(pw_ids)) {
    tryCatch({
      txt  <- readLines(paste0("https://rest.kegg.jp/get/", pw_ids[i]),
                        warn = FALSE)
      nl   <- txt[grepl("^NAME", txt)][1]
      if (!is.na(nl)) pw_names$Name[i] <- trimws(sub("^NAME\\s+", "", nl))
      Sys.sleep(0.3)
    }, error = function(e) {
      cat("  Warning: could not fetch", pw_ids[i], "\n")
    })
  }
  write.table(pw_names, OUT_PATHWAY_NAMES, sep = "\t",
              row.names = FALSE, quote = FALSE)
  saveRDS(pathway_genes,
          sub("\\.txt$", "_map.rds", OUT_PATHWAY_NAMES))
  cat("  Pathway names written:", nrow(pw_names), "pathways\n")
  cat("[", format(Sys.time()), "] build_maps done\n")

# ============================================================
# MODE 2: enrich
# Args: enrich <refgenes> <gene_kegg> <pathway_names> <out_result>
# ============================================================
} else if (MODE == "enrich") {
  if (length(args) < 5)
    stop("Usage: enrich <refgenes> <gene_kegg> <pathway_names> <out_result>")

  REFGENES_FILE      <- args[2]
  GENE_KEGG_FILE     <- args[3]
  PATHWAY_NAMES_FILE <- args[4]
  OUT_RESULT         <- args[5]

  dir.create(dirname(OUT_RESULT), recursive = TRUE, showWarnings = FALSE)

  cat("[", format(Sys.time()), "] KEGG enrichment for:", REFGENES_FILE, "\n")

  refgenes <- readLines(REFGENES_FILE)
  refgenes <- refgenes[nchar(trimws(refgenes)) > 0]

  # Strip isoform suffix to match annotation gene IDs
  refgenes <- strip_isoform(refgenes)
  cat("  Ref genes:", length(refgenes), "\n")

  # Load pathway→gene map
  rds_file <- sub("\\.txt$", "_map.rds", PATHWAY_NAMES_FILE)
  if (!file.exists(rds_file)) {
    cat("  WARNING: pathway map RDS not found — writing empty result\n")
    write.table(data.frame(), OUT_RESULT, sep = "\t",
                row.names = FALSE, quote = FALSE)
    quit(status = 0)
  }
  pathway_genes <- readRDS(rds_file)

  # Universe = all genes that appear in any pathway
  universe <- unique(unlist(pathway_genes))
  N <- length(universe)
  K <- sum(refgenes %in% universe)
  cat("  Universe:", N, "| Category genes in universe:", K, "\n")

  if (K == 0) {
    cat("  WARNING: No category genes in KEGG universe — writing empty result\n")
    write.table(data.frame(), OUT_RESULT, sep = "\t",
                row.names = FALSE, quote = FALSE)
    quit(status = 0)
  }

  results <- lapply(names(pathway_genes), function(pw) {
    pw_genes <- pathway_genes[[pw]]
    m  <- length(pw_genes)
    k  <- sum(refgenes %in% pw_genes)
    pv <- phyper(k - 1, m, N - m, K, lower.tail = FALSE)
    data.frame(Pathway = pw, Total = m,
               Expected = round(K * m / N, 2),
               Observed = k, pvalue = pv,
               stringsAsFactors = FALSE)
  })
  results_df     <- do.call(rbind, results)
  results_df$FDR <- p.adjust(results_df$pvalue, method = "BH")
  results_df     <- results_df[order(results_df$FDR), ]
  results_sig    <- results_df[results_df$FDR < KEGG_PVAL &
                                 results_df$Observed > 0, ]

  cat("  Significant pathways (FDR <", KEGG_PVAL, "):", nrow(results_sig), "\n")
  write.table(results_sig, OUT_RESULT, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("[", format(Sys.time()), "] enrich done\n")

} else {
  stop("Unknown mode '", MODE, "'. Use 'build_maps' or 'enrich'.")
}