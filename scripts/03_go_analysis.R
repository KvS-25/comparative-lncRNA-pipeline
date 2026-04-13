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
#       enrich <bed_file> <gene_go> <mstrg_map> \
#              <out_refgenes> <out_bp> <out_mf> <out_cc>
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
# GO term filter: exclude animal/human/disease-specific terms
#
# Applied after enrichment so topGO graph structure is not
# disturbed during scoring. Patterns matched case-insensitively.
# ============================================================
EXCLUDE_GO_PATTERNS <- c(
  # Vertebrate immune system
  "\\bimmunoglobulin\\b", "\\bimmune response\\b", "\\bimmunity\\b",
  "\\bleukocyte\\b", "\\blymphocyte\\b", "\\bB cell\\b", "\\bT cell\\b",
  "\\bnatural killer\\b", "\\bneutrophil\\b", "\\bmacrophage\\b",
  "\\bdendritic cell\\b", "\\bmast cell\\b", "\\bbasophil\\b",
  "\\beosinophil\\b", "\\bplasma cell\\b", "\\bisotype switching\\b",
  "\\bantibody\\b", "\\bantigen\\b", "\\bcytokine\\b", "\\binterleukin\\b",
  "\\btumor necrosis factor\\b", "\\binterferon\\b", "\\bcomplement\\b",
  "\\binnate immune\\b", "\\badaptive immune\\b",
  # Animal-specific cell types / tissues
  "\\bneuron\\b", "\\bneuronal\\b", "\\bneurotransmitter\\b",
  "\\bsynapse\\b", "\\bsynaptic\\b", "\\bnervous system\\b",
  "\\bmuscle cell\\b", "\\bcardiomyocyte\\b", "\\berythrocyte\\b",
  "\\bplatelet\\b", "\\bblood cell\\b", "\\bblood pressure\\b",
  "\\bvasodilation\\b", "\\bvasoconstriction\\b",
  "\\bbronchodilator\\b", "\\bbronchial\\b",
  # Mammalian-specific signalling / hormones
  "\\binsulin\\b", "\\bglucagon\\b", "\\bcortisol\\b",
  "\\badrenaline\\b", "\\bepinephrine\\b", "\\bnorepinephrine\\b",
  "\\bdopamine\\b", "\\bserotonin\\b", "\\bendorphin\\b",
  "\\btestosterone\\b", "\\bestrogen\\b", "\\bprogesterone\\b",
  "\\bandrogen\\b", "\\bsteroid hormone\\b",
  # Animal-specific processes
  "\\bToll[- ]Imd\\b", "\\bDrosophila\\b",
  "\\bsomatic diversification\\b", "\\bsomatic cell DNA recombination\\b",
  "\\bmulti-organism reproductive process\\b",
  "\\bfertilization\\b", "\\bhatching\\b"
)

is_excluded_go <- function(term_names) {
  pattern <- paste(EXCLUDE_GO_PATTERNS, collapse = "|")
  grepl(pattern, term_names, ignore.case = TRUE, perl = TRUE)
}

# ============================================================
# Shared helper: load eggNOG annotation (# bug-safe version)
# ============================================================
load_annotation <- function(annot_file) {
  cat("[", format(Sys.time()), "] Loading annotation:", annot_file, "\n")
  lines <- readLines(gzfile(annot_file))
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nchar(trimws(lines)) > 0]
  annot <- read.table(
    text = lines, sep = "\t", header = FALSE,
    stringsAsFactors = FALSE, fill = TRUE, quote = ""
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
# Helper: run topGO for one ontology, return filtered df
# ============================================================
run_topgo <- function(gene_list, gene_go_list, ontology,
                      algo, stat, top_n, pval_cutoff) {
  if (sum(gene_list == 1) == 0) {
    cat("  WARNING: No genes of interest — skipping", ontology, "\n")
    return(data.frame())
  }

  go_data <- tryCatch(
    new("topGOdata",
        ontology = ontology,
        allGenes = gene_list,
        annot    = annFUN.gene2GO,
        gene2GO  = gene_go_list,
        nodeSize = 5),
    error = function(e) {
      cat("  WARNING: topGOdata failed for", ontology, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(go_data)) return(data.frame())

  result   <- runTest(go_data, algorithm = algo, statistic = stat)
  go_table <- GenTable(go_data, pvalue = result,
                       topNodes = min(top_n, length(score(result))),
                       numChar  = 100)
  go_table$pvalue <- as.numeric(go_table$pvalue)
  go_table <- go_table[!is.na(go_table$pvalue), ]
  go_table$FDR    <- p.adjust(go_table$pvalue, method = "BH")
  go_table$domain <- ontology

  go_sig <- go_table[go_table$FDR < pval_cutoff, ]

  n_before <- nrow(go_sig)
  go_sig   <- go_sig[!is_excluded_go(go_sig$Term), ]
  if (n_before > nrow(go_sig))
    cat("  Filtered", n_before - nrow(go_sig),
        "non-plant GO terms from", ontology, "\n")

  cat("  Significant", ontology, "terms (after filter):", nrow(go_sig), "\n")
  go_sig
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
  cat("  gene->GO map:", nrow(gene_go_df), "entries written\n")

  write.table(data.frame(mstrg = character(), ref = character()),
              OUT_MSTRG, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  MSTRG map placeholder written\n")
  cat("[", format(Sys.time()), "] build_maps done\n")

# ============================================================
# MODE 2: enrich
# Args: enrich <bed> <gene_go> <mstrg_map>
#              <out_refgenes> <out_bp> <out_mf> <out_cc>
# ============================================================
} else if (MODE == "enrich") {
  if (length(args) < 8)
    stop("Usage: enrich <bed> <gene_go> <mstrg_map> <out_refgenes> <out_bp> <out_mf> <out_cc>")

  BED_FILE     <- args[2]
  GENE_GO_FILE <- args[3]
  MSTRG_FILE   <- args[4]
  OUT_REFGENES <- args[5]
  OUT_BP       <- args[6]
  OUT_MF       <- args[7]
  OUT_CC       <- args[8]

  for (p in c(OUT_REFGENES, OUT_BP, OUT_MF, OUT_CC))
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

  cat("[", format(Sys.time()), "] GO enrichment for:", BED_FILE, "\n")

  bed <- read.table(BED_FILE, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)

  if (all(grepl("^[0-9]+$", as.character(bed$V4)))) {
    cat("  Detected transcript-based BED (spruce-style)\n")
    all_ids   <- unique(as.character(bed$V1))
    mstrg_ids <- all_ids[grepl("^MSTRG", all_ids)]
    known_ids <- all_ids[!grepl("^MSTRG", all_ids)]
    resolved  <- unique(known_ids)
  } else {
    cat("  Detected chromosome-based BED (pine-style)\n")
    all_ids   <- unique(unlist(strsplit(as.character(bed$V4), ",")))
    mstrg_ids <- all_ids[grepl("^MSTRG", all_ids)]
    known_ids <- all_ids[!grepl("^MSTRG", all_ids)]

    mstrg_map <- data.frame(mstrg = character(), ref = character(),
                            stringsAsFactors = FALSE)
    for (col4 in bed$V4) {
      ids    <- strsplit(col4, ",")[[1]]
      m_here <- ids[grepl("^MSTRG", ids)]
      k_here <- ids[!grepl("^MSTRG", ids)]
      if (length(m_here) > 0 && length(k_here) > 0)
        for (m in m_here) for (k in k_here)
          mstrg_map <- rbind(mstrg_map,
                             data.frame(mstrg = m, ref = k,
                                        stringsAsFactors = FALSE))
    }
    mstrg_map <- unique(mstrg_map)
    resolved  <- unique(c(known_ids,
                          mstrg_map$ref[mstrg_map$mstrg %in% mstrg_ids]))
  }

  resolved <- strip_isoform(resolved)
  resolved <- resolved[nchar(resolved) > 0]
  write(resolved, OUT_REFGENES)
  cat("  Regions:", nrow(bed), "| MSTRG IDs:", length(mstrg_ids),
      "| Known IDs:", length(known_ids), "\n")
  cat("  Resolved ref genes:", length(resolved), "\n")

  gene_go_df   <- read.table(GENE_GO_FILE, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE)
  gene_go_list <- setNames(
    lapply(gene_go_df$GO, function(x) strsplit(x, ",")[[1]]),
    gene_go_df$gene
  )

  all_genes <- names(gene_go_list)
  gene_list <- factor(as.integer(all_genes %in% resolved))
  names(gene_list) <- all_genes

  if (sum(gene_list == 1) == 0) {
    cat("  WARNING: No genes of interest in GO annotation — writing empty results\n")
    for (p in c(OUT_BP, OUT_MF, OUT_CC))
      write.table(data.frame(), p, sep = "\t", row.names = FALSE, quote = FALSE)
    quit(status = 0)
  }

  for (domain_info in list(
    list(ont = "BP", out = OUT_BP),
    list(ont = "MF", out = OUT_MF),
    list(ont = "CC", out = OUT_CC)
  )) {
    result <- run_topgo(gene_list, gene_go_list,
                        ontology    = domain_info$ont,
                        algo        = ALGO,
                        stat        = STAT,
                        top_n       = TOP_N,
                        pval_cutoff = GO_PVAL)
    write.table(result, domain_info$out, sep = "\t",
                row.names = FALSE, quote = FALSE)
  }

  cat("[", format(Sys.time()), "] enrich done\n")

} else {
  stop("Unknown mode '", MODE, "'. Use 'build_maps' or 'enrich'.")
}