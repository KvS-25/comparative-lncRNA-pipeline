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
# Pathway filter: exclude animal/human/disease/xenobiotic
#
# EXCLUDED:
#   map05xxx  - Human diseases (cancer, infection, neurodegeneration)
#   map04xxx  - Animal signalling (immune, neural, endocrine, etc.)
#   map00980/982/983, map0152x - Xenobiotic/drug metabolism
#
# RETAINED:
#   map00xxx  - Core metabolism
#   map01xxx  - Global/overview metabolism (excl. drug resistance)
#   map03xxx  - Genetic information processing
#   map04016  - MAPK signaling - plant
#   map04075  - Plant hormone signal transduction
#   map04120  - Ubiquitin mediated proteolysis
#   map04130  - SNARE interactions
#   map04136/138/141 - Autophagy, ER protein processing
#   map04144/146  - Endocytosis, peroxisome
# ============================================================
EXCLUDE_PATHWAYS <- c(
  # Animal/human map04xxx
  "map04010","map04012","map04014","map04015","map04020","map04022","map04024",
  "map04060","map04061","map04062","map04064","map04066","map04068","map04070",
  "map04071","map04072",
  "map04110","map04113","map04114","map04115",
  "map04150","map04151","map04152",
  "map04210","map04212","map04213",
  "map04310","map04330","map04340","map04350","map04360","map04370","map04371",
  "map04380","map04390","map04392",
  "map04510","map04512","map04514","map04520","map04530","map04540","map04550",
  "map04610","map04611","map04612","map04614",
  "map04620","map04621","map04622","map04623","map04625","map04630","map04640",
  "map04650","map04657","map04658","map04659","map04660","map04662","map04664",
  "map04666","map04668","map04670","map04672",
  "map04710","map04711","map04713",
  "map04720","map04721","map04722","map04723","map04724","map04725","map04726",
  "map04727","map04728","map04730","map04740","map04742","map04744","map04750",
  "map04810",
  "map04910","map04911","map04912","map04913","map04914","map04915","map04916",
  "map04917","map04918","map04919","map04920","map04921","map04922","map04923",
  "map04924","map04925","map04926","map04927","map04928","map04929",
  "map04930","map04931","map04932","map04933","map04934",
  "map04960","map04961","map04962","map04964","map04966",
  "map04970","map04971","map04972","map04973","map04974","map04975","map04976",
  "map04977","map04978","map04979","map04980",
  # Xenobiotic and drug-resistance
  "map00980","map00982","map00983",
  "map01521","map01522","map01523","map01524","map01526"
)

is_excluded <- function(pathway_id) {
  pathway_id %in% EXCLUDE_PATHWAYS | grepl("^map05", pathway_id)
}

# ============================================================
# Shared helper: load eggNOG annotation
# ============================================================
load_annotation <- function(annot_file) {
  cat("[", format(Sys.time()), "] Loading annotation:", annot_file, "\n")
  lines <- readLines(gzfile(annot_file))
  lines <- lines[!grepl("^##", lines)]
  annot <- read.table(text = lines, sep = "\t", header = FALSE,
                      comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
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

  # pathway → gene map using KEGG_Pathway column (map-prefixed IDs)
  annot_path <- annot[!is.na(annot$KEGG_Pathway) &
                        annot$KEGG_Pathway != "" & annot$KEGG_Pathway != "-", ]
  pathway_genes <- list()
  for (i in seq_len(nrow(annot_path))) {
    gid      <- annot_path$gene_id[i]
    pathways <- trimws(strsplit(annot_path$KEGG_Pathway[i], ",")[[1]])
    pathways <- pathways[grepl("^map", pathways)]
    # Apply biological relevance filter at build time
    pathways <- pathways[!is_excluded(pathways)]
    for (p in pathways)
      pathway_genes[[p]] <- unique(c(pathway_genes[[p]], gid))
  }

  cat("  Pathways after biological relevance filter:", length(pathway_genes), "\n")

  # Fetch pathway names from KEGG REST API
  cat("[", format(Sys.time()), "] Fetching pathway names from KEGG API\n")
  pw_ids   <- names(pathway_genes)
  pw_names <- data.frame(Pathway = pw_ids, Name = NA_character_,
                         stringsAsFactors = FALSE)
  for (i in seq_along(pw_ids)) {
    tryCatch({
      url  <- paste0("https://rest.kegg.jp/get/", pw_ids[i])
      txt  <- readLines(url, warn = FALSE)
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
  cat("  Ref genes:", length(refgenes), "\n")

  # Load filtered pathway→gene map from RDS saved during build_maps
  rds_file <- sub("\\.txt$", "_map.rds", PATHWAY_NAMES_FILE)
  if (!file.exists(rds_file)) {
    cat("  WARNING: pathway map RDS not found — writing empty result\n")
    write.table(data.frame(), OUT_RESULT, sep = "\t",
                row.names = FALSE, quote = FALSE)
    quit(status = 0)
  }
  pathway_genes <- readRDS(rds_file)

  # Load pathway names for annotation
  pw_names <- read.table(PATHWAY_NAMES_FILE, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, quote = "")

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
    data.frame(pathway_id  = pw,
               name        = pw_names$Name[match(pw, pw_names$Pathway)],
               total       = m,
               expected    = round(K * m / N, 2),
               in_category = k,
               pvalue      = pv,
               stringsAsFactors = FALSE)
  })
  results_df     <- do.call(rbind, results)
  results_df$FDR <- p.adjust(results_df$pvalue, method = "BH")
  results_df     <- results_df[order(results_df$FDR), ]
  results_sig    <- results_df[results_df$FDR < KEGG_PVAL &
                                 results_df$in_category > 0, ]

  cat("  Significant pathways (FDR <", KEGG_PVAL, "):", nrow(results_sig), "\n")
  write.table(results_sig, OUT_RESULT, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("[", format(Sys.time()), "] enrich done\n")

} else {
  stop("Unknown mode '", MODE, "'. Use 'build_maps' or 'enrich'.")
}