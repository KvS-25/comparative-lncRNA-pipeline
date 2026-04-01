# ============================================================
# Script 04: KEGG pathway enrichment analysis
#
# Reads config/config.yaml. Iterates over all non-empty
# category refgene files produced by 03_go_analysis.R.
# Uses Fisher's exact test against all annotated genes.
#
# Usage:
#   micromamba activate goanalysis
#   Rscript scripts/04_kegg_analysis.R
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
})

# ---- Load config --------------------------------------------
config   <- yaml.load_file("config/config.yaml")
OUT_GO   <- config$output$go
OUT_KEGG <- config$output$kegg
ANNOT    <- config$genome$annotation
P_CUTOFF <- config$kegg$pvalue_cutoff

dir.create(OUT_KEGG, recursive = TRUE, showWarnings = FALSE)
dir.create(config$slurm$logs, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(config$slurm$logs, "kegg_analysis.log")
con <- file(log_file, open = "at")
msg <- function(...) {
  m <- paste0("[", format(Sys.time()), "] ", ...)
  message(m); cat(m, "\n", file = con)
}

# ---- Categories (match 03 output) ---------------------------
all_cats <- c(
  "needle_specific", "root_specific",
  "cold_specific",   "drought_specific",
  "somatic_embryo_specific", "zygotic_embryo_specific"
)

cats <- Filter(function(cat) {
  f <- file.path(OUT_GO, paste0(cat, "_refgenes.txt"))
  file.exists(f) && file.info(f)$size > 0
}, all_cats)

if (length(cats) == 0) {
  msg("No refgene files found from GO step — nothing to do.")
  quit(status = 0)
}
msg("Categories: ", paste(cats, collapse = ", "))

# ---- Load annotation ----------------------------------------
if (is.null(ANNOT) || nchar(ANNOT) == 0 || !file.exists(ANNOT)) {
  stop("Annotation file not set or not found.")
}
msg("Loading annotation: ", ANNOT)
annot <- read.table(gzfile(ANNOT), sep = "\t", header = FALSE,
                    comment.char = "#", quote = "", fill = TRUE)

# col 1 = gene, col 13 = KEGG pathway
gene2KEGG <- list()
for (i in seq_len(nrow(annot))) {
  gene  <- as.character(annot[i, 1])
  keggs <- as.character(annot[i, 13])
  if (nchar(keggs) > 0 && keggs != "-") {
    gene2KEGG[[gene]] <- unlist(strsplit(keggs, ","))
  }
}

# Build pathway -> gene list
pathway2gene <- list()
for (gene in names(gene2KEGG)) {
  for (pw in gene2KEGG[[gene]]) {
    pathway2gene[[pw]] <- c(pathway2gene[[pw]], gene)
  }
}

all_annotated <- unique(names(gene2KEGG))
N_all         <- length(all_annotated)

write.table(
  data.frame(gene    = names(gene2KEGG),
             pathway = sapply(gene2KEGG, paste, collapse = ",")),
  file = file.path(OUT_KEGG, "gene_to_KEGG.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
msg("gene_to_KEGG.txt: ", length(gene2KEGG), " genes, ",
    length(pathway2gene), " pathways")

# ---- Fetch KEGG pathway names (requires internet) -----------
pathway_names_file <- file.path(OUT_KEGG, "kegg_pathway_names.txt")
if (!file.exists(pathway_names_file)) {
  msg("Fetching KEGG pathway names...")
  tryCatch({
    url  <- "https://rest.kegg.jp/list/pathway"
    raw  <- readLines(url)
    df   <- do.call(rbind, strsplit(raw, "\t"))
    write.table(df, file = pathway_names_file,
                sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    msg("Pathway names saved: ", nrow(df), " pathways")
  }, error = function(e) {
    msg("WARNING: Could not fetch KEGG names (no internet?): ", e$message)
    file.create(pathway_names_file)
  })
}

pathway_names <- tryCatch({
  pn <- read.table(pathway_names_file, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
  setNames(pn[, 2], pn[, 1])
}, error = function(e) character(0))

# ---- Process each category ----------------------------------
for (cat in cats) {
  msg("Processing: ", cat)

  cat_genes <- readLines(file.path(OUT_GO, paste0(cat, "_refgenes.txt")))
  cat_genes <- cat_genes[nchar(cat_genes) > 0]
  n_cat     <- length(cat_genes)

  if (n_cat == 0) {
    msg("  No genes for ", cat, " — skipping")
    file.create(file.path(OUT_KEGG, paste0(cat, "_KEGG_enrichment.txt")))
    next
  }

  results <- lapply(names(pathway2gene), function(pw) {
    pw_genes  <- pathway2gene[[pw]]
    in_cat    <- sum(cat_genes   %in% pw_genes)
    not_cat   <- sum(all_annotated[!all_annotated %in% cat_genes] %in% pw_genes)
    mat <- matrix(c(in_cat,
                    n_cat - in_cat,
                    not_cat,
                    N_all - n_cat - not_cat),
                  nrow = 2)
    p <- fisher.test(mat, alternative = "greater")$p.value
    data.frame(
      pathway    = pw,
      name       = ifelse(pw %in% names(pathway_names), pathway_names[[pw]], pw),
      in_category = in_cat,
      in_pathway  = length(pw_genes),
      total_genes = N_all,
      pvalue      = p,
      stringsAsFactors = FALSE
    )
  })

  result_df  <- do.call(rbind, results)
  result_df  <- result_df[order(result_df$pvalue), ]
  sig_df     <- result_df[result_df$pvalue < P_CUTOFF, ]

  out_file <- file.path(OUT_KEGG, paste0(cat, "_KEGG_enrichment.txt"))
  write.table(sig_df, file = out_file,
              sep = "\t", row.names = FALSE, quote = FALSE)
  msg("  Significant pathways: ", nrow(sig_df), " -> ", out_file)
}

close(con)
message("KEGG analysis complete.")