# ============================================================
# Script 05: Generate all visualisation plots
#
# Reads config/config.yaml. Produces:
#   - UpSet plot (all samples)
#   - Region counts bar chart
#   - GO enrichment bar plots per category
#   - KEGG bubble plot per category
#
# Usage:
#   micromamba activate goanalysis
#   Rscript scripts/05_plots.R
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
  library(UpSetR)
})

# ---- Load config --------------------------------------------
config    <- yaml.load_file("config/config.yaml")
SAMPLES   <- config$samples$names
OUT_BASE  <- config$output$base
OUT_GO    <- config$output$go
OUT_KEGG  <- config$output$kegg
OUT_PLOTS <- config$output$plots

dir.create(OUT_PLOTS, recursive = TRUE, showWarnings = FALSE)
dir.create(config$slurm$logs, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(config$slurm$logs, "plots.log")
con <- file(log_file, open = "at")
msg <- function(...) {
  m <- paste0("[", format(Sys.time()), "] ", ...)
  message(m); cat(m, "\n", file = con)
}

# ---- Category colours (extended for embryo) -----------------
cat_colors <- c(
  "needle_specific"         = "#C6EFCE",
  "root_specific"           = "#FFEB9C",
  "cold_specific"           = "#BDD7EE",
  "drought_specific"        = "#FCE4D6",
  "somatic_embryo_specific" = "#E2CFEE",
  "zygotic_embryo_specific" = "#FFD6CC"
)

# Active categories (non-empty BED files)
all_cats <- names(cat_colors)
cats <- Filter(function(cat) {
  f <- file.path(OUT_BASE, paste0(cat, ".bed"))
  file.exists(f) && file.info(f)$size > 0
}, all_cats)

msg("Samples: ", paste(SAMPLES, collapse = ", "))
msg("Categories: ", paste(cats, collapse = ", "))

# ============================================================
# Plot 1: UpSet plot
# ============================================================
msg("Generating UpSet plot")

multiinter_file <- file.path(OUT_BASE, "multiinter_output.bed")
if (file.exists(multiinter_file) && file.info(multiinter_file)$size > 0) {
  mi <- read.table(multiinter_file, header = FALSE, sep = "\t")

  n_fixed <- 5   # chr, start, end, count, names
  n_samp  <- length(SAMPLES)

  if (ncol(mi) >= n_fixed + n_samp) {
    samp_cols <- mi[, (n_fixed + 1):(n_fixed + n_samp)]
    colnames(samp_cols) <- SAMPLES

    png(file.path(OUT_PLOTS, "upset_plot.png"),
        width = max(1400, 200 * n_samp), height = 900, res = 150)
    print(upset(
      samp_cols,
      sets     = SAMPLES,
      mb.ratio = c(0.65, 0.35),
      order.by = "freq",
      text.scale = 1.5,
      point.size = 3
    ))
    dev.off()
    msg("  upset_plot.png done")
  } else {
    msg("  WARNING: multiinter columns don't match sample count — skipping UpSet")
  }
} else {
  msg("  WARNING: multiinter_output.bed missing or empty — skipping UpSet")
}

# ============================================================
# Plot 2: Region counts bar chart
# ============================================================
msg("Generating region counts bar chart")

bed_files <- c(
  "conserved.bed",
  paste0(cats, ".bed")
)

counts <- sapply(bed_files, function(f) {
  full <- file.path(OUT_BASE, f)
  if (file.exists(full) && file.info(full)$size > 0) {
    nrow(read.table(full, header = FALSE, sep = "\t"))
  } else {
    0L
  }
})

count_df <- data.frame(
  category = sub("\\.bed$", "", names(counts)),
  count    = counts,
  stringsAsFactors = FALSE
)

count_df$color <- ifelse(
  count_df$category == "conserved",
  "#D9D9D9",
  cat_colors[count_df$category]
)
count_df$color[is.na(count_df$color)] <- "#CCCCCC"

p <- ggplot(count_df, aes(x = reorder(category, -count), y = count, fill = category)) +
  geom_col(fill = count_df$color, color = "grey40", width = 0.7) +
  geom_text(aes(label = count), vjust = -0.4, size = 3.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10)) +
  labs(title = paste("Region counts —", config$species),
       x = NULL, y = "Number of regions") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

ggsave(file.path(OUT_PLOTS, "region_counts_bar.png"),
       plot = p, width = max(8, 1.5 * length(bed_files)), height = 5, dpi = 150)
msg("  region_counts_bar.png done")

# ============================================================
# Plot 3: GO enrichment bar plots per category
# ============================================================
msg("Generating GO bar plots")

for (cat in cats) {
  go_file <- file.path(OUT_GO, paste0(cat, "_BP_GO_enrichment.txt"))
  if (!file.exists(go_file) || file.info(go_file)$size == 0) next

  go_df <- read.table(go_file, header = TRUE, sep = "\t", quote = "")
  if (nrow(go_df) == 0) next

  go_df$pvalue <- as.numeric(go_df$pvalue)
  go_df <- go_df[!is.na(go_df$pvalue), ]
  go_df <- head(go_df[order(go_df$pvalue), ], 20)
  go_df$Term <- factor(go_df$Term, levels = rev(go_df$Term))

  fill_col <- ifelse(!is.na(cat_colors[cat]), cat_colors[cat], "#CCCCCC")

  p_go <- ggplot(go_df, aes(x = Term, y = -log10(pvalue))) +
    geom_col(fill = fill_col, color = "grey40", width = 0.7) +
    coord_flip() +
    theme_classic() +
    labs(title = paste("GO enrichment (BP):", gsub("_", " ", cat)),
         subtitle = config$species,
         x = NULL, y = "-log10(p-value)") +
    theme(axis.text.y = element_text(size = 8))

  out_name <- paste0("GO_bar_", cat, ".png")
  ggsave(file.path(OUT_PLOTS, out_name),
         plot = p_go, width = 10, height = 6, dpi = 150)
  msg("  ", out_name, " done")
}

# ============================================================
# Plot 4: KEGG bubble plot per category
# ============================================================
msg("Generating KEGG bubble plots")

for (cat in cats) {
  kegg_file <- file.path(OUT_KEGG, paste0(cat, "_KEGG_enrichment.txt"))
  if (!file.exists(kegg_file) || file.info(kegg_file)$size == 0) next

  kegg_df <- read.table(kegg_file, header = TRUE, sep = "\t", quote = "")
  if (nrow(kegg_df) == 0) next

  kegg_df$pvalue <- as.numeric(kegg_df$pvalue)
  kegg_df <- kegg_df[!is.na(kegg_df$pvalue), ]
  kegg_df <- head(kegg_df[order(kegg_df$pvalue), ], 20)
  kegg_df$name <- factor(kegg_df$name, levels = rev(kegg_df$name))

  fill_col <- ifelse(!is.na(cat_colors[cat]), cat_colors[cat], "#CCCCCC")

  p_kegg <- ggplot(kegg_df,
                   aes(x = -log10(pvalue), y = name, size = in_category)) +
    geom_point(color = fill_col, alpha = 0.8) +
    theme_classic() +
    labs(title = paste("KEGG enrichment:", gsub("_", " ", cat)),
         subtitle = config$species,
         x = "-log10(p-value)", y = NULL,
         size = "Genes in category") +
    theme(axis.text.y = element_text(size = 8))

  out_name <- paste0("KEGG_bubble_", cat, ".png")
  ggsave(file.path(OUT_PLOTS, out_name),
         plot = p_kegg, width = 10, height = 6, dpi = 150)
  msg("  ", out_name, " done")
}

close(con)
message("Plots complete. Output in: ", OUT_PLOTS)