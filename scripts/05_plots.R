# ============================================================
# Script 05: Generate all visualisation plots
#
# Reads config/config.yaml. Produces:
#   - UpSet plot (all samples)
#   - Region counts bar chart
#   - GO enrichment bar plots per category (BP/MF/CC panels)
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
config       <- yaml.load_file("config/config.yaml")
SPECIES      <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(SPECIES) || SPECIES == "") SPECIES <- config$default_species
SP           <- config[[SPECIES]]

SAMPLES      <- SP$samples$names
OUT_BASE     <- SP$output$base
OUT_GO       <- SP$output$go
OUT_KEGG     <- SP$output$kegg
OUT_PLOTS    <- SP$output$plots
LOGS_DIR     <- SP$slurm$logs
SPECIES_NAME <- SP$species_name

dir.create(OUT_PLOTS, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR,  recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGS_DIR, "plots.log")
con <- file(log_file, open = "at")
msg <- function(...) {
  m <- paste0("[", format(Sys.time()), "] ", ...)
  message(m); cat(m, "\n", file = con)
}

# ---- Category colours ---------------------------------------
cat_colors <- c(
  "needle_specific"         = "#C6EFCE",
  "root_specific"           = "#FFEB9C",
  "cold_specific"           = "#BDD7EE",
  "drought_specific"        = "#FCE4D6",
  "somatic_embryo_specific" = "#E2CFEE",
  "zygotic_embryo_specific" = "#FFD6CC"
)

# Domain colours — Pine 1 style
domain_colors <- c(
  "Biological Process" = "#4575B4",
  "Molecular Function" = "#D73027",
  "Cellular Component" = "#1A9850"
)

# Active categories (non-empty BED files)
cats <- Filter(function(cat) {
  f <- file.path(OUT_BASE, paste0(cat, ".bed"))
  file.exists(f) && file.info(f)$size > 0
}, names(cat_colors))

msg("Species: ", SPECIES_NAME)
msg("Samples: ", paste(SAMPLES, collapse = ", "))
msg("Categories: ", paste(cats, collapse = ", "))

# ============================================================
# Plot 1: UpSet plot
# ============================================================
msg("Generating UpSet plot")

multiinter_file <- file.path(OUT_BASE, "multiinter_output.bed")
if (file.exists(multiinter_file) && file.info(multiinter_file)$size > 0) {
  mi     <- read.table(multiinter_file, header = FALSE, sep = "\t")
  n_samp <- length(SAMPLES)

  if (ncol(mi) >= 5 + n_samp) {
    samp_cols <- mi[, 6:(5 + n_samp)]
    colnames(samp_cols) <- SAMPLES

    png(file.path(OUT_PLOTS, "upset_plot.png"),
        width = max(1400, 200 * n_samp), height = 900, res = 150)
    print(upset(
      samp_cols,
      sets       = SAMPLES,
      mb.ratio   = c(0.65, 0.35),
      order.by   = "freq",
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

bed_files <- c("conserved.bed", paste0(cats, ".bed"))

counts <- sapply(bed_files, function(f) {
  full <- file.path(OUT_BASE, f)
  if (file.exists(full) && file.info(full)$size > 0)
    nrow(read.table(full, header = FALSE, sep = "\t"))
  else 0L
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

p <- ggplot(count_df, aes(x = reorder(category, -count), y = count)) +
  geom_col(fill = count_df$color, color = "grey40", width = 0.7) +
  geom_text(aes(label = count), vjust = -0.4, size = 3.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10)) +
  labs(title = paste("Region counts —", SPECIES_NAME),
       x = NULL, y = "Number of regions") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

ggsave(file.path(OUT_PLOTS, "region_counts_bar.png"),
       plot = p, width = max(8, 1.5 * length(bed_files)), height = 5, dpi = 150)
msg("  region_counts_bar.png done")

# ============================================================
# Plot 3: GO enrichment bar plots per category
#
# Style: faceted BP / MF / CC panels, top 5 terms each,
# coloured by domain. Matches Pine 1 publication style.
# ============================================================
msg("Generating GO bar plots")

TOP_GO_PER_DOMAIN <- 5

domain_map <- c(
  "BP" = "Biological Process",
  "MF" = "Molecular Function",
  "CC" = "Cellular Component"
)

for (cat in cats) {

  panels <- list()
  for (domain_code in c("BP", "MF", "CC")) {
    go_file <- file.path(OUT_GO, paste0(cat, "_", domain_code, "_GO_enrichment.txt"))
    if (!file.exists(go_file) || file.info(go_file)$size == 0) next

    df <- tryCatch(
      read.table(go_file, header = TRUE, sep = "\t", quote = ""),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) next

    df$pvalue <- as.numeric(df$pvalue)
    df <- df[!is.na(df$pvalue), ]
    if (nrow(df) == 0) next

    df <- head(df[order(df$pvalue), ], TOP_GO_PER_DOMAIN)
    df$Domain <- domain_map[[domain_code]]
    panels[[domain_code]] <- df
  }

  if (length(panels) == 0) {
    msg("  No GO results for ", cat, " — skipping")
    next
  }

  combined <- do.call(rbind, panels)
  combined$Domain <- factor(combined$Domain,
                            levels = c("Biological Process",
                                       "Molecular Function",
                                       "Cellular Component"))

  # Preserve within-domain p-value ordering (best at top of each panel)
  combined$Term <- factor(
    combined$Term,
    levels = rev(unlist(lapply(rev(levels(combined$Domain)), function(d) {
      sub_df <- combined[combined$Domain == d, ]
      as.character(sub_df$Term[order(sub_df$pvalue, decreasing = TRUE)])
    })))
  )

  # Gene count label
  if ("Significant" %in% colnames(combined)) {
    combined$count_label <- combined$Significant
  } else {
    combined$count_label <- ""
  }

  p_go <- ggplot(combined, aes(x = -log10(pvalue), y = Term, fill = Domain)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = count_label), hjust = -0.2, size = 3.5) +
    scale_fill_manual(values = domain_colors) +
    facet_wrap(~Domain, scales = "free_y", ncol = 1) +
    labs(
      title = paste0("GO Enrichment — ",
                     gsub("_", " ", tools::toTitleCase(cat))),
      x = "-log10(p-value)", y = NULL
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "#F0F0F0"),
      strip.text       = element_text(face = "bold", size = 11),
      axis.text.y      = element_text(size = 10),
      axis.text.x      = element_text(size = 10),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position  = "none"
    ) +
    xlim(0, max(-log10(combined$pvalue), na.rm = TRUE) * 1.15)

  out_file <- file.path(OUT_PLOTS, paste0("GO_bar_", cat, ".png"))
  ggsave(out_file, plot = p_go,
         width = 10, height = 4 * length(panels), dpi = 300)
  msg("  GO_bar_", cat, ".png done (", length(panels), " panel(s))")
}

# ============================================================
# Plot 4: KEGG bubble plots per category
# ============================================================
msg("Generating KEGG bubble plots")

for (cat in cats) {
  kegg_file <- file.path(OUT_KEGG, paste0(cat, "_KEGG_enrichment.txt"))
  if (!file.exists(kegg_file) || file.info(kegg_file)$size == 0) {
    msg("  Skipping ", cat, " KEGG — file missing or empty")
    next
  }

  kegg_df <- tryCatch(
    read.table(kegg_file, header = TRUE, sep = "\t", quote = ""),
    error = function(e) NULL
  )
  if (is.null(kegg_df) || nrow(kegg_df) == 0) {
    msg("  Skipping ", cat, " KEGG — no results")
    next
  }

  kegg_df$pvalue <- as.numeric(kegg_df$pvalue)
  kegg_df <- kegg_df[!is.na(kegg_df$pvalue) & kegg_df$pvalue > 0, ]
  if (nrow(kegg_df) == 0) next

  kegg_df <- head(kegg_df[order(kegg_df$pvalue), ], 20)
  kegg_df$name <- factor(kegg_df$name, levels = rev(kegg_df$name))

  p_kegg <- ggplot(kegg_df,
                   aes(x     = in_category / total,
                       y     = name,
                       size  = in_category,
                       color = -log10(pvalue))) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(
      colours = c("green3", "yellow2", "orange", "red3"),
      name    = "-log10(p-value)"
    ) +
    scale_size_continuous(name = "Gene Count", range = c(3, 12)) +
    labs(
      title = paste0("KEGG Pathway Enrichment — ",
                     gsub("_", " ", tools::toTitleCase(cat))),
      x = "Gene Ratio", y = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.y  = element_text(size = 10),
      axis.text.x  = element_text(size = 10),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "right"
    )

  out_file <- file.path(OUT_PLOTS, paste0("KEGG_bubble_", cat, ".png"))
  ggsave(out_file, plot = p_kegg,
         width = 10, height = max(5, 0.35 * nrow(kegg_df) + 2), dpi = 300)
  msg("  KEGG_bubble_", cat, ".png done")
}

close(con)
message("Plots complete. Output in: ", OUT_PLOTS)