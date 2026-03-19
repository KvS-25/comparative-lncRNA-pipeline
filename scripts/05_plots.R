# ============================================================
# Script 05: Generate all visualisation plots
#
# Usage: Rscript scripts/05_plots.R
# Requires: config/config.yaml, results from 03 and 04
# Environment: goanalysis (R, ggplot2, reshape2, UpSetR)
# ============================================================

library(yaml)
library(ggplot2)
library(reshape2)
library(UpSetR)

# ---- Load config ----
config <- yaml.load_file("config/config.yaml")

OUT_BASE <- config$output$base
OUT_GO   <- config$output$go
OUT_KEGG <- config$output$kegg
OUT_PLOTS <- config$output$plots

dir.create(OUT_PLOTS, recursive=TRUE, showWarnings=FALSE)

categories <- c("needle_specific", "root_specific",
                 "cold_specific", "drought_specific")

category_colors <- c(
  "needle_specific" = "#C6EFCE",
  "root_specific"   = "#FFEB9C",
  "cold_specific"   = "#BDD7EE",
  "drought_specific"= "#FCE4D6"
)

# ============================================================
# Plot 1: UpSet plot
# ============================================================
cat("[", format(Sys.time()), "] Generating UpSet plot\n")

multiinter <- read.table(file.path(OUT_BASE, "multiinter_output.bed"),
                          header=FALSE, sep="\t")
colnames(multiinter)[6:9] <- config$samples$names

png(file.path(OUT_PLOTS, "upset_plot.png"),
    width=1400, height=900, res=150, pointsize=12)
par(mar=c(6, 6, 4, 2))
upset(multiinter[,6:9],
      sets         = config$samples$names,
      mb.ratio     = c(0.65, 0.35),
      order.by     = "freq",
      main.bar.color = "steelblue",
      text.scale   = c(1.5, 1.2, 1.5, 1.3, 1.5, 0.9),
      point.size   = 3,
      line.size    = 1,
      shade.color  = "lightblue",
      shade.alpha  = 0.3,
      matrix.color = "gray30",
      mainbar.y.label = "Number of Genomic Regions",
      sets.x.label    = "Total Regions per Sample")
dev.off()

# ============================================================
# Plot 2: Region counts bar chart
# ============================================================
cat("[", format(Sys.time()), "] Generating region counts bar chart\n")

region_files <- c(
  "Conserved"        = "conserved.bed",
  "Needle-specific"  = "needle_specific.bed",
  "Root-specific"    = "root_specific.bed",
  "Cold-specific"    = "cold_specific.bed",
  "Drought-specific" = "drought_specific.bed"
)

counts_data <- data.frame(
  Category = names(region_files),
  Count    = sapply(region_files, function(f) {
    fp <- file.path(OUT_BASE, f)
    if (file.exists(fp)) as.integer(system(paste("wc -l <", fp), intern=TRUE)) else 0
  }),
  Group = c("Conserved", "Tissue", "Tissue", "Condition", "Condition")
)

counts_data$Category <- factor(counts_data$Category, levels=rev(counts_data$Category))

p <- ggplot(counts_data, aes(x=Count, y=Category, fill=Group)) +
  geom_bar(stat="identity", width=0.7) +
  geom_text(aes(label=format(Count, big.mark=",")), hjust=-0.1, size=3.5) +
  scale_fill_manual(values=c(
    "Conserved" = "#525252",
    "Tissue"    = "#2CA25F",
    "Condition" = "#2B8CBE"
  )) +
  scale_x_continuous(
    labels=function(x) format(x, big.mark=",", scientific=FALSE),
    expand=expansion(mult=c(0, 0.15))) +
  labs(title="Genomic Regions by Category",
       x="Number of Regions", y="", fill="Category Type") +
  theme_bw() +
  theme(axis.text.y=element_text(size=11),
        plot.title=element_text(hjust=0.5, face="bold", size=13))

ggsave(file.path(OUT_PLOTS, "region_counts_bar.png"),
       plot=p, width=10, height=7, dpi=300)

# ============================================================
# Plot 3: GO bar plots per category
# ============================================================
cat("[", format(Sys.time()), "] Generating GO bar plots\n")

make_go_bar <- function(category, title, filename) {
  bp <- read.table(file.path(OUT_GO, paste0(category, "_BP_GO_enrichment.txt")),
                   header=TRUE, sep="\t")
  mf <- read.table(file.path(OUT_GO, paste0(category, "_MF_GO_enrichment.txt")),
                   header=TRUE, sep="\t")
  cc <- read.table(file.path(OUT_GO, paste0(category, "_CC_GO_enrichment.txt")),
                   header=TRUE, sep="\t")

  bp <- head(bp[order(bp$weight01), ], 5)
  mf <- head(mf[order(mf$weight01), ], 5)
  cc <- head(cc[order(cc$weight01), ], 5)

  bp$Domain <- "Biological Process"
  mf$Domain <- "Molecular Function"
  cc$Domain <- "Cellular Component"

  combined <- rbind(
    data.frame(Term=bp$Term, pvalue=bp$weight01, Count=bp$Significant, Domain=bp$Domain),
    data.frame(Term=mf$Term, pvalue=mf$weight01, Count=mf$Significant, Domain=mf$Domain),
    data.frame(Term=cc$Term, pvalue=cc$weight01, Count=cc$Significant, Domain=cc$Domain)
  )

  combined$Term <- factor(combined$Term,
    levels=combined$Term[order(combined$Domain, combined$pvalue, decreasing=TRUE)])
  combined$Domain <- factor(combined$Domain,
    levels=c("Biological Process", "Molecular Function", "Cellular Component"))

  p <- ggplot(combined, aes(x=-log10(pvalue), y=Term, fill=Domain)) +
    geom_bar(stat="identity", width=0.7) +
    geom_text(aes(label=Count), hjust=-0.2, size=3.5) +
    scale_fill_manual(values=c(
      "Biological Process" = "#4575B4",
      "Molecular Function" = "#D73027",
      "Cellular Component" = "#1A9850"
    )) +
    facet_wrap(~Domain, scales="free_y", ncol=1) +
    labs(title=title, x="-log10(p-value)", y="") +
    theme_bw() +
    theme(strip.background=element_rect(fill="#F0F0F0"),
          strip.text=element_text(face="bold", size=11),
          axis.text.y=element_text(size=10),
          plot.title=element_text(hjust=0.5, face="bold", size=13),
          legend.position="none") +
    xlim(0, max(-log10(combined$pvalue)) * 1.15)

  ggsave(filename, plot=p, width=10, height=12, dpi=300)
}

make_go_bar("needle_specific", "GO Enrichment - Needle Specific",
            file.path(OUT_PLOTS, "GO_bar_needle.png"))
make_go_bar("root_specific",   "GO Enrichment - Root Specific",
            file.path(OUT_PLOTS, "GO_bar_root.png"))
make_go_bar("cold_specific",   "GO Enrichment - Cold Specific",
            file.path(OUT_PLOTS, "GO_bar_cold.png"))
make_go_bar("drought_specific","GO Enrichment - Drought Specific",
            file.path(OUT_PLOTS, "GO_bar_drought.png"))

cat("[", format(Sys.time()), "] All plots saved to", OUT_PLOTS, "\n")