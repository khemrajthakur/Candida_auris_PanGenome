################## PAN-GENOME PIE / DONUT CHART ##################

# Load required libraries
library(ggplot2)
library(dplyr)

cat("=============================================\n")
cat("PAN-GENOME COMPOSITION ANALYSIS\n")
cat("=============================================\n\n")

# Create output directory
dir.create("output", showWarnings = FALSE)

# ============================================================================
# READ AND PROCESS DATA
# ============================================================================

cat("Reading OrthoFinder gene count file...\n")

# Read Orthogroups.GeneCount.tsv (from OrthoFinder)
gene_counts <- read.delim(
  "input/Orthogroups.GeneCount.tsv",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

# Remove the 'Total' column if present
if ("Total" %in% colnames(gene_counts)) {
  gene_counts <- gene_counts[, -which(colnames(gene_counts) == "Total")]
}

# Number of genomes
n_genomes <- ncol(gene_counts)

# Presence / absence matrix
presence_absence <- gene_counts
presence_absence[presence_absence > 0] <- 1

# Number of genomes per gene family
genomes_per_family <- rowSums(presence_absence)

# ============================================================================
# CLASSIFY GENE FAMILIES
# ============================================================================

# Core genes: present in ≥95% genomes
core_cutoff <- ceiling(0.95 * n_genomes)

core_genes <- genomes_per_family >= core_cutoff
accessory_genes <- genomes_per_family > 1 & genomes_per_family < core_cutoff
singleton_genes <- genomes_per_family == 1

# Counts
n_core <- sum(core_genes)
n_accessory <- sum(accessory_genes)
n_singleton <- sum(singleton_genes)
n_total <- nrow(gene_counts)

# Print statistics
cat("=== PAN-GENOME STATISTICS ===\n")
cat("Total gene families:", format(n_total, big.mark = ","), "\n")
cat("Core genes:", format(n_core, big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * n_core / n_total))
cat("Accessory genes:", format(n_accessory, big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * n_accessory / n_total))
cat("Singleton genes:", format(n_singleton, big.mark = ","),
    sprintf("(%.1f%%)\n\n", 100 * n_singleton / n_total))

# ============================================================================
# PREPARE DATA FOR DONUT CHART
# ============================================================================

pie_data <- data.frame(
  category = factor(
    c("Core", "Accessory", "Singleton"),
    levels = c("Core", "Accessory", "Singleton")
  ),
  count = c(n_core, n_accessory, n_singleton),
  percentage = c(
    n_core / n_total * 100,
    n_accessory / n_total * 100,
    n_singleton / n_total * 100
  )
)

# Calculate positions
pie_data$fraction <- pie_data$count / sum(pie_data$count)
pie_data$ymax <- cumsum(pie_data$fraction)
pie_data$ymin <- c(0, head(pie_data$ymax, -1))
pie_data$labelPosition <- (pie_data$ymax + pie_data$ymin) / 2

# Labels
pie_data$label <- paste0(
  pie_data$category, "\n",
  format(pie_data$count, big.mark = ","), "\n",
  sprintf("(%.1f%%)", pie_data$percentage)
)

# Label positions
pie_data$label_x <- c(5.5, 5.5, 5.5)
pie_data$line_end_x <- c(5.0, 5.0, 5.0)

# Adjust spacing
pie_data$labelPosition[2] <- pie_data$labelPosition[2] + 0.05

# Colors
colors <- c(
  "Core" = "#7A5C33",
  "Accessory" = "#D4A574",
  "Singleton" = "#5F9EA0"
)

# ============================================================================
# CREATE DONUT CHART
# ============================================================================

cat("Creating pan-genome donut plot...\n")

p <- ggplot(
  pie_data,
  aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = category)
) +
  geom_rect(color = "white", size = 2) +
  scale_fill_manual(values = colors) +
  coord_polar(theta = "y") +
  xlim(c(0.5, 6.0)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5, size = 28, face = "bold",
      margin = margin(b = 5, t = 10)
    ),
    plot.subtitle = element_text(
      hjust = 0.5, size = 16,
      margin = margin(b = 25, t = 5)
    ),
    plot.background = element_rect(
      fill = "white", color = "grey30", size = 2
    )
  ) +
  ggtitle("Pan-genome Composition") +
  labs(
    subtitle = paste0(
      "Total Gene Families: ",
      format(n_total, big.mark = ",")
    )
  )

# Add straight connector lines
for (i in 1:nrow(pie_data)) {
  p <- p + geom_segment(
    x = 4, xend = pie_data$line_end_x[i],
    y = pie_data$labelPosition[i],
    yend = pie_data$labelPosition[i],
    color = "grey30", size = 0.8
  )
}

# Add labels
p <- p +
  geom_text(
    data = pie_data[1, ],
    aes(x = label_x, y = labelPosition, label = label),
    size = 6.5, fontface = "bold", hjust = 0
  ) +
  geom_text(
    data = pie_data[c(2, 3), ],
    aes(x = label_x, y = labelPosition, label = label),
    size = 6.5, fontface = "bold", hjust = 1
  )

# ============================================================================
# SAVE FIGURE
# ============================================================================

ggsave(
  "output/Pangenome_Composition_Final.png",
  plot = p,
  width = 13,
  height = 11,
  dpi = 300,
  bg = "white"
)

cat("Plot saved to output/Pangenome_Composition_Final.png\n")
cat("=============================================\n")
cat("PAN-GENOME ANALYSIS COMPLETE\n")
cat("=============================================\n")
