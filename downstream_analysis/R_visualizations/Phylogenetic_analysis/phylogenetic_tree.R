#!/usr/bin/env Rscript

##########################
# Phylogenetic tree plot
# Candida auris
##########################

library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(ggtreeExtra)

cat("=============================================\n")
cat("PHYLOGENETIC TREE VISUALIZATION\n")
cat("=============================================\n\n")

# Create output directory
dir.create("output", showWarnings = FALSE)

# ------------------------------------------------
# Load tree
# ------------------------------------------------
tree <- read.tree("input/SpeciesTree_rooted.txt")

# ------------------------------------------------
# Extract clade info and shorten labels
# ------------------------------------------------
clade_info <- data.frame(
  tip_label = tree$tip.label,
  ShortLabel = ifelse(
    grepl("NCCPF", tree$tip.label),
    paste0(
      "NCCPF_",
      substr(
        sub(".*NCCPF([0-9]+).*", "\\1", tree$tip.label),
        nchar(sub(".*NCCPF([0-9]+).*", "\\1", tree$tip.label)) - 2,
        nchar(sub(".*NCCPF([0-9]+).*", "\\1", tree$tip.label))
      )
    ),
    paste0(
      gsub("(CL\\d+).*", "\\1", tree$tip.label), "_",
      substr(
        gsub(".*([SDE]RR\\d+).*", "\\1", tree$tip.label),
        nchar(gsub(".*([SDE]RR\\d+).*", "\\1", tree$tip.label)) - 2,
        nchar(gsub(".*([SDE]RR\\d+).*", "\\1", tree$tip.label))
      )
    )
  ),
  Clade = ifelse(
    grepl("NCCPF", tree$tip.label),
    "NCCPF",
    gsub("CL(\\d+)_.*", "CL\\1", tree$tip.label)
  ),
  stringsAsFactors = FALSE
)

# Update tree tip labels
tree$tip.label <- clade_info$ShortLabel

# ------------------------------------------------
# Define clade colors and labels
# ------------------------------------------------
clade_colors <- c(
  "CL1" = "#3366CC",
  "CL2" = "#33A02C",
  "CL3" = "#E31A1C",
  "CL4" = "#FF7F00",
  "CL5" = "#6A3D9A",
  "CL6" = "#8B4513",
  "NCCPF" = "#FFD700"
)

clade_labels <- c(
  "CL1" = "Clade I - South Asia",
  "CL2" = "Clade II - East Asia",
  "CL3" = "Clade III - South Africa",
  "CL4" = "Clade IV - South America",
  "CL5" = "Clade V - Iran",
  "CL6" = "Clade VI - Singapore",
  "NCCPF" = "NCCPF - Reference"
)

# ------------------------------------------------
# Compatibility patch (ggplot2 ≥ 4.0)
# ------------------------------------------------
if (!exists("is.waive")) {
  is.waive <- function(x) inherits(x, "waiver")
}

# ------------------------------------------------
# Plot tree
# ------------------------------------------------
p <- ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tiplab(
    size = 2.8,
    align = TRUE,
    offset = 2.5,
    hjust = -0.1
  ) +
  geom_fruit(
    data = clade_info,
    geom = geom_tile,
    mapping = aes(y = ShortLabel, fill = Clade),
    width = 2.0,
    offset = 0.5
  ) +
  scale_fill_manual(
    values = clade_colors,
    labels = clade_labels,
    na.value = "gray"
  ) +
  theme_tree2() +
  labs(fill = "Clades") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(30, 30, 30, 30)
  ) +
  xlim(-2, NA)

# ------------------------------------------------
# Save plot
# ------------------------------------------------
ggsave(
  "output/phylogenetic_tree.png",
  plot = p,
  width = 14,
  height = 14,
  units = "in",
  dpi = 900,
  bg = "white"
)

cat("Tree saved to output/phylogenetic_tree.png\n")
cat("=============================================\n")
cat("PHYLOGENETIC ANALYSIS COMPLETE\n")
cat("=============================================\n")
