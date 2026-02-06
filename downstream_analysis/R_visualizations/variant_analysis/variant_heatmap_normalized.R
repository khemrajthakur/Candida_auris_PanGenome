#!/usr/bin/env Rscript

############################
# Variant heatmap (normalized)
# Candida auris
############################

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(grid)

cat("=============================================\n")
cat("NORMALIZED VARIANT HEATMAP\n")
cat("=============================================\n\n")

# Create output directory
dir.create("output", showWarnings = FALSE)

# -------------------------------------------------
# Target genes
# -------------------------------------------------
target_genes <- c(
  "ACC1", "IRA2", "NOC2", "SWI4",
  "HSF1", "HSP104", "HOG1", "PBS2", "SOD2",
  "FKS1", "CDR1", "ERG2", "YOR1",
  "PCK1", "ILV1",
  "CLN3", "CDC7",
  "FET3", "SMF1"
)

# -------------------------------------------------
# Clade sample sizes (used for normalization)
# -------------------------------------------------
clade_sample_sizes <- c(
  "CL1" = 50,
  "CL2" = 50,
  "CL3" = 50,
  "CL4" = 50,
  "CL5" = 14,
  "CL6" = 9,
  "NCCPF-470200" = 1
)

clades <- names(clade_sample_sizes)

# -------------------------------------------------
# Read variant count files
# -------------------------------------------------
all_data <- list()

for (clade in clades) {
  file_name <- paste0(
    "input/vcf_gene_symbols_with_variant_counts_",
    clade,
    "_filtered_target_genes.csv"
  )
  
  if (file.exists(file_name)) {
    df <- read.csv(file_name, stringsAsFactors = FALSE)
    df$Clade <- clade
    all_data[[clade]] <- df
  }
}

combined_data <- bind_rows(all_data)

# -------------------------------------------------
# Extract gene names
# -------------------------------------------------
combined_data <- combined_data %>%
  mutate(Gene_Name = case_when(
    str_detect(Product, "ACC1") ~ "ACC1",
    str_detect(Product, "IRA2") ~ "IRA2",
    str_detect(Product, "NOC2") ~ "NOC2",
    str_detect(Product, "SWI4") ~ "SWI4",
    str_detect(Product, "HSF1") ~ "HSF1",
    str_detect(Product, "HSP104") ~ "HSP104",
    str_detect(Product, "HOG1") ~ "HOG1",
    str_detect(Product, "PBS2") ~ "PBS2",
    str_detect(Product, "SOD2") ~ "SOD2",
    str_detect(Product, "FKS1") ~ "FKS1",
    str_detect(Product, "CDR1") ~ "CDR1",
    str_detect(Product, "ERG2") ~ "ERG2",
    str_detect(Product, "YOR1") ~ "YOR1",
    str_detect(Product, "PCK1") ~ "PCK1",
    str_detect(Product, "ILV1") ~ "ILV1",
    str_detect(Product, "CLN3") ~ "CLN3",
    str_detect(Product, "CDC7") ~ "CDC7",
    str_detect(Product, "FET3") ~ "FET3",
    str_detect(Product, "SMF1") ~ "SMF1",
    TRUE ~ NA_character_
  )) %>%
  filter(Gene_Name %in% target_genes)

# -------------------------------------------------
# Create normalized matrix
# -------------------------------------------------
heatmap_matrix <- combined_data %>%
  select(Gene_Name, Clade, Variant_Count) %>%
  pivot_wider(names_from = Clade, values_from = Variant_Count, values_fill = 0) %>%
  column_to_rownames("Gene_Name")

for (clade in clades) {
  if (clade %in% colnames(heatmap_matrix)) {
    heatmap_matrix[, clade] <-
      heatmap_matrix[, clade] / clade_sample_sizes[clade]
  }
}

# Order columns
colnames(heatmap_matrix) <- c(
  "CLADE I", "CLADE II", "CLADE III",
  "CLADE IV", "CLADE V", "CLADE VI",
  "NCCPF_470200"
)

heatmap_matrix <- heatmap_matrix[target_genes, ]

# Log10 transform
heatmap_log10 <- log10(heatmap_matrix + 1)

# -------------------------------------------------
# Plot normalized heatmap (manuscript figure)
# -------------------------------------------------
png(
  "output/variant_heatmap_normalized_log10.png",
  width = 16, height = 10, units = "in", res = 300
)

breaks_seq <- seq(0, max(heatmap_log10), length.out = 100)

pheatmap(
  heatmap_log10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("green", "black", "red"))(100),
  border_color = NA,
  fontsize = 13,
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 45,
  cellwidth = 55,
  cellheight = 30,
  breaks = breaks_seq,
  legend = TRUE
)

grid.text(
  "Log10 (variants per sample + 1)",
  x = unit(0.945, "npc"),
  y = unit(0.38, "npc"),
  gp = gpar(fontsize = 11, fontface = "bold")
)

dev.off()

# -------------------------------------------------
# Save matrix used for plotting
# -------------------------------------------------
write.csv(
  heatmap_log10,
  "output/variant_count_matrix_normalized_log10.csv",
  row.names = TRUE
)

cat("=============================================\n")
cat("Heatmap saved: output/variant_heatmap_normalized_log10.png\n")
cat("=============================================\n")
