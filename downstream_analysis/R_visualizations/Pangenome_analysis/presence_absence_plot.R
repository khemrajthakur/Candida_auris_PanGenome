#!/usr/bin/env Rscript

############################
# Presence–absence analysis
# Candida auris pan-genome
############################

library(tidyverse)
library(patchwork)
library(ggtext)

cat("=============================================\n")
cat("PRESENCE–ABSENCE ANALYSIS\n")
cat("=============================================\n\n")

# Create output directory
dir.create("output", showWarnings = FALSE)

# -------------------------------------------------
# Step 1: Read Orthogroups.tsv (OrthoFinder output)
# -------------------------------------------------
orthogroups <- read_tsv("input/Orthogroups.tsv")

# -------------------------------------------------
# Step 2: Reshape to long format
# -------------------------------------------------
ortholong <- orthogroups %>%
  pivot_longer(-Orthogroup, names_to = "Genome", values_to = "Gene") %>%
  mutate(Present = ifelse(is.na(Gene), 0, 1))

# -------------------------------------------------
# Step 3: Extract clade (CL1–CL6)
# -------------------------------------------------
ortholong <- ortholong %>%
  mutate(Clade = str_extract(Genome, "^CL[0-9]+")) %>%
  filter(!is.na(Clade))

# -------------------------------------------------
# Step 4: Presence summary per clade
# -------------------------------------------------
presence_summary <- ortholong %>%
  group_by(Orthogroup, Clade) %>%
  summarise(
    Presence = sum(Present),
    Total = n(),
    Percent = 100 * Presence / Total,
    .groups = "drop"
  )

# -------------------------------------------------
# Step 5: Color settings
# -------------------------------------------------
clade_order <- paste0("CL", 1:6)
clade_colors <- c(
  "CL1" = "#4C72B0", "CL2" = "#55A868", "CL3" = "#C44E52",
  "CL4" = "#8172B2", "CL5" = "#CCB974", "CL6" = "#64B5CD"
)
presence_summary$Clade <- factor(presence_summary$Clade, levels = clade_order)

# =================================================
# PANEL A: Clade-specific gene family gains
# =================================================
specific_gains <- presence_summary %>%
  group_by(Orthogroup) %>%
  mutate(num_clades_present = sum(Percent > 0)) %>%
  filter(num_clades_present == 1 & Percent > 0) %>%
  ungroup()

# Export table
write_tsv(specific_gains, "output/Clade_Specific_Gains.tsv")

pA <- ggplot(specific_gains, aes(x = Clade, y = Percent, color = Clade)) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.9) +
  scale_color_manual(values = clade_colors) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    title = "A) Clade-specific gene family gains",
    x = "Clade", y = "% of genomes per clade"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0)
  )

# =================================================
# PANEL B: Clade-specific gene family absences
# =================================================
specific_absences <- presence_summary %>%
  group_by(Orthogroup) %>%
  mutate(num_clades_absent = sum(Percent == 0)) %>%
  filter(num_clades_absent == 1) %>%
  mutate(Absent_in = Clade[Percent == 0][1]) %>%
  ungroup() %>%
  filter(!is.na(Absent_in)) %>%
  mutate(Absent_in = factor(Absent_in, levels = clade_order))

# Export table
write_tsv(specific_absences, "output/Clade_Specific_Losses.tsv")

facet_labels <- setNames(
  paste0("Absent in ", clade_order),
  clade_order
)

pB <- ggplot(specific_absences, aes(x = Clade, y = Percent, group = Orthogroup)) +
  geom_line(color = "grey80", linewidth = 0.3, alpha = 0.8) +
  geom_point(aes(color = Clade), size = 1.1) +
  scale_color_manual(values = clade_colors) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  facet_wrap(~Absent_in, ncol = 3, labeller = labeller(Absent_in = facet_labels)) +
  labs(
    title = "B) Clade-specific gene family absences",
    x = NULL, y = "% of genomes in a clade"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(face = "bold", hjust = 0)
  )

# =================================================
# Combine panels
# =================================================
final_plot <- pA + pB + plot_layout(ncol = 2, widths = c(1, 1.8))

ggsave(
  "output/Figure_Presence_Absence_CL1_CL6.png",
  final_plot,
  width = 14,
  height = 6.5,
  dpi = 500,
  bg = "white"
)

cat("Plot saved to output/Figure_Presence_Absence_CL1_CL6.png\n")
cat("=============================================\n")
cat("PRESENCE–ABSENCE ANALYSIS COMPLETE\n")
cat("=============================================\n")
