#!/usr/bin/env Rscript

# ============================================================
# DAPC ANALYSIS FOR CANDIDA AURIS
# ============================================================

library(adegenet)
library(vcfR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(RColorBrewer)

cat("=============================================================\n")
cat("DAPC ANALYSIS FOR CANDIDA AURIS\n")
cat("=============================================================\n\n")

# Output directory (GitHub-safe)
dir.create("output/dapc_results", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1: LOAD DATA
# ============================================================

cat("Step 1: Loading VCF and population data...\n")

vcf <- read.vcfR("input/all_clades_filtered.vcf.gz", verbose = FALSE)
gl <- vcfR2genlight(vcf)

pop_map <- read.table("input/population_map.txt",
                      header = FALSE,
                      col.names = c("Sample", "Clade"))

pop(gl) <- pop_map$Clade[match(indNames(gl), pop_map$Sample)]

cat("  ✓ Loaded:", nInd(gl), "samples,", nLoc(gl), "SNPs\n")
cat("  ✓ Clades:", paste(names(table(pop(gl))), collapse = ", "), "\n\n")

# ============================================================
# STEP 2: CROSS-VALIDATION FOR OPTIMAL PCS
# ============================================================

cat("Step 2: Finding optimal number of PCs...\n")

set.seed(123)
xval <- xvalDapc(tab(gl, NA.method = "mean"), pop(gl),
                 n.pca.max = 100,
                 training.set = 0.9,
                 result = "groupMean",
                 n.rep = 30,
                 xval.plot = FALSE)

optimal_pcs <- as.integer(xval$`Number of PCs Achieving Lowest MSE`)

if (is.na(optimal_pcs)) {
  optimal_pcs <- which.max(xval$`Mean Successful Assignment by Number of PCs`)
}

cat("  ✓ Optimal PCs:", optimal_pcs, "\n\n")

# ============================================================
# STEP 3: RUN DAPC
# ============================================================

cat("Step 3: Running DAPC...\n")

n_groups <- length(unique(pop(gl)))
n_da <- min(n_groups - 1, 5)

dapc_result <- dapc(gl, pop(gl),
                    n.pca = optimal_pcs,
                    n.da = n_da)

var_explained <- dapc_result$eig / sum(dapc_result$eig) * 100

cat(sprintf("  ✓ Variance explained (DA1+DA2): %.2f%%\n\n",
            sum(var_explained[1:2])))

# ============================================================
# STEP 4: ASSIGNMENT ACCURACY
# ============================================================

assignments <- apply(dapc_result$posterior, 1, which.max)
original <- as.numeric(factor(pop(gl)))

accuracy <- mean(assignments == original) * 100

cat(sprintf("  ✓ Overall accuracy: %.2f%%\n\n", accuracy))

# ============================================================
# STEP 5: FIND OPTIMAL K
# ============================================================

grp <- find.clusters(gl,
                     max.n.clust = 15,
                     n.pca = optimal_pcs,
                     choose.n.clust = FALSE,
                     criterion = "diffNgroup")

optimal_k <- length(unique(grp$grp))

cat(sprintf("  ✓ Optimal K: %d\n\n", optimal_k))

comparison <- table(Clade = pop(gl), Cluster = grp$grp)

# ============================================================
# STEP 6: SUMMARY TABLES
# ============================================================

summary_df <- data.frame(
  Sample = indNames(gl),
  Original = pop(gl),
  Assigned = levels(pop(gl))[assignments],
  Probability = apply(dapc_result$posterior, 1, max)
)

misassigned <- summary_df %>% filter(Original != Assigned)
ambiguous <- summary_df %>% filter(Probability < 0.8)

# ============================================================
# STEP 7: PLOTS
# ============================================================

clade_colors <- c(
  "CL1" = "#542788", "CL2" = "#D6604D", "CL3" = "#B8BA3C",
  "CL4" = "#F4A582", "CL5" = "#8C510A", "CL6" = "#4393C3"
)

dapc_coords <- as.data.frame(dapc_result$ind.coord)
dapc_coords$Clade <- pop(gl)

p_scatter <- ggplot(dapc_coords, aes(LD1, LD2, color = Clade)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = clade_colors) +
  theme_classic(base_size = 14)

ggsave("output/dapc_results/DAPC_Scatter.png",
       p_scatter, width = 8, height = 6, dpi = 600)

# ============================================================
# STEP 8: SAVE RESULTS
# ============================================================

write.csv(dapc_result$posterior,
          "output/dapc_results/posterior_probabilities.csv")
write.csv(dapc_result$ind.coord,
          "output/dapc_results/discriminant_scores.csv")
write.csv(comparison,
          "output/dapc_results/clade_vs_cluster.csv")

if (nrow(misassigned) > 0) {
  write.csv(misassigned,
            "output/dapc_results/misassigned_samples.csv",
            row.names = FALSE)
}

if (nrow(ambiguous) > 0) {
  write.csv(ambiguous,
            "output/dapc_results/ambiguous_samples.csv",
            row.names = FALSE)
}

saveRDS(list(dapc = dapc_result, grp = grp, gl = gl),
        "output/dapc_results/dapc_all_objects.rds")

cat("=============================================================\n")
cat("DAPC ANALYSIS COMPLETE\n")
cat("Results saved in output/dapc_results/\n")
cat("=============================================================\n")
