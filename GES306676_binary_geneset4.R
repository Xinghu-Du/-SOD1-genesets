library(dplyr)
library(ggplot2)
library(Seurat)
library(AUCell)
library(ggridges)
rm(list = ls())
gc()

options(future.globals.maxSize = 120 * 1024^3)
setwd("E:/R/Spinal core/GSE306676")

sobj <- readRDS("E:/R/Spinal core/GSE306676/GSE306676_all_integrated_annotated_final.rds")

save_path <- "E:/R/Spinal core/GSE306676/graph"
class_base_path <- file.path(save_path, "Classification")
geneset_path <- file.path(save_path, "GeneSets")

output_dir <- file.path(class_base_path, "Binary_ctl_vs_sod_GeneSet4")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot, path_base, width, height) {
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(path_base, ".tiff"), plot, width = width, height = height, dpi = 300, device = "tiff")
}

selected <- sobj$stage %in% c("ctl", "sod.early", "sod.end", "sod.mid")
sobj_subset <- sobj[, selected]
expr_matrix <- GetAssayData(sobj_subset, assay = "RNA", layer = "data")

genesets <- readRDS(file.path(geneset_path, "Correct_GeneSets_K15.rds"))
geneset4 <- genesets[["GeneSet4"]]

rankings <- AUCell_buildRankings(expr_matrix)
auc_score <- AUCell_calcAUC(list(GeneSet4 = geneset4), rankings, verbose = TRUE)
auc_vec <- as.numeric(getAUC(auc_score)[1, ])

auc_df <- data.frame(
  score = auc_vec,
  type_binary = ifelse(sobj_subset$stage == "ctl", "ctl", "sod")
)
auc_df$type_binary <- factor(auc_df$type_binary, levels = c("ctl", "sod"))

p_box <- ggplot(auc_df, aes(x = type_binary, y = score, fill = type_binary)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.4) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "GeneSet4 AUCell Score (ctl vs sod)", x = "Group", y = "AUCell Score")
save_plot(p_box, file.path(output_dir, "geneset4_boxplot"), width = 5, height = 5)

p_ridge <- ggplot(auc_df, aes(x = score, y = type_binary, fill = type_binary)) +
  geom_density_ridges(alpha = 0.6) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "GeneSet4 AUCell Score Distribution", x = "AUCell Score", y = "Group")
save_plot(p_ridge, file.path(output_dir, "geneset4_ridgeline"), width = 6, height = 3.5)

p_violin <- ggplot(auc_df, aes(x = type_binary, y = score, fill = type_binary)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "GeneSet4 AUCell Score (Violin)", x = "Group", y = "AUCell Score")
save_plot(p_violin, file.path(output_dir, "geneset4_violin"), width = 5, height = 5)

cat("GeneSet4 binary analysis complete. Outputs saved to:", output_dir, "\n")
