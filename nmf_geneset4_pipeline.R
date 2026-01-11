suppressPackageStartupMessages({
  library(Seurat)
  library(NMF)
  library(ggplot2)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggridges)
})

# ---- Config ----
correct_genesets_path <- "E:/R/Spinal core/GSE306676/graph/GeneSets/Correct_GeneSets_K15.rds"
sobj_subset_path <- "E:/R/Spinal core/GSE306676/sobj_subset.rds"
sobj_path <- "E:/R/Spinal core/GSE306676/GSE306676_all_integrated_annotated_final.rds"
output_dir <- "nmf_geneset4_outputs"

group_column <- "group" # update to metadata column that encodes ctl vs sod1
subset_stages <- c("ctl", "sod.early", "sod.end", "sod.mid")
k_range <- 2:10
nmf_runs <- 10
seed <- 42
max_dense_gb <- 4
max_cells <- 5000
max_genes <- 2000

# ---- Helpers ----
ensure_nonnegative <- function(mat) {
  if (min(mat) < 0) {
    warning("Negative values detected; shifting matrix to be non-negative.")
    mat <- mat - min(mat)
  }
  mat
}

select_elbow_k <- function(k_values, errors) {
  if (length(k_values) < 3) {
    return(k_values[1])
  }
  second_diff <- diff(errors, differences = 2)
  elbow_index <- which.max(abs(second_diff)) + 1
  k_values[elbow_index]
}

estimate_dense_gb <- function(n_rows, n_cols) {
  n_rows * n_cols * 8 / 1024^3
}

downsample_matrix <- function(mat, max_genes, max_cells, seed) {
  set.seed(seed)
  if (!is.null(max_genes) && nrow(mat) > max_genes) {
    mat <- mat[sample(rownames(mat), max_genes), , drop = FALSE]
  }
  if (!is.null(max_cells) && ncol(mat) > max_cells) {
    mat <- mat[, sample(colnames(mat), max_cells), drop = FALSE]
  }
  mat
}

# ---- Load inputs ----
set.seed(seed)

if (!file.exists(correct_genesets_path)) {
  stop(
    paste0(
      "Cannot find Correct_GeneSets_K15.rds at: ",
      normalizePath(correct_genesets_path, winslash = "/", mustWork = FALSE),
      ". Update correct_genesets_path to the directory where you saved it (e.g. ",
      "geneset_path/Correct_GeneSets_K15.rds)."
    ),
    call. = FALSE
  )
}

genesets <- readRDS(correct_genesets_path)
if (!"GeneSet4" %in% names(genesets)) {
  stop("GeneSet4 not found in Correct_GeneSets_K15.rds")
}

if (file.exists(sobj_subset_path)) {
  sobj_subset <- readRDS(sobj_subset_path)
} else if (file.exists(sobj_path)) {
  sobj <- readRDS(sobj_path)
  selected <- sobj$stage %in% subset_stages
  sobj_subset <- sobj[, selected]
} else {
  stop(
    paste0(
      "Cannot find sobj_subset.rds at: ",
      normalizePath(sobj_subset_path, winslash = "/", mustWork = FALSE),
      " and sobj_path is missing: ",
      normalizePath(sobj_path, winslash = "/", mustWork = FALSE),
      ". Update sobj_subset_path or sobj_path to a valid location."
    ),
    call. = FALSE
  )
}

# ---- Step 1-2: GeneSet4 submatrix ----
GeneSet4_genes <- unique(genesets$GeneSet4)
expr_data <- GetAssayData(sobj_subset, assay = "RNA", layer = "data")
expr_data <- expr_data[intersect(rownames(expr_data), GeneSet4_genes), , drop = FALSE]
expr_data <- downsample_matrix(expr_data, max_genes = max_genes, max_cells = max_cells, seed = seed)

if (anyNA(expr_data) || any(!is.finite(expr_data))) {
  stop("GeneSet4 matrix contains NA/Inf values; please clean the data before NMF.", call. = FALSE)
}

zero_gene <- rowSums(expr_data) == 0
if (any(zero_gene)) {
  expr_data <- expr_data[!zero_gene, , drop = FALSE]
}

zero_cell <- colSums(expr_data) == 0
if (any(zero_cell)) {
  expr_data <- expr_data[, !zero_cell, drop = FALSE]
}

if (nrow(expr_data) == 0 || ncol(expr_data) == 0) {
  stop("GeneSet4 matrix is empty after removing zero rows/columns.", call. = FALSE)
}

max_rank <- min(nrow(expr_data), ncol(expr_data))
if (max(k_range) > max_rank) {
  k_range <- k_range[k_range <= max_rank]
  if (length(k_range) == 0) {
    stop(
      paste0(
        "No valid k in k_range after filtering. max_rank=",
        max_rank,
        "; adjust k_range or increase max_cells/max_genes."
      ),
      call. = FALSE
    )
  }
}

estimated_gb <- estimate_dense_gb(nrow(expr_data), ncol(expr_data))
if (estimated_gb > max_dense_gb) {
  stop(
    paste0(
      "GeneSet4 matrix would require ~",
      round(estimated_gb, 2),
      " GB as a dense matrix. Reduce max_cells/max_genes or set a larger ",
      "max_dense_gb limit before running NMF."
    ),
    call. = FALSE
  )
}

expr_data <- as.matrix(expr_data)
expr_data <- ensure_nonnegative(expr_data)

if (nrow(expr_data) == 0) {
  stop("No overlapping genes between GeneSet4 and sobj_subset RNA data.")
}

# ---- Step 3: NMF across k ----
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

nmf_results <- vector("list", length(k_range))
recon_errors <- numeric(length(k_range))

for (i in seq_along(k_range)) {
  k <- k_range[i]
  fit <- nmf(expr_data, rank = k, method = "brunet", nrun = nmf_runs, seed = seed)
  nmf_results[[i]] <- fit
  recon_errors[i] <- rss(fit)
}

nmf_summary <- tibble(k = k_range, rss = recon_errors)
write.csv(nmf_summary, file = file.path(output_dir, "nmf_rss_summary.csv"), row.names = FALSE)

# ---- Step 4: Elbow curve ----
elbow_k <- select_elbow_k(k_range, recon_errors)

elbow_plot <- ggplot(nmf_summary, aes(x = k, y = rss)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = elbow_k, linetype = "dashed", color = "red") +
  labs(title = "NMF elbow curve", subtitle = paste("Elbow k:", elbow_k))

ggsave(file.path(output_dir, "nmf_elbow_curve.png"), elbow_plot, width = 6, height = 4)

# ---- Step 5: Submodule scores ----
fit_best <- nmf_results[[which(k_range == elbow_k)]]

W <- basis(fit_best)
H <- coef(fit_best)

write.csv(W, file = file.path(output_dir, "nmf_basis_W.csv"))
write.csv(H, file = file.path(output_dir, "nmf_coef_H.csv"))

# Define submodule gene sets from W loadings (top genes per module)
module_gene_sets <- apply(W, 2, function(weights) {
  names(sort(weights, decreasing = TRUE))[1:min(50, length(weights))]
})

sobj_subset <- AddModuleScore(
  sobj_subset,
  features = as.list(module_gene_sets),
  name = "NMF_Module"
)

module_score_cols <- grep("^NMF_Module", colnames(sobj_subset@meta.data), value = TRUE)
module_scores <- sobj_subset@meta.data[, module_score_cols, drop = FALSE]

write.csv(module_scores, file = file.path(output_dir, "nmf_module_scores.csv"))

# ---- Step 6: Visualization ----
score_matrix <- as.matrix(module_scores)

score_pca <- prcomp(score_matrix, center = TRUE, scale. = TRUE)
score_pca_df <- as.data.frame(score_pca$x[, 1:2])
score_pca_df[[group_column]] <- sobj_subset@meta.data[[group_column]]

pca_plot <- ggplot(score_pca_df, aes(x = PC1, y = PC2, color = .data[[group_column]])) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA on NMF module scores")

ggsave(file.path(output_dir, "nmf_module_score_pca.png"), pca_plot, width = 6, height = 4)

ridge_data <- module_scores %>%
  mutate(group = sobj_subset@meta.data[[group_column]]) %>%
  tidyr::pivot_longer(cols = all_of(module_score_cols), names_to = "module", values_to = "score")

ridge_plot <- ggplot(ridge_data, aes(x = score, y = module, fill = group)) +
  ggridges::geom_density_ridges(alpha = 0.6) +
  labs(title = "Module score distributions")

ggsave(file.path(output_dir, "nmf_module_score_ridges.png"), ridge_plot, width = 8, height = 5)

# ---- Step 7 (optional): Simple classifier ----
if (!is.null(sobj_subset@meta.data[[group_column]])) {
  classifier_data <- as.data.frame(score_matrix)
  classifier_data$group <- as.factor(sobj_subset@meta.data[[group_column]])

  set.seed(seed)
  fit_glm <- glm(group ~ ., data = classifier_data, family = binomial())
  saveRDS(fit_glm, file = file.path(output_dir, "nmf_module_glm.rds"))
}

message("Pipeline completed. Outputs saved to: ", output_dir)
