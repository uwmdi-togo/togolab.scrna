# Trajectory inference (Slingshot) setup helpers. These cores are dataset-
# agnostic; the ATTEMPT-specific pseudotime *plots* (PRE/POST x Placebo/Dapa)
# were intentionally not ported.

#' Build a SingleCellExperiment and PCA for trajectory inference
#'
#' Converts a Seurat object to a SingleCellExperiment, filters low-count genes,
#' applies full-quantile normalization, runs PCA, and returns an elbow plot.
#'
#' @param object A Seurat object.
#' @param title Title for the elbow plot (e.g. a cell-type label).
#' @param min_count,min_cells Gene filter: keep genes with at least `min_count`
#'   counts in at least `min_cells` cells.
#' @return A list: `sce`, `pca`, `var_explained`, `elbow_plot`.
#' @export
togo_slingshot_setup <- function(object, title = "", min_count = 3, min_cells = 10) {
  .togo_need(c("Seurat", "SingleCellExperiment", "ggplot2"))

  counts <- Seurat::GetAssayData(object, layer = "counts")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts), colData = object@meta.data)

  keep <- apply(SingleCellExperiment::counts(sce), 1,
                function(x) sum(x >= min_count) >= min_cells)
  sce <- sce[keep, ]

  fq_norm <- function(counts) {
    rk <- apply(counts, 2, rank, ties.method = "min")
    sorted <- apply(counts, 2, sort)
    refdist <- apply(sorted, 1, stats::median)
    norm <- apply(rk, 2, function(r) refdist[r])
    rownames(norm) <- rownames(counts)
    norm
  }
  SummarizedExperiment::assays(sce)$norm <- fq_norm(SingleCellExperiment::counts(sce))

  pca <- stats::prcomp(t(log1p(SummarizedExperiment::assays(sce)$norm)), scale. = FALSE)
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)

  n <- min(50, length(var_explained))
  elbow <- ggplot2::ggplot(
    data.frame(PC = seq_len(n), VarianceExplained = var_explained[seq_len(n)]),
    ggplot2::aes(x = .data$PC, y = .data$VarianceExplained)) +
    ggplot2::geom_point(size = 2) + ggplot2::geom_line() +
    ggplot2::labs(title = paste("Elbow Plot:", title),
                  x = "Principal Component", y = "Proportion of Variance Explained") +
    ggplot2::theme_bw()

  list(sce = sce, pca = pca, var_explained = var_explained, elbow_plot = elbow)
}

#' Run Slingshot trajectory inference
#'
#' Adds top-PC and UMAP reduced dims to a SingleCellExperiment and runs
#' [slingshot::slingshot()] with optional start/end clusters.
#'
#' @param sce A SingleCellExperiment (e.g. `$sce` from [togo_slingshot_setup()]).
#' @param pca_obj The PCA object (`$pca`).
#' @param n_pcs Number of principal components to use.
#' @param start_cluster,end_cluster Optional starting/ending cluster labels.
#' @param cluster_label Column in `colData(sce)` with cluster/cell-type labels.
#' @return A SingleCellExperiment with Slingshot results.
#' @export
togo_run_slingshot <- function(sce, pca_obj, n_pcs = 6, start_cluster = NULL,
                               end_cluster = NULL, cluster_label = "celltype") {
  .togo_need(c("slingshot", "uwot", "SingleCellExperiment", "S4Vectors"))
  rd1 <- pca_obj$x[, seq_len(n_pcs)]
  umap_mat <- uwot::umap(t(log1p(SummarizedExperiment::assays(sce)$norm)))
  colnames(umap_mat) <- c("UMAP1", "UMAP2")
  SingleCellExperiment::reducedDims(sce) <-
    S4Vectors::SimpleList(PCA = rd1, UMAP = umap_mat)
  slingshot::slingshot(sce, clusterLabels = cluster_label, reducedDim = "PCA",
                       start.clus = start_cluster, end.clus = end_cluster)
}
