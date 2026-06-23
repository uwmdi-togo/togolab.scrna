# Pooled-library size-factor offset for NEBULA / count models.

# Round a (sparse) counts matrix, optionally with stochastic rounding.
# Operates only on the non-zero entries (the @x slot) so sparsity is preserved
# and large matrices aren't densified.
.togo_round_counts <- function(counts, stochastic = FALSE, seed = 1) {
  .togo_need("Matrix")
  if (!methods::is(counts, "CsparseMatrix")) {
    counts <- methods::as(counts, "CsparseMatrix")
  }
  v <- counts@x
  if (isTRUE(stochastic)) {
    if (!is.null(seed)) set.seed(seed)
    # floor(x) + Bernoulli(x - floor(x)) -> integer, expectation = x
    counts@x <- floor(v) + stats::rbinom(length(v), size = 1, prob = v - floor(v))
  } else {
    counts@x <- round(v)
  }
  Matrix::drop0(counts)   # remove any entries rounded down to 0
}

#' Add a pooled-library size-factor offset to a Seurat object
#'
#' Computes scran pooled size factors (`scran::computeSumFactors()`) from the
#' counts matrix and stores them as a metadata column (default `pooled_offset`),
#' for use as the offset in NEBULA / negative-binomial models
#' (see [togo_run_nebula()]).
#'
#' Counts are rounded to integers first. With `stochastic = FALSE` (default)
#' this is ordinary rounding; with `stochastic = TRUE` it uses stochastic
#' rounding, `floor(x) + rbinom(n, 1, x - floor(x))`, seeded for reproducibility
#' — useful when counts are non-integer (e.g. ambient-corrected) and you want to
#' preserve expected totals.
#'
#' @param so A Seurat object.
#' @param layer Assay layer holding counts (Seurat v5). Default `"counts"`.
#' @param offset_col Metadata column to store the size factors in. Default
#'   `"pooled_offset"`.
#' @param stochastic If `TRUE`, use stochastic rounding (see Details). Default `FALSE`.
#' @param seed Seed used when `stochastic = TRUE`, for reproducibility. Default `1`.
#' @param workers Number of parallel workers for [scran::computeSumFactors()].
#'   `1` (default) runs serially; `>1` uses [BiocParallel::MulticoreParam()].
#' @param clusters Optional cluster assignment passed to
#'   [scran::computeSumFactors()] (e.g. from `scran::quickCluster()`); `NULL`
#'   (default) pools across all cells.
#' @param BPPARAM Optional `BiocParallelParam` to override `workers` entirely.
#' @return The Seurat object with the size-factor offset added to `meta.data`.
#' @export
#' @examples
#' \dontrun{
#' so <- togo_add_pooled_offset(so)                       # standard rounding
#' so <- togo_add_pooled_offset(so, stochastic = TRUE, seed = 42)
#' so <- togo_add_pooled_offset(so, workers = 63)         # parallel
#' fit <- togo_run_nebula(~ group, so, offset_col = "pooled_offset")
#' }
togo_add_pooled_offset <- function(so,
                                   layer      = "counts",
                                   offset_col = "pooled_offset",
                                   stochastic = FALSE,
                                   seed       = 1,
                                   workers    = 1,
                                   clusters   = NULL,
                                   BPPARAM    = NULL) {
  .togo_need(c("Seurat", "Matrix", "SingleCellExperiment", "scran", "BiocParallel"))

  counts <- Seurat::GetAssayData(so, layer = layer)
  counts <- .togo_round_counts(counts, stochastic = stochastic, seed = seed)

  if (is.null(BPPARAM)) {
    BPPARAM <- if (workers > 1) {
      BiocParallel::MulticoreParam(workers = workers)
    } else {
      BiocParallel::SerialParam()
    }
  }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce <- scran::computeSumFactors(sce, clusters = clusters, BPPARAM = BPPARAM)
  size_factors <- SingleCellExperiment::sizeFactors(sce)

  if (any(is.na(size_factors)) || any(size_factors <= 0)) {
    warning("Some size factors are non-positive or NA; consider supplying ",
            "`clusters` (e.g. from scran::quickCluster()).", call. = FALSE)
  }

  so[[offset_col]] <- size_factors
  so
}
