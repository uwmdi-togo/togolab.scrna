# Single-cell RNA-seq helpers: subsetting, QC, and NEBULA differential
# expression. Heavy packages (Seurat, nebula, DoubletFinder, doParallel) are
# Suggested and guarded at call time so the rest of togolab stays light.

#' Subset a Seurat object by cell-type groups
#'
#' Splits a Seurat object into a named list of subsets, one per group of cell
#' types. Generalizes across datasets via the `celltype_col` argument.
#'
#' @param object A Seurat object.
#' @param groups Named list mapping a group name to a vector of cell-type labels.
#' @param celltype_col Column in `object@meta.data` holding cell-type labels.
#' @param assign_global If `TRUE`, also assign each subset into the global
#'   environment as `<prefix><group>` (legacy behavior). Default `FALSE`.
#' @param prefix Object-name prefix used when `assign_global = TRUE`.
#' @return A list with `subsets` (named list of Seurat objects) and `summary`
#'   (data frame of group, n_cells).
#' @export
togo_make_subsets <- function(object, groups, celltype_col = "celltype",
                              assign_global = FALSE, prefix = "so_") {
  .togo_need("Seurat")
  stopifnot(celltype_col %in% colnames(object@meta.data))

  objs <- lapply(names(groups), function(grp) {
    types <- as.character(groups[[grp]])
    cells <- rownames(object@meta.data)[object@meta.data[[celltype_col]] %in% types]
    if (length(cells) == 0) return(NULL)
    subset(object, cells = cells)
  })
  names(objs) <- names(groups)
  objs <- objs[!vapply(objs, is.null, logical(1))]

  if (isTRUE(assign_global)) {
    for (nm in names(objs)) {
      obj_name <- paste0(prefix, tolower(gsub("[^A-Za-z0-9]+", "_", nm)))
      assign(obj_name, objs[[nm]], envir = .GlobalEnv)
    }
  }

  summary <- data.frame(
    group   = names(objs),
    n_cells = vapply(objs, function(o) as.integer(ncol(o)), integer(1)),
    row.names = NULL
  )
  list(subsets = objs, summary = summary)
}

#' Run DoubletFinder on a single sample
#'
#' Standard Seurat preprocessing followed by DoubletFinder, returning a data
#' frame of cell IDs and Singlet/Doublet classifications. Generalizable to any
#' 10x scRNA sample.
#'
#' @param object A Seurat object for a single sample.
#' @param multiplet_rate Expected multiplet rate. If `NULL`, estimated from the
#'   number of recovered cells using the 10x rate table.
#' @param sample_id_col Meta column identifying the sample (for messages).
#' @return A data frame with `row_names` (cell IDs) and `doublet_finder`.
#' @export
togo_run_doubletfinder <- function(object, multiplet_rate = NULL,
                                   sample_id_col = "SampleID") {
  .togo_need(c("Seurat", "DoubletFinder", "tibble", "dplyr"))

  if (sample_id_col %in% colnames(object@meta.data)) {
    message("Sample ", unique(object@meta.data[[sample_id_col]]))
  }

  if (is.null(multiplet_rate)) {
    message("multiplet_rate not provided; estimating from cell count")
    rates_10x <- data.frame(
      Multiplet_rate  = c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046,
                          0.054, 0.061, 0.069, 0.076),
      Recovered_cells = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000,
                          8000, 9000, 10000)
    )
    multiplet_rate <- rates_10x %>%
      dplyr::filter(.data$Recovered_cells < nrow(object@meta.data)) %>%
      dplyr::slice(which.max(.data$Recovered_cells)) %>%
      dplyr::pull(.data$Multiplet_rate)
    message("Setting multiplet rate to ", multiplet_rate)
  }

  s <- Seurat::NormalizeData(object)
  s <- Seurat::FindVariableFeatures(s)
  s <- Seurat::ScaleData(s)
  s <- Seurat::RunPCA(s, nfeatures.print = 10)

  stdv <- s[["pca"]]@stdev
  pct  <- stdv / sum(stdv) * 100
  cum  <- cumsum(pct)
  co1  <- which(cum > 90 & pct < 5)[1]
  co2  <- sort(which((pct[seq_len(length(pct) - 1)] - pct[-1]) > 0.1),
               decreasing = TRUE)[1] + 1
  min_pc <- min(co1, co2)

  s <- Seurat::RunUMAP(s, dims = 1:min_pc)
  s <- Seurat::FindNeighbors(s, dims = 1:min_pc)
  s <- Seurat::FindClusters(s, resolution = 0.1)

  sweep_list  <- DoubletFinder::paramSweep(s, PCs = 1:min_pc, sct = FALSE)
  sweep_stats <- DoubletFinder::summarizeSweep(sweep_list)
  bcmvn       <- DoubletFinder::find.pK(sweep_stats)
  optimal_pk  <- as.numeric(as.character(
    bcmvn$pK[which.max(bcmvn$BCmetric)]))

  homotypic <- DoubletFinder::modelHomotypic(s@meta.data$seurat_clusters)
  nexp      <- round(multiplet_rate * nrow(s@meta.data))
  nexp_adj  <- round(nexp * (1 - homotypic))

  s <- DoubletFinder::doubletFinder(s, PCs = 1:min_pc, pK = optimal_pk,
                                    nExp = nexp_adj)
  colnames(s@meta.data)[grepl("DF.classifications.*", colnames(s@meta.data))] <-
    "doublet_finder"

  res <- s@meta.data["doublet_finder"]
  tibble::rownames_to_column(res, "row_names")
}

#' NEBULA differential expression for a single model formula
#'
#' Fits a NEBULA negative-binomial mixed model across all genes for one model
#' formula. Works for any design: a treatment-by-visit interaction (ATTEMPT),
#' a disease-group main effect (PB90), or a continuous covariate.
#'
#' @param model A formula (e.g. `~ visit * treatment`, `~ group`) or a string
#'   passed to [stats::as.formula()].
#' @param object A Seurat object.
#' @param id_col Subject/grouping id column in `object@meta.data`.
#' @param offset_col Optional meta column with the per-cell offset (e.g.
#'   `"pooled_offset"`). `NULL` for no offset.
#' @param layer Assay layer holding counts (Seurat v5). Default `"counts"`.
#' @param ncore Cores passed to [nebula::nebula()].
#' @return A nebula result object, or `NULL` if too few cells/subjects.
#' @export
togo_run_nebula <- function(model, object, id_col = "record_id",
                            offset_col = NULL, layer = "counts", ncore = 4) {
  .togo_need(c("Seurat", "nebula", "dplyr"))
  meta <- object@meta.data
  form <- if (inherits(model, "formula")) model else stats::as.formula(model)
  vars_used <- all.vars(form)
  keep <- stats::complete.cases(meta[, vars_used, drop = FALSE])
  if (sum(keep) < 2 || dplyr::n_distinct(meta[[id_col]][keep]) < 2) {
    return(NULL)
  }
  meta_v   <- droplevels(meta[keep, , drop = FALSE])
  counts_v <- round(Seurat::GetAssayData(object, layer = layer))[, keep, drop = FALSE]
  offset_v <- if (!is.null(offset_col)) object@meta.data[[offset_col]][keep] else NULL
  pred     <- stats::model.matrix(form, data = meta_v)

  data_g <- nebula::group_cell(count = counts_v, id = meta_v[[id_col]],
                               pred = pred, offset = offset_v)
  if (is.null(data_g)) {
    data_g <- list(count = counts_v, id = meta_v[[id_col]],
                   pred = pred, offset = offset_v)
  }
  nebula::nebula(count = data_g$count, id = data_g$id, pred = data_g$pred,
                 offset = data_g$offset, ncore = ncore, reml = 1,
                 model = "NBLMM")
}

#' NEBULA differential expression in parallel, per gene
#'
#' Fits NEBULA per gene in parallel (one gene per task), returning a named list
#' of per-gene results plus run metadata. Optionally uploads the results to S3
#' via an aws.s3-style client. Design is fully formula-driven and generalizes
#' across datasets.
#'
#' @param object A Seurat object.
#' @param formula Model formula (default `~ group`).
#' @param id_col Subject/grouping id column. Default `"record_id"`.
#' @param offset_col Optional per-cell offset column. Default `"pooled_offset"`.
#' @param layer Counts layer. Default `"counts"`.
#' @param n_cores Number of parallel workers.
#' @param group Whether to call [nebula::group_cell()] (set `FALSE` if already grouped).
#' @param verbose Print per-gene warnings/errors and timing.
#' @return A list: `results` (named list of per-gene nebula fits),
#'   `n_genes_tested`, `n_genes_converged`, `nonconverged_percent`,
#'   `runtime_minutes`.
#' @export
togo_run_nebula_parallel <- function(object,
                                     formula = ~ group,
                                     id_col = "record_id",
                                     offset_col = "pooled_offset",
                                     layer = "counts",
                                     n_cores = max(parallel::detectCores() - 1, 1),
                                     group = TRUE,
                                     verbose = TRUE) {
  .togo_need(c("Seurat", "nebula", "doParallel", "foreach"))

  counts_mat <- round(Seurat::GetAssayData(object, layer = layer))
  genes_list <- rownames(counts_mat)

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  start_time <- Sys.time()
  g <- NULL  # for R CMD check (foreach iterator)
  `%dopar%` <- foreach::`%dopar%`
  res_list <- foreach::foreach(g = genes_list,
                               .packages = c("nebula", "Matrix"),
                               .errorhandling = "pass") %dopar% {
    warn <- err <- res <- NULL
    tryCatch({
      count_gene <- counts_mat[g, , drop = FALSE]
      meta_gene  <- subset(object, features = g)@meta.data
      pred_gene  <- stats::model.matrix(formula, data = meta_gene)
      data_g <- if (group) {
        nebula::group_cell(count = count_gene, id = meta_gene[[id_col]], pred = pred_gene)
      } else {
        list(count = count_gene, id = meta_gene[[id_col]], pred = pred_gene)
      }
      res <- withCallingHandlers(
        nebula::nebula(count = data_g$count, id = data_g$id, pred = data_g$pred,
                       ncore = 1, output_re = TRUE, covariance = TRUE, reml = 1,
                       model = "NBLMM",
                       offset = if (!is.null(offset_col)) meta_gene[[offset_col]] else NULL),
        warning = function(w) { warn <<- conditionMessage(w); invokeRestart("muffleWarning") }
      )
    }, error = function(e) { err <<- conditionMessage(e) })
    list(gene = g, result = res, warning = warn, error = err)
  }

  if (verbose) {
    for (x in res_list) {
      if (!is.null(x$warning)) message(sprintf("[WARN] %s: %s", x$gene, x$warning))
      if (!is.null(x$error))   message(sprintf("[ERR ] %s: %s", x$gene, x$error))
    }
  }

  names(res_list) <- vapply(res_list, function(x) x$gene, "")
  fits <- Filter(Negate(is.null), lapply(res_list, `[[`, "result"))

  runtime <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  nonconv <- (length(genes_list) - length(fits)) / length(genes_list)
  if (verbose) {
    message(sprintf("Runtime: %.2f min; %.2f%% of genes filtered (low expr/convergence).",
                    runtime, nonconv * 100))
  }

  list(results = fits,
       n_genes_tested = length(genes_list),
       n_genes_converged = length(fits),
       nonconverged_percent = nonconv,
       runtime_minutes = runtime)
}

#' Summarize a list of NEBULA results
#'
#' Filters to converged models, row-binds their summary tables, and adds an FDR
#' column for a chosen p-value column. The p-value column name depends on your
#' design (e.g. `"p_groupType 2 Diabetes"` or
#' `"p_treatmentDapagliflozin:visitPOST"`); inspect the summary names if unsure.
#'
#' @param nebula_list Named list of nebula result objects (e.g. the `results`
#'   element from [togo_run_nebula_parallel()]).
#' @param pval_col Name of the p-value column to FDR-adjust. If `NULL`, the
#'   first column starting with `"p_"` is used.
#' @param convergence_cut Minimum convergence code to keep. Default `-10`.
#' @return A list: `convergence` (df), `results` (combined summary with `fdr`),
#'   `overdispersion` (df).
#' @export
togo_process_nebula_results <- function(nebula_list, pval_col = NULL,
                                        convergence_cut = -10) {
  .togo_need(c("dplyr", "purrr"))

  convergence_df <- purrr::map_dfr(names(nebula_list), function(g) {
    data.frame(Gene = g, Convergence_Code = nebula_list[[g]]$convergence)
  })
  converged <- convergence_df$Gene[convergence_df$Convergence_Code >= convergence_cut]

  summary_df <- purrr::map_dfr(converged, function(g) {
    dplyr::mutate(nebula_list[[g]]$summary, Gene = g)
  })

  if (is.null(pval_col)) {
    pcols <- grep("^p_", names(summary_df), value = TRUE)
    pval_col <- if (length(pcols)) pcols[1] else NA_character_
    if (!is.na(pval_col)) message("Using p-value column: ", pval_col)
  }
  if (!is.na(pval_col) && pval_col %in% names(summary_df)) {
    summary_df$fdr <- stats::p.adjust(summary_df[[pval_col]], method = "fdr")
  } else {
    warning("p-value column '", pval_col, "' not found; FDR not computed.")
    summary_df$fdr <- NA_real_
  }

  overdisp_df <- purrr::map_dfr(names(nebula_list), function(g) {
    od <- nebula_list[[g]]$overdispersion
    od$Gene <- g
    od
  })

  list(convergence = convergence_df, results = summary_df, overdispersion = overdisp_df)
}
