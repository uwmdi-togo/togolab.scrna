# Export a Seurat object to per-participant GEO-style processed files on S3.

#' Export a Seurat object to per-participant processed files on S3 (GEO format)
#'
#' For each participant in a Seurat object, writes the counts matrix
#' (Matrix Market `.mtx`), cell barcodes (`.tsv`), and gene features (`.tsv`)
#' to the lab S3 store, computes MD5 checksums, and uploads a CSV manifest of
#' files + checksums. Everything is read from / written to S3 — no local
#' directories are created (temp files are used and cleaned up).
#'
#' Generalized from the lab's ATTEMPT GEO-export script: works on any Seurat
#' object, with a configurable participant-ID column and optional metadata
#' columns to attach to the manifest. All uploads route through
#' `togolab::s3write_using_region()` so `region` is always passed (Kopah
#' requires `region = ""`).
#'
#' Assumes S3 credentials are configured (run `togolab::togo_paths()` first).
#'
#' @param seurat_object A Seurat object.
#' @param project_path S3 key prefix (folder) to write under, e.g.
#'   `"P001/clean_data/GEO_ProcessedData"`.
#' @param bucket S3 bucket. Default `"togo.projects"`.
#' @param region S3 region. Default `""` (Kopah).
#' @param id_col Metadata column identifying participants. Default `"orig.ident"`.
#' @param layer Assay layer holding counts (Seurat v5). Default `"counts"`.
#' @param meta_cols Optional character vector of `meta.data` columns to attach
#'   to a second, metadata-joined manifest (e.g. `c("visit", "treatment")`).
#'   `NULL` (default) skips the metadata manifest.
#' @param suffix Filename suffix for the processed files. Default `"_processed"`.
#' @param write_manifest If `TRUE` (default), upload the file-info CSV(s).
#' @return Invisibly, the file-info `data.frame` (participant IDs, file names,
#'   MD5 checksums).
#' @export
#' @examples
#' \dontrun{
#' library(togolab); togo_paths()
#' so <- togolab::togo_s3_read_rds("P001/clean_data/example.rds",
#'                                 bucket = "togo.projects")
#' togo_export_seurat_geo(so,
#'                        project_path = "P001/clean_data/GEO_ProcessedData",
#'                        meta_cols = c("visit", "treatment"))
#' }
togo_export_seurat_geo <- function(seurat_object,
                                   project_path,
                                   bucket    = "togo.projects",
                                   region    = "",
                                   id_col    = "orig.ident",
                                   layer     = "counts",
                                   meta_cols = NULL,
                                   suffix    = "_processed",
                                   write_manifest = TRUE) {
  .togo_need(c("Seurat", "Matrix", "digest", "aws.s3", "dplyr"))

  md <- seurat_object@meta.data
  if (!id_col %in% colnames(md)) {
    stop("id_col '", id_col, "' not found in seurat_object@meta.data.", call. = FALSE)
  }
  if (!is.null(meta_cols)) {
    miss <- setdiff(meta_cols, colnames(md))
    if (length(miss)) {
      stop("meta_cols not found in meta.data: ", paste(miss, collapse = ", "),
           ".", call. = FALSE)
    }
  }

  participants <- unique(md[[id_col]])

  # Upload a local file to <bucket>/<key> via the region-aware Kopah helper.
  # Routing through togolab::s3write_using_region() guarantees `region` is
  # always passed (Kopah requires region = "").
  upload <- function(local, key) {
    togolab::s3write_using_region(
      FUN    = function(dest, src) file.copy(src, dest, overwrite = TRUE),
      src    = local,
      object = key,
      bucket = bucket,
      region = region
    )
  }
  write_tsv <- function(x, file) {
    utils::write.table(x, file = file, sep = "\t", col.names = FALSE,
                       row.names = FALSE, quote = FALSE)
  }

  info_list <- vector("list", length(participants))

  for (i in seq_along(participants)) {
    pid <- participants[[i]]
    pid_fixed <- gsub("-", "_", pid)

    cells <- rownames(md)[md[[id_col]] == pid]
    sub <- subset(seurat_object, cells = cells)
    counts <- Seurat::GetAssayData(sub, layer = layer)
    if (!methods::is(counts, "CsparseMatrix")) {
      counts <- methods::as(counts, "CsparseMatrix")
    }
    barcodes <- colnames(counts)
    features <- rownames(counts)

    # write to temp files (cleaned up at end of iteration)
    tmp_dir <- tempfile("geo_")
    dir.create(tmp_dir)
    mtx_local <- file.path(tmp_dir, paste0(pid_fixed, "_matrix",   suffix, ".mtx"))
    bc_local  <- file.path(tmp_dir, paste0(pid_fixed, "_barcodes", suffix, ".tsv"))
    ft_local  <- file.path(tmp_dir, paste0(pid_fixed, "_features", suffix, ".tsv"))

    Matrix::writeMM(counts, file = mtx_local)
    write_tsv(barcodes, bc_local)
    write_tsv(features, ft_local)

    md5 <- function(f) digest::digest(f, algo = "md5", file = TRUE)
    mtx_md5 <- md5(mtx_local); bc_md5 <- md5(bc_local); ft_md5 <- md5(ft_local)

    folder <- paste0(project_path, "/", pid_fixed)
    upload(mtx_local, paste0(folder, "/", basename(mtx_local)))
    upload(bc_local,  paste0(folder, "/", basename(bc_local)))
    upload(ft_local,  paste0(folder, "/", basename(ft_local)))

    info_list[[i]] <- data.frame(
      Participant_ID = pid_fixed,
      MTX_File       = basename(mtx_local), MTX_MD5      = mtx_md5,
      Barcodes_File  = basename(bc_local),  Barcodes_MD5 = bc_md5,
      Features_File  = basename(ft_local),  Features_MD5 = ft_md5,
      stringsAsFactors = FALSE
    )
    unlink(tmp_dir, recursive = TRUE)
    message(sprintf("Uploaded %s (%d/%d)", pid_fixed, i, length(participants)))
  }

  file_info <- do.call(rbind, info_list)

  if (isTRUE(write_manifest)) {
    upload_csv <- function(df, key) {
      tmp <- tempfile(fileext = ".csv")
      on.exit(unlink(tmp), add = TRUE)
      utils::write.csv(df, tmp, row.names = FALSE)
      upload(tmp, key)
    }
    upload_csv(file_info,
               paste0(project_path, "/processed_participant_file_info_with_md5.csv"))

    if (!is.null(meta_cols)) {
      meta_df <- md %>%
        dplyr::mutate(Participant_ID = gsub("-", "_", .data[[id_col]])) %>%
        dplyr::select("Participant_ID", dplyr::all_of(meta_cols)) %>%
        dplyr::right_join(file_info, by = "Participant_ID") %>%
        dplyr::distinct(.data$Participant_ID, .keep_all = TRUE)
      upload_csv(meta_df,
                 paste0(project_path,
                        "/processed_participant_file_info_with_md5_meta.csv"))
    }
    message("Manifest(s) uploaded to s3://", bucket, "/", project_path)
  }

  invisible(file_info)
}
