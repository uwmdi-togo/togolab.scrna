# Internal helpers and import declarations for togolab.scrna.

#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
NULL

# Internal: stop with a clear message if any required (Suggested) package is
# missing. Mirrors togolab's guard so each function fails gracefully.
.togo_need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("This function requires the package(s): ",
         paste(missing, collapse = ", "),
         ".\nInstall with install.packages(c(",
         paste(sprintf('\"%s\"', missing), collapse = ", "),
         ")) (or via Bioconductor where applicable).",
         call. = FALSE)
  }
  invisible(TRUE)
}
