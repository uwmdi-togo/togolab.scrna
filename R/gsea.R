# Pathway / GSEA helpers: GMT preparation, pathway-name cleaning, redundancy
# filtering, and an fgsea dot plot. All dataset-agnostic.

# Internal defaults for acronym capitalization in pathway names.
.togo_acronyms <- c("rna", "dna", "mtor", "foxo", "ppar", "nmd", "fgfr", "robo",
                    "bhl", "cov", "jak", "stat", "wnt", "hiv", "bcl", "mapk",
                    "pt", "tal", "pc", "ic", "ec")
.togo_special_mixed <- c("rrna", "mrna", "trna", "gtpase", "atpase", "robos", "slits")
.togo_special_repl  <- c("rRNA", "mRNA", "tRNA", "GTPase", "ATPase", "ROBOs", "SLITs")

.togo_replace_mixed_case <- function(text, from, to) {
  .togo_need("stringr")
  for (i in seq_along(from)) {
    text <- stringr::str_replace_all(text,
      stringr::regex(paste0("\\b", from[i], "\\b"), ignore_case = TRUE), to[i])
  }
  text
}

.togo_capitalize_acronyms <- function(text, terms) {
  .togo_need("stringr")
  for (term in terms) {
    text <- stringr::str_replace_all(text,
      stringr::regex(paste0("\\b", term, "\\b"), ignore_case = TRUE), toupper(term))
  }
  text
}

#' Convert a binary gene-set membership matrix to a list
#'
#' @param pws A matrix with genes as rows, gene sets as columns (0/1 entries).
#' @return A named list mapping each gene set to its member genes.
#' @export
togo_matrix_to_list <- function(pws) {
  lst <- list()
  for (pw in colnames(pws)) lst[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  lst
}

#' Prepare a GMT gene-set file for GSEA, subset to genes in your data
#'
#' Reads a GMT file, restricts gene sets to genes present in your data, and
#' drops sets with fewer than 5 annotated genes.
#'
#' @param gmt_file Path to a `.gmt` file.
#' @param genes_in_data Character vector of genes present in your dataset.
#' @param savefile If `TRUE`, save the subset list as an `.RData` file.
#' @return A named list of gene sets (the GSEA pathway list).
#' @export
togo_prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE) {
  .togo_need("fgsea")
  gmt <- fgsea::gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in seq_len(ncol(mat))) mat[, i] <- as.numeric(hidden %in% gmt[[i]])
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1, ]) > 5)]]
  final_list <- togo_matrix_to_list(mat)
  if (savefile) {
    saveRDS(final_list, file = paste0(gsub(".gmt", "", gmt_file), "_subset_",
                                      format(Sys.time(), "%d%m"), ".RData"))
  }
  final_list
}

#' Clean pathway names for display
#'
#' Strips database prefixes (REACTOME_, GOBP_, KEGG_, HALLMARK_), converts
#' underscores to spaces, applies title case, and restores common biological
#' acronyms (DNA, ATP, NAD, etc.) to uppercase.
#'
#' @param pathways Character vector of raw pathway names.
#' @return A character vector of cleaned names.
#' @export
togo_clean_pathway_names <- function(pathways) {
  cleaned <- gsub("^REACTOME_|^GOBP_|^GOMF_|^KEGG_|^HALLMARK_", "", pathways)
  cleaned <- gsub("_", " ", cleaned)
  cleaned <- tools::toTitleCase(tolower(cleaned))
  uppercase_words <- c("\\bI\\b", "\\bIi\\b", "\\bIii\\b", "\\bIv\\b", "\\bV\\b",
                       "\\bTca\\b", "\\bAtp\\b", "\\bAdp\\b", "\\bAmp\\b", "\\bGtp\\b",
                       "\\bNad\\b", "\\bNadh\\b", "\\bDna\\b", "\\bRna\\b", "\\bMrna\\b",
                       "\\bEr\\b", "\\bMhc\\b", "\\bTgf\\b", "\\bVegf\\b", "\\bRos\\b",
                       "\\bNos\\b", "\\bMapk\\b")
  for (pattern in uppercase_words) {
    word <- gsub("\\\\b", "", pattern)
    cleaned <- gsub(pattern, toupper(word), cleaned, ignore.case = TRUE)
  }
  cleaned <- gsub("\\bOf\\b", "of", cleaned)
  cleaned <- gsub("\\bBy\\b", "by", cleaned)
  cleaned <- gsub("\\bThe\\b", "the", cleaned)
  cleaned <- gsub("^the", "The", cleaned)
  cleaned
}

#' Drop redundant pathways by leading-edge overlap
#'
#' Collapses pathways whose leading-edge gene sets overlap above a Jaccard
#' threshold, keeping the most significant representative from each cluster.
#'
#' @param gsea_result Data frame with `pathway`, `leadingEdge` (list), `padj`.
#' @param overlap_pct Jaccard overlap threshold above which pathways are merged.
#' @return The filtered `gsea_result` data frame.
#' @export
togo_filter_redundant_pathways <- function(gsea_result, overlap_pct = 0.3) {
  .togo_need("igraph")
  if (!all(c("pathway", "leadingEdge", "padj") %in% colnames(gsea_result))) {
    stop("Input must have 'pathway', 'leadingEdge', and 'padj' columns.", call. = FALSE)
  }
  le <- gsea_result$leadingEdge
  names(le) <- gsea_result$pathway
  overlap <- sapply(le, function(x) sapply(le, function(y)
    length(intersect(x, y)) / length(union(x, y))))
  pairs <- which(overlap > overlap_pct & lower.tri(overlap), arr.ind = TRUE)
  if (nrow(pairs) == 0) {
    message("No redundant pathways above threshold.")
    return(gsea_result)
  }
  edges <- data.frame(from = rownames(overlap)[pairs[, 1]],
                      to   = colnames(overlap)[pairs[, 2]])
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  comp <- igraph::components(g)
  clusters <- split(names(comp$membership), comp$membership)
  reps <- sapply(clusters, function(terms) {
    sub <- gsea_result[gsea_result$pathway %in% terms, ]
    sub[which.min(sub$padj), "pathway"]
  })
  keep <- c(setdiff(gsea_result$pathway, unique(unlist(clusters))), reps)
  gsea_result[gsea_result$pathway %in% keep, ]
}

#' Dot plot of top enriched pathways from fgsea
#'
#' Plots the top pathways by p-value as a dot plot of |NES| with cleaned names,
#' colored by direction and significance.
#'
#' @param fgsea_res An fgsea results data frame (`pathway`, `pval`, `NES`, `size`).
#' @param top_n Number of pathways to show.
#' @param title Plot title.
#' @param xmin,xmax X-axis (|NES|) limits.
#' @param text_label,text_axis,text_title Sizing controls.
#' @return A ggplot object.
#' @export
togo_plot_fgsea <- function(fgsea_res, top_n = 30, title = "Top Enriched Pathways",
                            xmin = 0, xmax = 3,
                            text_label = 6.5, text_axis = 18, text_title = 20) {
  .togo_need(c("dplyr", "ggplot2", "stringr", "forcats"))
  fgsea_res <- fgsea_res %>%
    dplyr::arrange(.data$pval) %>% utils::head(top_n) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        .data$NES < 0 & .data$pval <= 0.05 ~ "Negative",
        .data$NES > 0 & .data$pval <= 0.05 ~ "Positive",
        .data$NES < 0 ~ "Negative p > 0.05",
        TRUE ~ "Positive p > 0.05"),
      face = ifelse(.data$pval <= 0.05, "bold", "plain"),
      pathway_clean = .togo_capitalize_acronyms(
        stringr::str_to_sentence(stringr::str_replace_all(togo_clean_pathway_names(.data$pathway), "_", " ")),
        .togo_acronyms),
      pathway_clean = paste0(.data$pathway_clean, " (", .data$size, ")"))

  fgsea_res$pathway_clean <- stats::reorder(fgsea_res$pathway_clean, -abs(fgsea_res$NES))

  ggplot2::ggplot(fgsea_res, ggplot2::aes(x = abs(.data$NES),
                                          y = forcats::fct_rev(.data$pathway_clean),
                                          label = .data$pathway_clean)) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(.data$pval), color = .data$direction), alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(color = .data$direction, fontface = .data$face),
                       hjust = 0, size = text_label, nudge_x = (xmax - xmin) / 100) +
    ggplot2::scale_color_manual(values = c("Positive" = "#c75146", "Negative" = "#2c7da0",
                                           "Positive p > 0.05" = "#e18c80",
                                           "Negative p > 0.05" = "#7ab6d1")) +
    ggplot2::scale_x_continuous(limits = c(xmin, xmax), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = "NES", y = "Pathways", color = "Direction",
                  size = "-log10(p)", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = text_axis),
                   plot.title = ggplot2::element_text(size = text_title))
}
