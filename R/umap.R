# UMAP feature plots and cell-type composition helpers.

#' Prepare UMAP metadata, embeddings, and expression for plotting
#'
#' Pulls a Seurat object's UMAP embeddings, metadata, and expression for chosen
#' genes into a tidy data frame, computes per-cell-type label centers, and
#' builds a cell-type color mapping. Generalizes via `celltype_col` and
#' `umap_reduction`.
#'
#' @param object A Seurat object.
#' @param genes Character vector of genes to fetch.
#' @param celltype_col Cell-type column in `object@meta.data`.
#' @param umap_reduction Name of the UMAP reduction (e.g. `"umap.harmony"`).
#' @param color_palette RColorBrewer palette name used when `custom_colors` is `NULL`.
#' @param custom_colors Optional named vector of colors overriding the palette.
#' @return A list: `metadata` (df with `umapharmony_1/2` + expression),
#'   `centers` (per-cell-type median x/y), `celltype_colors` (named vector).
#' @export
togo_prepare_umap_metadata <- function(object, genes,
                                       celltype_col = "KPMP_celltype",
                                       umap_reduction = "umap.harmony",
                                       color_palette = "Set3",
                                       custom_colors = NULL) {
  .togo_need(c("Seurat", "dplyr"))

  celltypes <- levels(object@meta.data[[celltype_col]])
  if (is.null(celltypes)) celltypes <- sort(unique(object@meta.data[[celltype_col]]))
  n_celltypes <- length(celltypes)

  if (!is.null(custom_colors)) {
    celltype_colors <- custom_colors
  } else {
    .togo_need("RColorBrewer")
    max_cols <- RColorBrewer::brewer.pal.info[color_palette, "maxcolors"]
    base_cols <- RColorBrewer::brewer.pal(min(max_cols, max(3, n_celltypes)), color_palette)
    umap_colors <- grDevices::colorRampPalette(base_cols)(n_celltypes)
    celltype_colors <- stats::setNames(umap_colors, celltypes)
  }

  expr_df <- as.data.frame(Seurat::FetchData(object, vars = genes))
  colnames(expr_df) <- genes

  embed <- object@reductions[[umap_reduction]]@cell.embeddings
  umap_embed_cols <- colnames(embed)
  metadata <- cbind(object@meta.data, embed, expr_df)
  colnames(metadata)[colnames(metadata) %in% umap_embed_cols] <-
    c("umapharmony_1", "umapharmony_2")

  centers <- metadata %>%
    dplyr::group_by(.data[[celltype_col]]) %>%
    dplyr::summarise(x = stats::median(.data$umapharmony_1),
                     y = stats::median(.data$umapharmony_2), .groups = "drop")

  list(metadata = metadata, centers = centers, celltype_colors = celltype_colors)
}

#' Feature UMAP plots highlighting expressing cell types
#'
#' For each gene, draws a UMAP where cells expressing the gene are colored by
#' cell type (others greyed out), labels cell types passing an expression
#' threshold, and adds a markdown caption of percent-expressed per cell type.
#'
#' @param genes Character vector of genes to plot.
#' @param metadata Data frame from [togo_prepare_umap_metadata()] (`$metadata`).
#' @param centers Label-center data frame (`$centers`) with `x`, `y`, and the
#'   cell-type column.
#' @param celltype_colors Named vector of cell-type colors (`$celltype_colors`).
#' @param celltype_col Cell-type column name (must match `centers`). Default `"KPMP_celltype"`.
#' @param pct_threshold Minimum percent expressed to highlight/label a cell type.
#' @param save_fun Optional function `(plot, gene) -> NULL` to save each plot
#'   (e.g. wrap [togo_s3_save_plot()] or `ggplot2::ggsave`). `NULL` to skip saving.
#' @param plot_height,plot_width Passed to `save_fun` via attributes (height/width).
#' @return Invisibly, a named list of ggplot objects (one per gene).
#' @export
togo_plot_feature_umap <- function(genes, metadata, centers, celltype_colors,
                                   celltype_col = "KPMP_celltype",
                                   pct_threshold = 10,
                                   save_fun = NULL,
                                   plot_height = 10, plot_width = 10) {
  .togo_need(c("dplyr", "ggplot2", "ggrepel", "ggtext"))
  plots <- list()

  for (gene in genes) {
    expr_stats <- metadata %>%
      dplyr::group_by(.data[[celltype_col]]) %>%
      dplyr::summarise(pct_expressed = round(mean(.data[[gene]] > 0) * 100, 1),
                       n_total = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(.data$pct_expressed > pct_threshold)

    pass <- expr_stats[[celltype_col]]
    centers_filtered <- centers %>%
      dplyr::mutate(passes = ifelse(.data[[celltype_col]] %in% pass, "Y", "N"))

    caption_text <- expr_stats %>%
      dplyr::arrange(dplyr::desc(.data$pct_expressed)) %>%
      dplyr::mutate(
        color = celltype_colors[as.character(.data[[celltype_col]])],
        label = paste0("<span style='color:", .data$color, "'>&#9679; </span>**",
                       .data[[celltype_col]], "**: ", .data$pct_expressed, "%"),
        grp = ceiling(dplyr::row_number() / 5)) %>%
      dplyr::group_by(.data$grp) %>%
      dplyr::summarise(line = paste(.data$label, collapse = " | "), .groups = "drop") %>%
      dplyr::pull(.data$line) %>%
      paste(collapse = "<br><br>")
    caption_text <- paste0("**% Expressed (>", pct_threshold, "%):**<br><br>", caption_text)

    feature_p <- metadata %>%
      dplyr::mutate(feature_color_logic = ifelse(.data[[gene]] > 0,
                                                 as.character(.data[[celltype_col]]), "No")) %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$umapharmony_1, y = .data$umapharmony_2,
                                   color = .data$feature_color_logic)) +
      ggplot2::geom_point(alpha = 0.4) +
      ggrepel::geom_text_repel(
        data = centers_filtered,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data[[celltype_col]]),
        inherit.aes = FALSE, size = 4,
        fontface = ifelse(centers_filtered$passes == "Y", "bold", "plain"),
        color = ifelse(centers_filtered$passes == "Y", "black", "#495057"),
        max.overlaps = Inf) +
      ggplot2::scale_color_manual(values = c(celltype_colors, "No" = "#edede9")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(size = 15),
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.caption = ggtext::element_markdown(size = 12, hjust = 0.5)) +
      togolab::theme_togo_transparent() +
      ggplot2::labs(x = "UMAP1", y = "UMAP2", title = gene, caption = caption_text)

    plots[[gene]] <- feature_p
    if (!is.null(save_fun)) {
      attr(feature_p, "height") <- plot_height
      attr(feature_p, "width")  <- plot_width
      save_fun(feature_p, gene)
    }
  }
  invisible(plots)
}

#' Cell-type proportions within bins of a continuous variable
#'
#' Bins a continuous variable (default pseudotime) and computes the percentage
#' of cells of each group within each bin. Useful for composition-vs-pseudotime
#' or composition-vs-any-axis plots.
#'
#' @param data A data frame of cells.
#' @param value_col Continuous column to bin (default `"slingPseudotime_1"`).
#' @param bin_width Width of each bin.
#' @param group_by Grouping column (e.g. cell type). Default `"KPMP_celltype"`.
#' @return A data frame of bins x groups with `n`, `total_in_bin`, `proportion` (percent).
#' @export
togo_celltype_proportions <- function(data, value_col = "slingPseudotime_1",
                                      bin_width = 10, group_by = "KPMP_celltype") {
  .togo_need(c("dplyr", "tidyr", "rlang"))
  grp <- rlang::sym(group_by)

  binned <- data %>%
    dplyr::mutate(bin_start = floor(.data[[value_col]] / bin_width) * bin_width,
                  bin_label = paste0(.data$bin_start, "-", .data$bin_start + bin_width))

  props <- binned %>%
    dplyr::group_by(.data$bin_label, .data$bin_start, !!grp) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(total_in_bin = sum(.data$n),
                  proportion = .data$n / .data$total_in_bin * 100) %>%
    dplyr::ungroup()

  props %>%
    tidyr::complete(bin_label = unique(props$bin_label),
                    !!grp := unique(data[[group_by]]),
                    fill = list(n = 0, proportion = 0)) %>%
    dplyr::group_by(.data$bin_label) %>%
    dplyr::mutate(bin_start = min(.data$bin_start, na.rm = TRUE),
                  total_in_bin = sum(.data$n)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$bin_start, !!grp)
}

#' Pie chart of group composition for one bin
#'
#' Draws a pie chart of group proportions for a single bin (e.g. one pseudotime
#' bin) with a colored markdown caption. Pairs with [togo_celltype_proportions()].
#'
#' @param data Output of [togo_celltype_proportions()].
#' @param bin_value `bin_start` value selecting the bin to plot.
#' @param color_palette Named vector of colors for groups.
#' @param group_by Grouping column. Default `"KPMP_celltype"`.
#' @param caption_groups Optional ordered subset of groups for the caption/legend.
#' @param digits Rounding for caption percentages.
#' @return A ggplot object.
#' @export
togo_celltype_pie <- function(data, bin_value, color_palette,
                              group_by = "KPMP_celltype",
                              caption_groups = NULL, digits = 1) {
  .togo_need(c("dplyr", "ggplot2", "rlang", "ggtext"))
  grp <- rlang::sym(group_by)

  bin_data <- dplyr::filter(data, .data$bin_start == bin_value)
  if (nrow(bin_data) == 0) stop(sprintf("No rows for bin_start == %s", bin_value))

  cap_df <- bin_data %>%
    dplyr::group_by(!!grp) %>%
    dplyr::summarise(proportion = sum(.data$proportion, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$proportion))

  if (!is.null(caption_groups)) {
    cap_df <- cap_df %>%
      dplyr::filter(!!grp %in% caption_groups) %>%
      dplyr::mutate(!!grp := factor(!!grp, levels = caption_groups)) %>%
      dplyr::arrange(!!grp)
    bin_data[[group_by]] <- factor(bin_data[[group_by]], levels = caption_groups)
  }

  caption_text <- paste(vapply(seq_len(nrow(cap_df)), function(i) {
    g   <- as.character(cap_df[[group_by]][i])
    pct <- round(cap_df$proportion[i], digits)
    col <- if (!is.null(color_palette[[g]])) color_palette[[g]] else "#000000"
    sprintf("<span style='color:%s'>%s: %s%%</span>", col, g, pct)
  }, character(1)), collapse = "<br>")

  present <- unique(as.character(bin_data[[group_by]]))
  pal_used <- color_palette[names(color_palette) %in% present]

  ggplot2::ggplot(bin_data, ggplot2::aes(x = "", y = .data$proportion, fill = !!grp)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::coord_polar("y") +
    ggplot2::theme_void() +
    ggplot2::labs(fill = NULL, caption = caption_text) +
    ggplot2::theme(plot.caption = ggtext::element_markdown(hjust = 0.5, face = "bold")) +
    ggplot2::scale_fill_manual(values = pal_used)
}
