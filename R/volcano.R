# Generalized volcano plot for differential-expression / association results.

#' Volcano plot of differential-expression results
#'
#' A general volcano plot: colors points by significance and direction, labels
#' the top genes in each direction, and draws direction arrows along the x-axis.
#' Generalized from the lab's ATTEMPT volcano (no treatment-specific defaults).
#'
#' @param data Data frame of results.
#' @param fc_col Column with effect size / log fold-change.
#' @param p_col Column with the p-value (or adjusted p-value) to plot on -log10 y.
#' @param label_col Column with point labels (e.g. gene symbols). Default `"Gene"`.
#' @param title Plot title.
#' @param x_lab,y_lab Axis labels.
#' @param p_thresh Significance threshold (horizontal dashed line). Default `0.05`.
#' @param positive_text,negative_text Arrow annotations for the up/down directions.
#' @param n_label Max genes to label per direction. Default `20`.
#' @param genes_to_label Optional explicit vector of labels to restrict labeling to.
#' @param fc_limit Drop points with `abs(fc) >= fc_limit` before plotting. Default `10`.
#' @param pos_color,neg_color,ns_color Point colors for up / down / non-significant.
#' @param text_size,point_size,label_size Sizing controls.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' togo_plot_volcano(res, fc_col = "logFC", p_col = "fdr",
#'                   x_lab = "log2 FC", y_lab = "-log10(FDR)")
#' }
togo_plot_volcano <- function(data, fc_col, p_col, label_col = "Gene",
                              title = NULL, x_lab = fc_col, y_lab = "-log10(p)",
                              p_thresh = 0.05,
                              positive_text = "Up", negative_text = "Down",
                              n_label = 20, genes_to_label = NULL,
                              fc_limit = 10,
                              pos_color = "#f28482", neg_color = "#457b9d",
                              ns_color = "#ced4da",
                              text_size = 15, point_size = 3, label_size = 3) {
  .togo_need(c("dplyr", "ggplot2", "ggrepel", "rlang"))
  fc <- rlang::sym(fc_col); pp <- rlang::sym(p_col); lab <- rlang::sym(label_col)

  data <- data %>%
    dplyr::mutate(.neg_log_p = -log10(.data[[p_col]] + 1e-300)) %>%
    dplyr::filter(abs(.data[[fc_col]]) < fc_limit)

  y_max <- max(data$.neg_log_p, na.rm = TRUE) * 1.1

  top_pos <- data %>% dplyr::filter(.data[[fc_col]] > 0 & .data[[p_col]] < p_thresh) %>%
    dplyr::arrange(.data[[p_col]])
  top_neg <- data %>% dplyr::filter(.data[[fc_col]] < 0 & .data[[p_col]] < p_thresh) %>%
    dplyr::arrange(.data[[p_col]])
  n_pos <- nrow(top_pos); n_neg <- nrow(top_neg)

  pick <- function(df) {
    if (!is.null(genes_to_label)) df <- df[df[[label_col]] %in% genes_to_label, ]
    utils::head(df, n_label)
  }
  label_genes <- c(pick(top_pos)[[label_col]], pick(top_neg)[[label_col]])

  data <- data %>%
    dplyr::mutate(
      .col = dplyr::case_when(.data[[label_col]] %in% top_pos[[label_col]] ~ pos_color,
                              .data[[label_col]] %in% top_neg[[label_col]] ~ neg_color,
                              TRUE ~ ns_color),
      .lab = ifelse(.data[[label_col]] %in% label_genes, as.character(.data[[label_col]]), ""))

  max_fc <- max(data[[fc_col]], na.rm = TRUE); min_fc <- min(data[[fc_col]], na.rm = TRUE)

  ggplot2::ggplot(data, ggplot2::aes(x = !!fc, y = .data$.neg_log_p)) +
    ggplot2::geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "darkgrey") +
    ggplot2::geom_point(ggplot2::aes(color = .data$.col), alpha = 0.5, size = point_size) +
    ggrepel::geom_label_repel(
      data = dplyr::filter(data, .data$.lab != ""),
      ggplot2::aes(label = .data$.lab, color = .data$.col),
      fontface = "bold", size = label_size, max.overlaps = Inf,
      segment.alpha = 0.3, label.size = 0, fill = ggplot2::alpha("white", 0.7)) +
    ggplot2::scale_color_identity() +
    ggplot2::annotate("segment", x = max_fc / 8, xend = max_fc * 7 / 8,
                      y = -y_max * 0.09, col = "darkgrey",
                      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"))) +
    ggplot2::annotate("text", x = mean(c(max_fc / 8, max_fc * 7 / 8)),
                      y = -y_max * 0.14, label = positive_text, color = "#343a40") +
    ggplot2::annotate("segment", x = min_fc / 8, xend = min_fc * 7 / 8,
                      y = -y_max * 0.09, col = "darkgrey",
                      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"))) +
    ggplot2::annotate("text", x = mean(c(min_fc / 8, min_fc * 7 / 8)),
                      y = -y_max * 0.14, label = negative_text, color = "#343a40") +
    ggplot2::coord_cartesian(ylim = c(0, y_max), clip = "off") +
    ggplot2::labs(title = title, x = x_lab, y = y_lab,
                  caption = paste0("Positive n = ", n_pos, " | Negative n = ", n_neg)) +
    ggplot2::theme(legend.position = "none", panel.grid = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = text_size),
                   panel.background = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(hjust = 0.5))
}
