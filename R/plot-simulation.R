# ==============================================================================
# plot.aucmat_simulation() + plot_correlation_heatmap()
# ==============================================================================

#' Plot simulated biomarker data
#'
#' Default plot for `aucmat_simulation` objects. Shows a pairs plot of
#' biomarker columns with class-coloured density panels.
#'
#' @param x An `aucmat_simulation` object.
#' @param ... Ignored.
#'
#' @return A `ggplot2` object (via GGally if available) or a base-R pairs plot.
#' @importFrom graphics pairs
#' @export
plot.aucmat_simulation <- function(x, ...) {
  data  <- x$data
  p     <- ncol(data) - 1L  # last column is 'truth'
  truth <- factor(data$truth, levels = c(0, 1), labels = c("neg", "pos"))

  if (requireNamespace("GGally", quietly = TRUE) && p <= 10L) {
    df_plot <- data[, seq_len(p), drop = FALSE]
    colnames(df_plot) <- paste0(colnames(df_plot), " (",
      formatC(x$target_aucs, 2, format = "f"), ")")
    df_plot$class <- truth
    p_out <- GGally::ggpairs(df_plot, columns = seq_len(p),
      ggplot2::aes(colour = class, alpha = 0.4),
      upper = list(continuous = "cor"),
      diag  = list(continuous = "densityDiag"),
      lower = list(continuous = "points")) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste0("Simulated biomarkers (n=", x$n,
        ", prevalence=", round(x$prevalence, 3), ")"))
    return(p_out)
  }

  # Fallback: base-R pairs
  grDevices::dev.new()
  pairs(data[, seq_len(p)], col = ifelse(truth == "pos", "#B2182B", "#2166AC"),
        main = paste0("Simulated biomarkers (n=", x$n,
          ", prevalence=", round(x$prevalence, 3), ")"),
        pch = 16, cex = 0.6, gap = 0)
  invisible(x)
}

#' Plot achieved vs requested correlation heatmap
#'
#' Side-by-side heatmaps comparing the requested and achieved correlation
#' matrices from a simulation object.
#'
#' @param sim An `aucmat_simulation` object.
#' @param difference If `TRUE`, show a third panel with the difference matrix
#'   (achieved - requested). Default `FALSE`.
#'
#' @return A `ggplot2` object.
#' @export
plot_correlation_heatmap <- function(sim, difference = FALSE) {
  p <- length(sim$target_aucs)
  if (p < 2L) stop("Need at least 2 biomarkers for a correlation heatmap.")

  req <- sim$requested_correlation
  ach <- sim$achieved_correlation
  nm  <- paste0("X", seq_len(p))

  .melt_mat <- function(mat, label) {
    m <- as.data.frame(as.table(mat))
    names(m) <- c("row", "col", "value")
    m$row <- factor(m$row, levels = nm)
    m$col <- factor(m$col, levels = nm)
    m$panel <- label
    m
  }

  df <- rbind(.melt_mat(req, "Requested"), .melt_mat(ach, "Achieved"))
  if (isTRUE(difference)) {
    diff_mat <- ach - req
    df <- rbind(df, .melt_mat(diff_mat, "Difference"))
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$col, y = .data$row,
                                     fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = if (isTRUE(difference)) 0 else 0.5,
      limits = if (isTRUE(difference)) c(-1, 1) else c(-1, 1),
      name = if (isTRUE(difference)) "Delta" else "r"
    ) +
    ggplot2::geom_text(ggplot2::aes(label = formatC(.data$value, 2, format = "f")),
                       size = 3) +
    ggplot2::facet_wrap(~ panel) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL,
      title = paste0("Correlation: requested vs achieved (n=", sim$n, ")"))
}
