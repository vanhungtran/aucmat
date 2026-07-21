# ==============================================================================
# Visualization Layer
#
# Five focused plotting functions returning ggplot2 objects.
# Every function labels only a limited subset of features in wide-matrix plots
# to avoid unreadable overplotting.
# ==============================================================================

#' Rank plot: ordered discrimination strengths
#'
#' @param fit An `aucmat_screen` object.
#' @param n_label Number of top biomarkers to label.  Default 20.
#' @param show_ci If `TRUE` and CIs are present, add error bars.
#'
#' @return A `ggplot2` object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_text theme_minimal labs
#' @importFrom utils head
plot_auc_rank <- function(fit, n_label = 20L, show_ci = TRUE) {
  df <- fit$results
  df <- df[order(df$rank), ]
  df$label <- ifelse(seq_len(nrow(df)) <= n_label, df$biomarker, "")
  df$col  <- ifelse(df$effect_direction == "higher_in_positive", "#2166AC",
             ifelse(df$effect_direction == "lower_in_positive", "#B2182B",
                    "grey50"))
  df$col[is.na(df$effect_direction)] <- "grey70"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$auc_strength)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$effect_direction),
                        size = 1.2, alpha = 0.7) +
    ggplot2::scale_colour_manual(
      values = c(higher_in_positive = "#2166AC", lower_in_positive = "#B2182B"),
      na.value = "grey70", name = "Direction"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Rank", y = "Discrimination strength (AUC strength)",
                  title = "Biomarker discrimination strength by rank")

  if (show_ci && "conf_low" %in% names(df) && any(!is.na(df$conf_low))) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
      alpha = 0.15, width = 0
    )
  }

  if (any(df$label != "")) {
    p <- p + ggrepel::geom_text_repel(
      data = df[df$label != "", ],
      ggplot2::aes(label = .data$label), size = 2.8, max.overlaps = 30
    )
  }

  p
}

#' Volcano plot: effect magnitude against statistical evidence
#'
#' @param fit An `aucmat_screen` object.
#' @param n_label Number of top biomarkers to label.  Default 20.
#' @param q_cutoff Highlight biomarkers with q-value below this threshold.
#'
#' @return A `ggplot2` object.
#' @export
plot_auc_volcano <- function(fit, n_label = 20L, q_cutoff = 0.05) {
  df <- fit$results
  df$neg_log10_q <- -log10(df$q_value)
  df$neg_log10_q[!is.finite(df$neg_log10_q)] <- NA_real_

  df$effect <- abs(df$auc_raw - 0.5)
  df$sig <- !is.na(df$q_value) & df$q_value < q_cutoff

  top_idx <- head(order(df$auc_strength, decreasing = TRUE), n_label)
  df$label <- ""
  df$label[top_idx] <- df$biomarker[top_idx]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$effect,
                                         y = .data$neg_log10_q)) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$sig), size = 1.5, alpha = 0.6
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "grey60", "TRUE" = "#B2182B"),
      labels = c("FALSE" = "NS", "TRUE" = paste0("q < ", q_cutoff)),
      name = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Effect size |AUC - 0.5|",
      y = expression(-log[10](q)),
      title = "Biomarker effect magnitude vs. statistical evidence"
    )

  if (any(df$label != "")) {
    p <- p + ggrepel::geom_text_repel(
      data = df[df$label != "", ],
      ggplot2::aes(label = .data$label), size = 2.8, max.overlaps = 30
    )
  }

  p
}

#' Forest plot: selected AUCs with confidence intervals
#'
#' @param fit An `aucmat_screen` object.
#' @param biomarkers Character vector of biomarker names to show, or `"top"`
#'   to take the top `n` by `auc_strength`.
#' @param n Number of biomarkers when `biomarkers = "top"`.  Default 20.
#'
#' @return A `ggplot2` object.
#' @export
plot_auc_forest <- function(fit, biomarkers = "top", n = 20L) {
  df <- fit$results
  if (identical(biomarkers, "top")) {
    df <- head(df[order(df$auc_strength, decreasing = TRUE), ], n)
  } else {
    df <- df[df$biomarker %in% biomarkers, ]
    if (nrow(df) == 0L) stop("No matching biomarkers found.")
  }

  df <- df[order(df$auc_raw), ]
  df$biomarker <- factor(df$biomarker, levels = df$biomarker)

  has_ci <- "conf_low" %in% names(df) && any(!is.na(df$conf_low))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$auc_raw, y = .data$biomarker)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                        colour = "grey50", alpha = 0.6) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "AUC (raw)", y = NULL,
                  title = "Selected biomarker AUCs")

  if (has_ci) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$conf_low, xmax = .data$conf_high),
      height = 0.2
    )
  }

  p
}

#' Stability plot: rank distributions and top-k probabilities
#'
#' @param stability An `aucmat_stability` object.
#' @param n_label Number of biomarkers to show.  Default 25.
#'
#' @return A `ggplot2` object.
#' @export
plot_auc_stability <- function(stability, n_label = 25L) {
  rs <- stability$rank_summary
  n_show <- min(n_label, nrow(rs))
  rs <- rs[seq_len(n_show), ]
  rs$biomarker <- factor(rs$biomarker, levels = rev(rs$biomarker))

  ggplot2::ggplot(rs, ggplot2::aes(x = .data$biomarker)) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$rank_median), size = 2
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$rank_q25, ymax = .data$rank_q75),
      width = 0.3
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = NULL, y = "Bootstrap rank (median and IQR)",
      title = "Rank stability of leading biomarkers"
    )
}

#' ROC curves for a small set of biomarkers
#'
#' Plots empirical ROC curves for deliberately selected biomarkers.
#' Designed for small sets; labelling is explicit.
#'
#' @param fit An `aucmat_screen` object (may not retain data).
#' @param X Numeric matrix (required if `fit` does not retain data).
#' @param y Binary outcome (required if `fit` does not retain data).
#' @param biomarkers Character vector of biomarker names.  Default: top 6
#'   by `auc_strength`.
#' @param show_auc Add AUC values to the legend.  Default `TRUE`.
#'
#' @return A `ggplot2` object.
#' @export
plot_roc_top <- function(fit, X = NULL, y = NULL,
                          biomarkers = NULL, show_auc = TRUE) {
  # Resolve data
  if (!is.null(fit$X) && !is.null(fit$y)) {
    X <- fit$X; y <- fit$y
  }
  if (is.null(X) || is.null(y))
    stop("X and y are required.  Re-run aucmat() with retain_data = TRUE or supply them.")

  if (is.null(biomarkers)) {
    biomarkers <- head(fit$results$biomarker, 6L)
  }
  biomarkers <- biomarkers[!is.na(biomarkers)]

  if (length(biomarkers) == 0L) stop("No biomarkers to plot.")
  if (length(biomarkers) > 12L)
    warning("More than 12 biomarkers; the plot may be hard to read.")

  X <- as.matrix(X)
  y_norm <- .normalize_binary_y(y)
  levs <- levels(y_norm)

  roc_list <- list()
  auc_vals <- numeric(length(biomarkers))
  names(auc_vals) <- biomarkers

  for (i in seq_along(biomarkers)) {
    bm <- biomarkers[i]
    if (!bm %in% colnames(X)) stop("Biomarker not found: ", bm)
    use <- !is.na(X[, bm]) & !is.na(y_norm)
    roc_list[[bm]] <- pROC::roc(
      response  = y_norm[use],
      predictor = X[use, bm],
      levels    = levs,
      direction = "auto",
      quiet     = TRUE
    )
    auc_vals[i] <- as.numeric(pROC::auc(roc_list[[bm]]))
  }

  labels <- biomarkers
  if (show_auc) {
    labels <- paste0(biomarkers, " (AUC=", formatC(auc_vals, 3, format = "f"), ")")
  }

  palette <- if (length(roc_list) <= 8L) {
    c("#2166AC", "#B2182B", "#4DAF4A", "#FF7F00",
      "#984EA3", "#A65628", "#F781BF", "#999999")[seq_along(roc_list)]
  } else {
    grDevices::rainbow(length(roc_list))
  }

  ggroc_obj <- pROC::ggroc(roc_list, legacy.axes = FALSE, size = 0.9)
  p <- ggroc_obj +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         alpha = 0.5, colour = "grey50") +
    ggplot2::labs(x = "1 - Specificity (FPR)", y = "Sensitivity",
                  title = "ROC curves for selected biomarkers") +
    ggplot2::scale_colour_manual(values = palette, labels = labels,
                                 name = "Biomarker") +
    ggplot2::theme(legend.position = "right")

  p
}
