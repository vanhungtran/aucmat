# ==============================================================================
# Precision-Recall Curve & PR-AUC
#
# PR curves are more informative than ROC when prevalence is low (<< 0.5).
# The PR-AUC (area under the precision-recall curve) complements the ROC AUC.
# ==============================================================================

#' Compute Precision-Recall AUC for a single biomarker
#'
#' Precision-Recall AUC (average precision) computed via the trapezoidal rule
#' on the empirical PR curve.  Uses linear interpolation between observed
#' thresholds.
#'
#' @param x Numeric biomarker vector.
#' @param pos,neg Logical vectors.
#'
#' @return A list with `prauc` (PR-AUC), `prevalence` (baseline).
#' @noRd
#' @keywords internal
.compute_single_prauc <- function(x, pos, neg) {
  use   <- !is.na(x)
  p_use <- pos & use; n_use <- neg & use
  np <- sum(p_use); nn <- sum(n_use)
  if (np < 2L || nn < 2L) return(list(prauc = NA_real_, prevalence = np / (np + nn)))

  thresholds <- sort(unique(x[use]), decreasing = TRUE)
  n_thresh   <- length(thresholds)
  prevalence <- np / (np + nn)

  # At each threshold, compute recall (TPR) and precision
  recall    <- numeric(n_thresh)
  precision <- numeric(n_thresh)

  for (i in seq_len(n_thresh)) {
    t <- thresholds[i]
    tp <- sum(x[p_use] >= t)
    fp <- sum(x[n_use] >= t)
    recall[i]    <- tp / np
    precision[i] <- if (tp + fp > 0) tp / (tp + fp) else 1
  }

  # Sort by recall ascending for integration
  ord <- order(recall)
  recall    <- recall[ord]
  precision <- precision[ord]

  # Add boundary points: recall=0 (precision intercept), recall=1 (precision=prevalence)
  recall    <- c(0, recall, 1)
  precision <- c(precision[1], precision, prevalence)

  # Trapezoidal integration
  prauc <- sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2)
  list(prauc = prauc, prevalence = prevalence)
}

#' Compute matrix-wide PR-AUC values
#'
#' @param X Numeric matrix.
#' @param pos,neg Logical vectors.
#' @return Data.frame with `prauc` and `prevalence` columns.
#' @noRd
.compute_matrix_prauc <- function(X, pos, neg) {
  X <- as.matrix(X)
  cn <- colnames(X)
  if (is.null(cn)) cn <- paste0("V", seq_len(ncol(X)))

  out <- lapply(seq_len(ncol(X)), function(j) {
    .compute_single_prauc(X[, j], pos, neg)
  })

  data.frame(
    biomarker  = cn,
    prauc      = vapply(out, `[[`, numeric(1L), "prauc"),
    prevalence = vapply(out, `[[`, numeric(1L), "prevalence"),
    stringsAsFactors = FALSE
  )
}

#' Plot Precision-Recall curves for selected biomarkers
#'
#' PR curves show the trade-off between precision (positive predictive value)
#' and recall (sensitivity).  More informative than ROC when the positive
#' class is rare (prevalence << 0.5).
#'
#' @param X Numeric matrix.
#' @param y Binary outcome.
#' @param biomarkers Character vector of biomarker names.  Default: top 6 by
#'   PR-AUC.
#' @param show_auc Add PR-AUC values to legend.  Default `TRUE`.
#' @param show_baseline Add horizontal line at prevalence (random classifier).
#'   Default `TRUE`.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_auc_matrix(n = 500, prevalence = 0.15,
#'   target_aucs = c(0.85, 0.75, 0.65),
#'   correlation = 0.3, structure = "exchangeable")
#' X <- as.matrix(sim$data[, 1:3])
#' y <- sim$data$truth
#' plot_auc_pr(X, y)
#' }
plot_auc_pr <- function(X, y, biomarkers = NULL, show_auc = TRUE,
                         show_baseline = TRUE) {
  X <- as.matrix(X)
  y_norm <- .normalize_binary_y(y)
  pos <- y_norm == levels(y_norm)[2L]
  neg <- !pos

  if (is.null(biomarkers)) {
    praucs <- .compute_matrix_prauc(X, pos, neg)
    praucs <- praucs[order(praucs$prauc, decreasing = TRUE, na.last = TRUE), ]
    biomarkers <- utils::head(praucs$biomarker, 6L)
  }
  biomarkers <- biomarkers[!is.na(biomarkers)]
  if (length(biomarkers) == 0L) stop("No biomarkers to plot.")

  prev <- sum(pos) / length(pos)

  # Build PR curve data for each biomarker
  curve_list <- lapply(biomarkers, function(bm) {
    x <- X[, bm]
    use <- !is.na(x)
    p_use <- pos & use; n_use <- neg & use
    np <- sum(p_use); nn <- sum(n_use)

    thresholds <- sort(unique(x[use]), decreasing = TRUE)
    recall <- precision <- numeric(length(thresholds))
    for (i in seq_along(thresholds)) {
      tp <- sum(x[p_use] >= thresholds[i])
      fp <- sum(x[n_use] >= thresholds[i])
      recall[i]    <- tp / np
      precision[i] <- if (tp + fp > 0) tp / (tp + fp) else 1
    }
    ord <- order(recall)
    data.frame(biomarker = bm, recall = c(0, recall[ord], 1),
               precision = c(precision[ord][1], precision[ord], prev))
  })

  df <- do.call(rbind, curve_list)

  # Add PR-AUC to legend labels
  if (show_auc) {
    praucs <- vapply(biomarkers, function(bm) {
      .compute_single_prauc(X[, bm], pos, neg)$prauc
    }, numeric(1L))
    df$biomarker <- factor(df$biomarker, levels = biomarkers)
    levels(df$biomarker) <- paste0(biomarkers, " (AUC-PR=",
      formatC(praucs, 3, format = "f"), ")")
  }

  palette <- c("#2166AC", "#B2182B", "#4DAF4A", "#FF7F00",
               "#984EA3", "#A65628", "#F781BF", "#999999")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$recall, y = .data$precision,
    colour = .data$biomarker)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = palette, name = "Biomarker") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Recall (Sensitivity)", y = "Precision (PPV)",
      title = "Precision-Recall curves for selected biomarkers")

  if (show_baseline) {
    p <- p + ggplot2::geom_hline(yintercept = prev, linetype = "dashed",
      colour = "grey50", alpha = 0.5) +
      ggplot2::annotate("text", x = 0.02, y = prev + 0.03,
        label = paste0("prevalence = ", formatC(prev, 3, format = "f")),
        hjust = 0, size = 3, colour = "grey50")
  }

  p
}

# ---- Speed-optimized matrix AUC (vectorized rank-sum) ----

#' Vectorized single-AUC computation via rank-sum
#'
#' Faster than the per-column loop for many biomarkers by using
#' `colSums()` on a precomputed rank matrix.
#'
#' @param X Numeric matrix (n x p), no missing values expected.
#' @param pos,neg Logical vectors.
#'
#' @return Numeric vector of length p with raw AUC values.
#' @noRd
#' @keywords internal
.matrix_auc_fast <- function(X, pos, neg) {
  n_pos <- sum(pos)
  n_neg <- sum(neg)
  n_total <- nrow(X)

  # Rank each column (ties = "average")
  # Use apply with rank — still per-column but cleaner
  # A truly vectorized approach uses colRanks() from matrixStats if available
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    R <- matrixStats::colRanks(X, ties.method = "average")
  } else {
    R <- apply(X, 2, rank, ties.method = "average")
  }

  # For each column: U = sum(ranks[pos]) - n_pos*(n_pos+1)/2
  # auc = U / (n_pos * n_neg)
  pos_idx <- which(pos)
  U <- colSums(R[pos_idx, , drop = FALSE]) - n_pos * (n_pos + 1) / 2
  as.numeric(U / (n_pos * n_neg))
}
