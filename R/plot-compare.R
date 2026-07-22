# ==============================================================================
# plot.aucmat_compare()  --  forest plot of pairwise AUC differences
# ==============================================================================

#' Forest plot of paired AUC differences
#'
#' Plots pairwise AUC differences with confidence intervals from a
#' [compare_auc()] result.
#'
#' @param x An `aucmat_compare` data.frame.
#' @param ... Ignored.
#'
#' @return A `ggplot2` object.
#' @examples
#' \donttest{
#' set.seed(42)
#' sim <- simulate_auc_matrix(n=100, prevalence=0.3,
#'   target_aucs=c(0.85,0.75,0.65), correlation=0.3, structure="exchangeable")
#' X <- as.matrix(sim$data[,1:3]); y <- sim$data$truth
#' fit <- aucmat(X, y, ci="none")
#' cmp <- compare_auc(fit, X, y, reference="X1")
#' plot(cmp)
#' }
#' @export
plot.aucmat_compare <- function(x, ...) {
  df <- as.data.frame(x)
  df <- df[!is.na(df$auc_diff), , drop = FALSE]
  if (nrow(df) == 0L) stop("No valid comparisons to plot.")

  df$pair <- paste0(df$biomarker_a, " - ", df$biomarker_b)
  has_ci <- "conf_low" %in% names(df) && any(!is.na(df$conf_low))

  # Order by auc_diff
  df$pair <- factor(df$pair, levels = df$pair[order(df$auc_diff)])

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$auc_diff, y = .data$pair)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50", alpha = 0.6) +
    ggplot2::geom_point(size = 2.5, colour = "#2c7da0") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "AUC difference", y = NULL,
      title = paste0("Paired AUC comparisons (", df$hypothesis[1],
        if (!is.na(df$margin[1]) && df$margin[1] > 0)
          paste0(", margin=", df$margin[1]) else "", ")"))

  if (has_ci) {
    df$conf_low[!is.finite(df$conf_low)] <- NA_real_
    df$conf_high[!is.finite(df$conf_high)] <- NA_real_
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$conf_low, xmax = .data$conf_high),
      height = 0.25, colour = "grey40", alpha = 0.7)
  }

  # Add non-inferiority/equivalence margin lines
  margin <- df$margin[1]
  if (!is.na(margin) && margin > 0) {
    hyp <- df$hypothesis[1]
    if (hyp == "noninferiority") {
      p <- p + ggplot2::geom_vline(xintercept = -margin, linetype = "dotted",
                                    colour = "#B2182B", alpha = 0.5)
    } else if (hyp == "equivalence") {
      p <- p + ggplot2::geom_vline(xintercept = c(-margin, margin),
        linetype = "dotted", colour = "#B2182B", alpha = 0.5)
    }
  }

  if (!is.null(attr(x, "selection_status")) &&
      attr(x, "selection_status") == "same_data") {
    p <- p + ggplot2::labs(caption = "NOTE: same-data selection  --  inference is exploratory")
  }

  p
}
