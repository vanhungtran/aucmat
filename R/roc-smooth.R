# ==============================================================================
# roc_test() + smoothed ROC in plot_roc_top()
# ==============================================================================

#' Single-biomarker ROC test against chance or another biomarker
#'
#' A focused wrapper around DeLong's test for one or two biomarkers on
#' common subjects.  Returns a clean summary suitable for reporting.
#'
#' @param x Numeric biomarker vector (or matrix column name if X supplied).
#' @param y Binary outcome.
#' @param X Optional numeric matrix.  If supplied, `x` and `x2` are column
#'   names or indices.
#' @param x2 Optional second biomarker for paired comparison.
#' @param positive Positive class label.
#' @param alternative `"two.sided"` (default), `"greater"`, or `"less"`.
#' @param method `"delong"` (default) or `"bootstrap"`.
#' @param boot_n Bootstrap replicates when `method = "bootstrap"`.
#' @param conf_level Confidence level.
#'
#' @return A list of class `aucmat_roc_test` with `auc`, `se`, `ci`,
#'   `p_value`, `n`, `method`.
#' @export
#'
#' @examples
#' set.seed(1)
#' y <- rbinom(100, 1, 0.3)
#' x <- rnorm(100) + y * 1.5
#' roc_test(x, y)
roc_test <- function(x, y, X = NULL, x2 = NULL, positive = NULL,
                      alternative = c("two.sided", "greater", "less"),
                      method = c("delong", "bootstrap"),
                      boot_n = 2000, conf_level = 0.95) {
  alternative <- match.arg(alternative)
  method      <- match.arg(method)

  y_norm <- .normalize_binary_y(y, positive)
  pos <- y_norm == levels(y_norm)[2L]
  neg <- !pos

  # Resolve x from X if needed
  if (!is.null(X)) {
    X <- as.matrix(X)
    if (is.character(x)) x <- X[, x] else x <- X[, x]
    if (!is.null(x2)) {
      if (is.character(x2)) x2 <- X[, x2] else x2 <- X[, x2]
    }
  }

  if (!is.numeric(x)) stop("x must be numeric.")

  use <- !is.na(x) & !is.na(y_norm)
  x <- x[use]; p_use <- pos[use]; n_use <- neg[use]
  np <- sum(p_use); nn <- sum(n_use)

  auc1 <- .compute_single_auc(x, p_use, n_use)

  if (!is.null(x2)) {
    # Paired comparison
    x2 <- x2[use]
    auc2 <- .compute_single_auc(x2, p_use, n_use)
    delta <- auc1$auc_raw - auc2$auc_raw

    if (method == "delong") {
      se <- .delong_diff_se(x, x2, p_use, n_use)
    } else {
      bs <- bootstrap_auc_distribution(x, p_use, n_use, boot_n)
      se <- bs$std_error
    }

    if (is.na(se) || se <= 0) {
      return(structure(list(
        auc_a = auc1$auc_raw, auc_b = auc2$auc_raw, delta = delta,
        se = NA_real_, ci = c(NA_real_, NA_real_), p_value = NA_real_,
        n = np + nn, n_pos = np, n_neg = nn,
        alternative = alternative, method = method,
        status = "inference_failed"
      ), class = "aucmat_roc_test"))
    }

    res <- .compare_superiority(delta, se, alternative, conf_level,
      "a", "b", auc1$auc_raw, auc2$auc_raw, np, nn)

    out <- list(
      auc_a = auc1$auc_raw, auc_b = auc2$auc_raw, delta = delta,
      se = se, ci = c(res$conf_low, res$conf_high),
      p_value = res$p_value,
      n = np + nn, n_pos = np, n_neg = nn,
      alternative = alternative, method = method,
      status = "ok"
    )
  } else {
    # One-sample: test against chance (0.5)
    auc_val <- auc1$auc_raw
    if (method == "delong") {
      se <- unname(.delong_variance(x, p_use, n_use)["std_error"])
    } else {
      bs <- .bootstrap_auc_se(x, p_use, n_use, boot_n)
      se <- unname(bs["std_error"])
    }

    if (is.na(se) || se <= 0) {
      return(structure(list(
        auc = auc_val, se = NA_real_,
        ci = c(NA_real_, NA_real_), p_value = NA_real_,
        n = np + nn, n_pos = np, n_neg = nn,
        alternative = alternative, method = method,
        status = "inference_failed"
      ), class = "aucmat_roc_test"))
    }

    z_stat <- (auc_val - 0.5) / se
    z_two  <- stats::qnorm((1 + conf_level) / 2)
    z_one  <- stats::qnorm(conf_level)

    if (alternative == "two.sided") {
      ci <- c(auc_val - z_two * se, auc_val + z_two * se)
      p  <- 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
    } else if (alternative == "greater") {
      ci <- c(auc_val - z_one * se, 1)
      p  <- stats::pnorm(z_stat, lower.tail = FALSE)
    } else {
      ci <- c(0, auc_val + z_one * se)
      p  <- stats::pnorm(z_stat, lower.tail = TRUE)
    }

    out <- list(
      auc = auc_val, se = se, ci = ci, p_value = p,
      n = np + nn, n_pos = np, n_neg = nn,
      alternative = alternative, method = method,
      status = "ok"
    )
  }

  class(out) <- "aucmat_roc_test"
  out
}

#' @export
print.aucmat_roc_test <- function(x, ...) {
  cat("<aucmat_roc_test>  ", x$method, "  |  ", x$alternative, "\n", sep = "")
  if (!is.null(x$auc_a)) {
    cat(sprintf("  AUC_a = %.4f, AUC_b = %.4f, delta = %.4f\n",
      x$auc_a, x$auc_b, x$delta))
  } else {
    cat(sprintf("  AUC = %.4f, SE = %.4f\n", x$auc, x$se))
  }
  cat(sprintf("  %.0f%% CI: [%.4f, %.4f]  p = %.4f\n",
    100 * if (is.null(x$alternative) || x$alternative == "two.sided") 0.95 else 0.90,
    x$ci[1], x$ci[2], x$p_value))
  cat(sprintf("  n = %d (%d + / %d -)\n", x$n, x$n_pos, x$n_neg))
  invisible(x)
}

# ==============================================================================
# Smoothed ROC curves in plot_roc_top()
# ==============================================================================

#' Plot smoothed ROC curves for selected biomarkers
#'
#' As [plot_roc_top()] but adds a binormal-smoothed ROC curve overlay
#' using [pROC::smooth()].  The smoothed curve reduces step-artifact
#' noise and is the standard for publication-quality ROC figures.
#'
#' @param fit An `aucmat_screen` object.
#' @param X Numeric matrix.
#' @param y Binary outcome.
#' @param biomarkers Character vector.  Default: top 6 by AUC.
#' @param smooth_method `"binormal"` (default) or `"density"`.
#' @param show_empirical Also plot the empirical (step) ROC curve.
#'   Default `TRUE`.
#' @param add_ci Add bootstrap CI ribbons.  Default `TRUE`.
#' @param boot_n Bootstrap replicates for CI ribbons.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
#'   target_aucs = c(0.85, 0.75, 0.65),
#'   correlation = 0.3, structure = "exchangeable")
#' fit <- aucmat(as.matrix(sim$data[, 1:3]), sim$data$truth, ci = "none")
#' plot_roc_smooth(fit, as.matrix(sim$data[, 1:3]), sim$data$truth)
#' }
plot_roc_smooth <- function(fit, X, y, biomarkers = NULL,
                             smooth_method = c("binormal", "density"),
                             show_empirical = TRUE,
                             add_ci = TRUE, boot_n = 500) {
  smooth_method <- match.arg(smooth_method)
  X <- as.matrix(X)
  y_norm <- .normalize_binary_y(y)

  if (is.null(biomarkers)) {
    biomarkers <- utils::head(fit$results$biomarker, 6L)
  }

  palette <- c("#2166AC", "#B2182B", "#4DAF4A", "#FF7F00",
               "#984EA3", "#A65628", "#F781BF", "#999999")

  # Build ROC objects (smoothed)
  roc_list <- list()
  labels   <- character(0)

  for (i in seq_along(biomarkers)) {
    bm <- biomarkers[i]
    use <- !is.na(X[, bm]) & !is.na(y_norm)
    roc_raw <- pROC::roc(y_norm[use], X[use, bm],
      levels = levels(y_norm), direction = "auto", quiet = TRUE)
    if (smooth_method == "binormal") {
      roc_sm <- try(pROC::smooth(roc_raw, method = "binormal"), silent = TRUE)
    } else {
      roc_sm <- try(pROC::smooth(roc_raw, method = "density"), silent = TRUE)
    }
    if (inherits(roc_sm, "try-error")) roc_sm <- roc_raw

    roc_list[[bm]] <- roc_sm
    auc_val <- as.numeric(pROC::auc(roc_raw))
    labels <- c(labels, paste0(bm, " (AUC=", formatC(auc_val, 3, format = "f"), ")"))
  }

  p <- pROC::ggroc(roc_list, legacy.axes = FALSE, size = 0.9)
  # pROC >= 1.19 returns S7 objects; strip class to ensure ggplot2 compatibility
  class(p) <- setdiff(class(p), c("ggroc", "ggroc_multiclass", "S7_object"))
  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
      alpha = 0.5, colour = "grey50") +
    ggplot2::labs(x = "1 - Specificity (FPR)", y = "Sensitivity",
      title = paste0("Smoothed ROC curves (", smooth_method, ")")) +
    ggplot2::scale_colour_manual(values = palette, labels = labels,
      name = "Biomarker")

  # Overlay empirical step curves
  if (show_empirical) {
    roc_emp <- list()
    for (bm in biomarkers) {
      use <- !is.na(X[, bm]) & !is.na(y_norm)
      roc_emp[[bm]] <- pROC::roc(y_norm[use], X[use, bm],
        levels = levels(y_norm), direction = "auto", quiet = TRUE)
    }
    p_emp <- try(pROC::ggroc(roc_emp, legacy.axes = FALSE, size = 0.3,
      alpha = 0.35, linetype = "dashed"), silent = TRUE)
	    if (!inherits(p_emp, "try-error")) p <- p + p_emp
  }

  # CI ribbons on the smoothed curves
  if (add_ci) {
    spec_grid <- seq(1, 0, length.out = 101)
    ci_list <- lapply(seq_along(roc_list), function(i) {
      ci <- try(pROC::ci.se(roc_list[[i]], specificities = spec_grid,
        boot.n = boot_n), silent = TRUE)
      if (inherits(ci, "try-error")) return(NULL)
      spec <- as.numeric(rownames(ci))
      data.frame(name = names(roc_list)[i], x = spec,
        lower = ci[, 1], upper = ci[, 3])
    })
    ci_df <- do.call(rbind, ci_list)
    if (!is.null(ci_df) && nrow(ci_df) > 0) {
      ci_df$name <- factor(ci_df$name, levels = names(roc_list))
      p <- p + ggplot2::geom_ribbon(data = ci_df,
        ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper,
          fill = .data$name), alpha = 0.15, inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = palette, guide = "none")
    }
  }

  p
}
