# ==============================================================================
# plot_roc_panel()  --  ROC curve for a combined panel vs. its components
#
# The non-deprecated, honestly-scored replacement for what
# plot_roc_with_combos() used to provide informally: a way to see whether
# combining biomarkers actually improved the ROC curve, not just the AUC
# number.
# ==============================================================================

#' ROC curve for a fitted biomarker panel vs. its individual components
#'
#' Plots the empirical ROC curve of a combined panel score -- the honest
#' out-of-fold score by default -- together with the ROC curves of its
#' leading individual biomarkers, for a direct visual read on whether
#' combining biomarkers improved discrimination.
#'
#' @param panel An `aucmat_panel` object from [fit_auc_panel()].
#' @param X Numeric matrix -- the same data used to fit `panel` (or a
#'   superset of its columns).
#' @param y Binary outcome -- the same outcome used to fit `panel`.
#' @param positive Optional positive class label (see [aucmat()]).
#' @param biomarkers Character vector of individual biomarkers to overlay.
#'   Default: the top `n_top` panel components by their own univariate AUC.
#' @param n_top Number of individual biomarkers to overlay when
#'   `biomarkers` is not supplied. Default 3.
#' @param use `"cv"` (default) plots the honest cross-validated panel
#'   score; `"train"` plots the in-sample (optimistic) score. Falls back
#'   to `"train"` with a warning if `panel` was fit with `n_folds < 2`.
#' @param add_ci Add bootstrap sensitivity confidence ribbons. Default
#'   `TRUE` -- set `FALSE` for a faster, ribbon-free plot.
#' @param boot_n Bootstrap replicates for CI ribbons. Default 500.
#' @param seed Optional seed for CI bootstrap reproducibility.
#'
#' @return A `ggplot2` object.
#' @examples
#' \donttest{
#' set.seed(42)
#' sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
#'   target_aucs = c(0.85, 0.75, 0.65, 0.60),
#'   correlation = 0.3, structure = "exchangeable")
#' X <- as.matrix(sim$data[, 1:4]); y <- sim$data$truth
#' panel <- fit_auc_panel(X, y, method = "ridge", n_folds = 5, seed = 1)
#' plot_roc_panel(panel, X, y)
#' }
#' @export
plot_roc_panel <- function(panel, X, y, positive = NULL,
                            biomarkers = NULL, n_top = 3L,
                            use = c("cv", "train"),
                            add_ci = TRUE, boot_n = 500, seed = NULL) {
  if (!inherits(panel, "aucmat_panel"))
    stop("panel must be an 'aucmat_panel' object from fit_auc_panel().")
  use <- match.arg(use)

  used_cols <- names(panel$coefficients)
  X <- as.matrix(X)
  if (!all(used_cols %in% colnames(X))) {
    stop("X is missing columns used to fit the panel: ",
         paste(setdiff(used_cols, colnames(X)), collapse = ", "))
  }

  y_norm <- .normalize_binary_y(y, positive)
  Xu <- X[, used_cols, drop = FALSE]
  keep <- stats::complete.cases(Xu) & !is.na(y_norm)
  Xu <- Xu[keep, , drop = FALSE]
  y_kept <- droplevels(y_norm[keep])

  score <- if (use == "cv") panel$predictions$score_cv else panel$predictions$score_train
  if (use == "cv" && all(is.na(score))) {
    warning("panel has no cross-validated score (it was fit with n_folds < 2); ",
            "falling back to the in-sample (optimistic) score.")
    score <- panel$predictions$score_train
    use <- "train"
  }
  if (length(score) != nrow(Xu)) {
    stop("X/y do not match the data used to fit this panel ",
         "(row count after complete-case filtering differs). ",
         "Pass the same X and y used in fit_auc_panel().")
  }

  levs <- levels(y_kept)

  if (is.null(biomarkers)) {
    pos <- y_kept == levs[2L]; neg <- !pos
    aucs <- vapply(used_cols, function(cn) {
      .compute_single_auc(Xu[, cn], pos, neg)$auc_strength
    }, numeric(1))
    n_top <- min(n_top, length(aucs))
    biomarkers <- names(sort(aucs, decreasing = TRUE))[seq_len(n_top)]
  }
  biomarkers <- intersect(biomarkers, used_cols)
  if (length(biomarkers) == 0L) stop("No matching biomarkers to overlay.")

  panel_label <- "Combined panel"

  # Individual biomarkers first, panel last so it draws on top.
  roc_list <- stats::setNames(
    lapply(biomarkers, function(bm) {
      pROC::roc(response = y_kept, predictor = Xu[, bm],
        levels = levs, direction = "auto", quiet = TRUE)
    }),
    biomarkers
  )
  roc_list[[panel_label]] <- pROC::roc(response = y_kept, predictor = score,
    levels = levs, direction = "auto", quiet = TRUE)

  auc_vals <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)), numeric(1))
  labels <- paste0(names(roc_list), " (AUC=", formatC(auc_vals, 3, format = "f"), ")")
  palette <- c(aucmat_colors$accent[seq_along(biomarkers)], aucmat_colors$dark)

  ggroc_obj <- pROC::ggroc(roc_list, legacy.axes = FALSE, size = 0.9)

  if (isTRUE(add_ci)) {
    ci_df <- with_seed(seed, {
      spec_grid <- seq(1, 0, length.out = 101)
      ci_list <- lapply(names(roc_list), function(nm) {
        ci <- try(pROC::ci.se(roc_list[[nm]], specificities = spec_grid,
          boot.n = boot_n), silent = TRUE)
        if (inherits(ci, "try-error")) return(NULL)
        spec <- as.numeric(rownames(ci))
        data.frame(name = nm, x = spec, lower = ci[, 1], upper = ci[, 3])
      })
      do.call(rbind, ci_list)
    })
    if (!is.null(ci_df) && nrow(ci_df) > 0) {
      ci_df$name <- factor(ci_df$name, levels = names(roc_list))
      ggroc_obj <- ggroc_obj +
        ggplot2::geom_ribbon(data = ci_df,
          ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper,
            fill = .data$name), alpha = 0.15, inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = palette, guide = "none")
    }
  }

  subtitle <- paste0("Panel score: ", panel$settings$method,
    if (use == "cv") " (honest, out-of-fold)" else " (in-sample -- optimistic)")

  ggroc_obj +
    theme_aucmat() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
      alpha = 0.5, colour = "grey50") +
    ggplot2::scale_colour_manual(values = palette, labels = labels, name = NULL) +
    ggplot2::labs(x = "1 - Specificity (FPR)", y = "Sensitivity",
      title = "ROC: combined panel vs. individual biomarkers",
      subtitle = subtitle) +
    ggplot2::theme(legend.position = "right")
}
