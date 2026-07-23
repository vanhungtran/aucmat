# ==============================================================================
# fit_auc_panel() — Build and validate multivariable biomarker panels
#
# Combines biomarkers into a single score via one of five methods:
#   ridge / lasso / elasticnet / logistic  -- penalized/plain logistic
#     regression (glmnet), with honest *nested* cross-validation: every
#     tuning step (penalty selection) and preprocessing step (centering,
#     scaling) is re-run on the training portion of each outer fold only.
#   su_liu     -- closed-form linear combination maximizing AUC under the
#     multivariate-normal (binormal) model (Su & Liu, 1993): w = Sigma^-1
#     (mu_pos - mu_neg), with shrinkage toward the diagonal for stability
#     when biomarkers are collinear or p approaches n.
#   unweighted -- equal-weight sum of direction-aligned, standardized
#     biomarkers; an assumption-light baseline with no fitted parameters.
# ==============================================================================

#' Fit and validate a multivariable biomarker panel
#'
#' Builds a combined biomarker score and reports both an in-sample AUC and
#' an honest, nested cross-validated AUC. "Nested" means every step that
#' looks at the outcome -- penalty selection for the penalized methods, and
#' centering/scaling for all methods -- is re-derived from the training
#' portion of each outer fold alone, never from data that includes the held
#' -out test rows. This prevents the optimistic bias that comes from tuning
#' on the full data and only refitting coefficients per fold.
#'
#' @param X Numeric matrix (n x p).
#' @param y Binary outcome.
#' @param positive Optional positive class label.
#' @param method `"ridge"` (default), `"lasso"`, `"elasticnet"`, `"logistic"`
#'   (unpenalized), `"su_liu"` (closed-form AUC-maximizing linear
#'   combination under the binormal model), or `"unweighted"` (equal-weight
#'   sum of direction-aligned standardized biomarkers; no fitting).
#' @param alpha Elastic-net mixing: 0 = ridge, 1 = lasso. Only used when
#'   `method = "elasticnet"`; defaults to 0.5 (a true ridge/lasso blend)
#'   when not supplied. Ignored for other methods.
#' @param shrinkage Shrinkage intensity in `[0, 1]` toward the diagonal
#'   covariance, used only by `method = "su_liu"`. `0` uses the raw pooled
#'   covariance (unstable or infeasible once `p` approaches `n`); `1`
#'   shrinks all the way to a diagonal (independence) covariance. Default
#'   0.1.
#' @param n_folds Outer CV folds for the honest AUC (0 = no CV, uses
#'   training AUC only). Default 5.
#' @param seed Optional seed for CV fold-assignment reproducibility.
#' @param standardize Center/scale predictors before fitting. Default
#'   `TRUE`. Centering and scaling are always computed from the relevant
#'   training data only (the full data for the final model; each outer
#'   training fold for the honest CV score), never from held-out rows.
#'
#' @return A list of class `aucmat_panel` with components:
#'   `coefficients`, `auc_train`, `auc_cv`, `fitted_model` (final-model
#'   fit object together with the centering/scaling used to predict on new
#'   data via [predict.aucmat_panel()]), `predictions`, `settings`,
#'   `sample_summary`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
#'   target_aucs = c(0.85, 0.75, 0.65, 0.60),
#'   correlation = 0.3, structure = "exchangeable")
#' X <- as.matrix(sim$data[, 1:4]); y <- sim$data$truth
#'
#' panel <- fit_auc_panel(X, y, method = "ridge", n_folds = 3, seed = 42)
#' print(panel)
#'
#' # Closed-form AUC-maximizing combination -- no glmnet dependency
#' panel_su <- fit_auc_panel(X, y, method = "su_liu", n_folds = 3, seed = 42)
#'
#' # Assumption-light equal-weight baseline
#' panel_uw <- fit_auc_panel(X, y, method = "unweighted", n_folds = 3, seed = 42)
#' }
fit_auc_panel <- function(X, y, positive = NULL,
                           method = c("ridge", "lasso", "elasticnet", "logistic",
                                      "su_liu", "unweighted"),
                           alpha = NULL, shrinkage = 0.1,
                           n_folds = 5L, seed = NULL,
                           standardize = TRUE) {
  method <- match.arg(method)
  if (method == "ridge")      alpha <- 0
  if (method == "lasso")      alpha <- 1
  if (method == "elasticnet") alpha <- if (is.null(alpha)) 0.5 else alpha
  if (is.null(alpha)) alpha <- 0  # unused by logistic / su_liu / unweighted

  if (method == "su_liu") {
    if (!is.numeric(shrinkage) || length(shrinkage) != 1L || is.na(shrinkage) ||
        shrinkage < 0 || shrinkage > 1) {
      stop("shrinkage must be a single number in [0, 1].")
    }
  }

  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  p <- ncol(X)
  if (p < 2L) stop("Need at least 2 biomarkers for a panel.")

  y_norm <- .normalize_binary_y(y, positive)
  y_num  <- .binary_factor_to_numeric(y_norm)

  use <- stats::complete.cases(X) & !is.na(y_num)
  X <- X[use, , drop = FALSE]
  y_num <- y_num[use]

  cn <- colnames(X)
  if (is.null(cn)) cn <- paste0("X", seq_len(p))
  colnames(X) <- cn

  penalized <- method %in% c("ridge", "lasso", "elasticnet")
  if (penalized && !requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for method = '", method, "'. ",
         "Install it with install.packages('glmnet') or use method = ",
         "'logistic', 'su_liu', or 'unweighted' instead.")
  }

  # ---- Fit coefficients + preprocessing on a (training) block. Used both
  # ---- for the final full-data model and for every honest outer CV fold,
  # ---- so nothing here ever sees the held-out test rows. ----
  .fit_block <- function(Xtr_raw, ytr) {
    if (isTRUE(standardize)) {
      m <- colMeans(Xtr_raw)
      s <- apply(Xtr_raw, 2, stats::sd)
      s[s == 0] <- 1
    } else {
      m <- rep(0, ncol(Xtr_raw)); s <- rep(1, ncol(Xtr_raw))
    }
    Xtr <- .scale_block(Xtr_raw, m, s)

    coefs <- switch(method,
      logistic = {
        df_tr <- as.data.frame(Xtr)
        df_tr$.y <- ytr
        fit_k <- stats::glm(.y ~ ., data = df_tr, family = stats::binomial())
        stats::coef(fit_k)[-1]
      },
      su_liu     = .fit_su_liu(Xtr, ytr, shrinkage = shrinkage),
      unweighted = .fit_unweighted(Xtr, ytr),
      {
        nf <- .safe_glmnet_nfolds(ytr)
        cv_k <- glmnet::cv.glmnet(Xtr, ytr, family = "binomial",
          alpha = alpha, standardize = FALSE, nfolds = nf)
        as.vector(stats::coef(cv_k, s = "lambda.min"))[-1]
      }
    )
    names(coefs) <- colnames(Xtr_raw)
    list(coefs = coefs, means = m, sds = s)
  }

  # ---- Final full-data fit + honest nested outer CV, sharing one seeded
  # ---- block. This matters even for the final fit: the penalized methods'
  # ---- own internal cv.glmnet() call for penalty selection draws random
  # ---- fold assignments, so without this the "final" coefficients would
  # ---- neither be reproducible under `seed` nor leave the caller's global
  # ---- RNG state untouched. ----
  final <- NULL
  auc_cv <- NA_real_
  oof_scores <- rep(NA_real_, length(y_num))

  with_seed(seed, {
    final <- .fit_block(X, y_num)

    if (n_folds >= 2L) {
      pos_idx <- which(y_num == 1)
      neg_idx <- which(y_num == 0)
      pos_folds <- sample(rep(seq_len(n_folds), length.out = length(pos_idx)))
      neg_folds <- sample(rep(seq_len(n_folds), length.out = length(neg_idx)))
      folds <- integer(length(y_num))
      folds[pos_idx] <- pos_folds
      folds[neg_idx] <- neg_folds

      for (k in seq_len(n_folds)) {
        test_idx  <- which(folds == k)
        train_idx <- which(folds != k)
        fold_fit <- .fit_block(X[train_idx, , drop = FALSE], y_num[train_idx])
        Xte <- .scale_block(X[test_idx, , drop = FALSE], fold_fit$means, fold_fit$sds)
        oof_scores[test_idx] <- as.vector(Xte %*% fold_fit$coefs)
      }
    }
  })

  coefs <- final$coefs
  X_scaled_full <- .scale_block(X, final$means, final$sds)
  score_train <- as.vector(X_scaled_full %*% coefs)
  auc_train <- .compute_single_auc(score_train, y_num == 1, y_num == 0)$auc_raw
  if (n_folds >= 2L) {
    auc_cv <- .compute_single_auc(oof_scores, y_num == 1, y_num == 0)$auc_raw
  }

  out <- list(
    coefficients  = coefs,
    auc_train     = auc_train,
    auc_cv        = auc_cv,
    fitted_model  = final,
    predictions   = data.frame(
      score_train = score_train,
      score_cv    = if (n_folds >= 2L) oof_scores else rep(NA_real_, length(y_num))
    ),
    settings = list(
      method      = method,
      alpha       = if (penalized) alpha else NA_real_,
      shrinkage   = if (method == "su_liu") shrinkage else NA_real_,
      n_folds     = n_folds,
      standardize = standardize,
      seed        = seed
    ),
    sample_summary = list(
      n_total    = length(y_num),
      n_positive = sum(y_num == 1),
      n_negative = sum(y_num == 0)
    )
  )
  class(out) <- "aucmat_panel"
  out
}

# ---- Internal: center/scale a block with supplied (training-derived)
# ---- means and sds -- never with the block's own statistics at test time ----
.scale_block <- function(M, means, sds) {
  out <- sweep(M, 2, means, "-")
  sweep(out, 2, sds, "/")
}

# ---- Internal: a defensive inner fold count for cv.glmnet, bounded by the
# ---- minority class size of the training block being tuned on ----
.safe_glmnet_nfolds <- function(ytr) {
  min_class <- min(sum(ytr == 1), sum(ytr == 0))
  as.integer(min(10L, max(3L, min_class)))
}

# ---- Su & Liu (1993): closed-form linear combination maximizing AUC under
# ---- the common-covariance multivariate normal (binormal) model.
# ---- w = Sigma^-1 (mu_pos - mu_neg); pooled covariance is shrunk toward
# ---- its diagonal for stability when biomarkers are collinear or p is
# ---- close to n. ----
.fit_su_liu <- function(Xtr, ytr, shrinkage = 0.1) {
  pos <- ytr == 1; neg <- ytr == 0
  np <- sum(pos); nn <- sum(neg)
  if (np < 2L || nn < 2L) {
    stop("Need at least 2 subjects per class to fit method = 'su_liu'.")
  }
  mu_pos <- colMeans(Xtr[pos, , drop = FALSE])
  mu_neg <- colMeans(Xtr[neg, , drop = FALSE])
  Sigma_pos <- stats::cov(Xtr[pos, , drop = FALSE])
  Sigma_neg <- stats::cov(Xtr[neg, , drop = FALSE])
  Sigma_pooled <- ((np - 1) * Sigma_pos + (nn - 1) * Sigma_neg) / (np + nn - 2)
  Sigma_pooled[!is.finite(Sigma_pooled)] <- 0

  target <- diag(diag(Sigma_pooled), nrow = nrow(Sigma_pooled))
  Sigma_reg <- (1 - shrinkage) * Sigma_pooled + shrinkage * target

  if (!is.finite(rcond(Sigma_reg)) || rcond(Sigma_reg) < 1e-10) {
    stop("Covariance matrix is numerically singular even after shrinkage = ",
         shrinkage, ". Increase 'shrinkage' (closer to 1) or reduce the ",
         "number of biomarkers in the panel.")
  }

  as.vector(solve(Sigma_reg, mu_pos - mu_neg))
}

# ---- Unweighted (equal-weight) direction-aligned z-score sum: a
# ---- fitting-free baseline. Direction per biomarker is the sign of its
# ---- own class-mean difference on the training block. ----
.fit_unweighted <- function(Xtr, ytr) {
  pos <- ytr == 1; neg <- ytr == 0
  p <- ncol(Xtr)
  direction <- vapply(seq_len(p), function(j) {
    d <- mean(Xtr[pos, j]) - mean(Xtr[neg, j])
    if (!is.finite(d) || d == 0) 1 else sign(d)
  }, numeric(1))
  direction / p
}

#' @export
print.aucmat_panel <- function(x, ...) {
  cat("<aucmat_panel>  ", x$settings$method, "\n", sep = "")
  if (!is.na(x$settings$shrinkage))
    cat("  Shrinkage:      ", x$settings$shrinkage, "\n", sep = "")
  cat("  AUC (training): ", round(x$auc_train, 4), "\n", sep = "")
  if (!is.na(x$auc_cv)) {
    cat("  AUC (CV):       ", round(x$auc_cv, 4), "\n", sep = "")
    cat("  Optimism:       ", round(x$auc_train - x$auc_cv, 4), "\n", sep = "")
  }
  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  coefs <- coefs[order(abs(coefs), decreasing = TRUE)]
  n_show <- min(10L, length(coefs))
  print(round(coefs[seq_len(n_show)], 4))
  if (length(coefs) > n_show)
    cat("... ", length(coefs) - n_show, " more\n", sep = "")
  invisible(x)
}

#' @export
plot.aucmat_panel <- function(x, ...) {
  coefs <- x$coefficients
  coefs <- coefs[order(coefs)]
  nm <- names(coefs)
  if (is.null(nm)) nm <- paste0("X", seq_along(coefs))

  df <- data.frame(
    biomarker = factor(nm, levels = nm),
    coefficient = unname(coefs),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$coefficient, y = .data$biomarker)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50", alpha = 0.5) +
    ggplot2::geom_point(size = 2.5,
      colour = ifelse(df$coefficient >= 0, "#2166AC", "#B2182B")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Coefficient", y = NULL,
      title = paste0("Panel coefficients (", x$settings$method, ")"),
      subtitle = paste0("AUC train=", round(x$auc_train, 3),
        if (!is.na(x$auc_cv)) paste0(", AUC CV=", round(x$auc_cv, 3)) else "")
    )
}

#' Predict panel scores for new subjects
#'
#' Applies a fitted panel's coefficients to new data, using the
#' centering/scaling stored from the original (training) fit.
#'
#' @param object An `aucmat_panel` object from [fit_auc_panel()].
#' @param newdata Numeric matrix or data.frame containing at least the
#'   columns used to fit `object`.
#' @param ... Ignored.
#'
#' @return A numeric vector of panel scores, one per row of `newdata`.
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
#'   target_aucs = c(0.85, 0.75, 0.65, 0.60),
#'   correlation = 0.3, structure = "exchangeable")
#' X <- as.matrix(sim$data[, 1:4]); y <- sim$data$truth
#' panel <- fit_auc_panel(X[1:200, ], y[1:200], method = "ridge", n_folds = 0)
#' predict(panel, X[201:300, ])
#' }
#' @export
predict.aucmat_panel <- function(object, newdata, ...) {
  newdata <- as.matrix(newdata)
  cn <- names(object$coefficients)
  if (!all(cn %in% colnames(newdata))) {
    stop("newdata is missing columns used to fit the panel: ",
         paste(setdiff(cn, colnames(newdata)), collapse = ", "))
  }
  newdata <- newdata[, cn, drop = FALSE]
  Xs <- .scale_block(newdata, object$fitted_model$means, object$fitted_model$sds)
  as.vector(Xs %*% object$coefficients)
}

# ==============================================================================
# compare_panel_auc() — Does the combined panel actually beat its parts?
# ==============================================================================

#' Test whether a combined panel improves on its individual biomarkers
#'
#' Treats the panel's honest cross-validated score as one more "biomarker"
#' and reuses [compare_auc()]'s DeLong or bootstrap paired-comparison engine
#' to test it against each individual biomarker used to build the panel, on
#' the same subjects. Answers "did combining help?" with a confidence
#' interval and p-value rather than a side-by-side comparison of two point
#' estimates.
#'
#' @param panel An `aucmat_panel` object from [fit_auc_panel()].
#' @param X,y The same data (and outcome) used to fit `panel`.
#' @param positive Optional positive class label.
#' @param biomarkers Character vector of individual biomarkers to compare
#'   against. Default: all biomarkers used in the panel.
#' @param use `"cv"` (default) compares the honest cross-validated panel
#'   score; `"train"` compares the in-sample (optimistic) score. Errors if
#'   `"cv"` is requested but `panel` was fit with `n_folds < 2`.
#' @param ... Passed to [compare_auc()] (`hypothesis`, `alternative`,
#'   `margin`, `conf_level`, `method`, `boot_n`, `seed`, `adjust`).
#'
#' @return A data.frame of class `aucmat_compare` (see [compare_auc()]).
#'   `biomarker_a` is always the panel score.
#' @examples
#' \donttest{
#' set.seed(42)
#' sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
#'   target_aucs = c(0.85, 0.75, 0.65, 0.60),
#'   correlation = 0.3, structure = "exchangeable")
#' X <- as.matrix(sim$data[, 1:4]); y <- sim$data$truth
#' panel <- fit_auc_panel(X, y, method = "ridge", n_folds = 5, seed = 1)
#' compare_panel_auc(panel, X, y)
#' }
#' @export
compare_panel_auc <- function(panel, X, y, positive = NULL,
                               biomarkers = NULL, use = c("cv", "train"),
                               ...) {
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
    stop("panel has no cross-validated score (it was fit with n_folds < 2). ",
         "Either refit with n_folds >= 2, or pass use = 'train' explicitly ",
         "(in-sample scores give optimistic, not honest, comparisons).")
  }
  if (length(score) != nrow(Xu)) {
    stop("X/y do not match the data used to fit this panel ",
         "(row count after complete-case filtering differs). ",
         "Pass the same X and y used in fit_auc_panel().")
  }

  biomarkers <- biomarkers %||% used_cols
  biomarkers <- intersect(biomarkers, used_cols)
  if (length(biomarkers) < 1L) stop("No matching biomarkers to compare against.")

  panel_name <- "PANEL_SCORE"
  while (panel_name %in% colnames(Xu)) panel_name <- paste0(".", panel_name)

  combined <- cbind(Xu, score)
  colnames(combined) <- c(colnames(Xu), panel_name)

  compare_auc(fit = NULL, X = combined, y = y_kept,
              reference = panel_name, biomarkers = biomarkers, ...)
}
