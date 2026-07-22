# ==============================================================================
# fit_auc_panel() — Build and validate multivariable biomarker panels
#
# Fits penalized logistic regression models (ridge, lasso, elastic-net) to
# combine biomarkers into a single score, with honest out-of-fold or
# holdout validation.
# ==============================================================================

#' Fit and validate a multivariable biomarker panel
#'
#' Builds a combined biomarker score via penalized logistic regression
#' (ridge, lasso, or elastic-net) with optional cross-validation for
#' honest performance estimation.  Returns the fitted model, coefficients,
#' and performance metrics on both training and validation data.
#'
#' @param X Numeric matrix (n × p).
#' @param y Binary outcome.
#' @param positive Optional positive class label.
#' @param method `"ridge"` (default), `"lasso"`, `"elasticnet"`, or
#'   `"logistic"` (unpenalized).
#' @param alpha Elastic-net mixing: 0 = ridge, 1 = lasso.  Default 0
#'   (ridge). Only used when `method = "elasticnet"`.
#' @param n_folds CV folds for honest AUC (0 = no CV, uses training AUC).
#'   Default 5.
#' @param seed Optional seed for CV reproducibility.
#' @param standardize Standardize predictors before fitting.  Default `TRUE`.
#'
#' @return A list of class `aucmat_panel` with components:
#'   `coefficients`, `auc_train`, `auc_cv`, `fitted_model`, `predictions`,
#'   `settings`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
#'   target_aucs = c(0.85, 0.75, 0.65, 0.60),
#'   correlation = 0.3, structure = "exchangeable")
#' panel <- fit_auc_panel(as.matrix(sim$data[, 1:4]), sim$data$truth,
#'   method = "ridge", n_folds = 3, seed = 42)
#' print(panel)
#' }
fit_auc_panel <- function(X, y, positive = NULL,
                           method = c("ridge", "lasso", "elasticnet", "logistic"),
                           alpha = 0, n_folds = 5L, seed = NULL,
                           standardize = TRUE) {
  method <- match.arg(method)
  if (method == "lasso")  alpha <- 1
  if (method == "ridge")  alpha <- 0

  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  p <- ncol(X)
  if (p < 2L) stop("Need at least 2 biomarkers for a panel.")

  y_norm <- .normalize_binary_y(y, positive)
  y_num  <- .binary_factor_to_numeric(y_norm)
  pos    <- y_norm == levels(y_norm)[2L]
  neg    <- !pos

  use <- stats::complete.cases(X) & !is.na(y_num)
  X <- X[use, , drop = FALSE]
  y_num <- y_num[use]

  if (isTRUE(standardize)) {
    X_means <- colMeans(X)
    X_sds   <- apply(X, 2, stats::sd)
    X_sds[X_sds == 0] <- 1
    X_scaled <- scale(X, center = X_means, scale = X_sds)
  } else {
    X_scaled <- X
    X_means <- rep(0, p)
    X_sds   <- rep(1, p)
  }

  cn <- colnames(X)
  if (is.null(cn)) cn <- paste0("X", seq_len(p))

  # Fit model
  if (method == "logistic") {
    df_fit <- as.data.frame(X_scaled)
    df_fit$y <- y_num
    fit <- stats::glm(y ~ ., data = df_fit, family = stats::binomial())
    coefs <- stats::coef(fit)[-1]
    names(coefs) <- cn
  } else {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package 'glmnet' is required for penalized methods. ",
           "Install it with install.packages('glmnet') or use method='logistic'.")
    }
    fit <- glmnet::cv.glmnet(X_scaled, y_num, family = "binomial",
      alpha = alpha, standardize = FALSE)
    coefs <- as.vector(stats::coef(fit, s = "lambda.min"))[-1]
    names(coefs) <- cn
  }

  # Training AUC
  score_train <- as.vector(X_scaled %*% coefs)
  auc_train <- .compute_single_auc(score_train,
    y_num == 1, y_num == 0)$auc_raw

  # Cross-validated AUC
  auc_cv <- NA_real_
  if (n_folds >= 2L) {
    if (!is.null(seed)) set.seed(seed)
    pos_idx <- which(y_num == 1)
    neg_idx <- which(y_num == 0)
    pos_folds <- sample(rep(seq_len(n_folds), length.out = length(pos_idx)))
    neg_folds <- sample(rep(seq_len(n_folds), length.out = length(neg_idx)))
    folds <- integer(length(y_num))
    folds[pos_idx] <- pos_folds
    folds[neg_idx] <- neg_folds

    oof_scores <- numeric(length(y_num))
    for (k in seq_len(n_folds)) {
      test_idx  <- which(folds == k)
      train_idx <- which(folds != k)
      if (method == "logistic") {
        df_tr <- as.data.frame(X_scaled[train_idx, , drop = FALSE])
        df_tr$y <- y_num[train_idx]
        fit_k <- stats::glm(y ~ ., data = df_tr, family = stats::binomial())
        coefs_k <- stats::coef(fit_k)[-1]
      } else {
        fit_k <- glmnet::glmnet(X_scaled[train_idx, , drop = FALSE],
          y_num[train_idx], family = "binomial", alpha = alpha,
          lambda = fit$lambda.min, standardize = FALSE)
        coefs_k <- as.vector(stats::coef(fit_k))[-1]
      }
      oof_scores[test_idx] <- as.vector(
        X_scaled[test_idx, , drop = FALSE] %*% coefs_k)
    }
    auc_cv <- .compute_single_auc(oof_scores, y_num == 1, y_num == 0)$auc_raw
  }

  out <- list(
    coefficients  = coefs,
    auc_train     = auc_train,
    auc_cv        = auc_cv,
    fitted_model  = fit,
    predictions   = data.frame(
      score_train = score_train,
      score_cv    = if (n_folds >= 2L) oof_scores else rep(NA_real_, length(y_num))
    ),
    settings = list(
      method      = method,
      alpha       = alpha,
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

#' @export
print.aucmat_panel <- function(x, ...) {
  cat("<aucmat_panel>  ", x$settings$method, "\n", sep = "")
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
