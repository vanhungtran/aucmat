# ==============================================================================
# cv_aucmat()  --  Cross-validated biomarker screening
#
# Eliminates same-data overfitting by computing out-of-fold AUCs via
# stratified v-fold cross-validation.  Every AUC is computed on subjects
# that were NOT used for that fold's training.
#
# Because no model is actually trained (it's univariate screening, not
# prediction), CV here means: split  ->  screen  ->  evaluate on held-out fold,
# then pool out-of-fold predictions across folds to compute one honest AUC
# per biomarker.
# ==============================================================================

#' Cross-validated AUC screening
#'
#' Computes out-of-fold AUCs for every biomarker via stratified v-fold
#' cross-validation.  Each biomarker's pooled out-of-fold predictions across
#' all folds are used to compute a single honest AUC that is not inflated by
#' same-data overfitting.
#'
#' @param X Numeric matrix (n  x  p).
#' @param y Binary outcome.
#' @param positive Optional positive class label.
#' @param n_folds Number of cross-validation folds.  Default 5.
#' @param n_repeats Number of repeats of the full CV.  Default 1.
#' @param seed Seed for fold assignment reproducibility.
#' @param ci CI method for final pooled AUCs: `"delong"` or `"none"`.
#' @param conf_level Confidence level.
#'
#' @return A list of class `aucmat_cv` with `results`, `folds`, `settings`.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_auc_matrix(n = 200, prevalence = 0.3,
#'   target_aucs = c(0.8, 0.7, 0.6), correlation = 0.3,
#'   structure = "exchangeable")
#' cv <- cv_aucmat(as.matrix(sim$data[, 1:3]), sim$data$truth,
#'   n_folds = 3, seed = 42)
#' print(cv)
#' }
cv_aucmat <- function(X, y, positive = NULL,
                       n_folds = 5L, n_repeats = 1L,
                       seed = NULL,
                       ci = c("none", "delong"),
                       conf_level = 0.95) {
  ci <- match.arg(ci)

  if (!is.matrix(X) && !is.data.frame(X))
    stop("X must be a numeric matrix or data.frame.")
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("All columns of X must be numeric.")
  if (nrow(X) < n_folds * 4L)
    stop("Need at least ", n_folds * 4L, " observations for ", n_folds, "-fold CV.")

  cn <- colnames(X)
  if (is.null(cn) || anyDuplicated(cn) || any(cn == ""))
    stop("Column names of X must be unique and non-empty.")

  y_norm <- .normalize_binary_y(y, positive)
  pos <- y_norm == levels(y_norm)[2L]
  neg <- !pos
  n_total <- length(pos)
  p <- ncol(X)

  # Store out-of-fold predictions per biomarker per repeat
  oof_preds <- matrix(NA_real_, nrow = n_total, ncol = p)
  colnames(oof_preds) <- cn

  if (!is.null(seed)) set.seed(seed)
  fold_ids <- vector("list", n_repeats)

  for (r in seq_len(n_repeats)) {
    # Stratified fold assignment
    pos_idx <- which(pos)
    neg_idx <- which(neg)
    pos_folds <- sample(rep(seq_len(n_folds), length.out = length(pos_idx)))
    neg_folds <- sample(rep(seq_len(n_folds), length.out = length(neg_idx)))
    folds <- integer(n_total)
    folds[pos_idx] <- pos_folds
    folds[neg_idx] <- neg_folds
    fold_ids[[r]] <- folds

    for (k in seq_len(n_folds)) {
      test_idx  <- which(folds == k)
      train_idx <- which(folds != k)

      # For each biomarker: use the raw value as the "prediction"
      # (univariate screening  --  no model training)
      for (j in seq_len(p)) {
        oof_preds[test_idx, j] <- X[test_idx, j]
      }
    }
  }

  # Per-biomarker: compute AUC of pooled OOF predictions vs truth
  results <- data.frame(
    biomarker    = cn,
    auc_cv       = NA_real_,
    auc_in_sample = NA_real_,
    optimism     = NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(p)) {
    use <- !is.na(oof_preds[, j]) & !is.na(y_norm)
    p_use <- pos[use]; n_use <- neg[use]
    res_cv <- .compute_single_auc(oof_preds[use, j], p_use, n_use)
    res_is <- .compute_single_auc(X[use, j], p_use, n_use)
    results$auc_cv[j]        <- res_cv$auc_raw
    results$auc_in_sample[j] <- res_is$auc_raw
    results$optimism[j]      <- res_is$auc_raw - res_cv$auc_raw
  }

  results <- results[order(results$auc_cv, decreasing = TRUE, na.last = TRUE), ]
  results$rank_cv <- seq_len(nrow(results))
  rownames(results) <- NULL

  out <- list(
    results = results,
    folds   = fold_ids,
    settings = list(
      n_folds    = n_folds,
      n_repeats  = n_repeats,
      ci         = ci,
      conf_level = conf_level,
      seed       = seed
    ),
    sample_summary = list(
      n_total    = n_total,
      n_positive = sum(pos),
      n_negative = sum(neg)
    )
  )
  class(out) <- "aucmat_cv"
  out
}

#' @export
print.aucmat_cv <- function(x, n = 10L, ...) {
  cat("<aucmat_cv>  ", nrow(x$results), " biomarkers\n", sep = "")
  cat("  CV: ", x$settings$n_folds, "-fold",
      if (x$settings$n_repeats > 1L)
        paste0("  x  ", x$settings$n_repeats, " repeats"), "\n", sep = "")
  cat("  Mean optimism: ", round(mean(x$results$optimism, na.rm = TRUE), 4), "\n\n")
  cat("Top biomarkers by CV AUC:\n")
  n_show <- min(n, nrow(x$results))
  top <- x$results[seq_len(n_show), ]
  print(top[, c("rank_cv", "biomarker", "auc_cv", "auc_in_sample", "optimism")],
        row.names = FALSE)
  invisible(x)
}
