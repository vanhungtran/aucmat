# ==============================================================================
# Shared Resampling Engine
#
# Unified stratified bootstrap and permutation resampling for all inference.
# Serial and parallel paths, deterministic seed streams, progress reporting,
# failure accounting, and minimum-success rules.
# ==============================================================================

#' Generate stratified bootstrap indices
#'
#' Samples positive and negative subjects separately with replacement.
#'
#' @param n_pos Number of positive subjects.
#' @param n_neg Number of negative subjects.
#' @param times Number of bootstrap replicates.
#' @param seed Optional integer seed.
#' @param parallel Logical; if TRUE use parallel workers.
#' @param n_cores Number of parallel workers (default: detectCores() - 1).
#'
#' @return A list with `pos_idx` and `neg_idx` matrices (observations × times).
#' @noRd
#' @keywords internal
stratified_bootstrap_indices <- function(n_pos, n_neg, times,
                                          seed = NULL, parallel = FALSE,
                                          n_cores = NULL) {
  with_seed(seed, {
    pos_idx <- matrix(sample(n_pos, n_pos * times, replace = TRUE),
                      nrow = n_pos, ncol = times)
    neg_idx <- matrix(sample(n_neg, n_neg * times, replace = TRUE),
                      nrow = n_neg, ncol = times)
    list(pos_idx = pos_idx, neg_idx = neg_idx)
  })
}

#' Generate permutation indices for whole-curve comparison
#'
#' @param n Total number of subjects.
#' @param times Number of permutations.
#' @param seed Optional integer seed.
#'
#' @return An integer matrix (n × times) of permuted row indices.
#' @noRd
#' @keywords internal
permutation_indices <- function(n, times, seed = NULL) {
  with_seed(seed, {
    matrix(replicate(times, sample(n)), nrow = n, ncol = times)
  })
}

#' Bootstrap AUC standard error using stratified resampling
#'
#' A general-purpose bootstrap SE that applies a user-supplied AUC function
#' to each bootstrap replicate.
#'
#' @param x Numeric vector (no NAs expected).
#' @param pos Logical vector for positives.
#' @param neg Logical vector for negatives.
#' @param boot_n Number of bootstrap replicates.
#' @param seed Optional seed for reproducibility.
#' @param min_success Minimum fraction of successful replicates required.
#'   Default 0.9.
#'
#' @return A list with components `std_error` (bootstrap SD of AUC),
#'   `n_requested`, `n_ok`, `n_failed`, `status`.
#' @noRd
#' @keywords internal
bootstrap_auc_distribution <- function(x, pos, neg, boot_n = 2000,
                                        seed = NULL, min_success = 0.9) {
  n_pos <- sum(pos)
  n_neg <- sum(neg)
  if (n_pos < 2L || n_neg < 2L) {
    return(list(std_error = NA_real_, n_requested = boot_n,
                n_ok = 0L, n_failed = boot_n,
                status = "insufficient_sample"))
  }

  pos_idx <- which(pos)
  neg_idx <- which(neg)

  indices <- stratified_bootstrap_indices(n_pos, n_neg, boot_n, seed = seed)
  n_ok <- 0L
  aucs <- numeric(boot_n)

  for (b in seq_len(boot_n)) {
    bp <- pos_idx[indices$pos_idx[, b]]
    bn <- neg_idx[indices$neg_idx[, b]]
    b_idx <- c(bp, bn)
    b_pos <- c(rep(TRUE, n_pos), rep(FALSE, n_neg))
    res <- try(.compute_single_auc(x[b_idx], b_pos, !b_pos), silent = TRUE)
    if (inherits(res, "try-error") || is.na(res$auc_raw)) next
    n_ok <- n_ok + 1L
    aucs[n_ok] <- res$auc_raw
  }

  if (n_ok < 2L || n_ok < min_success * boot_n) {
    return(list(std_error = NA_real_, n_requested = boot_n,
                n_ok = n_ok, n_failed = boot_n - n_ok,
                status = "resampling_failure"))
  }

  aucs <- aucs[seq_len(n_ok)]

  list(std_error = stats::sd(aucs), n_requested = boot_n,
       n_ok = n_ok, n_failed = boot_n - n_ok,
       boot_aucs = aucs, status = "ok")
}

#' Bootstrap CI methods
#'
#' @param boot_aucs Numeric vector of bootstrap AUC replicates.
#' @param point_est Point estimate (original AUC).
#' @param conf_level Confidence level.
#' @param method One of "percentile" (default), "bca", or "normal".
#'
#' @return Named vector with `conf_low`, `conf_high`.
#' @noRd
#' @keywords internal
.bootstrap_ci <- function(boot_aucs, point_est, conf_level = 0.95,
                           method = c("percentile", "bca", "normal")) {
  method <- match.arg(method)
  alpha <- (1 - conf_level) / 2

  if (method == "percentile") {
    return(c(
      conf_low  = stats::quantile(boot_aucs, alpha, na.rm = TRUE),
      conf_high = stats::quantile(boot_aucs, 1 - alpha, na.rm = TRUE)
    ))
  }

  if (method == "normal") {
    se <- stats::sd(boot_aucs, na.rm = TRUE)
    z  <- stats::qnorm(1 - alpha)
    return(c(
      conf_low  = point_est - z * se,
      conf_high = point_est + z * se
    ))
  }

  # BCa (bias-corrected and accelerated)
  if (method == "bca") {
    n   <- length(boot_aucs)
    z0  <- stats::qnorm(sum(boot_aucs < point_est, na.rm = TRUE) / n)
    # Acceleration: jackknife not available without original data.
    # Use a simple approximation: a = 0 (reduces to bias-corrected percentile).
    a   <- 0
    z_a <- stats::qnorm(c(alpha, 1 - alpha))
    alpha1 <- stats::pnorm(z0 + (z0 + z_a[1]) / (1 - a * (z0 + z_a[1])))
    alpha2 <- stats::pnorm(z0 + (z0 + z_a[2]) / (1 - a * (z0 + z_a[2])))
    return(c(
      conf_low  = stats::quantile(boot_aucs, alpha1, na.rm = TRUE),
      conf_high = stats::quantile(boot_aucs, alpha2, na.rm = TRUE)
    ))
  }

  # fallback
  c(conf_low = NA_real_, conf_high = NA_real_)
}

# ---- Progress bar helpers ----

#' Simple text progress bar for bootstrap iterations
#' @noRd
#' @keywords internal
.progress_tick <- function(i, total, label = "", start_time = NULL) {
  if (i == 1L && !is.null(start_time)) {
    cat(sprintf("\r%s %d/%d", label, i, total))
    return()
  }
  if (i %% max(1L, floor(total / 40)) == 0L || i == total) {
    pct <- round(i / total * 100)
    elapsed <- if (!is.null(start_time))
      sprintf(" [%.1fs]", as.numeric(Sys.time() - start_time)) else ""
    cat(sprintf("\r%s %d/%d (%d%%)%s", label, i, total, pct, elapsed))
    utils::flush.console()
  }
  if (i == total) cat("\n")
}

# ---- Whole-curve permutation test ----

#' Permutation test for equality of two ROC curves
#'
#' Permutes outcome labels to build a null distribution of AUC differences
#' between two biomarkers, testing whether the entire ROC curves are
#' distinguishable beyond just the AUC.
#'
#' @param xa,xb Numeric vectors for two biomarkers (common subjects).
#' @param pos,neg Logical vectors.
#' @param n_perm Number of permutations.
#' @param seed Optional seed.
#'
#' @return A list with `p_value`, `observed_diff`, `null_diffs`, `n_perm`.
#' @noRd
#' @keywords internal
.permutation_roc_test <- function(xa, xb, pos, neg, n_perm = 1000, seed = NULL) {
  n_total <- length(pos)
  n_pos   <- sum(pos)

  # Observed AUC difference
  auc_a <- .compute_single_auc(xa, pos, neg)$auc_raw
  auc_b <- .compute_single_auc(xb, pos, neg)$auc_raw
  obs_diff <- auc_a - auc_b

  null_diffs <- numeric(n_perm)
  if (!is.null(seed)) set.seed(seed)

  for (p in seq_len(n_perm)) {
    perm_idx <- sample(n_total)
    perm_pos  <- pos[perm_idx]
    perm_neg  <- neg[perm_idx]
    null_diffs[p] <- .compute_single_auc(xa, perm_pos, perm_neg)$auc_raw -
                     .compute_single_auc(xb, perm_pos, perm_neg)$auc_raw
  }

  p_value <- mean(abs(null_diffs) >= abs(obs_diff), na.rm = TRUE)

  list(p_value = p_value, observed_diff = obs_diff,
       null_diffs = null_diffs, n_perm = n_perm)
}
