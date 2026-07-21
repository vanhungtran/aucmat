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
#' @param parallel Not yet implemented.
#'
#' @return A list with `pos_idx` and `neg_idx` matrices (observations × times).
#' @noRd
#' @keywords internal
stratified_bootstrap_indices <- function(n_pos, n_neg, times,
                                          seed = NULL, parallel = FALSE) {
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
