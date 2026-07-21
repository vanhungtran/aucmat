# ==============================================================================
# Inference Engine
#
# DeLong and stratified-bootstrap variance estimation, confidence intervals,
# and two-sided p-values for departure from AUC = 0.5.  Not exported.
# ==============================================================================

#' DeLong variance for a single biomarker
#'
#' Implements the DeLong et al. (1988) variance estimator for the area under
#' the empirical ROC curve.
#'
#' @param x Numeric vector (no missing values expected).
#' @param pos Logical; TRUE for positive-class observations.
#' @param neg Logical; TRUE for negative-class observations.
#'
#' @return Named vector with `std_error`.  Returns `NA` when variance is
#'   not estimable.
#' @noRd
#' @keywords internal
.delong_variance <- function(x, pos, neg) {
  x_pos <- x[pos]
  x_neg <- x[neg]
  np <- length(x_pos)
  nn <- length(x_neg)

  if (np < 2L || nn < 2L) return(c(std_error = NA_real_))

  # placement values for positives
  V10 <- vapply(x_pos, function(xi) {
    mean((x_neg < xi) + 0.5 * (x_neg == xi))
  }, numeric(1L))

  # placement values for negatives (fraction of positives *above* them)
  V01 <- vapply(x_neg, function(xj) {
    mean((x_pos > xj) + 0.5 * (x_pos == xj))
  }, numeric(1L))

  S10 <- stats::var(V10)  # sample variance (divides by np-1)
  S01 <- stats::var(V01)

  se <- sqrt(S10 / np + S01 / nn)
  c(std_error = se)
}

#' Stratified bootstrap AUC standard error
#'
#' Resamples positive and negative subjects independently with replacement,
#' recomputes AUC for each replicate, and returns the standard deviation of
#' the bootstrap AUC distribution.
#'
#' @param x Numeric vector (no missing values).
#' @param pos Logical; TRUE for positives.
#' @param neg Logical; TRUE for negatives.
#' @param boot_n Number of bootstrap replicates.
#' @param seed Optional seed for this biomarker (applied per call).
#'
#' @return Named vector with `std_error` (bootstrap SD of AUC).
#' @noRd
#' @keywords internal
.bootstrap_auc_se <- function(x, pos, neg, boot_n = 2000, seed = NULL) {
  np <- sum(pos)
  nn <- sum(neg)
  if (np < 2L || nn < 2L) return(c(std_error = NA_real_))

  if (!is.null(seed)) {
    old <- globalenv()$.Random.seed
    on.exit({
      if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
      else rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
  }

  pos_idx <- which(pos)
  neg_idx <- which(neg)

  aucs <- vapply(seq_len(boot_n), function(b) {
    bp <- sample(pos_idx, np, replace = TRUE)
    bn <- sample(neg_idx, nn, replace = TRUE)
    b_use <- c(bp, bn)
    b_pos <- c(rep(TRUE, np), rep(FALSE, nn))
    res <- .compute_single_auc(x[b_use], b_pos, !b_pos)
    res$auc_raw
  }, numeric(1L))

  c(std_error = stats::sd(aucs, na.rm = TRUE))
}

#' Apply inference to a matrix-wide AUC result
#'
#' Adds `std_error`, `conf_low`, `conf_high`, `p_value` columns to the
#' results data.frame.
#'
#' @param results Data.frame from `.compute_matrix_auc()`.
#' @param X Numeric matrix.
#' @param pos,neg Logical vectors.
#' @param ci Method: `"delong"`, `"bootstrap"`, or `"none"`.
#' @param conf_level Confidence level in (0, 1).
#' @param boot_n Bootstrap replicates.
#'
#' @return The results data.frame with additional columns.
#' @noRd
#' @keywords internal
.apply_inference <- function(results, X, pos, neg,
                              ci = c("delong", "bootstrap", "none"),
                              conf_level = 0.95, boot_n = 2000) {
  ci <- match.arg(ci)
  X <- as.matrix(X)

  z <- stats::qnorm((1 + conf_level) / 2)

  results$std_error <- NA_real_
  results$conf_low  <- NA_real_
  results$conf_high <- NA_real_
  results$p_value   <- NA_real_

  ok <- results$status == "ok"
  if (!any(ok)) return(results)

  for (j in which(ok)) {
    x <- X[, results$biomarker[j]]

    if (ci == "delong") {
      se_val <- unname(.delong_variance(x, pos, neg)["std_error"])
    } else if (ci == "bootstrap") {
      se_val <- unname(.bootstrap_auc_se(x, pos, neg, boot_n = boot_n)["std_error"])
    } else {
      next  # "none"
    }

    if (is.na(se_val) || se_val <= 0) {
      results$status[j] <- "inference_failed"
      next
    }

    auc_raw <- results$auc_raw[j]
    results$std_error[j] <- se_val
    results$conf_low[j]  <- auc_raw - z * se_val
    results$conf_high[j] <- auc_raw + z * se_val
    results$p_value[j]   <- 2 * stats::pnorm(abs(auc_raw - 0.5) / se_val,
                                             lower.tail = FALSE)
  }

  results
}
