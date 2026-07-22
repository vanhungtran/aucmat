# ==============================================================================
# Partial AUC: estimation, inference, and comparison
# ==============================================================================

#' Compute partial AUC (pAUC) for a single biomarker
#'
#' Partial AUC over a restricted specificity or sensitivity range.
#' Uses the trapezoidal rule on the empirical ROC curve.
#'
#' @param x Numeric biomarker vector.
#' @param pos,neg Logical vectors.
#' @param focus `"specificity"` (default) or `"sensitivity"`.
#' @param interval Numeric length-2 vector: range of the focus measure in (0, 1).
#' @param standardize If `TRUE`, standardize to the interval
#'   \eqn{[0.5, 1]} via McClish correction.
#'
#' @return A list with `pauc_raw`, `pauc_standardized`, `focus`, `interval`.
#' @importFrom stats approx
#' @noRd
#' @keywords internal
.compute_partial_auc <- function(x, pos, neg,
                                  focus = c("specificity", "sensitivity"),
                                  interval = c(0.8, 1.0),
                                  standardize = TRUE) {
  focus <- match.arg(focus)
  if (length(interval) != 2L || any(interval < 0 | interval > 1) ||
      interval[1] >= interval[2])
    stop("interval must be a length-2 vector with 0 <= lower < upper <= 1.")

  # Build empirical ROC
  thresholds <- sort(unique(x), decreasing = TRUE)
  n_thresh <- length(thresholds)

  tpr <- fpr <- numeric(n_thresh + 2L)
  tpr[1] <- 0; fpr[1] <- 0
  tpr[n_thresh + 2L] <- 1; fpr[n_thresh + 2L] <- 1

  for (i in seq_len(n_thresh)) {
    t <- thresholds[i]
    tpr[i + 1L] <- mean(x[pos] >= t)
    fpr[i + 1L] <- mean(x[neg] >= t)
  }

  if (focus == "specificity") {
    # Restrict to fpr in [1 - interval[2], 1 - interval[1]]
    fpr_lo <- 1 - interval[2]
    fpr_hi <- 1 - interval[1]
    ord <- order(fpr)
    fpr_sorted <- fpr[ord]
    tpr_sorted <- tpr[ord]

    # Clip to interval
    keep <- fpr_sorted >= fpr_lo & fpr_sorted <= fpr_hi
    if (sum(keep) < 2L) {
      # Add boundary points
      fpr_use <- c(fpr_lo, fpr_sorted[keep], fpr_hi)
      tpr_use <- c(
        approx(fpr_sorted, tpr_sorted, fpr_lo, rule = 2)$y,
        tpr_sorted[keep],
        approx(fpr_sorted, tpr_sorted, fpr_hi, rule = 2)$y
      )
    } else {
      fpr_use <- c(fpr_lo, fpr_sorted[keep], fpr_hi)
      tpr_use <- c(
        approx(fpr_sorted, tpr_sorted, fpr_lo, rule = 2)$y,
        tpr_sorted[keep],
        approx(fpr_sorted, tpr_sorted, fpr_hi, rule = 2)$y
      )
    }
    pauc_raw <- .trapezoid(fpr_use, tpr_use)
    # McClish standardization to [0.5, 1]
    if (standardize) {
      fpr_width <- fpr_hi - fpr_lo
      pauc_min  <- 0.5 * fpr_width^2
      pauc_max  <- fpr_width
      pauc_std  <- 0.5 + 0.5 * (pauc_raw - pauc_min) / (pauc_max - pauc_min)
    } else {
      pauc_std <- NA_real_
    }
  } else {
    # focus == "sensitivity"
    tpr_lo <- interval[1]
    tpr_hi <- interval[2]
    ord <- order(tpr)
    fpr_sorted <- fpr[ord]
    tpr_sorted <- tpr[ord]

    keep <- tpr_sorted >= tpr_lo & tpr_sorted <= tpr_hi
    if (sum(keep) < 2L) {
      tpr_use <- c(tpr_lo, tpr_sorted[keep], tpr_hi)
      fpr_use <- c(
        approx(tpr_sorted, fpr_sorted, tpr_lo, rule = 2)$y,
        fpr_sorted[keep],
        approx(tpr_sorted, fpr_sorted, tpr_hi, rule = 2)$y
      )
    } else {
      tpr_use <- c(tpr_lo, tpr_sorted[keep], tpr_hi)
      fpr_use <- c(
        approx(tpr_sorted, fpr_sorted, tpr_lo, rule = 2)$y,
        fpr_sorted[keep],
        approx(tpr_sorted, fpr_sorted, tpr_hi, rule = 2)$y
      )
    }
    pauc_raw <- .trapezoid(tpr_use, fpr_use)
    if (standardize) {
      tpr_width <- tpr_hi - tpr_lo
      pauc_min  <- 0
      pauc_max  <- tpr_width
      pauc_std  <- 0.5 * pauc_raw / (0.5 * pauc_max)
    } else {
      pauc_std <- NA_real_
    }
  }

  list(
    pauc_raw          = pauc_raw,
    pauc_standardized  = pauc_std,
    focus             = focus,
    interval          = interval
  )
}

#' Trapezoidal integration
#' @noRd
.trapezoid <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(0)
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  sum(diff(x) * (y[-1] + y[-n])) / 2
}

#' Matrix-wide partial AUC computation
#'
#' @param X Numeric matrix.
#' @param pos,neg Logical vectors.
#' @param focus,interval,standardize Passed to `.compute_partial_auc()`.
#' @return Data.frame like `.compute_matrix_auc()` but with pAUC columns.
#' @noRd
.compute_matrix_partial_auc <- function(X, pos, neg,
                                         focus = "specificity",
                                         interval = c(0.8, 1.0),
                                         standardize = TRUE) {
  X <- as.matrix(X)
  cn <- colnames(X)
  if (is.null(cn)) cn <- paste0("V", seq_len(ncol(X)))

  out <- lapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    use <- !is.na(x)
    p_use <- pos & use; n_use <- neg & use
    np <- sum(p_use); nn <- sum(n_use)
    if (np < 2L || nn < 2L) {
      return(list(pauc_raw = NA_real_, pauc_std = NA_real_,
                  status = "insufficient_sample"))
    }
    pa <- .compute_partial_auc(x[use], p_use, n_use,
      focus = focus, interval = interval, standardize = standardize)
    pa$status <- "ok"
    pa
  })

  data.frame(
    biomarker         = cn,
    pauc_raw          = vapply(out, `[[`, numeric(1L), "pauc_raw"),
    pauc_standardized  = vapply(out, `[[`, numeric(1L), "pauc_standardized"),
    focus             = focus,
    interval_lower    = interval[1],
    interval_upper    = interval[2],
    status            = vapply(out, `[[`, character(1L), "status"),
    stringsAsFactors  = FALSE
  )
}
