# ==============================================================================
# S3 Methods for aucmat_screen and related classes
# ==============================================================================

# ---- aucmat_screen ----------------------------------------------------------

#' @export
print.aucmat_screen <- function(x, n = 10L, ...) {
  cat("<aucmat_screen>  ", nrow(x$results), " biomarkers\n", sep = "")
  cat("  Outcome: ", x$sample_summary$n_positive, " positive / ",
      x$sample_summary$n_negative, " negative",
      "  (positive = ", x$sample_summary$positive_class, ")\n", sep = "")
  cat("  CI: ", x$settings$ci, "  |  adjust: ", x$settings$adjust,
      "  |  na_action: ", x$settings$na_action, "\n", sep = "")

  n_failed <- sum(x$results$status != "ok")
  if (n_failed > 0L) {
    cat("  ", n_failed, " biomarker(s) with non-ok status\n", sep = "")
  }

  cat("\nTop biomarkers by discrimination strength:\n")
  n_show <- min(n, nrow(x$results))
  top <- x$results[seq_len(n_show), ]
  cols <- c("rank", "biomarker", "auc_raw", "auc_strength",
            "effect_direction", "p_value", "q_value", "status")
  show_cols <- intersect(cols, names(top))
  print(top[, show_cols], row.names = FALSE)

  invisible(x)
}

#' @export
summary.aucmat_screen <- function(object, ...) {
  res <- object$results
  ss  <- object$sample_summary

  cat("aucmat screening summary\n")
  cat("========================\n")
  cat("Total samples:      ", ss$n_total, "\n")
  cat("Usable samples:     ", ss$n_used, "\n")
  cat("Excluded (NA y):    ", ss$n_excluded, "\n")
  cat("Positive class:     ", ss$positive_class, " (n = ", ss$n_positive, ")\n", sep = "")
  cat("Negative class:     ", ss$negative_class, " (n = ", ss$n_negative, ")\n", sep = "")
  cat("\n")

  cat("Biomarkers screened:", nrow(res), "\n")

  status_counts <- table(res$status)
  cat("Status counts:\n")
  print(status_counts)

  warn_counts <- table(res$warning, useNA = "no")
  if (length(warn_counts) > 0L) {
    cat("\nWarning counts:\n")
    print(warn_counts)
  }

  cat("\nMissingness:\n")
  cat("  Mean missing fraction: ", round(mean(res$missing_fraction), 4), "\n")
  cat("  Max  missing fraction: ", round(max(res$missing_fraction), 4), "\n")
  cat("  Fully observed:        ", sum(res$missing_fraction == 0), " biomarkers\n")

  cat("\nMultiplicity (", object$settings$adjust, "):\n", sep = "")
  n_sig <- sum(res$q_value < 0.05, na.rm = TRUE)
  cat("  q < 0.05: ", n_sig, " biomarkers\n")
  cat("  q < 0.01: ", sum(res$q_value < 0.01, na.rm = TRUE), " biomarkers\n")

  invisible(object)
}

#' @export
as.data.frame.aucmat_screen <- function(x, ...) {
  x$results
}

#' @export
plot.aucmat_screen <- function(x, ...) {
  plot_auc_rank(x, ...)
}

#' @export
subset.aucmat_screen <- function(x, subset, ...) {
  r <- x$results
  e <- substitute(subset)
  idx <- eval(e, r, parent.frame())
  if (!is.logical(idx)) idx <- seq_len(nrow(r)) %in% idx
  r[idx, , drop = FALSE]
}

# ---- aucmat_compare ---------------------------------------------------------

#' @export
print.aucmat_compare <- function(x, n = 10L, ...) {
  cat("<aucmat_compare>  ", nrow(x), " pairwise comparisons\n", sep = "")
  n_show <- min(n, nrow(x))
  # strip class to avoid recursive dispatch
  print.data.frame(x[seq_len(n_show), , drop = FALSE], row.names = FALSE)
  if (nrow(x) > n_show) cat("... ", nrow(x) - n_show, " more rows\n", sep = "")
  invisible(x)
}

# ---- aucmat_stability -------------------------------------------------------

#' @export
print.aucmat_stability <- function(x, n = 10L, ...) {
  cat("<aucmat_stability>  ", x$settings$n_ok, "/", x$settings$times,
      " successful replicates\n", sep = "")
  cat("Top biomarkers by median rank:\n")
  n_show <- min(n, nrow(x$rank_summary))
  top <- x$rank_summary[seq_len(n_show), ]
  print(top, row.names = FALSE)
  invisible(x)
}
