# ==============================================================================
# Hurdle-AUC: Two-stage model for zero-inflated biomarker data
#
# Stage 1: P(X > 0 | y) = logistic(beta0 + beta1 * y)  — zero hurdle
# Stage 2: AUC on X | X > 0                                  — magnitude discrimination
#
# The hurdle model separates the binary decision (expressed vs not)
# from the continuous decision (expression level when expressed).
# This preserves zero-inflation AND class-discriminative signal.
# ==============================================================================

#' Hurdle-AUC for zero-inflated biomarkers
#'
#' Computes a two-stage Hurdle-AUC for data where a large fraction of
#' values are exactly zero (e.g., scRNA-seq, microbiome, mass-spec with
#' detection limits).  Stage 1 models zero vs non-zero via logistic
#' regression; Stage 2 computes standard AUC on the non-zero values only.
#'
#' @param X Numeric matrix (n x p).  Zeros are treated as the hurdle.
#' @param y Binary outcome.
#' @param positive Optional positive class label.
#' @param zero_threshold Values at or below this are treated as zero.
#'   Default 0.
#'
#' @return A list of class `aucmat_hurdle` with components:
#'   `results` (data.frame with per-biomarker hurdle diagnostics),
#'   `stage1_models` (list of logistic regression fits),
#'   `settings`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_hurdle_auc(n = 300, prevalence = 0.3,
#'   target_hurdle_aucs = c(0.85, 0.72, 0.55),
#'   zero_rate_neg = c(0.55, 0.30, 0.80),
#'   zero_rate_pos = c(0.25, 0.10, 0.70))
#' X <- as.matrix(sim$data[, 1:3])
#' y <- sim$data$truth
#' res <- hurdle_auc(X, y)
#' print(res)
#' }
hurdle_auc <- function(X, y, positive = NULL, zero_threshold = 0) {
  X <- as.matrix(X)
  p <- ncol(X)
  cn <- colnames(X)
  if (is.null(cn)) cn <- paste0("X", seq_len(p))

  y_norm <- .normalize_binary_y(y, positive)
  pos <- y_norm == levels(y_norm)[2L]
  neg <- !pos
  y_num <- .binary_factor_to_numeric(y_norm)
  n_total <- length(y_num)

  results <- data.frame(
    biomarker       = cn,
    zero_rate_total = NA_real_,
    zero_rate_pos   = NA_real_,
    zero_rate_neg   = NA_real_,
    zero_odds_ratio = NA_real_,
    zero_auc        = NA_real_,    # AUC of zero/nonzero indicator
    nonzero_auc     = NA_real_,    # AUC on X > 0 only
    hurdle_auc      = NA_real_,    # combined two-stage AUC
    n_nonzero       = NA_integer_,
    n_pos           = NA_integer_,
    n_neg           = NA_integer_,
    stage1_p_value  = NA_real_,    # p-value for beta1 in logistic model
    stringsAsFactors = FALSE
  )

  stage1_models <- vector("list", p)
  names(stage1_models) <- cn

  for (j in seq_len(p)) {
    x <- X[, j]
    use <- !is.na(x)
    xu <- x[use]; yu <- y_num[use]; pu <- pos[use]; nu <- neg[use]
    np <- sum(pu); nn <- sum(nu)
    if (np < 5L || nn < 5L) next

    is_zero <- xu <= zero_threshold
    nz_pos  <- sum(is_zero[pu])
    nz_neg  <- sum(is_zero[nu])

    results$zero_rate_pos[j]   <- nz_pos / np
    results$zero_rate_neg[j]   <- nz_neg / nn
    results$zero_rate_total[j] <- mean(is_zero)
    results$n_pos[j] <- np
    results$n_neg[j] <- nn

    # Odds ratio for zero in positive vs negative
    if (nz_pos > 0 && nz_neg > 0 && nz_pos < np && nz_neg < nn) {
      odds_pos <- nz_pos / (np - nz_pos)
      odds_neg <- nz_neg / (nn - nz_neg)
      results$zero_odds_ratio[j] <- odds_pos / odds_neg
    }

    # ---- Stage 1: Logistic regression zero ~ class ----
    stage1_fit <- try(
      stats::glm(is_zero ~ yu, family = stats::binomial()), silent = TRUE)
    if (!inherits(stage1_fit, "try-error")) {
      stage1_models[[j]] <- stage1_fit
      results$stage1_p_value[j] <- summary(stage1_fit)$coefficients[2, 4]
    }

    # ---- Stage 1 AUC: zero/nonzero indicator as classifier ----
    if (np >= 2L && nn >= 2L) {
      zero_score <- as.numeric(!is_zero)
      z_auc <- .compute_single_auc(zero_score, pu, nu)
      results$zero_auc[j] <- z_auc$auc_raw
    }

    # ---- Stage 2: AUC on non-zero values only ----
    nonzero_idx <- !is_zero
    if (sum(nonzero_idx & pu) >= 2L && sum(nonzero_idx & nu) >= 2L) {
      results$nonzero_auc[j] <- .compute_single_auc(
        xu[nonzero_idx], pu[nonzero_idx], nu[nonzero_idx])$auc_raw
      results$n_nonzero[j] <- sum(nonzero_idx)
    }

    # ---- Combined Hurdle-AUC ----
    # Use predicted probability from Stage 1 for zeros/nonzeros,
    # combine with Stage 2 score for nonzeros to produce one score.
    # For zero observations: score = predicted P(zero | positive class)
    #   → flip sign so higher = more likely positive
    # For nonzero observations: use the observed value as score
    # Then compute AUC on this combined score.

    if (!inherits(stage1_fit, "try-error") &&
        sum(nonzero_idx & pu) >= 2L && sum(nonzero_idx & nu) >= 2L) {

      # Build composite score
      score <- numeric(length(xu))

      # For zeros: predicted probability of being nonzero in positive class
      # Use logistic model: P(nonzero | y=1) vs P(nonzero | y=0)
      # score = predicted probability of being positive class
      p_pred <- stats::predict(stage1_fit, type = "response")

      # For nonzero obs: use the actual value (normalized to [0,1] within expressed)
      x_nonzero <- xu
      x_nonzero[is_zero] <- NA
      x_range <- range(x_nonzero, na.rm = TRUE)
      if (x_range[2] > x_range[1]) {
        x_norm <- (x_nonzero - x_range[1]) / (x_range[2] - x_range[1])
      } else {
        x_norm <- rep(0.5, length(xu))
      }

      # Composite: for zeros use p_pred (probability of being nonzero, which
      # correlates with class), for nonzeros blend p_pred + normalized value
      score[is_zero]  <- p_pred[is_zero]
      score[!is_zero] <- 0.5 * p_pred[!is_zero] + 0.5 * x_norm[!is_zero]

      if (sum(pu) >= 2L && sum(nu) >= 2L) {
        results$hurdle_auc[j] <- .compute_single_auc(score, pu, nu)$auc_raw
      }
    }
  }

  out <- list(
    results        = results,
    stage1_models  = stage1_models,
    settings       = list(
      zero_threshold = zero_threshold,
      n_total        = n_total
    )
  )
  class(out) <- "aucmat_hurdle"
  out
}

#' @export
print.aucmat_hurdle <- function(x, n = 10L, ...) {
  cat("<aucmat_hurdle>  ", nrow(x$results), " biomarkers\n", sep = "")
  cat("  Zero threshold:", x$settings$zero_threshold, "\n\n")

  # Sort by hurdle_auc descending, fall back to nonzero_auc
  df <- x$results
  df <- df[order(df$hurdle_auc, df$nonzero_auc,
    decreasing = TRUE, na.last = TRUE), ]
  n_show <- min(n, nrow(df))

  cols <- c("biomarker", "zero_rate_total", "zero_auc", "nonzero_auc",
    "hurdle_auc", "n_nonzero")
  show_cols <- intersect(cols, names(df))
  print(df[seq_len(n_show), show_cols], row.names = FALSE)

  if (nrow(df) > n_show)
    cat("... ", nrow(df) - n_show, " more\n", sep = "")
  invisible(x)
}

#' @export
summary.aucmat_hurdle <- function(object, ...) {
  df <- object$results
  cat("Hurdle-AUC Summary\n")
  cat("==================\n")
  cat("Biomarkers:         ", nrow(df), "\n")
  cat("Zero threshold:     ", object$settings$zero_threshold, "\n\n")

  cat("Zero inflation:\n")
  cat("  Mean zero rate:   ", round(mean(df$zero_rate_total, na.rm = TRUE), 4), "\n")
  cat("  Mean zero pos:    ", round(mean(df$zero_rate_pos, na.rm = TRUE), 4), "\n")
  cat("  Mean zero neg:    ", round(mean(df$zero_rate_neg, na.rm = TRUE), 4), "\n\n")

  cat("AUC components:\n")
  cat("  Mean zero AUC:    ", round(mean(df$zero_auc, na.rm = TRUE), 4),
    " (discrimination via expression probability)\n")
  cat("  Mean nonzero AUC: ", round(mean(df$nonzero_auc, na.rm = TRUE), 4),
    " (discrimination in expression level)\n")
  cat("  Mean hurdle AUC:  ", round(mean(df$hurdle_auc, na.rm = TRUE), 4),
    " (combined two-stage)\n")
  invisible(object)
}

#' @export
plot.aucmat_hurdle <- function(x, ...) {
  df <- x$results
  df <- df[!is.na(df$hurdle_auc), , drop = FALSE]
  if (nrow(df) == 0L) stop("No valid results to plot.")

  df$biomarker <- factor(df$biomarker,
    levels = df$biomarker[order(df$hurdle_auc)])

  # Reshape for grouped bar/point plot
  long <- data.frame(
    biomarker = rep(df$biomarker, 3),
    auc = c(df$zero_auc, df$nonzero_auc, df$hurdle_auc),
    stage = rep(c("Zero AUC", "Nonzero AUC", "Hurdle AUC"), each = nrow(df))
  )

  ggplot2::ggplot(long, ggplot2::aes(x = .data$auc, y = .data$biomarker,
    colour = .data$stage, shape = .data$stage)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(
      values = c("Zero AUC" = "#FDB863", "Nonzero AUC" = "#B2ABD2",
        "Hurdle AUC" = "#2c7da0"),
      name = NULL) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
      colour = "grey50", alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "AUC", y = NULL,
      title = "Hurdle-AUC decomposition: zero vs nonzero discrimination",
      subtitle = "Hurdle AUC = combined two-stage score")
}
