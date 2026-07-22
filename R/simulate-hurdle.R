# ==============================================================================
# simulate_hurdle_auc() — Zero-inflated biomarker generator
#
# Two-stage generation:
#   1. Bernoulli hurdle: P(X=0 | y) from logistic model
#   2. Magnitude: X | X>0, y ~ shifted normal (binormal AUC calibration)
#
# This preserves realistic zero-inflation AND class-discriminative signal.
# ==============================================================================

#' Simulate zero-inflated biomarker data with controlled Hurdle-AUC
#'
#' Generates biomarkers where a large fraction of values are exactly zero
#' (e.g., scRNA-seq, microbiome, detection-limited assays).  Each biomarker
#' is generated via a two-stage hurdle model:
#'
#' 1. **Zero hurdle**: P(X = 0 | y) = logistic(beta0 + beta1 * y)
#' 2. **Magnitude**: X | X > 0 follows shifted log-normal with
#'    binormal-calibrated mean separation.
#'
#' Values are clipped to [0, Inf).
#'
#' @param n Number of observations.
#' @param prevalence Proportion of positive class.
#' @param target_hurdle_aucs Target combined Hurdle-AUC values in (0, 1).
#' @param zero_rate_neg Vector of zero rates in the negative class.
#' @param zero_rate_pos Vector of zero rates in the positive class.
#' @param nonzero_target_aucs Optional vector of target AUCs on nonzero
#'   values only.  If NULL, derived from `target_hurdle_aucs`.
#' @param cv_between Coefficient of variation for nonzero values.
#'   Default 0.5.
#' @param seed Optional seed.
#'
#' @return An object of class `aucmat_simulation` with `data`, `target_*`,
#'   `achieved_*`, and hurdle-specific metadata.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_hurdle_auc(
#'   n = 500, prevalence = 0.3,
#'   target_hurdle_aucs = c(0.85, 0.72, 0.55),
#'   zero_rate_neg = c(0.55, 0.30, 0.80),
#'   zero_rate_pos = c(0.25, 0.10, 0.70)
#' )
#' # Compare standard vs hurdle AUC
#' fit_std <- aucmat(as.matrix(sim$data[, 1:3]), sim$data$truth, ci = "none")
#' fit_hur <- hurdle_auc(as.matrix(sim$data[, 1:3]), sim$data$truth)
#' data.frame(
#'   biomarker = paste0("X", 1:3),
#'   standard_AUC = round(fit_std$results$auc_strength, 3),
#'   hurdle_AUC   = round(fit_hur$results$hurdle_auc, 3)
#' )
#' }
simulate_hurdle_auc <- function(n, prevalence,
                                 target_hurdle_aucs,
                                 zero_rate_neg,
                                 zero_rate_pos,
                                 nonzero_target_aucs = NULL,
                                 cv_between = 0.5,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  p <- length(target_hurdle_aucs)
  if (length(zero_rate_neg) != p || length(zero_rate_pos) != p)
    stop("zero_rate_neg and zero_rate_pos must have same length as target_hurdle_aucs.")

  if (any(zero_rate_neg <= 0 | zero_rate_neg >= 1))
    stop("zero_rate_neg must be in (0, 1).")
  if (any(zero_rate_pos <= 0 | zero_rate_pos >= 1))
    stop("zero_rate_pos must be in (0, 1).")
  if (any(target_hurdle_aucs <= 0 | target_hurdle_aucs >= 1))
    stop("target_hurdle_aucs must be in (0, 1).")

  # Derive nonzero AUCs from hurdle AUCs if not supplied
  # Rule: Hurdle-AUC ≈ weighted combination of zero-AUC and nonzero-AUC
  # zero-AUC is driven by difference in zero rates between classes
  if (is.null(nonzero_target_aucs)) {
    nonzero_target_aucs <- numeric(p)
    for (j in seq_len(p)) {
      # Approximate zero-AUC from the difference in zero rates
      dz <- zero_rate_neg[j] - zero_rate_pos[j]
      zero_auc_approx <- 0.5 + 0.5 * abs(dz)  # rough approximation

      # Target nonzero AUC to achieve target hurdle AUC
      # Hurdle AUC ≈ w * zero_auc + (1-w) * nonzero_auc
      # with w proportional to zero rate
      w <- (zero_rate_neg[j] + zero_rate_pos[j]) / 2
      nonzero_target_aucs[j] <- (target_hurdle_aucs[j] - w * zero_auc_approx) / (1 - w)
      nonzero_target_aucs[j] <- pmax(0.51, pmin(0.99, nonzero_target_aucs[j]))
    }
  }

  # Generate outcome
  n_pos <- max(1L, round(n * prevalence))
  n_neg <- n - n_pos
  if (n_pos < 2L) n_pos <- 2L
  if (n_neg < 2L) n_neg <- 2L
  n <- n_pos + n_neg

  truth <- c(rep(1L, n_pos), rep(0L, n_neg))
  truth <- sample(truth)  # randomize order

  X_mat <- matrix(0, nrow = n, ncol = p)
  colnames(X_mat) <- paste0("X", seq_len(p))

  achieved_zero_pos <- numeric(p)
  achieved_zero_neg <- numeric(p)
  achieved_nonzero_aucs <- numeric(p)

  for (j in seq_len(p)) {
    # ---- Stage 1: Zero hurdle via logistic model ----
    # Solve for logistic parameters:
    # P(zero | y=0) = zero_rate_neg[j] => logit(p0) = beta0
    # P(zero | y=1) = zero_rate_pos[j] => logit(p1) = beta0 + beta1
    logit_p0 <- log(zero_rate_neg[j] / (1 - zero_rate_neg[j]))
    logit_p1 <- log(zero_rate_pos[j] / (1 - zero_rate_pos[j]))
    beta0 <- logit_p0
    beta1 <- logit_p1 - logit_p0

    # Generate zero indicators
    eta <- beta0 + beta1 * truth
    p_zero <- 1 / (1 + exp(-eta))
    is_zero <- rbinom(n, 1, p_zero) == 1

    # ---- Stage 2: Nonzero magnitude with binormal AUC ----
    n_nonzero <- sum(!is_zero)
    pos_nonzero <- which(!is_zero & truth == 1)
    neg_nonzero <- which(!is_zero & truth == 0)

    n_pos_nz <- length(pos_nonzero)
    n_neg_nz <- length(neg_nonzero)

    if (n_pos_nz >= 2L && n_neg_nz >= 2L) {
      delta <- sqrt(2) * stats::qnorm(nonzero_target_aucs[j])

      # Means: centered with n-weighted grand mean = 0
      mu0 <- delta * (-n_pos_nz / n_nonzero)
      mu1 <- delta * ( n_neg_nz / n_nonzero)

      # Log-normal to ensure positivity
      sigma <- cv_between * abs(delta)
      if (sigma < 0.1) sigma <- 0.1

      # Generate on log scale, then exponentiate
      vals <- numeric(n)
      vals[neg_nonzero] <- stats::rlnorm(n_neg_nz,
        meanlog = log(max(0.01, mu0 + 0.5)),
        sdlog = sigma)
      vals[pos_nonzero] <- stats::rlnorm(n_pos_nz,
        meanlog = log(max(0.01, mu1 + 0.5)),
        sdlog = sigma)

      X_mat[pos_nonzero, j] <- vals[pos_nonzero]
      X_mat[neg_nonzero, j] <- vals[neg_nonzero]

      # Compute achieved nonzero AUC
      if (n_pos_nz >= 2L && n_neg_nz >= 2L) {
        achieved_nonzero_aucs[j] <- .compute_single_auc(
          X_mat[!is_zero, j],
          truth[!is_zero] == 1,
          truth[!is_zero] == 0
        )$auc_raw
      }
    }

    # Record achieved zero rates
    achieved_zero_pos[j] <- mean(is_zero[truth == 1])
    achieved_zero_neg[j] <- mean(is_zero[truth == 0])
  }

  data_df <- as.data.frame(X_mat)
  data_df$truth <- truth

  out <- list(
    data                    = data_df,
    target_hurdle_aucs      = target_hurdle_aucs,
    target_zero_rate_neg    = zero_rate_neg,
    target_zero_rate_pos    = zero_rate_pos,
    target_nonzero_aucs     = nonzero_target_aucs,
    achieved_zero_rate_neg  = achieved_zero_neg,
    achieved_zero_rate_pos  = achieved_zero_pos,
    achieved_nonzero_aucs   = achieved_nonzero_aucs,
    n = n, prevalence = n_pos / n,
    seed = seed
  )
  class(out) <- c("aucmat_hurdle_simulation", "aucmat_simulation")
  out
}

#' Plot Hurdle-AUC diagnostics
#'
#' Visualizes the zero-inflation pattern and AUC decomposition for
#' hurdle-model results.
#'
#' @param fit An `aucmat_hurdle` object from [hurdle_auc()].
#' @param n_label Number of biomarkers to label.  Default 15.
#'
#' @return A `ggplot2` object.
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_hurdle_auc(n=100, prevalence=0.3,
#'   target_hurdle_aucs=c(0.8,0.7), zero_rate_neg=c(0.5,0.3),
#'   zero_rate_pos=c(0.2,0.1))
#' fit <- hurdle_auc(as.matrix(sim$data[,1:2]), sim$data$truth)
#' plot_hurdle_diagnostics(fit)
#' }
#' @export
plot_hurdle_diagnostics <- function(fit, n_label = 15L) {
  df <- fit$results
  df <- df[!is.na(df$hurdle_auc), , drop = FALSE]
  if (nrow(df) == 0L) stop("No valid results.")

  df$label <- ""
  top_idx <- head(order(df$hurdle_auc, decreasing = TRUE), n_label)
  df$label[top_idx] <- df$biomarker[top_idx]

  # Scatter: zero AUC vs nonzero AUC, coloured by hurdle AUC
  ggplot2::ggplot(df, ggplot2::aes(x = .data$zero_auc, y = .data$nonzero_auc,
    colour = .data$hurdle_auc, size = .data$zero_rate_total)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC", mid = "grey70", high = "#B2182B",
      midpoint = 0.7, name = "Hurdle AUC") +
    ggplot2::scale_size_continuous(range = c(1, 5), name = "Zero rate") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
      colour = "grey50", alpha = 0.4) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
      colour = "grey50", alpha = 0.4) +
    ggrepel::geom_text_repel(
      data = df[df$label != "", ],
      ggplot2::aes(label = .data$label), size = 2.8, max.overlaps = 30,
      show.legend = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Zero AUC (expression probability discrimination)",
      y = "Nonzero AUC (magnitude discrimination)",
      title = "Hurdle-AUC diagnostics: zero inflation vs magnitude discrimination")
}
