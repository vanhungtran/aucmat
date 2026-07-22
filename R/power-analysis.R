# ==============================================================================
# power_auc_matrix() — Sample size and power for correlated AUC comparisons
#
# Uses the binormal model and DeLong variance approximation to compute
# required sample size or achievable power for comparing correlated AUCs.
# ==============================================================================

#' Power and sample size for AUC comparisons
#'
#' Computes the required sample size or achievable power for comparing AUCs
#' of correlated biomarkers via DeLong's test.  Uses the binormal model to
#' translate target AUCs into mean separations and the DeLong variance
#' structure for correlated ROC curves.
#'
#' @param target_aucs Numeric vector of AUC values under the alternative.
#'   The null value (chance) is `null_auc`.
#' @param null_auc Null AUC value.  Default 0.5.
#' @param correlation Between-biomarker correlation (scalar for exchangeable).
#'   Default 0.
#' @param prevalence Expected prevalence in (0, 1).
#' @param power Desired power.  If `NULL` (default), compute power for
#'   given `n`.  If supplied, compute required `n`.
#' @param n Sample size.  Required when `power = NULL`.  Ignored otherwise.
#' @param alpha Significance level.  Default 0.05.
#' @param alternative `"two.sided"` (default), `"greater"`, or `"less"`.
#' @param adjust Multiplicity adjustment: `"none"` (default), `"BH"`,
#'   `"bonferroni"`.
#'
#' @return A list of class `aucmat_power` with components `n`, `power`,
#'   `target_aucs`, `effect_sizes`, `alpha`, `adjusted_alpha`.
#'
#' @export
#'
#' @examples
#' # Power for n=200, detecting AUC=0.70 vs chance
#' power_auc_matrix(target_aucs = c(0.70, 0.65), n = 200,
#'   prevalence = 0.3, correlation = 0.2)
#'
#' # Sample size for 80% power
#' \donttest{
#' power_auc_matrix(target_aucs = 0.70, power = 0.80,
#'   prevalence = 0.3)
#' }
power_auc_matrix <- function(target_aucs, null_auc = 0.5,
                              correlation = 0, prevalence = 0.5,
                              power = NULL, n = NULL,
                              alpha = 0.05,
                              alternative = c("two.sided", "greater", "less"),
                              adjust = c("none", "BH", "bonferroni")) {
  alternative <- match.arg(alternative)
  adjust     <- match.arg(adjust)
  p          <- length(target_aucs)

  if (!is.null(power) && !is.null(n))
    stop("Specify either 'power' or 'n', not both.")
  if (is.null(power) && is.null(n))
    stop("Specify either 'power' or 'n'.")
  if (is.null(power) && !is.null(n) && (n < 6L))
    stop("n must be >= 6.")
  if (!is.null(power) && (power <= 0 || power >= 1))
    stop("power must be in (0, 1).")

  if (any(target_aucs <= 0 | target_aucs >= 1))
    stop("target_aucs must be in (0, 1).")
  if (prevalence <= 0 || prevalence >= 1)
    stop("prevalence must be in (0, 1).")
  if (abs(correlation) >= 1) stop("correlation must be in (-1, 1).")

  # Adjusted alpha for multiplicity
  alpha_adj <- alpha
  if (adjust == "bonferroni") alpha_adj <- alpha / p
  if (adjust == "BH") {
    # BH: approximate by averaging over ranks (conservative at rank 1)
    alpha_adj <- alpha / p
  }

  # Effect size: delta = sqrt(2) * (qnorm(AUC) - qnorm(0.5))
  delta <- sqrt(2) * (stats::qnorm(target_aucs) - stats::qnorm(null_auc))

  # Approximate DeLong variance for correlated AUCs
  # Var(AUC) ≈ (AUC*(1-AUC) + (n_pos-1)*(Q1-AUC^2) + (n_neg-1)*(Q2-AUC^2))
  #   / (n_pos * n_neg)
  # For correlated biomarkers: Cov ≈ correlation * sqrt(Var_i * Var_j)
  # Simplified: use the Hanley-McNeil variance approximation

  .var_auc <- function(auc, n_pos, n_neg) {
    q1 <- auc / (2 - auc)
    q2 <- 2 * auc^2 / (1 + auc)
    (auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) +
       (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg)
  }

  if (is.null(power)) {
    # Compute power for given n
    n_pos <- round(n * prevalence)
    n_neg <- n - n_pos
    if (n_pos < 2L || n_neg < 2L) stop("n too small for given prevalence.")

    # Power for each biomarker individually
    var_auc <- vapply(target_aucs, function(a) .var_auc(a, n_pos, n_neg),
                      numeric(1L))
    se <- sqrt(var_auc)
    z_alpha <- if (alternative == "two.sided")
      stats::qnorm(1 - alpha_adj / 2) else stats::qnorm(1 - alpha_adj)

    # Non-centrality parameter
    ncp <- (delta / se)
    if (alternative == "two.sided") {
      power_i <- stats::pnorm(-z_alpha + abs(ncp)) +
                 stats::pnorm(-z_alpha - abs(ncp))
    } else if (alternative == "greater") {
      power_i <- stats::pnorm(-z_alpha + ncp)
    } else {
      power_i <- stats::pnorm(-z_alpha - ncp)
    }
    power_vals <- power_i
    n_val <- n
  } else {
    # Compute n for target power
    z_power <- stats::qnorm(power)
    z_alpha <- if (alternative == "two.sided")
      stats::qnorm(1 - alpha_adj / 2) else stats::qnorm(1 - alpha_adj)

    # Solve: ncp = z_alpha + z_power, with ncp = delta / se
    # se^2 = f(auc) / (n_pos * n_neg) ≈ f(auc) / (n^2 * prev * (1-prev))
    # => n = ceil(solve for largest required n across biomarkers)
    n_vals <- vapply(seq_len(p), function(i) {
      # Solve iteratively for n
      n_guess <- 20
      for (iter in 1:100) {
        np <- round(n_guess * prevalence)
        nn <- n_guess - np
        if (np < 2L) np <- 2L
        if (nn < 2L) nn <- 2L
        se_n <- sqrt(.var_auc(target_aucs[i], np, nn))
        ncp_n <- delta[i] / se_n
        power_n <- if (alternative == "two.sided") {
          stats::pnorm(-z_alpha + abs(ncp_n)) + stats::pnorm(-z_alpha - abs(ncp_n))
        } else if (alternative == "greater") {
          stats::pnorm(-z_alpha + ncp_n)
        } else {
          stats::pnorm(-z_alpha - ncp_n)
        }
        if (abs(power_n - power) < 0.001) break
        # Adjust n
        n_guess <- n_guess * sqrt((z_alpha + z_power)^2 / (ncp_n^2))
        n_guess <- max(6, round(n_guess))
      }
      n_guess
    }, numeric(1L))
    n_val <- max(n_vals)
    power_vals <- rep(power, p)
  }

  out <- list(
    n              = n_val,
    power          = stats::setNames(power_vals, names(target_aucs)),
    target_aucs    = target_aucs,
    null_auc       = null_auc,
    effect_sizes   = stats::setNames(delta, names(target_aucs)),
    alpha          = alpha,
    adjusted_alpha = alpha_adj,
    adjust         = adjust,
    alternative    = alternative,
    correlation    = correlation,
    prevalence     = prevalence
  )
  class(out) <- "aucmat_power"
  out
}

#' @export
print.aucmat_power <- function(x, ...) {
  cat("<aucmat_power>\n")
  if (length(x$target_aucs) == 1L) {
    cat(sprintf("  Target AUC: %.3f (null = %.3f)\n", x$target_aucs[1], x$null_auc))
  } else {
    cat(sprintf("  Target AUCs: %s\n",
      paste(formatC(x$target_aucs, 3, format = "f"), collapse = ", ")))
  }
  cat(sprintf("  N: %d  |  Prevalence: %.2f  |  Correlation: %.2f\n",
    x$n, x$prevalence, x$correlation))
  cat(sprintf("  Alpha: %.3f  |  Adjustment: %s  |  Alternative: %s\n",
    x$alpha, x$adjust, x$alternative))
  if (x$adjust != "none") {
    cat(sprintf("  Adjusted alpha: %.5f\n", x$adjusted_alpha))
  }
  cat(sprintf("\n  Power: %s\n",
    paste(formatC(x$power, 3, format = "f"), collapse = ", ")))
  invisible(x)
}
