# ==============================================================================
# simulate_auc_copula() — Two-phase Copula + AUC perturbation generator
#
# Phase 1 — Gaussian Copula (class-conditional):
#   Rank-transform → probit → MVN(0, R) per class → probit⁻¹ → scale
#   Preserves correlation structure and produces normal-like marginals.
#
# Phase 2 — Iterative AUC perturbation:
#   For each feature where achieved AUC deviates from target:
#     Add minimal class-conditional shift scaled by feature variance.
#     Small steps (decreasing step size) preserve correlations.
#     Repeat until convergence or max iterations.
#
# This cleanly separates:
#   - Correlation structure  → Phase 1 (copula owns it)
#   - Discriminative signal  → Phase 2 (perturbation owns it)
#
# Compared to generate_data_analytical() (sequential binormal):
#   - No AUC-vs-correlation tradeoff
#   - No sequential signal erosion
#   - More accurate correlation preservation
# ==============================================================================

#' Simulate correlated biomarkers with two-phase Copula-AUC generation
#'
#' Generates biomarkers in two phases: (1) a class-conditional Gaussian copula
#' captures the full correlation structure, then (2) iterative perturbation
#' fine-tunes each biomarker's AUC against the outcome while minimally
#' disturbing correlations.
#'
#' @param n Number of observations.
#' @param prevalence Proportion of positive class, in (0, 1).
#' @param target_aucs Numeric vector of target AUC values, each in (0, 1).
#' @param corr_matrix A p x p target correlation matrix.
#' @param n_iterations Number of perturbation iterations.  Default 5.
#' @param step_decay Rate at which perturbation step size decreases.
#'   Default 0.5.
#' @param convergence_tol Stop perturbing features where
#'   |achieved - target| < tol.  Default 0.02.
#' @param verify Logical.  If TRUE (default), verify achieved AUCs and
#'   correlations.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class `aucmat_simulation` with components `data`,
#'   `target_aucs`, `achieved_aucs`, `requested_correlation`,
#'   `achieved_correlation`, `phase1_aucs`, `n`, `prevalence`.
#'
#' @details
#' **Comparison to other aucmat simulators:**
#'
#' \describe{
#'   \item{`generate_data_analytical()`}{Sequential binormal decomposition.
#'     Good AUC control, weak correlation control.}
#'   \item{`simulate_auc_matrix()`}{Class-conditional MVN.  Good correlation
#'     control via parametrized structures, normal marginals only.}
#'   \item{`generate_data_probit()`}{Latent probit model.  Simultaneously
#'     controls AUC and correlation via (p+1)x(p+1) matrix.}
#'   \item{`simulate_auc_copula()`}{**This function.**  Gaussian copula +
#'     iterative AUC perturbation.  Strong correlation preservation
#'     without the AUC-correlation tradeoff.}
#' }
#'
#' @export
#' @importFrom stats qnorm rnorm cor
#' @examples
#' \donttest{
#' set.seed(42)
#' sim <- simulate_auc_copula(
#'   n = 500, prevalence = 0.3,
#'   target_aucs = c(0.85, 0.75, 0.65),
#'   corr_matrix = matrix(c(1.0, 0.4, 0.2, 0.4, 1.0, 0.3, 0.2, 0.3, 1.0), 3, 3)
#' )
#' data.frame(
#'   target   = sim$target_aucs,
#'   phase1   = round(sim$phase1_aucs, 3),
#'   achieved = round(sim$achieved_aucs, 3)
#' )
#' }
simulate_auc_copula <- function(n,
                                 prevalence,
                                 target_aucs,
                                 corr_matrix,
                                 n_iterations    = 5,
                                 step_decay      = 0.5,
                                 convergence_tol = 0.02,
                                 verify          = TRUE,
                                 seed            = NULL) {

  p <- length(target_aucs)

  # ---- Seed handling ----
  with_seed(seed, {
    .simulate_auc_copula_body(n, prevalence, target_aucs, corr_matrix,
      n_iterations, step_decay, convergence_tol, verify)
  })
}

.simulate_auc_copula_body <- function(n, prevalence, target_aucs, corr_matrix,
                                       n_iterations, step_decay,
                                       convergence_tol, verify) {
  p <- length(target_aucs)

  # ---- Validate ----
  if (n < 6L) stop("n must be at least 6.")
  if (prevalence <= 0 || prevalence >= 1)
    stop("prevalence must be in (0, 1).")
  if (any(target_aucs <= 0 | target_aucs >= 1))
    stop("target_aucs must each be in (0, 1).")
  if (!is.matrix(corr_matrix) || nrow(corr_matrix) != p || ncol(corr_matrix) != p)
    stop("corr_matrix must be p x p.")
  if (n_iterations < 0L) stop("n_iterations must be >= 0.")

  # ---- 1. Generate outcome ----
  n_pos <- round(n * prevalence)
  n_pos <- max(1L, min(n - 1L, n_pos))
  n_neg <- n - n_pos
  truth <- sample(c(rep(1L, n_pos), rep(0L, n_neg)))
  pos <- truth == 1L
  neg <- !pos

  # Clamp target AUCs for qnorm stability
  target_aucs_clamped <- pmin(pmax(target_aucs, 0.5001), 0.9999)

  # ---- 2. Phase 1: Class-conditional Gaussian Copula ----
  # For each class independently:
  #   a. Build (p+1)x(p+1) correlation matrix (biomarkers + latent outcome)
  #   b. Draw from MVN(0, R)
  #   c. Rank-transform to uniform (copula step)
  #   d. Probit-transform to normal scores
  #   e. Apply class-conditional mean shift for AUC signal

  X_phase1 <- matrix(NA_real_, n, p)

  # Binormal delta per biomarker
  delta <- sqrt(2) * stats::qnorm(target_aucs_clamped)

  # Approximate latent correlation from AUC
  latent_rhos_corr <- delta / sqrt(1 + delta^2)

  for (cl in c(0, 1)) {
    cl_mask <- if (cl == 1) pos else neg
    n_cl <- sum(cl_mask)

    # Build (p+1) x (p+1) correlation matrix
    R_work <- diag(p + 1L)
    R_work[1:p, 1:p] <- corr_matrix
    R_work[1:p, p + 1L] <- latent_rhos_corr
    R_work[p + 1L, 1:p] <- latent_rhos_corr

    # Ensure positive definite (reuse shared PD projection)
    R_work <- .nearest_pd_correlation(R_work)$matrix

    # Draw from multivariate normal
    Z <- mvtnorm::rmvnorm(n_cl, rep(0, p + 1L), R_work)

    # Copula step: rank-transform to uniform, then probit to normal scores
    # This preserves the rank correlation structure exactly
    U <- apply(Z[, 1:p, drop = FALSE], 2, function(z) {
      rank(z, ties.method = "average") / (length(z) + 1)
    })
    Z_norm <- apply(U, 2, stats::qnorm)
    Z_norm <- pmin(pmax(Z_norm, -5), 5)  # clamp extremes

    # Apply class-conditional mean shift for AUC signal
    class_shift <- if (cl == 1) delta * (n_neg / n) else delta * (-n_pos / n)
    Z_norm <- sweep(Z_norm, 2, class_shift, "+")

    X_phase1[cl_mask, ] <- Z_norm
  }

  colnames(X_phase1) <- paste0("X", seq_len(p))

  # ---- 3. Verify Phase 1 AUCs ----
  phase1_aucs <- vapply(seq_len(p), function(j) {
    .compute_single_auc(X_phase1[, j], pos, neg)$auc_raw
  }, numeric(1L))

  # ---- 4. Phase 2: Iterative AUC perturbation ----
  X_final <- X_phase1

  if (n_iterations > 0) {
    for (iter in seq_len(n_iterations)) {
      current_aucs <- vapply(seq_len(p), function(j) {
        .compute_single_auc(X_final[, j], pos, neg)$auc_raw
      }, numeric(1L))

      auc_gap <- target_aucs - current_aucs
      adjust_mask <- abs(auc_gap) > convergence_tol
      if (!any(adjust_mask)) break

      step <- step_decay^iter
      feature_sd <- apply(X_final, 2, stats::sd)

      # Order by largest gap relative to scale
      importance <- abs(auc_gap) / pmax(feature_sd, 0.01)
      order_j <- order(importance, decreasing = TRUE)

      for (j in order_j) {
        if (!adjust_mask[j]) next

        target_j <- target_aucs_clamped[j]
        current_j <- pmin(pmax(current_aucs[j], 0.5001), 0.9999)

        delta_target  <- sqrt(2) * stats::qnorm(target_j)
        delta_current <- sqrt(2) * stats::qnorm(current_j)
        delta_diff    <- delta_target - delta_current

        shift <- delta_diff * feature_sd[j] * step

        X_final[pos, j] <- X_final[pos, j] + shift * (n_neg / n)
        X_final[neg, j] <- X_final[neg, j] - shift * (n_pos / n)
      }
    }
  }

  # ---- 5. Standardise (affine — preserves AUC and correlation) ----
  X_final <- scale(X_final)
  colnames(X_final) <- paste0("X", seq_len(p))

  # ---- 6. Assemble output ----
  data_df <- as.data.frame(X_final)
  data_df$truth <- as.integer(pos)

  # ---- 7. Verify final ----
  if (isTRUE(verify)) {
    achieved_aucs <- vapply(seq_len(p), function(j) {
      .compute_single_auc(data_df[[j]], pos, neg)$auc_raw
    }, numeric(1L))

    X_centered <- X_final
    X_centered[pos, ] <- sweep(X_final[pos, , drop = FALSE], 2,
                               colMeans(X_final[pos, , drop = FALSE]), "-")
    X_centered[neg, ] <- sweep(X_final[neg, , drop = FALSE], 2,
                               colMeans(X_final[neg, , drop = FALSE]), "-")
    achieved_correlations <- stats::cor(X_centered)
  } else {
    achieved_aucs <- rep(NA_real_, p)
    achieved_correlations <- matrix(NA_real_, p, p)
  }

  out <- list(
    data                  = data_df,
    target_aucs           = target_aucs,
    phase1_aucs           = stats::setNames(phase1_aucs, colnames(X_phase1)),
    achieved_aucs         = stats::setNames(achieved_aucs, colnames(X_final)),
    requested_correlation = corr_matrix,
    achieved_correlation  = achieved_correlations,
    n_iterations_used     = n_iterations,
    n = n,
    prevalence            = n_pos / n
  )
  class(out) <- "aucmat_simulation"
  out
}
