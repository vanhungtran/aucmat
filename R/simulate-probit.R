# ==============================================================================
# Latent Probit Simulator
#
# Generates correlated biomarkers with specified AUCs and correlation structure
# via a latent multivariate normal model with a thresholded binary outcome.
# This decouples AUC control from correlation control — both are free parameters
# within the positive-definiteness constraints of the (p+1)×(p+1) matrix.
# ==============================================================================

#' Calibrate latent correlation for a target AUC under the probit model
#'
#' Given a binary outcome with prevalence `pi`, the latent probit model relates
#' each biomarker $X_k$ to a latent continuous $Y^*$ via
#' $\\text{Cor}(X_k, Y^*) = \\rho_k$.  The resulting AUC is a monotone function
#' of $\\rho_k$.  This function numerically finds $\\rho_k$ for a target AUC.
#'
#' @param target_auc Target area under the ROC curve in (0, 1).
#' @param prevalence Proportion of positive class, in (0, 1).
#' @param n_cal Number of observations used for empirical calibration.
#'   Default 100000 gives very stable estimates.
#' @param tol Tolerance for the root finder.  Default 1e-6.
#'
#' @return The latent correlation $\\rho$ that yields the target AUC.
#' @noRd
#' @keywords internal
.calibrate_rho <- function(target_auc, prevalence, n_cal = 100000, tol = 1e-6) {
  if (target_auc <= 0.5) return(0)
  if (target_auc >= 1)   return(1)

  tau <- stats::qnorm(1 - prevalence)

  # Empirical AUC for a given rho, using a large calibration sample
  auc_for_rho <- function(rho) {
    # Bound rho away from ±1 for numerical stability
    rho <- max(-0.99999, min(0.99999, rho))

    # Generate bivariate normal (X, Y*)
    Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    z <- mvtnorm::rmvnorm(n_cal, c(0, 0), Sigma)

    x <- z[, 1L]
    y_star <- z[, 2L]
    y <- as.integer(y_star > tau)

    # Compute empirical AUC via the Mann-Whitney rank method
    pos <- y == 1L
    neg <- !pos
    np <- sum(pos)
    nn <- sum(neg)

    if (np < 2L || nn < 2L) return(NA_real_)

    r <- rank(x, ties.method = "average")
    U <- sum(r[pos]) - np * (np + 1L) / 2L
    as.numeric(U / (np * nn))
  }

  # Find rho via bisection on monotone auc(rho)
  lo <- -0.99999
  hi <-  0.99999

  auc_lo <- auc_for_rho(lo)
  auc_hi <- auc_for_rho(hi)

  if (target_auc <= auc_lo) return(lo)
  if (target_auc >= auc_hi) return(hi)

  for (i in seq_len(60)) {
    mid <- (lo + hi) / 2
    auc_mid <- auc_for_rho(mid)

    if (abs(auc_mid - target_auc) < tol) return(mid)
    if (auc_mid < target_auc) lo <- mid else hi <- mid
  }

  (lo + hi) / 2
}

#' Generate correlated biomarkers with specified AUCs via latent probit model
#'
#' Models the binary outcome as arising from a latent continuous variable
#' Y* ~ N(0,1) with threshold tau = qnorm(1 - pi).
#' Biomarkers (X_1, ..., X_p, Y*) are jointly multivariate normal with
#' a full $(p+1) \\times (p+1)$ correlation matrix $\\Sigma$.
#'
#' This approach simultaneously controls **both** the between-biomarker
#' correlations **and** each biomarker's AUC against the outcome, subject
#' only to the positive-definiteness constraint on $\\Sigma$.
#'
#' @param n Number of observations.
#' @param target_aucs Numeric vector of target AUC values, each in (0, 1).
#' @param corr_matrix A $p \\times p$ target correlation matrix for the
#'   biomarkers.  Must be symmetric and positive definite.
#' @param prevalence Proportion of positive class, in (0, 1).
#' @param n_cal Calibration sample size for numerical root-finding of latent
#'   correlations.  Larger values give more precise AUC targeting at the
#'   cost of calibration time.  Default 100000.
#' @param verify Logical.  If `TRUE` (default), empirically verify the
#'   achieved AUCs and correlations.
#'
#' @return A list with components:
#'   \item{data}{A data.frame with columns \code{X1, ..., Xp} and \code{truth}
#'     (binary outcome, 0/1).}
#'   \item{achieved_aucs}{Empirical AUCs of each biomarker against the outcome.}
#'   \item{achieved_correlations}{Empirical correlation matrix of the biomarkers.}
#'   \item{latent_rhos}{The calibrated latent correlations used for generation.}
#'   \item{Sigma}{The full $(p+1) \\times (p+1)$ correlation matrix used.}
#'
#' @details
#' **Theoretical basis**:
#'
#' Under the latent probit model, a biomarker X_k and the binary outcome
#' Y are linked through rho_k = Cor(X_k, Y*), the latent correlation.
#' For a given prevalence pi, the resulting AUC is a strictly increasing
#' function of rho_k.
#'
#' \deqn{AUC(rho_k, pi) = f(rho_k, pi)}{
#' AUC(rho_k, pi) = f(rho_k, pi)}
#'
#' The function f is not available in closed form (it involves the bivariate
#' normal CDF), so we calibrate rho_k numerically via bisection with a
#' large calibration sample.
#'
#' **Relationship to the binormal model**:
#'
#' The standard binormal model used in `generate_data_analytical()` assumes
#' X_k | Y=0 ~ N(mu_0, sigma^2) and X_k | Y=1 ~ N(mu_1, sigma^2)
#' (equal variance).  This links AUC and Cor(X_k, Y) through a single
#' parameter delta = (mu_1 - mu_0) / sigma.
#'
#' The latent probit model **relaxes this constraint** by allowing the
#' correlation between biomarkers to be specified independently of the
#' biomarker-outcome AUC, through the full (p+1) x (p+1) matrix.
#'
#' **Positive definiteness**:
#'
#' The matrix $\\Sigma$ must be positive definite.  If the user-supplied
#' `corr_matrix` is incompatible with the calibrated latent correlations
#' $\\rho_k$, the function uses the nearest positive-definite matrix
#' (via Higham's algorithm) and issues a warning.
#'
#' @importFrom stats qnorm rbinom
#' @export
#'
#' @examples
#' \donttest{
#' # Three biomarkers with AUCs 0.9, 0.8, 0.7, moderately correlated
#' set.seed(42)
#' sim <- generate_data_probit(
#'   n = 500,
#'   target_aucs = c(0.9, 0.8, 0.7),
#'   corr_matrix = matrix(c(1.0, 0.3, 0.1, 0.3, 1.0, 0.2, 0.1, 0.2, 1.0), 3, 3),
#'   prevalence = 0.3
#' )
#'
#' # Compare targets to achieved values
#' data.frame(
#'   Target  = sim$target_aucs,
#'   Achieved = round(sim$achieved_aucs, 4)
#' )
#' print(round(sim$achieved_correlations, 3))
#' }
generate_data_probit <- function(n,
                                  target_aucs,
                                  corr_matrix,
                                  prevalence,
                                  n_cal   = 100000,
                                  verify  = TRUE) {

  p <- length(target_aucs)

  # ---- Validate inputs ----
  if (n < 1L) stop("n must be a positive integer.")
  if (prevalence <= 0 || prevalence >= 1)
    stop("prevalence must be in (0, 1).")
  if (any(target_aucs <= 0 | target_aucs >= 1))
    stop("target_aucs must each be in (0, 1).")
  if (!is.matrix(corr_matrix) || nrow(corr_matrix) != p || ncol(corr_matrix) != p)
    stop("corr_matrix must be a p x p matrix (p = length(target_aucs)).")
  if (!all(diag(corr_matrix) == 1))
    stop("corr_matrix must have ones on the diagonal.")

  # ---- 1. Calibrate latent correlations ----
  latent_rhos <- vapply(target_aucs, function(auc) {
    .calibrate_rho(auc, prevalence, n_cal = n_cal)
  }, numeric(1L))

  # ---- 2. Build full (p+1) x (p+1) correlation matrix ----
  Sigma <- diag(p + 1L)
  Sigma[1:p, 1:p] <- corr_matrix
  Sigma[1:p, p + 1L] <- latent_rhos
  Sigma[p + 1L, 1:p] <- latent_rhos

  # Check positive definiteness and repair if needed
  eig <- eigen(Sigma, symmetric = TRUE)
  if (any(eig$values <= 0)) {
    warning(
      "The combined (biomarkers + outcome) correlation matrix is not ",
      "positive definite.  Using the nearest positive-definite approximation.  ",
      "This may affect achieved AUCs and correlations."
    )
    eig$values[eig$values < 1e-8] <- 1e-8
    Sigma <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    # Re-normalize to correlation matrix
    d <- sqrt(diag(Sigma))
    Sigma <- Sigma / outer(d, d)
  }

  # ---- 3. Generate from multivariate normal ----
  Z <- mvtnorm::rmvnorm(n, rep(0, p + 1L), Sigma)

  X_mat <- Z[, 1:p, drop = FALSE]
  colnames(X_mat) <- paste0("X", seq_len(p))

  y_star <- Z[, p + 1L]
  tau <- stats::qnorm(1 - prevalence)
  truth <- as.integer(y_star > tau)

  # ---- 4. Assemble output ----
  data_df <- as.data.frame(X_mat)
  data_df$truth <- truth

  # ---- 5. Verify ----
  if (isTRUE(verify)) {
    achieved_aucs <- vapply(seq_len(p), function(j) {
      .compute_single_auc(data_df[[j]], data_df$truth == 1L,
                          data_df$truth == 0L)$auc_raw
    }, numeric(1L))
    achieved_correlations <- stats::cor(X_mat)
  } else {
    achieved_aucs <- rep(NA_real_, p)
    achieved_correlations <- matrix(NA_real_, p, p)
  }

  list(
    data                  = data_df,
    target_aucs           = target_aucs,
    achieved_aucs         = achieved_aucs,
    achieved_correlations = achieved_correlations,
    latent_rhos           = latent_rhos,
    Sigma                 = Sigma
  )
}
