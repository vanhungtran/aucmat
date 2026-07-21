# ==============================================================================
# Nearest Positive-Semidefinite Correlation Projection
#
# Shared Higham-style eigenvalue-clipping projection used when a requested
# latent correlation matrix is not positive definite.  Reports adjustment
# diagnostics so callers can report requested vs. achieved matrices rather
# than silently repairing them.
# ==============================================================================

#' Project a matrix to the nearest positive-semidefinite correlation matrix
#'
#' Uses eigenvalue clipping (a simple Higham-style projection) to find a
#' positive-semidefinite matrix with unit diagonal close to the input.  This
#' is not the full alternating-projections Higham (2002) algorithm, but a
#' single-step eigenvalue floor followed by re-normalisation to a correlation
#' matrix; adequate for the small perturbations expected from latent-AUC
#' calibration infeasibility.
#'
#' @param Sigma A symmetric numeric matrix (assumed square).
#' @param eig_floor Minimum eigenvalue after clipping.  Default `1e-8`.
#'
#' @return A list with:
#'   \item{matrix}{The projected positive-semidefinite correlation matrix.}
#'   \item{adjusted}{`TRUE` if any eigenvalue was below `eig_floor`.}
#'   \item{frobenius_adjustment}{Frobenius norm of the difference between
#'     the input and the projected matrix.}
#'   \item{max_abs_adjustment}{Maximum absolute element-wise difference.}
#'   \item{min_eigenvalue}{The smallest eigenvalue of the input matrix,
#'     before projection.}
#' @noRd
#' @keywords internal
.nearest_pd_correlation <- function(Sigma, eig_floor = 1e-8) {
  Sigma <- as.matrix(Sigma)
  eig <- eigen(Sigma, symmetric = TRUE)
  min_eig <- min(eig$values)

  if (min_eig > eig_floor) {
    return(list(
      matrix = Sigma, adjusted = FALSE,
      frobenius_adjustment = 0, max_abs_adjustment = 0,
      min_eigenvalue = min_eig
    ))
  }

  vals <- eig$values
  vals[vals < eig_floor] <- eig_floor
  proj <- eig$vectors %*% diag(vals) %*% t(eig$vectors)

  # Re-normalise to a correlation matrix (unit diagonal)
  d <- sqrt(diag(proj))
  proj <- proj / outer(d, d)
  proj <- (proj + t(proj)) / 2  # enforce exact symmetry
  diag(proj) <- 1

  list(
    matrix = proj, adjusted = TRUE,
    frobenius_adjustment = sqrt(sum((proj - Sigma)^2)),
    max_abs_adjustment = max(abs(proj - Sigma)),
    min_eigenvalue = min_eig
  )
}
