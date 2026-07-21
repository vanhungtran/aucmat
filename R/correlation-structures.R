# ==============================================================================
# Correlation Structure Builders
#
# Constructs a p x p target correlation matrix from a named structure and
# a small set of parameters, so users do not have to hand-build a full
# matrix for common dependence patterns.
# ==============================================================================

#' Build a p x p correlation matrix from a named structure
#'
#' @param p Number of biomarkers.
#' @param structure One of `"user"`, `"exchangeable"`, `"ar1"`, `"block"`.
#' @param correlation For `"user"`: the full p x p matrix.  For
#'   `"exchangeable"` and `"ar1"`: a single numeric correlation in `(-1, 1)`.
#'   Ignored for `"block"`.
#' @param block_sizes Integer vector of block sizes summing to `p`.  Required
#'   for `structure = "block"`.
#' @param rho_within Numeric within-block correlation.  Required for
#'   `structure = "block"`.
#' @param rho_between Numeric between-block correlation.  Required for
#'   `structure = "block"`.
#'
#' @return A symmetric p x p matrix with unit diagonal.
#' @noRd
#' @keywords internal
.build_correlation_matrix <- function(p, structure,
                                       correlation = NULL,
                                       block_sizes = NULL,
                                       rho_within = NULL,
                                       rho_between = NULL) {
  structure <- match.arg(structure,
    c("user", "exchangeable", "ar1", "block"))

  if (structure == "user") {
    if (is.null(correlation) || !is.matrix(correlation))
      stop("structure = 'user' requires `correlation` to be a p x p matrix.")
    if (nrow(correlation) != p || ncol(correlation) != p)
      stop(sprintf("`correlation` must be %d x %d (length(target_aucs) = %d).",
                   p, p, p))
    if (!isTRUE(all.equal(correlation, t(correlation), tolerance = 1e-8)))
      stop("`correlation` must be symmetric.")
    if (!isTRUE(all.equal(diag(correlation), rep(1, p), tolerance = 1e-8)))
      stop("`correlation` must have a unit diagonal.")
    return(correlation)
  }

  if (structure == "exchangeable") {
    if (is.null(correlation) || !is.numeric(correlation) ||
        length(correlation) != 1L || is.na(correlation) ||
        correlation <= -1 || correlation >= 1) {
      stop("structure = 'exchangeable' requires `correlation` to be a ",
           "single value in (-1, 1).")
    }
    R <- matrix(correlation, p, p)
    diag(R) <- 1
    return(R)
  }

  if (structure == "ar1") {
    if (is.null(correlation) || !is.numeric(correlation) ||
        length(correlation) != 1L || is.na(correlation) ||
        correlation <= -1 || correlation >= 1) {
      stop("structure = 'ar1' requires `correlation` to be a single value ",
           "in (-1, 1).")
    }
    idx <- seq_len(p)
    R <- correlation ^ abs(outer(idx, idx, "-"))
    return(R)
  }

  # structure == "block"
  if (is.null(block_sizes) || !is.numeric(block_sizes) ||
      any(block_sizes < 1) || sum(block_sizes) != p) {
    stop("structure = 'block' requires `block_sizes` to be positive ",
         "integers summing to length(target_aucs) (", p, ").")
  }
  if (is.null(rho_within) || !is.numeric(rho_within) ||
      length(rho_within) != 1L || is.na(rho_within) ||
      rho_within <= -1 || rho_within >= 1) {
    stop("structure = 'block' requires `rho_within` to be a single value ",
         "in (-1, 1).")
  }
  if (is.null(rho_between) || !is.numeric(rho_between) ||
      length(rho_between) != 1L || is.na(rho_between) ||
      rho_between <= -1 || rho_between >= 1) {
    stop("structure = 'block' requires `rho_between` to be a single value ",
         "in (-1, 1).")
  }

  block_id <- rep(seq_along(block_sizes), times = block_sizes)
  R <- ifelse(outer(block_id, block_id, "=="), rho_within, rho_between)
  diag(R) <- 1
  R
}
