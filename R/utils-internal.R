# ==============================================================================
# Internal helpers (not exported)
# ==============================================================================

#' Normalize a binary response vector to a 2-level factor
#'
#' Accepts logical, numeric 0/1, factor, or character vectors and returns a
#' factor with exactly 2 levels. The first level is always the negative class.
#'
#' @param y A binary outcome vector (logical, numeric 0/1, factor, or character).
#' @param positive Optional positive class label. If NULL, the second level of
#'   the factor is treated as positive.
#'
#' @return A factor with 2 levels: c(negative, positive).
#' @noRd
#' @keywords internal
.normalize_binary_y <- function(y, positive = NULL) {
  if (is.logical(y)) {
    return(factor(y, levels = c(FALSE, TRUE), labels = c("neg", "pos")))
  }

  if (is.numeric(y)) {
    u <- unique(y)
    u <- u[!is.na(u)]
    if (length(u) == 2L && all(u %in% c(0, 1))) {
      return(factor(y, levels = c(0, 1), labels = c("neg", "pos")))
    }
  }

  if (!is.factor(y)) y <- factor(y)

  if (nlevels(y) != 2L) {
    stop("y must have exactly 2 unique non-missing values/classes.")
  }

  lev <- levels(y)

  if (is.null(positive)) {
    return(y)
  }

  if (!(positive %in% lev)) {
    stop("`positive` must be one of levels(y).")
  }

  neg <- setdiff(lev, positive)
  if (length(neg) != 1L) {
    stop("Could not infer the negative class from y.")
  }

  factor(y, levels = c(neg, positive))
}

#' Convert a 2-level factor to numeric 0/1
#'
#' @param y A 2-level factor.
#' @return Integer vector with 1 for the second (positive) level, 0 otherwise.
#' @noRd
#' @keywords internal
.binary_factor_to_numeric <- function(y) {
  lev <- levels(y)
  as.integer(y == lev[2L])
}

#' NULL-coalescing infix operator
#'
#' Defined locally rather than relying on base R's `%||%` (only available
#' from R 4.4.0) since the package supports R >= 4.0.0.
#'
#' @param x,y Values; `y` is returned when `x` is `NULL`.
#' @noRd
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
