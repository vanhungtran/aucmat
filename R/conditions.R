# ==============================================================================
# Structured Condition Classes
# ==============================================================================

#' Invalid outcome
#' @noRd
#' @keywords internal
aucmat_invalid_outcome <- function(message, call = sys.call(-1L)) {
  structure(list(message = message, call = call),
            class = c("aucmat_invalid_outcome", "error", "condition"))
}

#' Invalid correlation specification
#' @noRd
#' @keywords internal
aucmat_invalid_correlation <- function(message, call = sys.call(-1L)) {
  structure(list(message = message, call = call),
            class = c("aucmat_invalid_correlation", "error", "condition"))
}

#' Infeasible simulation targets
#' @noRd
#' @keywords internal
aucmat_infeasible_targets <- function(message, call = sys.call(-1L)) {
  structure(list(message = message, call = call),
            class = c("aucmat_infeasible_targets", "error", "condition"))
}

#' Insufficient sample
#' @noRd
#' @keywords internal
aucmat_insufficient_sample <- function(message, call = sys.call(-1L)) {
  structure(list(message = message, call = call),
            class = c("aucmat_insufficient_sample", "error", "condition"))
}

#' Resampling failure
#' @noRd
#' @keywords internal
aucmat_resampling_failure <- function(message, call = sys.call(-1L)) {
  structure(list(message = message, call = call),
            class = c("aucmat_resampling_failure", "error", "condition"))
}

# ---- RNG helpers ----

#' Set seed locally, restoring global state on exit
#' @param seed Integer seed or NULL.
#' @param expr Expression to evaluate.
#' @return Result of expr.
#' @noRd
#' @keywords internal
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old <- globalenv()$.Random.seed
  on.exit({
    if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv()))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  expr
}
