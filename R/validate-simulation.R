# ==============================================================================
# validate_simulation() — Repeated-draw calibration check
#
# Confirms that simulate_auc_matrix() actually hits its targets across many
# independent draws, rather than trusting a single sample.  Reports bias,
# RMSE, Monte Carlo standard error, and target-interval hit rates for both
# AUCs and pairwise correlations.
# ==============================================================================

#' Validate a `simulate_auc_matrix()` specification by repeated simulation
#'
#' Repeats a [simulate_auc_matrix()] call `times` times with independent
#' seeds and reports how closely the achieved AUCs and correlations track
#' their targets, with Monte Carlo-aware uncertainty.  A single simulated
#' draw is never sufficient evidence that a simulator is calibrated; this
#' function is the check.
#'
#' @param ... Arguments forwarded to [simulate_auc_matrix()] on every
#'   replicate (e.g. `n`, `prevalence`, `target_aucs`, `correlation`,
#'   `structure`, `feasibility`).  The replicate stream is controlled by
#'   `validate_simulation()`'s own `seed` argument, below.
#' @param times Number of independent replicates.  Default 100.
#' @param tolerance Named list with `auc` and `correlation` tolerances used
#'   for the target-interval hit rate.  Default
#'   `list(auc = 0.02, correlation = 0.05)`.
#' @param seed Optional integer seed controlling the replicate stream.
#'   Each replicate uses seed `seed + replicate_index - 1` when supplied.
#'
#' @return An object of class `aucmat_simulation_validation`, a list with
#'   components `auc_bias`, `auc_rmse`, `auc_mc_se`, `auc_hit_rate`,
#'   `corr_bias`, `corr_rmse`, `corr_mc_se`, `corr_hit_rate`, `n_requested`,
#'   `n_ok`, `n_failed`, `settings`.
#'
#' @export
#' @examples
#' \donttest{
#' val <- validate_simulation(
#'   n = 200, prevalence = 0.3,
#'   target_aucs = c(0.8, 0.7), correlation = 0.3, structure = "exchangeable",
#'   times = 50, seed = 1
#' )
#' print(val)
#' }
validate_simulation <- function(..., times = 100,
                                 tolerance = list(auc = 0.02, correlation = 0.05),
                                 seed = NULL) {

  if (!is.numeric(times) || length(times) != 1L || times < 2L) {
    stop("times must be a single integer >= 2.")
  }
  times <- as.integer(times)

  dots <- list(...)
  target_aucs <- dots$target_aucs
  if (is.null(target_aucs)) stop("target_aucs must be supplied.")
  p <- length(target_aucs)

  auc_mat  <- matrix(NA_real_, nrow = times, ncol = p)
  corr_arr <- vector("list", times)
  n_ok <- 0L
  n_failed <- 0L

  for (i in seq_len(times)) {
    rep_seed <- if (is.null(seed)) NULL else seed + i - 1L
    args <- c(dots, list(seed = rep_seed, verify = TRUE))

    res <- try(do.call(simulate_auc_matrix, args), silent = TRUE)

    if (inherits(res, "try-error")) {
      n_failed <- n_failed + 1L
      next
    }

    n_ok <- n_ok + 1L
    auc_mat[i, ] <- res$achieved_aucs
    corr_arr[[i]] <- res$achieved_correlation
  }

  if (n_ok < 2L) {
    stop("Fewer than 2 successful replicates; cannot validate. ",
         "n_failed = ", n_failed, " of ", times, " requested.")
  }

  ok_rows <- !apply(auc_mat, 1, function(r) all(is.na(r)))
  auc_mat <- auc_mat[ok_rows, , drop = FALSE]

  # ---- AUC bias / RMSE / MC SE / hit rate ----
  auc_bias <- colMeans(auc_mat) - target_aucs
  auc_rmse <- sqrt(colMeans(sweep(auc_mat, 2, target_aucs, "-")^2))
  auc_mc_se <- apply(auc_mat, 2, stats::sd) / sqrt(n_ok)
  auc_hit_rate <- colMeans(
    abs(sweep(auc_mat, 2, target_aucs, "-")) < tolerance$auc
  )
  names(auc_bias) <- names(auc_rmse) <- names(auc_mc_se) <-
    names(auc_hit_rate) <- paste0("X", seq_len(p))

  # ---- Correlation bias / RMSE / MC SE / hit rate (upper triangle) ----
  corr_result <- NULL
  if (p >= 2L) {
    ok_corr <- corr_arr[!vapply(corr_arr, is.null, logical(1))]
    ut <- upper.tri(ok_corr[[1]])
    target_R <- .build_correlation_matrix(
      p, dots$structure %||% "user",
      correlation = dots$correlation,
      block_sizes = dots$block_sizes,
      rho_within = dots$rho_within, rho_between = dots$rho_between
    )
    target_vec <- target_R[ut]

    achieved_mat <- do.call(rbind, lapply(ok_corr, function(m) m[ut]))
    corr_bias <- colMeans(achieved_mat) - target_vec
    corr_rmse <- sqrt(colMeans(sweep(achieved_mat, 2, target_vec, "-")^2))
    corr_mc_se <- apply(achieved_mat, 2, stats::sd) / sqrt(nrow(achieved_mat))
    corr_hit_rate <- colMeans(
      abs(sweep(achieved_mat, 2, target_vec, "-")) < tolerance$correlation
    )

    corr_result <- list(
      bias = corr_bias, rmse = corr_rmse,
      mc_se = corr_mc_se, hit_rate = corr_hit_rate
    )
  }

  out <- list(
    auc_bias = auc_bias, auc_rmse = auc_rmse,
    auc_mc_se = auc_mc_se, auc_hit_rate = auc_hit_rate,
    correlation = corr_result,
    n_requested = times, n_ok = n_ok, n_failed = n_failed,
    settings = list(tolerance = tolerance, seed = seed)
  )
  class(out) <- "aucmat_simulation_validation"
  out
}
