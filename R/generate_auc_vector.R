#' Generate a continuous score vector for a target AUC
#'
#' Construct a numeric score vector for a binary outcome so that the empirical
#' ROC AUC matches a requested value exactly when that value is attainable for
#' the observed numbers of positive and negative samples. When the requested AUC
#' is not attainable on the finite-sample grid, the function can either use the
#' nearest attainable value or stop with an error.
#'
#' @param y Binary outcome vector. Accepted inputs are logical, numeric 0/1,
#'   factor, or character with exactly two non-missing classes.
#' @param target_auc Target empirical AUC in `[0, 1]`.
#' @param positive Optional positive class label when `y` is a factor or
#'   character vector. By default, the second factor level is treated as the
#'   positive class.
#' @param unattainable One of `"nearest"` or `"error"`. If the requested AUC is
#'   not exactly attainable for the current numbers of cases and controls,
#'   `"nearest"` uses the closest attainable AUC and warns, while `"error"`
#'   stops.
#' @param seed Optional seed used only when shuffling scores within the positive
#'   and negative classes.
#' @param shuffle_within_class Logical. If `TRUE`, randomly permute the assigned
#'   scores within each class while preserving the achieved AUC. Default `TRUE`.
#'
#' @return A list with components:
#'   \item{x}{Numeric score vector with the same length as `y`. Missing values in
#'   `y` are propagated as `NA` in `x`.}
#'   \item{requested_auc}{The requested AUC.}
#'   \item{achieved_auc}{The empirical AUC achieved by `x` on the non-missing
#'   observations.}
#'   \item{auc_step}{The smallest attainable AUC increment,
#'   `1 / (n_pos * n_neg)`.}
#'   \item{n_pos}{Number of positive samples used.}
#'   \item{n_neg}{Number of negative samples used.}
#'
#' @examples
#' set.seed(1)
#' y <- rbinom(1000, size = 1, prob = 0.3)
#' out <- generate_auc_vector(y, target_auc = 0.8)
#' head(out$x)
#' out$achieved_auc
#'
#' @export
generate_auc_vector <- function(
    y,
    target_auc,
    positive = NULL,
    unattainable = c("nearest", "error"),
    seed = NULL,
    shuffle_within_class = TRUE
) {
  unattainable <- match.arg(unattainable)

  if (!is.numeric(target_auc) || length(target_auc) != 1L || is.na(target_auc) ||
      target_auc < 0 || target_auc > 1) {
    stop("target_auc must be a single numeric value in [0, 1].")
  }

  if (!is.logical(shuffle_within_class) || length(shuffle_within_class) != 1L) {
    stop("shuffle_within_class must be TRUE or FALSE.")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("seed must be NULL or a single non-missing numeric value.")
    }
    set.seed(seed)
  }

  keep <- !is.na(y)
  if (!any(keep)) stop("y contains only NA values.")

  y_used <- .normalize_binary_y(y[keep], positive = positive)
  lev <- levels(y_used)
  pos_idx_used <- which(y_used == lev[2L])
  neg_idx_used <- which(y_used == lev[1L])
  n_pos <- length(pos_idx_used)
  n_neg <- length(neg_idx_used)

  if (n_pos == 0L || n_neg == 0L) {
    stop("y must contain at least one positive and one negative observation.")
  }

  total_pairs <- n_pos * n_neg
  auc_step <- 1 / total_pairs
  target_u_raw <- target_auc * total_pairs
  target_u_rounded <- round(target_u_raw)

  if (unattainable == "error" && abs(target_u_raw - target_u_rounded) > 1e-10) {
    stop(
      sprintf(
        paste0(
          "Requested AUC %.10f is not exactly attainable with n_pos = %d and ",
          "n_neg = %d. The AUC grid step is %.10f."
        ),
        target_auc, n_pos, n_neg, auc_step
      )
    )
  }

  target_u <- as.integer(max(0, min(total_pairs, target_u_rounded)))
  achieved_auc <- target_u / total_pairs

  if (unattainable == "nearest" && abs(target_auc - achieved_auc) > 1e-10) {
    warning(
      sprintf(
        paste0(
          "Requested AUC %.10f is not exactly attainable with n_pos = %d and ",
          "n_neg = %d. Using nearest attainable AUC %.10f instead."
        ),
        target_auc, n_pos, n_neg, achieved_auc
      ),
      call. = FALSE
    )
  }

  min_rank_sum <- n_pos * (n_pos + 1L) / 2
  target_rank_sum <- min_rank_sum + target_u
  pos_ranks <- .build_rank_subset(
    n = length(y_used),
    k = n_pos,
    target_sum = target_rank_sum
  )
  neg_ranks <- setdiff(seq_len(length(y_used)), pos_ranks)

  ordered_scores <- stats::qnorm((seq_len(length(y_used)) - 0.5) / length(y_used))
  x_used <- numeric(length(y_used))

  if (shuffle_within_class) {
    pos_assign <- sample(pos_idx_used)
    neg_assign <- sample(neg_idx_used)
  } else {
    pos_assign <- pos_idx_used
    neg_assign <- neg_idx_used
  }

  x_used[pos_assign] <- ordered_scores[pos_ranks]
  x_used[neg_assign] <- ordered_scores[neg_ranks]

  x <- rep(NA_real_, length(y))
  x[keep] <- x_used
  if (!is.null(names(y))) names(x) <- names(y)

  structure(
    list(
      x = x,
      requested_auc = target_auc,
      achieved_auc = .empirical_auc_from_scores(y_used, x_used),
      auc_step = auc_step,
      n_pos = n_pos,
      n_neg = n_neg
    ),
    class = "aucmat_auc_vector"
  )
}

#' Generate a continuous score vector for target AUC and correlation
#'
#' Construct a numeric score vector for a binary outcome with an empirical AUC
#' matched exactly when attainable, then tune a strictly increasing Box-Cox
#' transform so the Pearson correlation between `y` and `x` matches a requested
#' value as closely as possible without changing the AUC.
#'
#' Theory:
#' - The empirical AUC depends only on the rank ordering of the score vector.
#' - Pearson correlation depends on the numeric spacing of the score values.
#' - A strictly increasing transform preserves ranks, so it preserves AUC.
#' - Therefore the method first constructs a base vector with the requested AUC,
#'   then tunes a monotone Box-Cox transform to move the Pearson correlation
#'   toward the requested target without changing the AUC.
#' - Not every `(AUC, correlation)` pair is achievable for a fixed `y`, because
#'   the monotone transform can change spacing but cannot change ordering.
#'
#' @param y Binary outcome vector. Accepted inputs are logical, numeric 0/1,
#'   factor, or character with exactly two non-missing classes.
#' @param target_auc Target empirical AUC in `[0, 1]`.
#' @param target_cor Target Pearson correlation between `y` and the generated
#'   score vector.
#' @param positive Optional positive class label when `y` is a factor or
#'   character vector. By default, the second factor level is treated as the
#'   positive class.
#' @param auc_unattainable One of `"nearest"` or `"error"`. Passed to
#'   [generate_auc_vector()] when the requested AUC is not exactly attainable on
#'   the finite-sample grid.
#' @param cor_unattainable One of `"nearest"` or `"error"`. If the requested
#'   correlation is outside the achievable range of the monotone transform
#'   family, `"nearest"` returns the closest match and warns, while `"error"`
#'   stops.
#' @param seed Optional seed used only when constructing the base AUC-matched
#'   score vector.
#' @param shuffle_within_class Logical. Passed to [generate_auc_vector()].
#' @param lambda_bounds Numeric length-2 vector giving the search interval for
#'   the Box-Cox transform parameter. Default `c(-8, 8)`.
#' @param lambda_grid_n Number of grid points used to bracket the correlation
#'   search before local refinement. Default `401`.
#' @param cor_tol Numeric tolerance used to decide whether the requested
#'   correlation was achieved. Default `1e-6`.
#'
#' @return A list with components:
#'   \item{x}{Numeric score vector with the same length as `y`. Missing values in
#'   `y` are propagated as `NA` in `x`.}
#'   \item{requested_auc}{The requested AUC.}
#'   \item{requested_cor}{The requested correlation.}
#'   \item{achieved_auc}{The empirical AUC achieved by `x`.}
#'   \item{achieved_cor}{The Pearson correlation achieved by `x`.}
#'   \item{lambda}{The fitted Box-Cox transform parameter.}
#'   \item{cor_range}{Length-2 numeric vector with the approximate minimum and
#'   maximum correlations achievable over `lambda_bounds`.}
#'   \item{base_vector}{The object returned by [generate_auc_vector()] before the
#'   monotone transform was applied.}
#'
#' @examples
#' set.seed(1)
#' y <- rbinom(200, size = 1, prob = 0.3)
#' out <- generate_auc_cor_vector(y, target_auc = 0.8, target_cor = 0.45)
#' c(out$achieved_auc, out$achieved_cor)
#'
#' @export
generate_auc_cor_vector <- function(
    y,
    target_auc,
    target_cor,
    positive = NULL,
    auc_unattainable = c("nearest", "error"),
    cor_unattainable = c("nearest", "error"),
    seed = NULL,
    shuffle_within_class = TRUE,
    lambda_bounds = c(-8, 8),
    lambda_grid_n = 401,
    cor_tol = 1e-6
) {
  auc_unattainable <- match.arg(auc_unattainable)
  cor_unattainable <- match.arg(cor_unattainable)

  if (!is.numeric(target_cor) || length(target_cor) != 1L || is.na(target_cor) ||
      target_cor < -1 || target_cor > 1) {
    stop("target_cor must be a single numeric value in [-1, 1].")
  }

  if (!is.numeric(lambda_bounds) || length(lambda_bounds) != 2L ||
      any(is.na(lambda_bounds)) || lambda_bounds[1] >= lambda_bounds[2]) {
    stop("lambda_bounds must be a numeric length-2 vector with lower < upper.")
  }

  if (!is.numeric(lambda_grid_n) || length(lambda_grid_n) != 1L || is.na(lambda_grid_n) ||
      lambda_grid_n < 3) {
    stop("lambda_grid_n must be a single integer >= 3.")
  }
  lambda_grid_n <- as.integer(lambda_grid_n)

  if (!is.numeric(cor_tol) || length(cor_tol) != 1L || is.na(cor_tol) || cor_tol < 0) {
    stop("cor_tol must be a single non-negative numeric value.")
  }

  base_vector <- generate_auc_vector(
    y = y,
    target_auc = target_auc,
    positive = positive,
    unattainable = auc_unattainable,
    seed = seed,
    shuffle_within_class = shuffle_within_class
  )

  keep <- !is.na(base_vector$x)
  y_used <- .normalize_binary_y(y[keep], positive = positive)
  y_num <- .binary_factor_to_numeric(y_used)
  x_base <- base_vector$x[keep]
  x_shifted <- x_base - min(x_base) + 1

  cor_at_lambda <- function(lambda) {
    x_lambda <- .box_cox_increasing(x_shifted, lambda)
    stats::cor(y_num, x_lambda)
  }

  lambda_grid <- seq(lambda_bounds[1], lambda_bounds[2], length.out = lambda_grid_n)
  cor_grid <- vapply(lambda_grid, cor_at_lambda, numeric(1))
  cor_range <- range(cor_grid)

  if (target_cor < cor_range[1] - cor_tol || target_cor > cor_range[2] + cor_tol) {
    message_text <- sprintf(
      paste0(
        "Requested correlation %.10f is outside the approximate achievable ",
        "range [%.10f, %.10f] for this AUC and transform family."
      ),
      target_cor, cor_range[1], cor_range[2]
    )
    if (cor_unattainable == "error") stop(message_text)
    warning(message_text, call. = FALSE)
  }

  closest_idx <- which.min(abs(cor_grid - target_cor))
  lower_idx <- max(1L, closest_idx - 1L)
  upper_idx <- min(length(lambda_grid), closest_idx + 1L)

  objective <- function(lambda) {
    (cor_at_lambda(lambda) - target_cor)^2
  }

  local_fit <- stats::optimize(
    f = objective,
    interval = c(lambda_grid[lower_idx], lambda_grid[upper_idx])
  )
  lambda_hat <- local_fit$minimum
  achieved_cor <- cor_at_lambda(lambda_hat)

  if (abs(achieved_cor - target_cor) > cor_tol && cor_unattainable == "error") {
    stop(
      sprintf(
        paste0(
          "Could not achieve correlation %.10f within tolerance %.3g. ",
          "Closest value found was %.10f."
        ),
        target_cor, cor_tol, achieved_cor
      )
    )
  }

  if (abs(achieved_cor - target_cor) > cor_tol && cor_unattainable == "nearest") {
    warning(
      sprintf(
        paste0(
          "Requested correlation %.10f was not achieved exactly. ",
          "Using nearest value %.10f instead."
        ),
        target_cor, achieved_cor
      ),
      call. = FALSE
    )
  }

  x_transformed_used <- .box_cox_increasing(x_shifted, lambda_hat)
  x_out <- rep(NA_real_, length(y))
  x_out[keep] <- x_transformed_used
  if (!is.null(names(y))) names(x_out) <- names(y)

  structure(
    list(
      x = x_out,
      requested_auc = target_auc,
      requested_cor = target_cor,
      achieved_auc = .empirical_auc_from_scores(y_used, x_transformed_used),
      achieved_cor = achieved_cor,
      lambda = lambda_hat,
      cor_range = cor_range,
      base_vector = base_vector
    ),
    class = "aucmat_auc_cor_vector"
  )
}

#' Simulate repeated X vectors for a fixed y and record AUC and correlation
#'
#' For a fixed binary outcome vector `y`, repeatedly simulate continuous
#' predictors `X` from two normal distributions with equal variance and a mean
#' separation chosen so that the theoretical AUC equals `target_auc`. For each
#' simulation, the function calculates the empirical AUC and the correlation
#' between `y` and `X`.
#'
#' @param y Binary outcome vector. Accepted inputs are logical, numeric 0/1,
#'   factor, or character with exactly two non-missing classes.
#' @param target_auc One target AUC or a numeric vector of target AUC values in
#'   `[0, 1]`.
#' @param n_sim Number of simulations. Default 10000.
#' @param positive Optional positive class label when `y` is a factor or
#'   character vector. By default, the second factor level is treated as the
#'   positive class.
#' @param seed Optional random seed.
#' @param cor_method Correlation method passed to [stats::cor()]. One of
#'   `"pearson"`, `"spearman"`, or `"kendall"`. Default `"pearson"`.
#' @param keep_x Logical. If `TRUE`, also return the simulated `X` values as a
#'   matrix with one column per simulation. Default `FALSE`.
#'
#' @return A list with components:
#' * `results` — A data.frame with one row per simulation per requested
#'   target and columns `target_auc`, `sim`, `auc`, and `correlation`.
#' * `target_auc` — The requested target AUC value(s).
#' * `mean_shift` — A data.frame with columns `target_auc` and
#'   `mean_shift`. For `target_auc = 1`, `mean_shift` is reported as `Inf` because
#'   perfect separation is handled as a special case rather than a finite normal
#'   shift.
#' * `n_sim` — Number of simulations performed per target AUC.
#' * `x_matrix` — Optional matrix of simulated predictors if `keep_x = TRUE`
#'   and `length(target_auc) == 1`.
#' * `x_matrix_list` — Optional named list of simulated predictor matrices if
#'   `keep_x = TRUE` and multiple target AUC values are supplied.
#'
#' @examples
#' set.seed(1)
#' y <- rbinom(1000, size = 1, prob = 0.3)
#' sim <- simulate_auc_correlation(y, target_auc = c(0.8, 0.9, 1), n_sim = 100)
#' head(sim$results)
#' aggregate(auc ~ target_auc, data = sim$results, mean)
#'
#' @export
simulate_auc_correlation <- function(
    y,
    target_auc,
    n_sim = 10000,
    positive = NULL,
    seed = NULL,
    cor_method = c("pearson", "spearman", "kendall"),
    keep_x = FALSE
) {
  cor_method <- match.arg(cor_method)

  if (!is.numeric(target_auc) || any(is.na(target_auc)) ||
      any(target_auc < 0) || any(target_auc > 1)) {
    stop("target_auc must contain numeric value(s) in [0, 1].")
  }

  if (!is.numeric(n_sim) || length(n_sim) != 1L || is.na(n_sim) || n_sim < 1) {
    stop("n_sim must be a single positive integer.")
  }
  n_sim <- as.integer(n_sim)

  if (!is.logical(keep_x) || length(keep_x) != 1L) {
    stop("keep_x must be TRUE or FALSE.")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("seed must be NULL or a single non-missing numeric value.")
    }
    set.seed(seed)
  }

  keep <- !is.na(y)
  if (!any(keep)) stop("y contains only NA values.")

  y_used <- .normalize_binary_y(y[keep], positive = positive)
  y_num <- .binary_factor_to_numeric(y_used)
  pos <- y_num == 1
  neg <- !pos

  n_pos <- sum(pos)
  n_neg <- sum(neg)
  if (n_pos == 0L || n_neg == 0L) {
    stop("y must contain at least one positive and one negative observation.")
  }

  sim_list <- lapply(target_auc, function(one_target) {
    .simulate_auc_correlation_single(
      y = y,
      y_used = y_used,
      y_num = y_num,
      keep = keep,
      pos = pos,
      neg = neg,
      n_pos = n_pos,
      n_neg = n_neg,
      target_auc = one_target,
      n_sim = n_sim,
      cor_method = cor_method,
      keep_x = keep_x
    )
  })

  results <- do.call(rbind, lapply(sim_list, function(item) item$results))
  rownames(results) <- NULL

  out <- list(
    results = results,
    target_auc = target_auc,
    mean_shift = do.call(rbind, lapply(sim_list, function(item) {
      data.frame(target_auc = item$target_auc, mean_shift = item$mean_shift)
    })),
    n_sim = n_sim
  )

  if (keep_x) {
    if (length(sim_list) == 1L) {
      out$x_matrix <- sim_list[[1L]]$x_matrix
    } else {
      x_matrix_list <- lapply(sim_list, function(item) item$x_matrix)
      names(x_matrix_list) <- format(target_auc, trim = TRUE, scientific = FALSE)
      out$x_matrix_list <- x_matrix_list
    }
  }

  out
}

#' @noRd
.simulate_auc_correlation_single <- function(
    y,
    y_used,
    y_num,
    keep,
    pos,
    neg,
    n_pos,
    n_neg,
    target_auc,
    n_sim,
    cor_method,
    keep_x
) {
  mean_shift <- if (target_auc < 1) sqrt(2) * stats::qnorm(target_auc) else Inf
  auc_values <- numeric(n_sim)
  cor_values <- numeric(n_sim)
  x_matrix <- if (keep_x) matrix(NA_real_, nrow = length(y), ncol = n_sim) else NULL

  for (sim_idx in seq_len(n_sim)) {
    x_used <- numeric(length(y_used))

    if (target_auc < 1) {
      x_used[neg] <- stats::rnorm(n_neg, mean = 0, sd = 1)
      x_used[pos] <- stats::rnorm(n_pos, mean = mean_shift, sd = 1)
    } else {
      # Perfect separation case: simulate within-class variation, then shift the
      # positive class just enough to lie strictly above all negatives.
      x_used[neg] <- stats::rnorm(n_neg, mean = 0, sd = 1)
      x_pos <- stats::rnorm(n_pos, mean = 0, sd = 1)
      gap <- max(x_used[neg]) - min(x_pos) + 1e-8
      if (gap < 0) gap <- 1e-8
      x_used[pos] <- x_pos + gap
    }

    auc_values[sim_idx] <- .empirical_auc_from_scores(y_used, x_used)
    cor_values[sim_idx] <- stats::cor(y_num, x_used, method = cor_method)

    if (keep_x) {
      x_full <- rep(NA_real_, length(y))
      x_full[keep] <- x_used
      x_matrix[, sim_idx] <- x_full
    }
  }

  results <- data.frame(
    target_auc = rep(target_auc, n_sim),
    sim = seq_len(n_sim),
    auc = auc_values,
    correlation = cor_values
  )

  if (keep_x) {
    colnames(x_matrix) <- paste0("sim_", seq_len(n_sim))
  }

  list(
    results = results,
    target_auc = target_auc,
    mean_shift = mean_shift,
    x_matrix = x_matrix
  )
}

#' @noRd
.box_cox_increasing <- function(x, lambda) {
  if (any(x <= 0, na.rm = TRUE)) {
    stop("Box-Cox transform requires strictly positive x values.")
  }

  if (abs(lambda) < 1e-8) {
    return(log(x))
  }

  (x^lambda - 1) / lambda
}

#' @noRd
.build_rank_subset <- function(n, k, target_sum) {
  min_sum <- k * (k + 1L) / 2
  max_sum <- k * (2L * n - k + 1L) / 2

  if (target_sum < min_sum || target_sum > max_sum) {
    stop("target_sum is outside the attainable range for the requested subset size.")
  }

  ranks <- seq_len(k)
  remaining <- target_sum - min_sum

  for (i in k:1) {
    max_rank_i <- n - (k - i)
    max_increment <- max_rank_i - ranks[i]
    increment <- min(remaining, max_increment)
    ranks[i] <- ranks[i] + increment
    remaining <- remaining - increment
  }

  if (remaining != 0) stop("Failed to construct ranks for the requested AUC.")
  ranks
}

#' @noRd
.empirical_auc_from_scores <- function(y, x) {
  lev <- levels(y)
  pos <- y == lev[2L]
  n_pos <- sum(pos)
  n_neg <- sum(!pos)
  x_ranks <- rank(x, ties.method = "average")
  u_stat <- sum(x_ranks[pos]) - n_pos * (n_pos + 1L) / 2
  as.numeric(u_stat / (n_pos * n_neg))
}