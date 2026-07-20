test_that("generate_auc_vector reaches an attainable AUC exactly", {
  y <- c(rep(0, 20), rep(1, 20))

  out <- generate_auc_vector(
    y = y,
    target_auc = 0.75,
    shuffle_within_class = FALSE
  )

  expect_length(out$x, length(y))
  expect_equal(out$achieved_auc, 0.75)
  expect_equal(out$auc_step, 1 / (20 * 20))
})

test_that("generate_auc_vector falls back to the nearest attainable AUC", {
  y <- c(0, 0, 0, 1, 1)

  expect_warning(
    out <- generate_auc_vector(y = y, target_auc = 0.7, shuffle_within_class = FALSE),
    "not exactly attainable"
  )

  expect_equal(out$achieved_auc, 4 / 6)
})

test_that("generate_auc_vector can preserve NA positions", {
  y <- c(0, 1, NA, 0, 1, NA)

  out <- generate_auc_vector(y = y, target_auc = 0.75, shuffle_within_class = FALSE)

  expect_true(all(is.na(out$x[is.na(y)])))
  expect_equal(sum(!is.na(out$x)), sum(!is.na(y)))
})

test_that("simulate_auc_correlation returns per-simulation AUC and correlation", {
  y <- c(rep(0, 30), rep(1, 20))

  sim <- simulate_auc_correlation(
    y = y,
    target_auc = 0.8,
    n_sim = 50,
    seed = 123
  )

  expect_true(is.list(sim))
  expect_true("results" %in% names(sim))
  expect_equal(nrow(sim$results), 50)
  expect_true(all(c("target_auc", "sim", "auc", "correlation") %in% names(sim$results)))
  expect_true(all(sim$results$auc >= 0 & sim$results$auc <= 1))
  expect_true(all(is.finite(sim$results$correlation)))
})

test_that("simulate_auc_correlation can optionally keep simulated X values", {
  y <- c(0, 1, NA, 0, 1)

  sim <- simulate_auc_correlation(
    y = y,
    target_auc = 0.75,
    n_sim = 5,
    seed = 1,
    keep_x = TRUE
  )

  expect_true("x_matrix" %in% names(sim))
  expect_equal(dim(sim$x_matrix), c(length(y), 5))
  expect_true(all(is.na(sim$x_matrix[3, ])))
})

test_that("simulate_auc_correlation accepts a vector of target AUC values", {
  y <- c(rep(0, 30), rep(1, 20))

  sim <- simulate_auc_correlation(
    y = y,
    target_auc = c(0.7, 0.8, 1),
    n_sim = 10,
    seed = 123
  )

  expect_equal(nrow(sim$results), 30)
  expect_equal(sort(unique(sim$results$target_auc)), c(0.7, 0.8, 1.0))
  expect_true(all(c("target_auc", "mean_shift") %in% names(sim$mean_shift)))
  expect_equal(sim$mean_shift$mean_shift[sim$mean_shift$target_auc == 1], Inf)
  expect_true(all(sim$results$auc[sim$results$target_auc == 1] == 1))
})

test_that("generate_auc_cor_vector preserves AUC while matching a feasible correlation", {
  y <- c(rep(0, 20), rep(1, 20))
  base <- generate_auc_vector(y, target_auc = 0.75, shuffle_within_class = FALSE)
  target_cor <- cor(y, base$x)

  out <- generate_auc_cor_vector(
    y = y,
    target_auc = 0.75,
    target_cor = target_cor,
    shuffle_within_class = FALSE,
    cor_tol = 1e-8
  )

  expect_equal(out$achieved_auc, 0.75)
  expect_equal(out$achieved_cor, target_cor, tolerance = 1e-6)
})

test_that("generate_auc_cor_vector warns and returns nearest match when correlation is infeasible", {
  y <- c(rep(0, 30), rep(1, 20))

  expect_warning(
    out <- generate_auc_cor_vector(
      y = y,
      target_auc = 0.6,
      target_cor = 0.95,
      shuffle_within_class = FALSE,
      cor_unattainable = "nearest"
    ),
    "outside the approximate achievable range|not achieved exactly"
  )

  expect_equal(out$achieved_auc, 0.6)
  expect_true(out$achieved_cor < 0.95)
})