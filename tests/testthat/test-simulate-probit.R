# Tests for generate_data_probit() and its .calibrate_rho() helper

test_that(".calibrate_rho does not overflow near prevalence = 0.5", {
  # Regression test: np * nn as 32-bit integers overflows
  # .Machine$integer.max (~2.1e9) once both class counts exceed ~46341 at
  # the default n_cal = 100000, which happens for any prevalence in the
  # broad middle of (0, 1). This silently produced NA and crashed the
  # bisection root-finder.
  rho <- .calibrate_rho(target_auc = 0.75, prevalence = 0.5, n_cal = 100000)
  expect_true(is.finite(rho))
  expect_true(rho > 0 && rho < 1)
})

test_that(".calibrate_rho is finite across a range of prevalences", {
  for (prev in c(0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9)) {
    rho <- .calibrate_rho(target_auc = 0.8, prevalence = prev, n_cal = 100000)
    expect_true(is.finite(rho), info = paste("prevalence =", prev))
  }
})

test_that("generate_data_probit works at prevalence = 0.5", {
  set.seed(1)
  sim <- generate_data_probit(
    n = 300, target_aucs = c(0.8, 0.7),
    corr_matrix = matrix(c(1, 0.3, 0.3, 1), 2, 2),
    prevalence = 0.5
  )
  expect_true(all(is.finite(sim$achieved_aucs)))
})
