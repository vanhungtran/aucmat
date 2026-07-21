# Tests for validate_simulation()

test_that("validate_simulation runs and reports calibration close to zero bias", {
  val <- validate_simulation(
    n = 300, prevalence = 0.3,
    target_aucs = c(0.8, 0.7), correlation = 0.3, structure = "exchangeable",
    times = 40, seed = 100
  )
  expect_s3_class(val, "aucmat_simulation_validation")
  expect_equal(val$n_requested, 40)
  expect_equal(val$n_ok, 40)
  expect_equal(val$n_failed, 0)

  # Monte Carlo-aware tolerance: bias should be small relative to MC SE
  expect_true(all(abs(val$auc_bias) < 5 * val$auc_mc_se + 0.01))
})

test_that("validate_simulation reports correlation calibration for p >= 2", {
  val <- validate_simulation(
    n = 300, prevalence = 0.3,
    target_aucs = c(0.8, 0.7, 0.6), correlation = 0.3,
    structure = "exchangeable", times = 30, seed = 1
  )
  expect_false(is.null(val$correlation))
  expect_length(val$correlation$bias, 3)  # 3 choose 2 = 3 pairs
})

test_that("validate_simulation omits correlation block for p = 1", {
  val <- validate_simulation(
    n = 200, prevalence = 0.3,
    target_aucs = 0.8, correlation = matrix(1, 1, 1), structure = "user",
    times = 20, seed = 1
  )
  expect_null(val$correlation)
})

test_that("seed binds to validate_simulation's own argument, not ...", {
  # R's argument matching binds a same-named formal before `...`, so `seed`
  # always controls the replicate stream regardless of call style.
  val <- validate_simulation(
    n = 100, prevalence = 0.3, target_aucs = 0.8,
    correlation = matrix(1, 1, 1), structure = "user",
    seed = 1, times = 5
  )
  expect_s3_class(val, "aucmat_simulation_validation")
})

test_that("validate_simulation reproducible replicate stream with same seed", {
  args <- list(n = 150, prevalence = 0.3, target_aucs = c(0.8, 0.7),
               correlation = 0.3, structure = "exchangeable",
               times = 15, seed = 55)
  val1 <- do.call(validate_simulation, args)
  val2 <- do.call(validate_simulation, args)
  expect_identical(val1$auc_bias, val2$auc_bias)
})

test_that("validate_simulation requires times >= 2", {
  expect_error(
    validate_simulation(n = 100, prevalence = 0.3, target_aucs = 0.8,
                        correlation = matrix(1, 1, 1), structure = "user",
                        times = 1),
    "times must be"
  )
})

test_that("print.aucmat_simulation_validation runs without error", {
  val <- validate_simulation(
    n = 150, prevalence = 0.3, target_aucs = c(0.8, 0.7),
    correlation = 0.3, structure = "exchangeable", times = 15, seed = 1
  )
  expect_output(print(val), "aucmat_simulation_validation")
})
