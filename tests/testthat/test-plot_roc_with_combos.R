test_that("plot_roc_with_combos returns plot and AUCs", {
  skip_on_cran()
  skip_if_not_installed("pROC")
  skip_if_not_installed("mvtnorm")

  set.seed(42)
  sim <- generate_data_analytical(
    n = 400,
    prevalence = 0.3,
    target_aucs = c(0.85, 0.75, 0.65),
    corr_matrix = diag(3)
  )

  preds <- setdiff(names(sim$data), "truth")

  res <- plot_roc_with_combos(
    data = sim$data,
    outcome = "truth",
    predictors = preds,
    combo_sizes = 1:2,
    add_ci = FALSE,
    boot_n = 50,
    seed = 123
  )

  expect_true(is.list(res))
  expect_true("plot" %in% names(res))
  expect_true("auc" %in% names(res))
  expect_true(is.numeric(res$auc))
  expect_true(length(res$auc) >= 1)
})
