# Tests for deprecated functions

test_that("tableroc() is deprecated but still works", {
  set.seed(123)
  X <- data.frame(M1 = c(rnorm(95), rep(NA, 5)), M2 = rnorm(100))
  y <- factor(rbinom(100, 1, 0.5), labels = c("neg", "pos"))

  expect_warning(
    result <- tableroc(X, y, na_impute = "median"),
    "deprecated"
  )
  expect_true("auc_raw" %in% names(result$results) ||
              "auc" %in% names(as.data.frame(result)))
})

test_that("plot_roc_with_combos() is deprecated but still works", {
  skip_if_not_installed("pROC")
  set.seed(42)
  data <- data.frame(
    M1 = rnorm(200, 0, 1),
    M2 = rnorm(200, 5, 2),
    Outcome = sample(c(0, 1), 200, replace = TRUE)
  )

  expect_warning(
    res <- plot_roc_with_combos(
      data = data, outcome = "Outcome",
      predictors = c("M1", "M2"),
      combo_sizes = 1, add_ci = FALSE
    ),
    "deprecated"
  )
  expect_true(is.list(res))
  expect_true("plot" %in% names(res))
})
