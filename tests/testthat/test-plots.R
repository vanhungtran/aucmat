# Tests for plot functions

test_that("plot_auc_rank returns a ggplot object", {
  set.seed(42)
  X <- matrix(rnorm(100 * 10), nrow = 100)
  colnames(X) <- paste0("bm", 1:10)
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")
  p <- plot_auc_rank(fit, n_label = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_auc_volcano returns a ggplot object", {
  set.seed(42)
  X <- matrix(rnorm(100 * 10), nrow = 100)
  colnames(X) <- paste0("bm", 1:10)
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")
  p <- plot_auc_volcano(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_auc_forest returns a ggplot object", {
  set.seed(42)
  X <- matrix(rnorm(100 * 10), nrow = 100)
  colnames(X) <- paste0("bm", 1:10)
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")
  p <- plot_auc_forest(fit, n = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_auc_stability returns a ggplot object", {
  set.seed(42)
  X <- matrix(rnorm(100 * 20), nrow = 100)
  colnames(X) <- paste0("bm", 1:20)
  y <- rep(c(0, 1), each = 50)

  stab <- auc_stability(X, y, times = 30, seed = 1)
  p <- plot_auc_stability(stab)
  expect_s3_class(p, "ggplot")
})

test_that("plot_roc_top works with explicit X and y", {
  skip_if_not_installed("pROC")
  set.seed(42)
  X <- cbind(
    bm1 = c(rnorm(50, 0, 1), rnorm(50, 2, 1)),
    bm2 = c(rnorm(50, 0, 1), rnorm(50, 1, 1))
  )
  colnames(X) <- c("bm1", "bm2")
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")
  p <- plot_roc_top(fit, X = X, y = y, biomarkers = c("bm1", "bm2"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_roc_top with retain_data uses stored data", {
  skip_if_not_installed("pROC")
  set.seed(42)
  X <- cbind(
    bm1 = c(rnorm(50, 0, 1), rnorm(50, 2, 1)),
    bm2 = c(rnorm(50, 0, 1), rnorm(50, 1, 1))
  )
  colnames(X) <- c("bm1", "bm2")
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none", retain_data = TRUE)
  p <- plot_roc_top(fit, biomarkers = c("bm1", "bm2"))
  expect_s3_class(p, "ggplot")
})
