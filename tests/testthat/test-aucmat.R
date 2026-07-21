# Tests for aucmat() main screening function

test_that("aucmat() accepts a numeric matrix and binary outcome", {
  set.seed(42)
  X <- matrix(rnorm(100 * 10), nrow = 100)
  colnames(X) <- paste0("bm", 1:10)
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")
  expect_s3_class(fit, "aucmat_screen")
  expect_equal(nrow(fit$results), 10)
  expect_true(all(c("biomarker", "auc_raw", "auc_strength",
    "effect_direction", "n_used", "p_value", "rank") %in% names(fit$results)))
})

test_that("aucmat() handles a data.frame input", {
  set.seed(42)
  X <- as.data.frame(matrix(rnorm(50 * 5), nrow = 50))
  colnames(X) <- paste0("g", 1:5)
  y <- c(rep("ctrl", 25), rep("case", 25))

  fit <- aucmat(X, y, positive = "case", ci = "none")
  expect_s3_class(fit, "aucmat_screen")
  expect_equal(nrow(fit$results), 5)
})

test_that("aucmat() rejects bad inputs", {
  X <- matrix(rnorm(20), nrow = 10)
  colnames(X) <- c("a", "b")

  expect_error(aucmat(X, c(0, 1)), "must match nrow")
  expect_error(aucmat(X, rep(1, 10)), "exactly 2 unique")
  # column without a name (empty string or NA)
  X_bad <- matrix(rnorm(20), nrow = 10)
  expect_error(aucmat(X_bad, rep(c(0, 1), each = 5)),
               "unique and non-empty")
})

test_that("aucmat() reports direction-preserving AUC", {
  # biomarker higher in positives
  y <- rep(c(0, 1), each = 50)
  X <- cbind(
    higher = c(rnorm(50, 0, 1), rnorm(50, 2, 1)),
    lower  = c(rnorm(50, 2, 1), rnorm(50, 0, 1))
  )
  colnames(X) <- c("higher", "lower")

  fit <- aucmat(X, y, ci = "none")
  expect_gt(fit$results$auc_raw[fit$results$biomarker == "higher"], 0.5)
  expect_lt(fit$results$auc_raw[fit$results$biomarker == "lower"], 0.5)
  expect_equal(fit$results$effect_direction[fit$results$biomarker == "higher"],
               "higher_in_positive")
  expect_equal(fit$results$effect_direction[fit$results$biomarker == "lower"],
               "lower_in_positive")
})

test_that("aucmat() handles missing values with featurewise action", {
  X <- matrix(rnorm(100 * 5), nrow = 100)
  X[1:5, 1] <- NA  # 5 missing in bm1
  colnames(X) <- paste0("bm", 1:5)
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none", na_action = "featurewise")
  expect_equal(fit$results$n_missing[fit$results$biomarker == "bm1"], 5)
  expect_equal(fit$results$n_used[fit$results$biomarker == "bm1"], 95)
  expect_equal(fit$results$n_used[fit$results$biomarker == "bm2"], 100)
})

test_that("aucmat() handles complete-case action", {
  X <- matrix(rnorm(30 * 4), nrow = 30)
  X[1, 1] <- NA
  colnames(X) <- paste0("bm", 1:4)
  y <- rep(c(0, 1), length.out = 30)

  fit <- aucmat(X, y, ci = "none", na_action = "complete")
  # all biomarkers should have the same n_used (row 1 excluded)
  expect_equal(length(unique(fit$results$n_used)), 1)
})

test_that("aucmat() handles constant and near-constant biomarkers", {
  X <- cbind(
    ok     = rnorm(50),
    const  = rep(1, 50),
    nearly = c(rep(1, 49), 2)
  )
  colnames(X) <- c("ok", "const", "nearly")
  y <- rep(c(0, 1), each = 25)

  fit <- aucmat(X, y, ci = "none")
  expect_equal(fit$results$status[fit$results$biomarker == "const"], "constant")
  expect_equal(fit$results$status[fit$results$biomarker == "ok"], "ok")
})

test_that("aucmat() results are sorted by auc_strength descending", {
  set.seed(123)
  X <- matrix(rnorm(100 * 20), nrow = 100)
  colnames(X) <- paste0("f", 1:20)
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")
  strengths <- fit$results$auc_strength
  expect_true(all(diff(strengths) <= 0) || all(is.na(diff(strengths))))
})

test_that("S3 methods work on aucmat_screen", {
  set.seed(42)
  X <- matrix(rnorm(100 * 10), nrow = 100)
  colnames(X) <- paste0("bm", 1:10)
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")

  expect_output(print(fit), "aucmat_screen")
  expect_output(summary(fit), "screening summary")
  expect_s3_class(as.data.frame(fit), "data.frame")
  expect_equal(nrow(as.data.frame(fit)), 10)

  sub <- subset(fit, q_value < 1)
  expect_s3_class(sub, "data.frame")
})

test_that("aucmat() with DeLong CI matches pROC approximately", {
  skip_if_not_installed("pROC")
  set.seed(42)
  y <- rep(c(0, 1), each = 100)
  x <- c(rnorm(100, 0, 1), rnorm(100, 1.5, 1))
  X <- matrix(x, ncol = 1)
  colnames(X) <- "bm1"

  fit <- aucmat(X, y, ci = "delong")

  roc_obj <- pROC::roc(y, x, quiet = TRUE)
  auc_pROC <- as.numeric(pROC::auc(roc_obj))
  ci_pROC  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))

  expect_equal(fit$results$auc_raw[1], auc_pROC, tolerance = 1e-10)
  expect_equal(fit$results$std_error[1],
    (ci_pROC[3] - ci_pROC[1]) / (2 * qnorm(0.975)),
    tolerance = 1e-6)
})
