# Tests for compare_auc()

test_that("compare_auc() computes paired comparisons", {
  set.seed(42)
  X <- cbind(
    bm1 = c(rnorm(50, 0, 1), rnorm(50, 2, 1)),
    bm2 = c(rnorm(50, 0, 1), rnorm(50, 1, 1)),
    bm3 = c(rnorm(50, 0, 1), rnorm(50, 0.5, 1))
  )
  colnames(X) <- c("bm1", "bm2", "bm3")
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")
  cmp <- compare_auc(fit, X, y, biomarkers = c("bm1", "bm2", "bm3"))

  expect_s3_class(cmp, "aucmat_compare")
  expect_equal(nrow(cmp), 3)  # 3 choose 2 = 3 pairs
  expect_true(all(c("biomarker_a", "biomarker_b", "auc_diff",
    "std_error", "p_value") %in% names(cmp)))
})

test_that("compare_auc() with reference compares all others to reference", {
  set.seed(42)
  X <- cbind(
    bm1 = c(rnorm(50, 0, 1), rnorm(50, 2, 1)),
    bm2 = c(rnorm(50, 0, 1), rnorm(50, 1, 1)),
    bm3 = c(rnorm(50, 0, 1), rnorm(50, 0.5, 1))
  )
  colnames(X) <- c("bm1", "bm2", "bm3")
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")
  cmp <- compare_auc(fit, X, y, reference = "bm1")

  expect_equal(unique(cmp$biomarker_a), "bm1")
  expect_equal(nrow(cmp), 2)
})

test_that("compare_auc() enforces max_pairs safety limit", {
  set.seed(42)
  X <- matrix(rnorm(50 * 15), nrow = 50)
  colnames(X) <- paste0("f", 1:15)
  y <- rep(c(0, 1), each = 25)

  fit <- aucmat(X, y, ci = "none")
  expect_error(
    compare_auc(fit, X, y, top_n = 15, max_pairs = 5),
    "max_pairs"
  )
})

test_that("compare_auc() handles adjust parameter", {
  set.seed(42)
  X <- cbind(
    bm1 = c(rnorm(50, 0, 1), rnorm(50, 2, 1)),
    bm2 = c(rnorm(50, 0, 1), rnorm(50, 1, 1)),
    bm3 = c(rnorm(50, 0, 1), rnorm(50, 0.5, 1))
  )
  colnames(X) <- c("bm1", "bm2", "bm3")
  y <- rep(c(0, 1), each = 50)

  fit <- aucmat(X, y, ci = "none")
  cmp_bh <- compare_auc(fit, X, y, biomarkers = c("bm1", "bm2", "bm3"),
                         adjust = "BH")

  expect_true(all(c("q_value") %in% names(cmp_bh)))
  expect_true(any(!is.na(cmp_bh$q_value)))
})
