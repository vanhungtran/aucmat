# Tests for auc_stability()

test_that("auc_stability() returns correct structure", {
  set.seed(42)
  X <- matrix(rnorm(100 * 20), nrow = 100)
  colnames(X) <- paste0("bm", 1:20)
  y <- rep(c(0, 1), each = 50)

  stab <- auc_stability(X, y, times = 50, top_k = c(5, 10), seed = 123)

  expect_s3_class(stab, "aucmat_stability")
  expect_true(all(c("rank_summary", "top_k_probs", "coselection", "settings")
    %in% names(stab)))
  expect_equal(nrow(stab$rank_summary), 20)
  expect_true("rank_median" %in% names(stab$rank_summary))
  expect_true("top1_freq" %in% names(stab$rank_summary))
})

test_that("auc_stability() is reproducible with a fixed seed", {
  X <- matrix(rnorm(60 * 8), nrow = 60)
  colnames(X) <- paste0("g", 1:8)
  y <- rep(c(0, 1), each = 30)

  s1 <- auc_stability(X, y, times = 30, seed = 42)
  s2 <- auc_stability(X, y, times = 30, seed = 42)

  expect_equal(s1$rank_summary$rank_median, s2$rank_summary$rank_median)
  expect_equal(s1$rank_summary$top1_freq, s2$rank_summary$top1_freq)
})

test_that("auc_stability() top_k_probs has correct structure", {
  set.seed(42)
  X <- matrix(rnorm(80 * 12), nrow = 80)
  colnames(X) <- paste0("m", 1:12)
  y <- rep(c(0, 1), each = 40)

  stab <- auc_stability(X, y, times = 40, top_k = c(3, 5), seed = 1)

  expect_true(all(c("biomarker", "k", "prob") %in% names(stab$top_k_probs)))
  expect_equal(sort(unique(stab$top_k_probs$k)), c(3, 5))
  # probabilities should be in [0, 1]
  expect_true(all(stab$top_k_probs$prob >= 0 & stab$top_k_probs$prob <= 1))
})

test_that("auc_stability() respects max_pairs_for_coselection", {
  X <- matrix(rnorm(60 * 30), nrow = 60)
  colnames(X) <- paste0("f", 1:30)
  y <- rep(c(0, 1), each = 30)

  stab <- auc_stability(X, y, times = 20, top_k = c(5),
    max_pairs_for_coselection = 5, seed = 1)

  expect_equal(length(stab$coselection$biomarkers), 5)
  expect_equal(dim(stab$coselection$matrix), c(5, 5))
})
