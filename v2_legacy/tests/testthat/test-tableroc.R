# File: tests/testthat/test-tableroc.R

test_that("basic functionality: numeric columns only, structure, and values", {
  set.seed(1)
  X <- data.frame(a = rnorm(200), b = rnorm(200), c = letters[1:200])
  y <- factor(rbinom(200, 1, 0.5), labels = c("neg", "pos"))

  tb <- tableroc(X, y)

  expect_true(is.data.frame(tb))
  expect_equal(nrow(tb), 2L)  # only numeric columns a, b
  expect_true(all(tb$biomarker %in% c("a", "b")))
  expect_true(all(c("biomarker","auc","ci_low","ci_high","n_pos","n_neg","na_handling","n_na","n_imputed") %in% names(tb)))
  # Sample sizes used are positive and sum to n
  expect_true(all(tb$n_pos + tb$n_neg == 200))
})

test_that("ranking works and constant columns yield NA AUC", {
  set.seed(2)
  X <- data.frame(a = rnorm(100), b = 1)  # constant column b
  y <- factor(rbinom(100, 1, 0.6), labels = c("neg", "pos"))

  tb <- tableroc(X, y, rank = TRUE)
  expect_equal(nrow(tb), 2L)
  expect_true("rank" %in% names(tb))
  expect_equal(tb$rank, seq_len(nrow(tb)))
  expect_true(is.na(tb$auc[tb$biomarker == "b"]))
})

test_that("y types: logical and numeric 0/1 are supported", {
  set.seed(3)
  X <- data.frame(a = rnorm(120), b = rnorm(120))
  y_log <- sample(c(TRUE, FALSE), 120, TRUE)
  y_num <- as.numeric(y_log)

  tb1 <- tableroc(X, y_log)
  tb2 <- tableroc(X, y_num)

  expect_equal(nrow(tb1), 2L)
  expect_equal(nrow(tb2), 2L)
  expect_true(all(!is.na(tb1$auc)))
  expect_true(all(!is.na(tb2$auc)))
})

test_that("NA handling: na_impute='none' with na_rm TRUE vs FALSE affects per-column n", {
  set.seed(10)
  X <- data.frame(a = rnorm(150), b = rnorm(150))
  # Inject NAs in 'a' only
  idx_na <- sample(1:150, 15)
  X$a[idx_na] <- NA
  y <- factor(rbinom(150, 1, 0.5), labels = c("neg","pos"))

  # Global row drop (drop any row with NA in numeric X)
  tb_drop <- tableroc(X, y, na_impute = "none", na_rm = TRUE)
  expected_n_drop <- sum(complete.cases(X))  # complete cases across a and b
  expect_true(all(tb_drop$n_pos + tb_drop$n_neg == expected_n_drop))

  # Per-column omission (no global drop)
  tb_per <- tableroc(X, y, na_impute = "none", na_rm = FALSE)
  n_used_a <- sum(!is.na(X$a))
  n_used_b <- sum(!is.na(X$b))
  expect_equal((tb_per$n_pos + tb_per$n_neg)[tb_per$biomarker == "a"], n_used_a)
  expect_equal((tb_per$n_pos + tb_per$n_neg)[tb_per$biomarker == "b"], n_used_b)
})

test_that("imputation methods produce results and diagnostics", {
  set.seed(42)
  X <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
  # Add NAs in a (10) and c (5)
  na_a <- sample(1:100, 10); na_c <- sample(setdiff(1:100, na_a), 5)
  X$a[na_a] <- NA; X$c[na_c] <- NA
  y <- factor(rbinom(100, 1, 0.5), labels = c("neg", "pos"))

  tb_mean   <- tableroc(X, y, na_impute = "mean")
  tb_median <- tableroc(X, y, na_impute = "median")
  tb_const  <- tableroc(X, y, na_impute = "constant", na_constant = -999)
  tb_quant  <- tableroc(X, y, na_impute = "quantile", na_quantile = 0.2)
  tb_half   <- tableroc(X, y, na_impute = "halfmin")
  tb_bycls  <- tableroc(X, y, na_impute = "median_by_class")

  for (tb in list(tb_mean, tb_median, tb_const, tb_quant, tb_half, tb_bycls)) {
    expect_true(all(tb$biomarker %in% c("a","b","c")))
    expect_true(all(!is.na(tb$auc)))
    expect_true(all(c("n_na","n_imputed","na_handling") %in% names(tb)))
  }
  # Check imputed counts track NA counts for 'mean'
  imp_a <- tb_mean$n_imputed[tb_mean$biomarker == "a"]
  imp_b <- tb_mean$n_imputed[tb_mean$biomarker == "b"]
  imp_c <- tb_mean$n_imputed[tb_mean$biomarker == "c"]
  expect_equal(imp_a, length(na_a))
  expect_equal(imp_b, 0L)
  expect_equal(imp_c, length(na_c))
})

test_that("KNN imputation works and enforces max_n error", {
  set.seed(123)
  X <- data.frame(a = rnorm(120), b = rnorm(120))
  X$a[sample(1:120, 7)] <- NA
  y <- factor(rbinom(120, 1, 0.5), labels = c("neg", "pos"))
  tb_knn <- tableroc(X, y, na_impute = "knn", knn_k = 3)
  expect_true(all(!is.na(tb_knn$auc)))

  # Create a dataset with rows > knn_max_n to trigger early stop (no big allocation)
  n_big <- 2001
  Xbig <- as.data.frame(matrix(rnorm(n_big * 3L), ncol = 3))
  Xbig[1:10, 1] <- NA
  ybig <- factor(sample(c("neg","pos"), n_big, TRUE))
  expect_error(tableroc(Xbig, y = ybig, na_impute = "knn", knn_max_n = 1000))
})

test_that("mice single imputation runs if installed", {
  skip_if_not_installed("mice")
  set.seed(2020)
  X <- data.frame(a = rnorm(80), b = rnorm(80), c = rnorm(80))
  X$a[sample(1:80, 9)] <- NA
  X$b[sample(1:80, 5)] <- NA
  y <- factor(rbinom(80, 1, 0.5), labels = c("neg", "pos"))
  tb_mice <- tableroc(X, y, na_impute = "mice", mice_m = 1, mice_maxit = 2)
  expect_true(all(!is.na(tb_mice$auc)))
})

test_that("missForest imputation runs if installed and data are usable", {
  skip_if_not_installed("missForest")
  set.seed(3030)
  X <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
  # Insert some NAs but keep columns non-constant and not all-NA
  X$a[sample(1:100, 8)] <- NA
  X$b[sample(1:100, 6)] <- NA
  y <- factor(rbinom(100, 1, 0.5), labels = c("neg", "pos"))

  # Ensure at least two non-constant columns remain after removing NAs
  stopifnot(any(apply(X, 2, function(v) !all(is.na(v)) && length(unique(v[!is.na(v)])) > 1)))

  tb_mf <- tableroc(X, y, na_impute = "missForest", missforest_maxiter = 1, missforest_ntree = 50)
  expect_true(all(!is.na(tb_mf$auc)))
})
