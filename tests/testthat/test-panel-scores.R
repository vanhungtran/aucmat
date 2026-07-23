# Tests for fit_auc_panel(), predict.aucmat_panel(), plot_roc_panel(),
# and compare_panel_auc()

.panel_sim <- function(seed = 1, n = 240) {
  set.seed(seed)
  sim <- simulate_auc_matrix(
    n = n, prevalence = 0.3,
    target_aucs = c(0.85, 0.75, 0.65, 0.60),
    correlation = 0.3, structure = "exchangeable", seed = seed
  )
  list(X = as.matrix(sim$data[, 1:4]), y = sim$data$truth)
}

# ---- Fitting: all methods ---------------------------------------------

test_that("fit_auc_panel() fits all six methods and returns sane objects", {
  skip_if_not_installed("glmnet")
  d <- .panel_sim()
  for (m in c("ridge", "lasso", "elasticnet", "logistic", "su_liu", "unweighted")) {
    panel <- fit_auc_panel(d$X, d$y, method = m, n_folds = 3, seed = 1)
    expect_s3_class(panel, "aucmat_panel")
    expect_equal(panel$settings$method, m)
    expect_true(is.finite(panel$auc_train))
    expect_true(panel$auc_train >= 0.5 && panel$auc_train <= 1)
    expect_true(is.finite(panel$auc_cv))
    expect_true(panel$auc_cv >= 0.3 && panel$auc_cv <= 1)
    expect_named(panel$coefficients, colnames(d$X))
    expect_equal(nrow(panel$predictions), nrow(d$X))
    expect_false(anyNA(panel$predictions$score_cv))
  }
})

test_that("fit_auc_panel() requires at least 2 biomarkers", {
  d <- .panel_sim()
  expect_error(fit_auc_panel(d$X[, 1, drop = FALSE], d$y), "at least 2 biomarkers")
})

test_that("fit_auc_panel() rejects non-numeric X", {
  d <- .panel_sim()
  Xc <- d$X
  storage.mode(Xc) <- "character"
  expect_error(fit_auc_panel(Xc, d$y), "numeric")
})

test_that("fit_auc_panel() n_folds = 0 skips CV cleanly", {
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 0)
  expect_true(is.na(panel$auc_cv))
  expect_true(all(is.na(panel$predictions$score_cv)))
  expect_false(is.na(panel$auc_train))
})

# ---- Bug fix regression: elastic-net default alpha ---------------------

test_that("method = 'elasticnet' defaults to alpha = 0.5, not ridge's 0", {
  skip_if_not_installed("glmnet")
  d <- .panel_sim()
  panel_default <- fit_auc_panel(d$X, d$y, method = "elasticnet", n_folds = 0)
  expect_equal(panel_default$settings$alpha, 0.5)

  panel_explicit <- fit_auc_panel(d$X, d$y, method = "elasticnet", alpha = 0.2,
    n_folds = 0)
  expect_equal(panel_explicit$settings$alpha, 0.2)

  panel_ridge <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 0)
  expect_equal(panel_ridge$settings$alpha, 0)
})

# ---- Bug fix regression: honest nested CV (no leakage) ------------------

test_that("fit_auc_panel() honest CV AUC stays near chance on pure noise", {
  skip_if_not_installed("glmnet")
  set.seed(99)
  n <- 160; p <- 8
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("noise", seq_len(p))))
  y <- rep(c(0, 1), each = n / 2)  # independent of X by construction

  for (m in c("ridge", "lasso", "elasticnet", "su_liu", "unweighted")) {
    panel <- fit_auc_panel(X, y, method = m, n_folds = 5, seed = 7)
    # A leaky procedure that tunes (penalty selection, centering/scaling)
    # on data overlapping the held-out fold can fit noise and inflate this
    # well above chance; an honestly nested procedure should stay close to
    # 0.5 give or take sampling noise.
    expect_true(panel$auc_cv > 0.25 && panel$auc_cv < 0.75,
      info = paste0("method=", m, " auc_cv=", panel$auc_cv))
  }
})

test_that("fit_auc_panel() CV AUC does not wildly exceed train AUC (no gross leakage)", {
  skip_if_not_installed("glmnet")
  d <- .panel_sim()
  for (m in c("ridge", "lasso", "elasticnet", "logistic")) {
    panel <- fit_auc_panel(d$X, d$y, method = m, n_folds = 5, seed = 3)
    # Honest out-of-fold AUC should not exceed in-sample AUC by a
    # meaningful margin -- that pattern is a leakage signature.
    expect_true(panel$auc_cv <= panel$auc_train + 0.02,
      info = paste0("method=", m, " train=", panel$auc_train, " cv=", panel$auc_cv))
  }
})

# ---- su_liu ---------------------------------------------------------------

test_that("fit_auc_panel(method = 'su_liu') validates shrinkage range", {
  d <- .panel_sim()
  expect_error(fit_auc_panel(d$X, d$y, method = "su_liu", shrinkage = -0.1),
    "shrinkage")
  expect_error(fit_auc_panel(d$X, d$y, method = "su_liu", shrinkage = 1.5),
    "shrinkage")
})

test_that("fit_auc_panel(method = 'su_liu') errors on near-singular covariance at shrinkage = 0", {
  d <- .panel_sim()
  Xc <- cbind(d$X, X1_dup = d$X[, "X1"] + rnorm(nrow(d$X), sd = 1e-9))
  expect_error(
    fit_auc_panel(Xc, d$y, method = "su_liu", shrinkage = 0, n_folds = 0),
    "singular"
  )
  # Shrinkage rescues it
  panel <- fit_auc_panel(Xc, d$y, method = "su_liu", shrinkage = 0.3, n_folds = 0)
  expect_true(is.finite(panel$auc_train))
})

test_that("fit_auc_panel(method = 'su_liu') combined score direction is higher-in-positive", {
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "su_liu", n_folds = 0)
  pos <- d$y == 1  # simulate_auc_matrix() encodes truth as integer 0/1
  expect_true(mean(panel$predictions$score_train[pos]) >
              mean(panel$predictions$score_train[!pos]))
})

# ---- unweighted -------------------------------------------------------

test_that("fit_auc_panel(method = 'unweighted') gives equal-magnitude, direction-aligned weights", {
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "unweighted", n_folds = 0)
  expect_equal(unname(abs(panel$coefficients)), rep(1 / 4, 4))
  # All target AUCs were simulated "higher in positive", so all signs
  # should agree
  expect_true(all(panel$coefficients > 0))
})

# ---- Reproducibility & RNG safety --------------------------------------

test_that("fit_auc_panel() is fully reproducible given the same seed", {
  skip_if_not_installed("glmnet")
  d <- .panel_sim()
  p1 <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 4, seed = 11)
  p2 <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 4, seed = 11)
  expect_identical(p1$coefficients, p2$coefficients)
  expect_identical(p1$predictions$score_cv, p2$predictions$score_cv)
})

test_that("fit_auc_panel() restores the caller's global RNG state", {
  skip_if_not_installed("glmnet")
  d <- .panel_sim()
  set.seed(555)
  before <- .Random.seed
  invisible(fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 4, seed = 321))
  expect_identical(.Random.seed, before)
})

# ---- print / plot / predict --------------------------------------------

test_that("print.aucmat_panel() and plot.aucmat_panel() run without error", {
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "su_liu", n_folds = 3, seed = 1)
  expect_output(print(panel), "aucmat_panel")
  p <- plot(panel)
  expect_s3_class(p, "ggplot")
})

test_that("predict.aucmat_panel() reproduces score_train on the training data", {
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 0)
  pred <- predict(panel, d$X)
  expect_equal(as.vector(pred), panel$predictions$score_train, tolerance = 1e-8)
})

test_that("predict.aucmat_panel() errors when newdata is missing panel columns", {
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 0)
  expect_error(predict(panel, d$X[, 1:2, drop = FALSE]), "missing columns")
})

# ---- plot_roc_panel() --------------------------------------------------

test_that("plot_roc_panel() returns a ggplot object with default and explicit biomarkers", {
  skip_if_not_installed("pROC")
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 3, seed = 1)

  p1 <- plot_roc_panel(panel, d$X, d$y)
  expect_s3_class(p1, "ggplot")

  p2 <- plot_roc_panel(panel, d$X, d$y, biomarkers = c("X1", "X2"))
  expect_s3_class(p2, "ggplot")
})

test_that("plot_roc_panel() falls back to train score with a warning when CV is unavailable", {
  skip_if_not_installed("pROC")
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 0)
  expect_warning(plot_roc_panel(panel, d$X, d$y), "n_folds < 2")
})

test_that("plot_roc_panel() errors when X is missing panel columns", {
  skip_if_not_installed("pROC")
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 3, seed = 1)
  expect_error(plot_roc_panel(panel, d$X[, 1:2, drop = FALSE], d$y), "missing columns")
})

test_that("plot_roc_panel() add_ci = TRUE runs without error", {
  skip_if_not_installed("pROC")
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 3, seed = 1)
  p <- plot_roc_panel(panel, d$X, d$y, add_ci = TRUE, boot_n = 50, seed = 2)
  expect_s3_class(p, "ggplot")
})

# ---- compare_panel_auc() -----------------------------------------------

test_that("compare_panel_auc() compares the panel against each individual biomarker", {
  skip_if_not_installed("glmnet")
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 5, seed = 4)
  cmp <- compare_panel_auc(panel, d$X, d$y)

  expect_s3_class(cmp, "aucmat_compare")
  expect_equal(nrow(cmp), 4)
  expect_true(all(cmp$biomarker_a == "PANEL_SCORE"))
  expect_setequal(cmp$biomarker_b, colnames(d$X))
  expect_true(all(c("auc_diff", "p_value", "conf_low", "conf_high") %in% names(cmp)))
})

test_that("compare_panel_auc() detects that the panel beats a weak biomarker more clearly than a strong one", {
  skip_if_not_installed("glmnet")
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 5, seed = 4)
  cmp <- compare_panel_auc(panel, d$X, d$y)

  p_vs_strong <- cmp$p_value[cmp$biomarker_b == "X1"]  # target AUC 0.85, close to panel
  p_vs_weak   <- cmp$p_value[cmp$biomarker_b == "X4"]  # target AUC 0.60, far from panel
  expect_true(p_vs_weak < p_vs_strong)
})

test_that("compare_panel_auc() errors on use = 'cv' without CV, works with use = 'train'", {
  skip_if_not_installed("glmnet")
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 0)
  expect_error(compare_panel_auc(panel, d$X, d$y), "n_folds < 2")
  cmp <- compare_panel_auc(panel, d$X, d$y, use = "train")
  expect_s3_class(cmp, "aucmat_compare")
})

test_that("compare_panel_auc() respects a biomarkers subset", {
  skip_if_not_installed("glmnet")
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 5, seed = 4)
  cmp <- compare_panel_auc(panel, d$X, d$y, biomarkers = c("X1", "X3"))
  expect_equal(nrow(cmp), 2)
  expect_setequal(cmp$biomarker_b, c("X1", "X3"))
})

test_that("compare_panel_auc() errors when X is missing panel columns", {
  d <- .panel_sim()
  panel <- fit_auc_panel(d$X, d$y, method = "ridge", n_folds = 3, seed = 1)
  expect_error(compare_panel_auc(panel, d$X[, 1:2, drop = FALSE], d$y), "missing columns")
})
