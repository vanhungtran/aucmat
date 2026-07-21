# Tests for simulate_auc_matrix() and its internal helpers

# ---- .build_correlation_matrix() ----

test_that(".build_correlation_matrix builds exchangeable correctly", {
  R <- .build_correlation_matrix(4, "exchangeable", correlation = 0.3)
  expect_equal(dim(R), c(4, 4))
  expect_equal(diag(R), rep(1, 4))
  expect_true(all(R[upper.tri(R)] == 0.3))
  expect_true(isSymmetric(R))
})

test_that(".build_correlation_matrix builds ar1 correctly", {
  R <- .build_correlation_matrix(4, "ar1", correlation = 0.5)
  expect_equal(R[1, 2], 0.5)
  expect_equal(R[1, 3], 0.25)
  expect_equal(R[1, 4], 0.125)
  expect_equal(R[2, 4], 0.25)
  expect_equal(diag(R), rep(1, 4))
})

test_that(".build_correlation_matrix builds block correctly", {
  R <- .build_correlation_matrix(4, "block", block_sizes = c(2, 2),
                                  rho_within = 0.6, rho_between = 0.1)
  expect_equal(R[1, 2], 0.6)   # within block 1
  expect_equal(R[3, 4], 0.6)   # within block 2
  expect_equal(R[1, 3], 0.1)   # between blocks
  expect_equal(R[2, 4], 0.1)   # between blocks
  expect_equal(diag(R), rep(1, 4))
})

test_that(".build_correlation_matrix passes through user matrix", {
  R_in <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  R_out <- .build_correlation_matrix(2, "user", correlation = R_in)
  expect_identical(R_out, R_in)
})

test_that(".build_correlation_matrix rejects malformed inputs", {
  expect_error(.build_correlation_matrix(3, "user", correlation = matrix(1, 2, 2)),
               "must be")
  expect_error(.build_correlation_matrix(2, "user",
                 correlation = matrix(c(1, 0.5, 0.9, 1), 2, 2)),
               "symmetric")
  expect_error(.build_correlation_matrix(3, "exchangeable", correlation = 1.5),
               "single value in")
  expect_error(.build_correlation_matrix(4, "block", block_sizes = c(2, 3),
                 rho_within = 0.5, rho_between = 0.1),
               "summing to")
})

# ---- .nearest_pd_correlation() ----

test_that(".nearest_pd_correlation leaves a valid PD matrix unchanged", {
  R <- .build_correlation_matrix(3, "exchangeable", correlation = 0.3)
  proj <- .nearest_pd_correlation(R)
  expect_false(proj$adjusted)
  expect_equal(proj$frobenius_adjustment, 0)
  expect_equal(proj$matrix, R)
})

test_that(".nearest_pd_correlation repairs an invalid matrix", {
  bad <- matrix(c(1, 0.9, 0.9, 0.9, 1, 0.9, 0.9, 0.9, 1), 3, 3)
  # Force clearly infeasible: high pairwise correlation but push toward
  # non-PD via a manually broken entry
  bad[1, 3] <- -0.95
  bad[3, 1] <- -0.95
  proj <- .nearest_pd_correlation(bad)
  eig_after <- eigen(proj$matrix, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig_after > -1e-6))
  expect_equal(diag(proj$matrix), rep(1, 3), tolerance = 1e-6)
})

# ---- simulate_auc_matrix(): outcome modes ----

test_that("simulate_auc_matrix generates outcome from n/prevalence", {
  set.seed(1)
  sim <- simulate_auc_matrix(
    n = 300, prevalence = 0.3,
    target_aucs = c(0.8, 0.7),
    correlation = 0.3, structure = "exchangeable"
  )
  expect_s3_class(sim, "aucmat_simulation")
  expect_equal(sim$n, 300)
  expect_equal(nrow(sim$data), 300)
  expect_equal(sum(sim$data$truth), round(300 * 0.3))  # exact mode
})

test_that("simulate_auc_matrix accepts a supplied outcome", {
  y <- rbinom(250, 1, 0.4)
  sim <- simulate_auc_matrix(
    y = y, target_aucs = c(0.75, 0.65),
    correlation = 0.2, structure = "exchangeable", seed = 1
  )
  expect_equal(sim$n, length(y))
  expect_equal(sim$prevalence, mean(y))
})

test_that("simulate_auc_matrix rejects both y and n/prevalence", {
  expect_error(
    simulate_auc_matrix(y = rbinom(50, 1, 0.3), n = 50, prevalence = 0.3,
                        target_aucs = 0.8, correlation = matrix(1, 1, 1)),
    "not both"
  )
})

test_that("simulate_auc_matrix rejects neither y nor n/prevalence", {
  expect_error(
    simulate_auc_matrix(target_aucs = 0.8, correlation = matrix(1, 1, 1)),
    "Supply either"
  )
})

test_that("bernoulli outcome_mode produces random (not exact) counts", {
  set.seed(1)
  sim <- simulate_auc_matrix(
    n = 200, prevalence = 0.3, outcome_mode = "bernoulli",
    target_aucs = c(0.8), correlation = matrix(1, 1, 1), structure = "user"
  )
  expect_equal(sim$n, 200)
  expect_true(sim$prevalence > 0 && sim$prevalence < 1)
})

# ---- simulate_auc_matrix(): target validation ----

test_that("simulate_auc_matrix rejects invalid target_aucs", {
  expect_error(
    simulate_auc_matrix(n = 100, prevalence = 0.3,
                        target_aucs = c(0.8, 1.5), correlation = 0.3,
                        structure = "exchangeable"),
    "\\(0, 1\\)"
  )
  expect_error(
    simulate_auc_matrix(n = 100, prevalence = 0.3,
                        target_aucs = numeric(0), correlation = 0.3,
                        structure = "exchangeable"),
    "in \\(0, 1\\)"
  )
})

# ---- simulate_auc_matrix(): achieved values are close to targets ----

test_that("achieved AUCs are reasonably close to targets", {
  set.seed(42)
  sim <- simulate_auc_matrix(
    n = 2000, prevalence = 0.3,
    target_aucs = c(0.9, 0.75),
    correlation = 0.2, structure = "exchangeable"
  )
  expect_equal(unname(sim$achieved_aucs), sim$target_aucs, tolerance = 0.03)
})

test_that("achieved correlation matches the requested structure", {
  set.seed(42)
  sim <- simulate_auc_matrix(
    n = 3000, prevalence = 0.3,
    target_aucs = c(0.85, 0.75, 0.65),
    correlation = 0.4, structure = "exchangeable"
  )
  off_diag <- sim$achieved_correlation[upper.tri(sim$achieved_correlation)]
  expect_true(all(abs(off_diag - 0.4) < 0.08))
})

# ---- feasibility ----
#
# Infeasibility now comes purely from the p x p correlation matrix R being
# non-PD (biomarker means are an unconstrained mean-shift, so near-perfect
# AUCs alone never cause infeasibility -- unlike the retired latent-probit
# construction, where AUC and correlation shared a single bounded
# parameter).  An exchangeable matrix with p biomarkers is PD iff
# correlation > -1/(p-1); for p = 4, correlation = -0.5 is well below the
# -1/3 boundary and is reliably non-PD.

test_that("feasibility = 'error' stops with a structured condition", {
  expect_error(
    simulate_auc_matrix(
      n = 200, prevalence = 0.3,
      target_aucs = c(0.8, 0.7, 0.6, 0.55),
      correlation = -0.5, structure = "exchangeable",
      feasibility = "error"
    ),
    class = "aucmat_infeasible_targets"
  )
})

test_that("feasibility = 'nearest' succeeds and reports adjustment", {
  sim <- simulate_auc_matrix(
    n = 200, prevalence = 0.3,
    target_aucs = c(0.8, 0.7, 0.6, 0.55),
    correlation = -0.5, structure = "exchangeable",
    feasibility = "nearest", seed = 1
  )
  expect_equal(sim$feasibility$status, "adjusted")
  expect_true(sim$feasibility$adjusted)
  expect_gt(sim$feasibility$frobenius_adjustment, 0)
})

test_that("feasible targets report status 'ok' with zero adjustment", {
  sim <- simulate_auc_matrix(
    n = 200, prevalence = 0.3,
    target_aucs = c(0.7, 0.6),
    correlation = 0.2, structure = "exchangeable",
    feasibility = "error", seed = 1
  )
  expect_equal(sim$feasibility$status, "ok")
  expect_false(sim$feasibility$adjusted)
})

# ---- reproducibility ----

test_that("same seed reproduces identical results", {
  args <- list(n = 200, prevalence = 0.3, target_aucs = c(0.8, 0.7),
               correlation = 0.3, structure = "exchangeable", seed = 42)
  sim1 <- do.call(simulate_auc_matrix, args)
  sim2 <- do.call(simulate_auc_matrix, args)
  expect_identical(sim1$achieved_aucs, sim2$achieved_aucs)
  expect_identical(sim1$data, sim2$data)
})

test_that("seeded call does not leak global RNG state", {
  set.seed(123)
  before <- .Random.seed
  simulate_auc_matrix(
    n = 100, prevalence = 0.3, target_aucs = c(0.8),
    correlation = matrix(1, 1, 1), structure = "user", seed = 999
  )
  after <- .Random.seed
  expect_identical(before, after)
})

# ---- S3 print method ----

test_that("print.aucmat_simulation runs without error", {
  sim <- simulate_auc_matrix(
    n = 100, prevalence = 0.3, target_aucs = c(0.8, 0.7),
    correlation = 0.3, structure = "exchangeable", seed = 1
  )
  expect_output(print(sim), "aucmat_simulation")
})
