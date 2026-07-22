# ==============================================================================
# Simulation study for aucmat manuscript
# ==============================================================================
devtools::load_all(".", quiet = TRUE)

cat("=== aucmat Manuscript Simulation Study ===\n")
cat("Package version:", as.character(packageVersion("aucmat")), "\n\n")

# ---- 1. Coverage and Power (representative conditions) ----
cat("1. Coverage and Power\n--------------------\n")

design <- expand.grid(
  n          = c(100, 300, 500),
  prevalence = c(0.2, 0.3, 0.5),
  auc        = c(0.65, 0.75, 0.85),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
n_sim <- 200L

results <- list()
for (d in seq_len(nrow(design))) {
  dd <- design[d, ]
  covered <- 0L; rejected <- 0L; n_valid <- 0L
  auc_ests <- numeric(n_sim)

  for (r in seq_len(n_sim)) {
    seed_val <- d * 10000L + r
    sim <- tryCatch(
      simulate_auc_matrix(
        n = dd$n, prevalence = dd$prevalence,
        target_aucs = rep(dd$auc, 2),
        correlation = 0.2, structure = "exchangeable",
        seed = seed_val, verify = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(sim)) next

    fit <- aucmat(as.matrix(sim$data[, 1:2]), sim$data$truth,
      ci = "delong", adjust = "none")

    n_valid <- n_valid + 1L
    ci_lo <- fit$results$conf_low[1]
    ci_hi <- fit$results$conf_high[1]
    if (!is.na(ci_lo) && !is.na(ci_hi) && dd$auc >= ci_lo && dd$auc <= ci_hi)
      covered <- covered + 1L
    if (!is.na(fit$results$p_value[1]) && fit$results$p_value[1] < 0.05)
      rejected <- rejected + 1L
    auc_ests[r] <- fit$results$auc_raw[1]
  }

  if (n_valid < 10L) next

  results[[d]] <- data.frame(
    n = dd$n, prevalence = dd$prevalence, auc = dd$auc,
    n_valid = n_valid,
    coverage = round(covered / n_valid, 4),
    power    = round(rejected / n_valid, 4),
    bias     = round(mean(auc_ests[1:n_valid]) - dd$auc, 5),
    rmse     = round(sqrt(mean((auc_ests[1:n_valid] - dd$auc)^2)), 5)
  )
  cat(sprintf("  n=%d prev=%.1f auc=%.2f: cov=%.3f pow=%.3f bias=%.5f\n",
    dd$n, dd$prevalence, dd$auc, tail(results, 1)[[1]]$coverage,
    tail(results, 1)[[1]]$power, tail(results, 1)[[1]]$bias))
}

res_df <- do.call(rbind, results)
cat("\nCoverage by n and AUC:\n")
print(aggregate(coverage ~ n + auc, data = res_df, mean))
cat("\nPower by n and AUC:\n")
print(aggregate(power ~ n + auc, data = res_df, mean))

# ---- 2. Type I Error under Null ----
cat("\n2. Type I Error (Null)\n----------------------\n")
set.seed(42)
null_rej <- 0L; null_valid <- 0L
for (r in seq_len(500)) {
  y <- factor(sample(c("neg", "pos"), 200, replace = TRUE))
  X <- matrix(rnorm(200 * 5), 200, 5, dimnames = list(NULL, paste0("X", 1:5)))
  fit <- aucmat(X, y, ci = "delong", adjust = "none")
  if (!is.na(fit$results$p_value[1]) && fit$results$p_value[1] < 0.05)
    null_rej <- null_rej + 1L
  null_valid <- null_valid + 1L
}
cat(sprintf("  Type I error: %.4f (nominal 0.05, n_sim=%d)\n",
  null_rej / null_valid, null_valid))

# ---- 3. Benchmarks ----
cat("\n3. Benchmarks\n-------------\n")

bench <- list(
  c(name = "Screen p=100",    p = 100,  n = 500),
  c(name = "Screen p=1000",   p = 1000, n = 500),
  c(name = "Screen p=10000",  p = 10000,n = 500)
)

for (b in bench) {
  p <- as.integer(b["p"]); n <- as.integer(b["n"])
  Xb <- matrix(rnorm(n * p), n, p,
    dimnames = list(NULL, paste0("X", seq_len(p))))
  yb <- sample(c(0, 1), n, replace = TRUE)
  t0 <- Sys.time()
  fit <- aucmat(Xb, yb, ci = "delong")
  t1 <- Sys.time()
  cat(sprintf("  %-20s: %.2f s\n", b["name"],
    as.numeric(difftime(t1, t0, units = "secs"))))
}

# Bootstrap stability
Xb <- matrix(rnorm(500 * 50), 500, 50,
  dimnames = list(NULL, paste0("X", seq_len(50))))
yb <- sample(c(0, 1), 500, replace = TRUE)
t0 <- Sys.time()
stab <- auc_stability(Xb, yb, times = 200, seed = 42)
t1 <- Sys.time()
cat(sprintf("  Stability B=200,p=50 : %.2f s\n",
  as.numeric(difftime(t1, t0, units = "secs"))))

# Validate
t0 <- Sys.time()
val <- validate_simulation(n = 300, prevalence = 0.3,
  target_aucs = c(0.85, 0.75, 0.65),
  correlation = 0.3, structure = "exchangeable",
  times = 50, seed = 1)
t1 <- Sys.time()
cat(sprintf("  Validate reps=50,p=3 : %.2f s\n",
  as.numeric(difftime(t1, t0, units = "secs"))))

saveRDS(res_df, "manuscript/simulation_results.rds")
cat("\nDone. Results saved to manuscript/simulation_results.rds\n")
