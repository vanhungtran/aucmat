# ==============================================================================
# Generate all README plot images
# Run from the package root: source("generate_readme_plots.R")
# ==============================================================================
devtools::load_all(".", quiet = TRUE)
set.seed(42)

fig <- "man/figures"
if (!dir.exists(fig)) dir.create(fig, recursive = TRUE)

# ---- Common data for all plots ----
sim <- generate_data_probit(
  n = 500,
  target_aucs = c(0.90, 0.80, 0.70, 0.65),
  corr_matrix = matrix(c(
    1.0, 0.4, 0.2, 0.1,
    0.4, 1.0, 0.3, 0.2,
    0.2, 0.3, 1.0, 0.3,
    0.1, 0.2, 0.3, 1.0
  ), 4, 4),
  prevalence = 0.3
)

X <- as.matrix(sim$data[, 1:4])
y <- sim$data$truth
fit <- aucmat(X, y, ci = "delong", adjust = "BH")
stab <- auc_stability(X, y, times = 200, top_k = c(3, 5), seed = 42)

# ---- 1. Rank plot ----
png(file.path(fig, "README-rank.png"), width = 900, height = 480, res = 110)
print(plot_auc_rank(fit, n_label = 4))
dev.off()

# ---- 2. Volcano plot ----
png(file.path(fig, "README-volcano.png"), width = 900, height = 480, res = 110)
print(plot_auc_volcano(fit, q_cutoff = 0.05))
dev.off()

# ---- 3. Forest plot ----
png(file.path(fig, "README-forest.png"), width = 900, height = 400, res = 110)
print(plot_auc_forest(fit, n = 4))
dev.off()

# ---- 4. ROC curves with CI ribbons ----
png(file.path(fig, "README-roc.png"), width = 600, height = 550, res = 110)
print(plot_roc_top(fit, X = X, y = y,
  biomarkers = head(fit$results$biomarker, 3),
  add_ci = TRUE, boot_n = 500))
dev.off()

# ---- 5. Stability ----
png(file.path(fig, "README-stability.png"), width = 900, height = 450, res = 110)
print(plot_auc_stability(stab, n_label = 8))
dev.off()

# ---- 6. Simulator comparison: achieved vs target ----
# Use simulate_auc_matrix with all 4 structures to show flexibility
set.seed(1)
sim_ex  <- simulate_auc_matrix(n = 500, prevalence = 0.3,
  target_aucs = c(0.9, 0.8, 0.7, 0.65),
  correlation = 0.4, structure = "exchangeable")
sim_ar1 <- simulate_auc_matrix(n = 500, prevalence = 0.3,
  target_aucs = c(0.9, 0.8, 0.7, 0.65),
  correlation = 0.6, structure = "ar1")
sim_blk <- simulate_auc_matrix(n = 500, prevalence = 0.3,
  target_aucs = c(0.9, 0.8, 0.7, 0.65),
  structure = "block", block_sizes = c(2, 2),
  rho_within = 0.6, rho_between = 0.1)

# ---- 7. Validation plot: AUC bias across replicates ----
val <- validate_simulation(
  n = 300, prevalence = 0.3,
  target_aucs = c(0.85, 0.75, 0.65),
  correlation = 0.3, structure = "exchangeable",
  times = 100, seed = 1
)

# ---- 8. Non-inferiority comparison example ----
cmp_ni <- compare_auc(fit, X, y, reference = "X1",
  hypothesis = "noninferiority", margin = 0.05)

# ---- 9. Equivalence comparison example ----
cmp_eq <- compare_auc(fit, X, y, reference = "X1",
  hypothesis = "equivalence", margin = 0.15)

# ---- 10. Global test on more biomarkers ----
sim_big <- simulate_auc_matrix(n = 500, prevalence = 0.3,
  target_aucs = c(0.85, 0.80, 0.75, 0.70, 0.65),
  correlation = 0.3, structure = "exchangeable", seed = 3)
X_big <- as.matrix(sim_big$data[, 1:5])
y_big <- sim_big$data$truth
fit_big <- aucmat(X_big, y_big, ci = "delong", adjust = "BH")
global_test <- compare_auc_global(fit_big, X_big, y_big)

cat("All README plots saved to man/figures/\n")
cat("\n=== Summary of generated data ===\n")
cat("Exchangeable sim AUCs:", round(sim_ex$achieved_aucs, 4), "\n")
cat("AR1 sim AUCs:", round(sim_ar1$achieved_aucs, 4), "\n")
cat("Block sim AUCs:", round(sim_blk$achieved_aucs, 4), "\n")
cat("\n=== Validation ===\n")
print(val)
cat("\n=== Non-inferiority ===\n")
print(cmp_ni)
cat("\n=== Equivalence ===\n")
print(cmp_eq)
cat("\n=== Global test ===\n")
print(global_test)
