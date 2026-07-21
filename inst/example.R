
# ==============================================================================
# aucmat: Example usage
# ==============================================================================

library(aucmat)

# --- Define Simulation Parameters ---
set.seed(42)

n_rows <- 5000
p_variables <- 5
prevalence <- 0.25

# Generate target AUCs and correlation matrix
target_aucs <- round(seq(from = 0.90, to = 0.70, length.out = p_variables), 2)

rho <- 0.6
corr_matrix <- matrix(nrow = p_variables, ncol = p_variables)
for (i in 1:p_variables) {
  for (j in 1:p_variables) {
    corr_matrix[i, j] <- rho^abs(i - j)
  }
}

# --- Run Simulation ---
simulated_dataset <- generate_data_analytical(
  n = n_rows,
  prevalence = prevalence,
  target_aucs = target_aucs,
  corr_matrix = corr_matrix
)

# --- Check Results ---
cat("--- TARGETS ---\n")
cat("Target AUCs:\n")
print(target_aucs)
cat("\nTarget Correlation Matrix:\n")
print(round(corr_matrix, 3))

cat("\n\n--- ACHIEVED RESULTS ---\n")
cat("Achieved AUCs:\n")
print(round(simulated_dataset$achieved_aucs, 3))
cat("\nAchieved Correlation Matrix:\n")
print(round(simulated_dataset$achieved_correlations, 3))


# ==============================================================================
# plot_roc_with_combos example
# ==============================================================================

set.seed(42)

# Generate data with known AUC properties
sim_data <- generate_data_analytical(
  n = 200,
  prevalence = 0.3,
  target_aucs = c(0.85, 0.75, 0.65),
  corr_matrix = diag(3)
)

preds <- setdiff(names(sim_data$data), "truth")

res <- plot_roc_with_combos(
  data = sim_data$data,
  outcome = "truth",
  predictors = preds,
  combo_sizes = 1:2,
  add_ci = TRUE,
  boot_n = 1500,
  title = "ROC curves for single biomarkers and pairs"
)

print(res$plot)
cat("\nAUC values:\n")
print(res$auc)


# ==============================================================================
# tableroc example
# ==============================================================================

set.seed(123)
X <- data.frame(
  M1 = c(rnorm(95), rep(NA, 5)),
  M2 = rnorm(100),
  M3 = c(rep(NA, 3), rnorm(97))
)
y <- factor(rbinom(100, 1, 0.5), labels = c("neg", "pos"))

auc_table <- tableroc(X = X, y = y, na_impute = "median", rank = TRUE)
print(auc_table)


# ==============================================================================
# generate_auc_vector example
# ==============================================================================

set.seed(1)
y <- rbinom(200, size = 1, prob = 0.3)
out <- generate_auc_vector(y, target_auc = 0.8)
cat("\ngenerate_auc_vector:\n")
cat("Achieved AUC:", out$achieved_auc, "\n")
