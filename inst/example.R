
# --- Define Simulation Parameters ---
set.seed(42) # For reproducibility

# Set the desired number of rows (n) and variables (p)
n_rows <- 5000
p_variables <- 5

# Set the prevalence of the positive class
prevalence <- 0.25

# --- Generate Targets for p Variables ---

# Create a vector of p target AUCs
target_aucs <- round(seq(from = 0.90, to = 0.70, length.out = p_variables), 2)

# Create a p x p target correlation matrix
# Here, we use a common structure where correlation decays exponentially
rho <- 0.6
corr_matrix <- matrix(nrow = p_variables, ncol = p_variables)
for (i in 1:p_variables) {
  for (j in 1:p_variables) {
    corr_matrix[i, j] <- rho^abs(i - j)
  }
}

# --- Run the Simulation ---
simulated_dataset <- generate_data_analytical(
  n = n_rows,
  prevalence = prevalence,
  target_aucs = target_aucs,
  corr_matrix = corr_matrix
)

# --- Check the Results ---
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

cat("\n\n--- SAMPLE OF GENERATED DATA ---\n")
print(head(simulated_dataset$data))














######     plot_roc_with_combos ####################################





# Set seed for reproducibility
set.seed(42)

# Number of samples
n <- 200

# Generate synthetic data
data <- data.frame(
  M1 = rnorm(n, mean = 0, sd = 1),
  M2 = rnorm(n, mean = 5, sd = 2),
  M3 = rnorm(n, mean = 10, sd = 3),
  M4 = rnorm(n, mean = 15, sd = 4),
  Outcome = sample(c(0, 1), n, replace = TRUE)
)

# Split into training (80%) and validation (20%) sets
train_indices <- sample(1:n, size = 0.8 * n)
train_data <- data[train_indices, ]
valid_data <- data[-train_indices, ]

# View first few rows
head(train_data)
head(valid_data)


res <- plot_roc_with_combos(
  train = train_data,
  valid = valid_data,
  outcome = "Outcome",
  positive = 1,            # <- adjust to your dataset
  predictors = c("M1", "M2", "M3", "M4"),
  use_predictor_pool = 10,      # mimic your vars[1:10] pool
  combo_sizes = 1:4,            # singles through 5-way combos
  unify_valid_rows = FALSE,     # set TRUE for identical VALID rows across curves
  add_ci = TRUE,
  boot_n = 1500,
  ci_points = 201,
  skip_ci_for_auc1 = TRUE,
  smooth = FALSE,
  title = "Endo: ROC (VALID) for single biomarkers and combinations"
)
























