library(aucmat)
library(ggplot2)

cat("=== Testing ROC plot with legend and CI ribbons ===\n\n")

# Generate test data
set.seed(123)
target_aucs <- c(0.9, 0.8, 0.7)
corr_matrix <- matrix(c(
  1.0, 0.3, 0.1,
  0.3, 1.0, 0.2,
  0.1, 0.2, 1.0
), nrow = 3)

cat("Generating synthetic data...\n")
sim_data <- generate_data_analytical(
  n = 500,
  prevalence = 0.3,
  target_aucs = target_aucs,
  corr_matrix = corr_matrix
)

# Get predictor names
preds <- setdiff(names(sim_data$data), "truth")
cat("Predictors:", paste(preds, collapse = ", "), "\n\n")

# Create plot with legends and CI ribbons
cat("Creating ROC plot with CI ribbons and AUC ± SD legend...\n")
res_plot <- plot_roc_with_combos(
  data = sim_data$data,
  outcome = "truth",
  predictors = preds,
  combo_sizes = 1:3,
  add_ci = TRUE,
  conf_level = 0.95,
  boot_n = 500,
  seed = 123,
  auc_digits = 4,
  show_auc_labels = TRUE,
  show_auc_sd = TRUE,  # Show AUC ± SD in legend (default)
  # sd_digits defaults to auc_digits (4), so both AUC and SD will have 4 decimal places
  title = "ROC curves for predictor combinations (simulated data)"
)

# Check the plot object
cat("\n=== Result structure ===\n")
cat("Components:", paste(names(res_plot), collapse = ", "), "\n")
cat("Plot class:", class(res_plot$plot), "\n")
cat("Number of ROC curves:", length(res_plot$roc_list), "\n")
cat("\nAUC values:\n")
print(res_plot$auc)
cat("\nAUC standard deviations:\n")
print(res_plot$auc_sd)

# Check CI data
if (!is.null(res_plot$ci_df) && nrow(res_plot$ci_df) > 0) {
  cat("\nCI data frame structure:\n")
  cat("Columns:", paste(names(res_plot$ci_df), collapse = ", "), "\n")
  cat("Number of CI points:", nrow(res_plot$ci_df), "\n")
  cat("Models with CI:", length(unique(res_plot$ci_df$name)), "\n")
}

# Print the plot
cat("\n=== Displaying plot ===\n")
cat("The plot should show:\n")
cat("1. ROC curves for all predictor combinations\n")
cat("2. A legend on the right with format: 'Model (AUC=0.xxxx ± 0.xxxx)'\n")
cat("3. AUC text labels on the curves\n")
cat("4. CI ribbons OVERLAPPING the ROC curves (same direction)\n\n")

print(res_plot$plot)

# Save to file to inspect
ggsave("test_plot.png", res_plot$plot, width = 8, height = 6, dpi = 150)
cat("\nPlot saved to test_plot.png\n")
cat("\n=== Test completed! ===\n")
cat("Check that CI ribbons overlap with ROC curves in the saved PNG file.\n")
