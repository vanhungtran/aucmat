#!/usr/bin/env Rscript
# Run this script to regenerate documentation and test the legend fix
# Usage: Rscript run_documentation_and_test.R
# or source("run_documentation_and_test.R") from R console

cat("=== Regenerating R package documentation ===\n")
library(roxygen2)
roxygenise()

cat("\n=== Running legend test ===\n")
library(aucmat)
library(ggplot2)

# Generate test data
set.seed(123)
target_aucs <- c(0.9, 0.8, 0.7)
corr_matrix <- matrix(c(
  1.0, 0.3, 0.1,
  0.3, 1.0, 0.2,
  0.1, 0.2, 1.0
), nrow = 3)

sim_data <- generate_data_analytical(
  n = 500,
  prevalence = 0.3,
  target_aucs = target_aucs,
  corr_matrix = corr_matrix
)

# Get predictor names
preds <- setdiff(names(sim_data$data), "truth")

# Create plot with legends
cat("\n=== Creating ROC plot with legend ===\n")
res_plot <- plot_roc_with_combos(
  data = sim_data$data,
  outcome = "truth",
  predictors = preds,
  combo_sizes = 1:3,
  add_ci = TRUE,
  conf_level = 0.95,
  boot_n = 500,
  seed = 123,
  auc_digits = 7,
  show_auc_labels = TRUE,
  title = "ROC curves for predictor combinations (simulated data)"
)

# Check the structure
cat("\n=== Result structure ===\n")
cat("Components:", paste(names(res_plot), collapse = ", "), "\n")
cat("Number of ROC curves:", length(res_plot$roc_list), "\n")
cat("\nAUC values:\n")
print(res_plot$auc)

# Display the plot
cat("\n=== Displaying plot ===\n")
cat("The plot should show:\n")
cat("1. ROC curves for all predictor combinations\n")
cat("2. A legend on the right side showing model names with AUC values\n")
cat("3. AUC text labels on the curves (if show_auc_labels = TRUE)\n")
cat("4. Confidence interval ribbons around the curves\n\n")

print(res_plot$plot)

# Save to file
output_file <- "roc_plot_with_legend.png"
ggsave(output_file, res_plot$plot, width = 8, height = 6, dpi = 150)
cat("\nPlot saved to:", output_file, "\n")

cat("\n=== Test completed successfully! ===\n")
cat("If the legend is still not visible, please check:\n")
cat("1. You're accessing res_plot$plot (not just res_plot)\n")
cat("2. Your graphics device supports legends\n")
cat("3. The saved PNG file shows the legend correctly\n")
