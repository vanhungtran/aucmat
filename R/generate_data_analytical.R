

# Load necessary packages
# Ensure you have these installed: install.packages(c("pROC", "mvtnorm"))
library(pROC)
library(mvtnorm)

# This function analytically generates a dataset with specified properties.
# n: number of rows (observations)
# p: number of variables (determined by the inputs)
# prevalence: the proportion of the positive class for the binary outcome
# target_aucs: a vector of p target AUC values
# corr_matrix: a p x p target correlation matrix

generate_data_analytical <- function(n, prevalence, target_aucs, corr_matrix) {

  # Sort target_aucs and reorder corr_matrix to process from highest AUC to lowest
  auc_order <- order(target_aucs, decreasing = TRUE)
  sorted_aucs <- target_aucs[auc_order]
  sorted_corr_matrix <- corr_matrix[auc_order, auc_order]
  p <- length(sorted_aucs)

  # Initialize data and the binary outcome variable 'truth'
  final_data <- matrix(NA, nrow = n, ncol = p)
  truth <- rbinom(n = n, size = 1, prob = prevalence)
  n_pos <- sum(truth == 1)
  n_neg <- sum(truth == 0)

  # --- Step 1: Simulate the first variable ---
  # This variable has no preceding variables to correlate with.
  mu_1 <- sqrt(2) * qnorm(sorted_aucs[1])
  X1 <- rnorm(n)

  target_mu_0 <- mu_1 * (-n_pos / n)
  target_mu_1 <- mu_1 * (n_neg / n)

  X1[truth == 0] <- X1[truth == 0] - mean(X1[truth == 0]) + target_mu_0
  X1[truth == 1] <- X1[truth == 1] - mean(X1[truth == 1]) + target_mu_1
  final_data[, 1] <- as.vector(scale(X1))

  # --- Step 2: Sequentially simulate subsequent variables ---
  for (i in 2:p) {
    prev_data <- final_data[, 1:(i-1), drop = FALSE]
    target_cor_vec <- sorted_corr_matrix[i, 1:(i-1)]

    # 2a. Calculate the part of the new variable explained by previous variables
    prev_cor_matrix <- cor(prev_data)
    betas <- solve(prev_cor_matrix) %*% target_cor_vec
    xi_correlated_part <- scale(prev_data) %*% betas

    # 2b. Calculate the properties required for the residual (unexplained) part
    residual_variance <- as.numeric(1 - t(betas) %*% prev_cor_matrix %*% betas)
    if (residual_variance < 0) residual_variance <- 0 # Handle potential floating point errors
    residual_sd <- sqrt(residual_variance)

    # Determine the mean separation needed from the residual to hit the target AUC
    target_total_sep <- sqrt(2) * qnorm(sorted_aucs[i])
    sep_from_correlated_part <- mean(xi_correlated_part[truth == 1]) - mean(xi_correlated_part[truth == 0])
    required_sep_from_residual <- target_total_sep - sep_from_correlated_part

    # 2c. Generate a residual with these exact properties
    xi_residual <- rnorm(n, mean = 0, sd = residual_sd)

    res_target_mu_0 <- required_sep_from_residual * (-n_pos / n)
    res_target_mu_1 <- required_sep_from_residual * (n_neg / n)

    xi_residual[truth == 0] <- xi_residual[truth == 0] - mean(xi_residual[truth == 0]) + res_target_mu_0
    xi_residual[truth == 1] <- xi_residual[truth == 1] - mean(xi_residual[truth == 1]) + res_target_mu_1

    # 2d. Combine the correlated and residual parts to create the final variable
    Xi <- xi_correlated_part + xi_residual
    final_data[, i] <- as.vector(scale(Xi))
  }

  # --- Step 3: Format and return the final dataset ---
  data_df <- as.data.frame(final_data)

  # Return columns to their original, user-specified order
  original_order_colnames <- paste0("X", order(auc_order))
  colnames(data_df) <- original_order_colnames
  data_df <- data_df[, paste0("X", 1:p)]

  data_df$truth <- truth

  # --- Step 4: Final verification and return ---
  achieved_aucs <- sapply(1:p, function(j) pROC::auc(pROC::roc(data_df$truth, data_df[, j], quiet = TRUE)))
  achieved_correlations <- cor(data_df[, 1:p])

  return(list(
    "data" = data_df,
    "achieved_aucs" = achieved_aucs,
    "achieved_correlations" = achieved_correlations
  ))
}


