source("c:/Users/tranh/OneDrive/aucmat/R/generate_auc_vector.R")

set.seed(123)
AUC <- c(seq(0.5, 0.99, by = 0.01), 1)
y <- rbinom(1000, size = 1, prob = 0.3)

sim <- simulate_auc_correlation(
  y = y,
  target_auc = AUC,
  n_sim = 10000,
  seed = 123
)

results <- sim$results

summary_df <- do.call(rbind, lapply(split(results, results$target_auc), function(df) {
  data.frame(
    target_auc = df$target_auc[1],
    auc_mean = mean(df$auc),
    auc_sd = sd(df$auc),
    auc_q025 = unname(quantile(df$auc, 0.025)),
    auc_q50 = unname(quantile(df$auc, 0.5)),
    auc_q975 = unname(quantile(df$auc, 0.975)),
    cor_mean = mean(df$correlation),
    cor_sd = sd(df$correlation),
    cor_q025 = unname(quantile(df$correlation, 0.025)),
    cor_q50 = unname(quantile(df$correlation, 0.5)),
    cor_q975 = unname(quantile(df$correlation, 0.975))
  )
}))
rownames(summary_df) <- NULL

write.csv(summary_df, "c:/Users/tranh/OneDrive/aucmat/outputs/auc_grid_summary.csv", row.names = FALSE)

png("c:/Users/tranh/OneDrive/aucmat/outputs/auc_grid_summary_plot.png", width = 1800, height = 900, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(summary_df$target_auc, summary_df$auc_mean,
     type = "l", lwd = 2, col = "#1b6ca8", ylim = c(0.45, 1.02),
     xlab = "Target AUC", ylab = "Empirical AUC", main = "AUC Mean and 95% Interval")
lines(summary_df$target_auc, summary_df$auc_q025, lty = 2, col = "#4c78a8")
lines(summary_df$target_auc, summary_df$auc_q975, lty = 2, col = "#4c78a8")
abline(0, 1, lty = 3, col = "grey50")
legend("topleft",
       legend = c("Mean empirical AUC", "2.5% and 97.5% quantiles", "Target = empirical"),
       lty = c(1, 2, 3), col = c("#1b6ca8", "#4c78a8", "grey50"), bty = "n", cex = 0.85)
plot(summary_df$target_auc, summary_df$cor_mean,
     type = "l", lwd = 2, col = "#c84c09",
     xlab = "Target AUC", ylab = "Correlation(y, X)", main = "Correlation Mean and 95% Interval")
lines(summary_df$target_auc, summary_df$cor_q025, lty = 2, col = "#dd8452")
lines(summary_df$target_auc, summary_df$cor_q975, lty = 2, col = "#dd8452")
par(op)
dev.off()

print(head(summary_df, 10))
print(tail(summary_df, 10))
cat(sprintf("rows=%d targets=%d\n", nrow(results), length(unique(results$target_auc))))
