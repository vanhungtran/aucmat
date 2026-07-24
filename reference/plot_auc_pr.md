# Plot Precision-Recall curves for selected biomarkers

PR curves show the trade-off between precision (positive predictive
value) and recall (sensitivity). More informative than ROC when the
positive class is rare (prevalence \<\< 0.5).

## Usage

``` r
plot_auc_pr(X, y, biomarkers = NULL, show_auc = TRUE, show_baseline = TRUE)
```

## Arguments

- X:

  Numeric matrix.

- y:

  Binary outcome.

- biomarkers:

  Character vector of biomarker names. Default: top 6 by PR-AUC.

- show_auc:

  Add PR-AUC values to legend. Default `TRUE`.

- show_baseline:

  Add horizontal line at prevalence (random classifier). Default `TRUE`.

## Value

A `ggplot2` object.

## Examples

``` r
# \donttest{
set.seed(1)
sim <- simulate_auc_matrix(n = 500, prevalence = 0.15,
  target_aucs = c(0.85, 0.75, 0.65),
  correlation = 0.3, structure = "exchangeable")
X <- as.matrix(sim$data[, 1:3])
y <- sim$data$truth
plot_auc_pr(X, y)

# }
```
