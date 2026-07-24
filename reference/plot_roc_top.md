# ROC curves for a small set of biomarkers

Plots empirical ROC curves for deliberately selected biomarkers.
Designed for small sets; labelling is explicit.

## Usage

``` r
plot_roc_top(
  fit,
  X = NULL,
  y = NULL,
  biomarkers = NULL,
  show_auc = TRUE,
  add_ci = TRUE,
  boot_n = 500
)
```

## Arguments

- fit:

  An `aucmat_screen` object (may not retain data).

- X:

  Numeric matrix (required if `fit` does not retain data).

- y:

  Binary outcome (required if `fit` does not retain data).

- biomarkers:

  Character vector of biomarker names. Default: top 6 by `auc_strength`.

- show_auc:

  Add AUC values to the legend. Default `TRUE`.

- add_ci:

  Add confidence ribbons to ROC curves via bootstrap. Default `TRUE` –
  set `FALSE` for a faster, ribbon-free plot.

- boot_n:

  Number of bootstrap replicates for CI ribbons. Default 500.

## Value

A `ggplot2` object.

## Examples

``` r
# \donttest{
set.seed(42)
sim <- simulate_auc_matrix(n=100, prevalence=0.3,
  target_aucs=c(0.85,0.75), correlation=0.3, structure="exchangeable")
X <- as.matrix(sim$data[,1:2]); y <- sim$data$truth
fit <- aucmat(X, y, ci="none")
plot_roc_top(fit, X, y, biomarkers=c("X1","X2"))

# }
```
