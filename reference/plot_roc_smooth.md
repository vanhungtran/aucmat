# Plot smoothed ROC curves for selected biomarkers

As
[`plot_roc_top()`](https://vanhungtran.github.io/aucmat/reference/plot_roc_top.md)
but adds a binormal-smoothed ROC curve overlay using `pROC::smooth()`.
The smoothed curve reduces step-artifact noise and is the standard for
publication-quality ROC figures.

## Usage

``` r
plot_roc_smooth(
  fit,
  X,
  y,
  biomarkers = NULL,
  smooth_method = c("binormal", "density"),
  show_empirical = TRUE,
  add_ci = TRUE,
  boot_n = 500
)
```

## Arguments

- fit:

  An `aucmat_screen` object.

- X:

  Numeric matrix.

- y:

  Binary outcome.

- biomarkers:

  Character vector. Default: top 6 by AUC.

- smooth_method:

  `"binormal"` (default) or `"density"`.

- show_empirical:

  Also plot the empirical (step) ROC curve. Default `TRUE`.

- add_ci:

  Add bootstrap CI ribbons. Default `TRUE`.

- boot_n:

  Bootstrap replicates for CI ribbons.

## Value

A `ggplot2` object.

## Examples

``` r
# \donttest{
set.seed(1)
sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
  target_aucs = c(0.85, 0.75, 0.65),
  correlation = 0.3, structure = "exchangeable")
fit <- aucmat(as.matrix(sim$data[, 1:3]), sim$data$truth, ci = "none")
plot_roc_smooth(fit, as.matrix(sim$data[, 1:3]), sim$data$truth)
#> NULL
# }
```
