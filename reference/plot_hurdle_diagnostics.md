# Plot Hurdle-AUC diagnostics

Visualizes the zero-inflation pattern and AUC decomposition for
hurdle-model results.

## Usage

``` r
plot_hurdle_diagnostics(fit, n_label = 15L)
```

## Arguments

- fit:

  An `aucmat_hurdle` object from
  [`hurdle_auc()`](https://vanhungtran.github.io/aucmat/reference/hurdle_auc.md).

- n_label:

  Number of biomarkers to label. Default 15.

## Value

A `ggplot2` object.

## Examples

``` r
# \donttest{
set.seed(1)
sim <- simulate_hurdle_auc(n=100, prevalence=0.3,
  target_hurdle_aucs=c(0.8,0.7), zero_rate_neg=c(0.5,0.3),
  zero_rate_pos=c(0.2,0.1))
fit <- hurdle_auc(as.matrix(sim$data[,1:2]), sim$data$truth)
plot_hurdle_diagnostics(fit)

# }
```
