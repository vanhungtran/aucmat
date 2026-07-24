# Rank plot: ordered discrimination strengths

Rank plot: ordered discrimination strengths

## Usage

``` r
plot_auc_rank(fit, n_label = 20L, show_ci = TRUE)
```

## Arguments

- fit:

  An `aucmat_screen` object.

- n_label:

  Number of top biomarkers to label. Default 20.

- show_ci:

  If `TRUE` and CIs are present, add error bars.

## Value

A `ggplot2` object.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(100*20), 100, 20, dimnames=list(NULL, paste0("bm",1:20)))
fit <- aucmat(X, rep(0:1, each=50), ci="none")
plot_auc_rank(fit, n_label=5)
```
