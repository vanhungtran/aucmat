# Volcano plot: effect magnitude against statistical evidence

Volcano plot: effect magnitude against statistical evidence

## Usage

``` r
plot_auc_volcano(fit, n_label = 20L, q_cutoff = 0.05)
```

## Arguments

- fit:

  An `aucmat_screen` object.

- n_label:

  Number of top biomarkers to label. Default 20.

- q_cutoff:

  Highlight biomarkers with q-value below this threshold.

## Value

A `ggplot2` object.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(100*20), 100, 20, dimnames=list(NULL, paste0("bm",1:20)))
fit <- aucmat(X, rep(0:1, each=50), ci="none")
plot_auc_volcano(fit, q_cutoff=0.1)
#> Warning: Removed 20 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 20 rows containing missing values or values outside the scale range
#> (`geom_text_repel()`).
```
