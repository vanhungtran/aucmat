# Forest plot: selected AUCs with confidence intervals

Forest plot: selected AUCs with confidence intervals

## Usage

``` r
plot_auc_forest(fit, biomarkers = "top", n = 20L)
```

## Arguments

- fit:

  An `aucmat_screen` object.

- biomarkers:

  Character vector of biomarker names to show, or `"top"` to take the
  top `n` by `auc_strength`.

- n:

  Number of biomarkers when `biomarkers = "top"`. Default 20.

## Value

A `ggplot2` object.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(100*20), 100, 20, dimnames=list(NULL, paste0("bm",1:20)))
fit <- aucmat(X, rep(0:1, each=50), ci="none")
plot_auc_forest(fit, n=5)
```
