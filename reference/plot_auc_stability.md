# Stability plot: rank distributions and top-k probabilities

Stability plot: rank distributions and top-k probabilities

## Usage

``` r
plot_auc_stability(stability, n_label = 25L)
```

## Arguments

- stability:

  An `aucmat_stability` object.

- n_label:

  Number of biomarkers to show. Default 25.

## Value

A `ggplot2` object.

## Examples

``` r
# \donttest{
set.seed(42)
X <- matrix(rnorm(200*10), 200, 10, dimnames=list(NULL, paste0("bm",1:10)))
y <- rep(c(0,1), each=100)
stab <- auc_stability(X, y, times=50, seed=1)
plot_auc_stability(stab, n_label=8)

# }
```
