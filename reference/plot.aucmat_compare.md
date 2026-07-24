# Forest plot of paired AUC differences

Plots pairwise AUC differences with confidence intervals from a
[`compare_auc()`](https://vanhungtran.github.io/aucmat/reference/compare_auc.md)
result.

## Usage

``` r
# S3 method for class 'aucmat_compare'
plot(x, ...)
```

## Arguments

- x:

  An `aucmat_compare` data.frame.

- ...:

  Ignored.

## Value

A `ggplot2` object.

## Examples

``` r
# \donttest{
set.seed(42)
sim <- simulate_auc_matrix(n=100, prevalence=0.3,
  target_aucs=c(0.85,0.75,0.65), correlation=0.3, structure="exchangeable")
X <- as.matrix(sim$data[,1:3]); y <- sim$data$truth
fit <- aucmat(X, y, ci="none")
cmp <- compare_auc(fit, X, y, reference="X1")
plot(cmp)
#> `height` was translated to `width`.

# }
```
