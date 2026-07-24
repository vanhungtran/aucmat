# Plot simulated biomarker data

Default plot for `aucmat_simulation` objects. Shows a pairs plot of
biomarker columns with class-coloured density panels.

## Usage

``` r
# S3 method for class 'aucmat_simulation'
plot(x, ...)
```

## Arguments

- x:

  An `aucmat_simulation` object.

- ...:

  Ignored.

## Value

A `ggplot2` object (via GGally if available) or a base-R pairs plot.

## Examples

``` r
# \donttest{
set.seed(42)
sim <- simulate_auc_matrix(n=50, prevalence=0.3,
  target_aucs=c(0.8,0.7,0.6), correlation=0.3, structure="exchangeable")
plot(sim)

# }
```
