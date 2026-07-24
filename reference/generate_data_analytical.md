# Generate Synthetic Data with Specific AUC and Correlation

This function analytically generates a dataset with specified
properties.

## Usage

``` r
generate_data_analytical(n, prevalence, target_aucs, corr_matrix)
```

## Arguments

- n:

  number of rows (observations)

- prevalence:

  the proportion of the positive class for the binary outcome

- target_aucs:

  a vector of p target AUC values

- corr_matrix:

  a p x p target correlation matrix

## Examples

``` r
# \donttest{
set.seed(42)
sim <- generate_data_analytical(n=100, prevalence=0.3,
  target_aucs=c(0.8,0.7), corr_matrix=matrix(c(1,0.3,0.3,1),2,2))
round(sim$achieved_aucs, 3)
#> [1] 0.82 0.68
# }
```
