# Global Wald test for equality of multiple correlated AUCs

Tests H0: AUC_1 = AUC_2 = ... = AUC_p on a common subject set using the
joint DeLong covariance matrix.

## Usage

``` r
compare_auc_global(fit = NULL, X, y, biomarkers = NULL, max_biomarkers = 100)
```

## Arguments

- fit:

  An `aucmat_screen` object (optional, for outcome encoding).

- X:

  Numeric matrix.

- y:

  Binary outcome.

- biomarkers:

  Character vector of biomarker names to include. Default: all.

- max_biomarkers:

  Safety limit. Default 100.

## Value

A list of class `aucmat_global_test` with components: `statistic`, `df`,
`p_value`, `aucs`, `covariance`, `contrasts`, `n_common`, `n_pos`,
`n_neg`, `status`.

## Examples

``` r
# \donttest{
set.seed(42)
sim <- simulate_auc_matrix(n=100, prevalence=0.3,
  target_aucs=c(0.85,0.75,0.65), correlation=0.3, structure="exchangeable")
X <- as.matrix(sim$data[,1:3]); y <- sim$data$truth
compare_auc_global(X=X, y=y)
#> <aucmat_global_test>
#>   H0: all AUCs equal
#>   n =100 (30+ / 70-)
#>   Wald chi2(2) = 10.412, p = 0.0055
#> 
#> AUC estimates:
#>     X1     X2     X3 
#> 0.8338 0.7848 0.6167 
# }
```
