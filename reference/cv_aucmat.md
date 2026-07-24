# Cross-validated AUC screening

Computes out-of-fold AUCs for every biomarker via stratified v-fold
cross-validation. Each biomarker's pooled out-of-fold predictions across
all folds are used to compute a single honest AUC that is not inflated
by same-data overfitting.

## Usage

``` r
cv_aucmat(
  X,
  y,
  positive = NULL,
  n_folds = 5L,
  n_repeats = 1L,
  seed = NULL,
  ci = c("none", "delong"),
  conf_level = 0.95
)
```

## Arguments

- X:

  Numeric matrix (n x p).

- y:

  Binary outcome.

- positive:

  Optional positive class label.

- n_folds:

  Number of cross-validation folds. Default 5.

- n_repeats:

  Number of repeats of the full CV. Default 1.

- seed:

  Seed for fold assignment reproducibility.

- ci:

  CI method for final pooled AUCs: `"delong"` or `"none"`.

- conf_level:

  Confidence level.

## Value

A list of class `aucmat_cv` with `results`, `folds`, `settings`.

## Examples

``` r
# \donttest{
set.seed(1)
sim <- simulate_auc_matrix(n = 200, prevalence = 0.3,
  target_aucs = c(0.8, 0.7, 0.6), correlation = 0.3,
  structure = "exchangeable")
cv <- cv_aucmat(as.matrix(sim$data[, 1:3]), sim$data$truth,
  n_folds = 3, seed = 42)
print(cv)
#> <aucmat_cv>  3 biomarkers
#>   CV: 3-fold
#>   Mean optimism:  0 
#> 
#> Top biomarkers by CV AUC:
#>  rank_cv biomarker    auc_cv auc_in_sample optimism
#>        1        X1 0.8464286     0.8464286        0
#>        2        X2 0.6970238     0.6970238        0
#>        3        X3 0.6398810     0.6398810        0
# }
```
