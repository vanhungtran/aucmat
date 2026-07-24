# Test whether a combined panel improves on its individual biomarkers

Treats the panel's honest cross-validated score as one more "biomarker"
and reuses
[`compare_auc()`](https://vanhungtran.github.io/aucmat/reference/compare_auc.md)'s
DeLong or bootstrap paired-comparison engine to test it against each
individual biomarker used to build the panel, on the same subjects.
Answers "did combining help?" with a confidence interval and p-value
rather than a side-by-side comparison of two point estimates.

## Usage

``` r
compare_panel_auc(
  panel,
  X,
  y,
  positive = NULL,
  biomarkers = NULL,
  use = c("cv", "train"),
  ...
)
```

## Arguments

- panel:

  An `aucmat_panel` object from
  [`fit_auc_panel()`](https://vanhungtran.github.io/aucmat/reference/fit_auc_panel.md).

- X, y:

  The same data (and outcome) used to fit `panel`.

- positive:

  Optional positive class label.

- biomarkers:

  Character vector of individual biomarkers to compare against. Default:
  all biomarkers used in the panel.

- use:

  `"cv"` (default) compares the honest cross-validated panel score;
  `"train"` compares the in-sample (optimistic) score. Errors if `"cv"`
  is requested but `panel` was fit with `n_folds < 2`.

- ...:

  Passed to
  [`compare_auc()`](https://vanhungtran.github.io/aucmat/reference/compare_auc.md)
  (`hypothesis`, `alternative`, `margin`, `conf_level`, `method`,
  `boot_n`, `seed`, `adjust`).

## Value

A data.frame of class `aucmat_compare` (see
[`compare_auc()`](https://vanhungtran.github.io/aucmat/reference/compare_auc.md)).
`biomarker_a` is always the panel score.

## Examples

``` r
# \donttest{
set.seed(42)
sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
  target_aucs = c(0.85, 0.75, 0.65, 0.60),
  correlation = 0.3, structure = "exchangeable")
X <- as.matrix(sim$data[, 1:4]); y <- sim$data$truth
panel <- fit_auc_panel(X, y, method = "ridge", n_folds = 5, seed = 1)
compare_panel_auc(panel, X, y)
#> <aucmat_compare>  4 pairwise comparisons
#>  biomarker_a biomarker_b     auc_a     auc_b     auc_diff  std_error
#>  PANEL_SCORE          X1 0.8066138 0.8119048 -0.005291005 0.00793003
#>  PANEL_SCORE          X2 0.8066138 0.6686772  0.137936508 0.03249139
#>  PANEL_SCORE          X3 0.8066138 0.5491005  0.257513228 0.04262894
#>  PANEL_SCORE          X4 0.8066138 0.5843915  0.222222222 0.04062107
#>     conf_low  conf_high      p_value q_value n_common n_pos n_neg  hypothesis
#>  -0.02083358 0.01025157 5.046372e-01      NA      300    90   210 superiority
#>   0.07425455 0.20161846 2.182769e-05      NA      300    90   210 superiority
#>   0.17396204 0.34106441 1.533444e-09      NA      300    90   210 superiority
#>   0.14260638 0.30183806 4.484780e-08      NA      300    90   210 superiority
#>  margin
#>       0
#>       0
#>       0
#>       0
# }
```
