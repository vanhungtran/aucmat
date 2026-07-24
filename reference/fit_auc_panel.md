# Fit and validate a multivariable biomarker panel

Builds a combined biomarker score and reports both an in-sample AUC and
an honest, nested cross-validated AUC. "Nested" means every step that
looks at the outcome – penalty selection for the penalized methods, and
centering/scaling for all methods – is re-derived from the training
portion of each outer fold alone, never from data that includes the held
-out test rows. This prevents the optimistic bias that comes from tuning
on the full data and only refitting coefficients per fold.

## Usage

``` r
fit_auc_panel(
  X,
  y,
  positive = NULL,
  method = c("ridge", "lasso", "elasticnet", "logistic", "su_liu", "unweighted"),
  alpha = NULL,
  shrinkage = 0.1,
  n_folds = 5L,
  seed = NULL,
  standardize = TRUE
)
```

## Arguments

- X:

  Numeric matrix (n x p).

- y:

  Binary outcome.

- positive:

  Optional positive class label.

- method:

  `"ridge"` (default), `"lasso"`, `"elasticnet"`, `"logistic"`
  (unpenalized), `"su_liu"` (closed-form AUC-maximizing linear
  combination under the binormal model), or `"unweighted"` (equal-weight
  sum of direction-aligned standardized biomarkers; no fitting).

- alpha:

  Elastic-net mixing: 0 = ridge, 1 = lasso. Only used when
  `method = "elasticnet"`; defaults to 0.5 (a true ridge/lasso blend)
  when not supplied. Ignored for other methods.

- shrinkage:

  Shrinkage intensity in `[0, 1]` toward the diagonal covariance, used
  only by `method = "su_liu"`. `0` uses the raw pooled covariance
  (unstable or infeasible once `p` approaches `n`); `1` shrinks all the
  way to a diagonal (independence) covariance. Default 0.1.

- n_folds:

  Outer CV folds for the honest AUC (0 = no CV, uses training AUC only).
  Default 5.

- seed:

  Optional seed for CV fold-assignment reproducibility.

- standardize:

  Center/scale predictors before fitting. Default `TRUE`. Centering and
  scaling are always computed from the relevant training data only (the
  full data for the final model; each outer training fold for the honest
  CV score), never from held-out rows.

## Value

A list of class `aucmat_panel` with components: `coefficients`,
`auc_train`, `auc_cv`, `fitted_model` (final-model fit object together
with the centering/scaling used to predict on new data via
[`predict.aucmat_panel()`](https://vanhungtran.github.io/aucmat/reference/predict.aucmat_panel.md)),
`predictions`, `settings`, `sample_summary`.

## Examples

``` r
# \donttest{
set.seed(1)
sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
  target_aucs = c(0.85, 0.75, 0.65, 0.60),
  correlation = 0.3, structure = "exchangeable")
X <- as.matrix(sim$data[, 1:4]); y <- sim$data$truth

panel <- fit_auc_panel(X, y, method = "ridge", n_folds = 3, seed = 42)
print(panel)
#> <aucmat_panel>  ridge
#>   AUC (training): 0.8824
#>   AUC (CV):       0.8756
#>   Optimism:       0.0068
#> 
#> Coefficients:
#>      X1      X2      X4      X3 
#>  1.3201  0.4382 -0.2257  0.1952 

# Closed-form AUC-maximizing combination -- no glmnet dependency
panel_su <- fit_auc_panel(X, y, method = "su_liu", n_folds = 3, seed = 42)

# Assumption-light equal-weight baseline
panel_uw <- fit_auc_panel(X, y, method = "unweighted", n_folds = 3, seed = 42)
# }
```
