# Generate a continuous score vector for target AUC and correlation

Construct a numeric score vector for a binary outcome with an empirical
AUC matched exactly when attainable, then tune a strictly increasing
Box-Cox transform so the Pearson correlation between `y` and `x` matches
a requested value as closely as possible without changing the AUC.

## Usage

``` r
generate_auc_cor_vector(
  y,
  target_auc,
  target_cor,
  positive = NULL,
  auc_unattainable = c("nearest", "error"),
  cor_unattainable = c("nearest", "error"),
  seed = NULL,
  shuffle_within_class = TRUE,
  lambda_bounds = c(-8, 8),
  lambda_grid_n = 401,
  cor_tol = 1e-06
)
```

## Arguments

- y:

  Binary outcome vector. Accepted inputs are logical, numeric 0/1,
  factor, or character with exactly two non-missing classes.

- target_auc:

  Target empirical AUC in `[0, 1]`.

- target_cor:

  Target Pearson correlation between `y` and the generated score vector.

- positive:

  Optional positive class label when `y` is a factor or character
  vector. By default, the second factor level is treated as the positive
  class.

- auc_unattainable:

  One of `"nearest"` or `"error"`. Passed to
  [`generate_auc_vector()`](https://vanhungtran.github.io/aucmat/reference/generate_auc_vector.md)
  when the requested AUC is not exactly attainable on the finite-sample
  grid.

- cor_unattainable:

  One of `"nearest"` or `"error"`. If the requested correlation is
  outside the achievable range of the monotone transform family,
  `"nearest"` returns the closest match and warns, while `"error"`
  stops.

- seed:

  Optional seed used only when constructing the base AUC-matched score
  vector.

- shuffle_within_class:

  Logical. Passed to
  [`generate_auc_vector()`](https://vanhungtran.github.io/aucmat/reference/generate_auc_vector.md).

- lambda_bounds:

  Numeric length-2 vector giving the search interval for the Box-Cox
  transform parameter. Default `c(-8, 8)`.

- lambda_grid_n:

  Number of grid points used to bracket the correlation search before
  local refinement. Default `401`.

- cor_tol:

  Numeric tolerance used to decide whether the requested correlation was
  achieved. Default `1e-6`.

## Value

A list with components:

- x:

  Numeric score vector with the same length as `y`. Missing values in
  `y` are propagated as `NA` in `x`.

- requested_auc:

  The requested AUC.

- requested_cor:

  The requested correlation.

- achieved_auc:

  The empirical AUC achieved by `x`.

- achieved_cor:

  The Pearson correlation achieved by `x`.

- lambda:

  The fitted Box-Cox transform parameter.

- cor_range:

  Length-2 numeric vector with the approximate minimum and maximum
  correlations achievable over `lambda_bounds`.

- base_vector:

  The object returned by
  [`generate_auc_vector()`](https://vanhungtran.github.io/aucmat/reference/generate_auc_vector.md)
  before the monotone transform was applied.

## Details

Theory:

- The empirical AUC depends only on the rank ordering of the score
  vector.

- Pearson correlation depends on the numeric spacing of the score
  values.

- A strictly increasing transform preserves ranks, so it preserves AUC.

- Therefore the method first constructs a base vector with the requested
  AUC, then tunes a monotone Box-Cox transform to move the Pearson
  correlation toward the requested target without changing the AUC.

- Not every `(AUC, correlation)` pair is achievable for a fixed `y`,
  because the monotone transform can change spacing but cannot change
  ordering.

## Examples

``` r
set.seed(1)
y <- rbinom(200, size = 1, prob = 0.3)
out <- generate_auc_cor_vector(y, target_auc = 0.8, target_cor = 0.45)
#> Warning: Requested AUC 0.8000000000 is not exactly attainable with n_pos = 62 and n_neg = 138. Using nearest attainable AUC 0.8000233754 instead.
#> Warning: Requested correlation 0.4500000000 was not achieved exactly. Using nearest value 0.4500019938 instead.
c(out$achieved_auc, out$achieved_cor)
#> [1] 0.8000234 0.4500020
```
