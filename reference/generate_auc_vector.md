# Generate a continuous score vector for a target AUC

Construct a numeric score vector for a binary outcome so that the
empirical ROC AUC matches a requested value exactly when that value is
attainable for the observed numbers of positive and negative samples.
When the requested AUC is not attainable on the finite-sample grid, the
function can either use the nearest attainable value or stop with an
error.

## Usage

``` r
generate_auc_vector(
  y,
  target_auc,
  positive = NULL,
  unattainable = c("nearest", "error"),
  seed = NULL,
  shuffle_within_class = TRUE
)
```

## Arguments

- y:

  Binary outcome vector. Accepted inputs are logical, numeric 0/1,
  factor, or character with exactly two non-missing classes.

- target_auc:

  Target empirical AUC in `[0, 1]`.

- positive:

  Optional positive class label when `y` is a factor or character
  vector. By default, the second factor level is treated as the positive
  class.

- unattainable:

  One of `"nearest"` or `"error"`. If the requested AUC is not exactly
  attainable for the current numbers of cases and controls, `"nearest"`
  uses the closest attainable AUC and warns, while `"error"` stops.

- seed:

  Optional seed used only when shuffling scores within the positive and
  negative classes.

- shuffle_within_class:

  Logical. If `TRUE`, randomly permute the assigned scores within each
  class while preserving the achieved AUC. Default `TRUE`.

## Value

A list with components:

- x:

  Numeric score vector with the same length as `y`. Missing values in
  `y` are propagated as `NA` in `x`.

- requested_auc:

  The requested AUC.

- achieved_auc:

  The empirical AUC achieved by `x` on the non-missing observations.

- auc_step:

  The smallest attainable AUC increment, `1 / (n_pos * n_neg)`.

- n_pos:

  Number of positive samples used.

- n_neg:

  Number of negative samples used.

## Examples

``` r
set.seed(1)
y <- rbinom(1000, size = 1, prob = 0.3)
out <- generate_auc_vector(y, target_auc = 0.8)
#> Warning: Requested AUC 0.8000000000 is not exactly attainable with n_pos = 304 and n_neg = 696. Using nearest attainable AUC 0.7999990547 instead.
head(out$x)
#> [1] -1.45020988  0.39749835 -0.04388007  1.35631175 -0.07401305  0.74379584
out$achieved_auc
#> [1] 0.7999991
```
