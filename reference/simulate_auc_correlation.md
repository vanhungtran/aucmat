# Simulate repeated X vectors for a fixed y and record AUC and correlation

For a fixed binary outcome vector `y`, repeatedly simulate continuous
predictors `X` from two normal distributions with equal variance and a
mean separation chosen so that the theoretical AUC equals `target_auc`.
For each simulation, the function calculates the empirical AUC and the
correlation between `y` and `X`.

## Usage

``` r
simulate_auc_correlation(
  y,
  target_auc,
  n_sim = 10000,
  positive = NULL,
  seed = NULL,
  cor_method = c("pearson", "spearman", "kendall"),
  keep_x = FALSE
)
```

## Arguments

- y:

  Binary outcome vector. Accepted inputs are logical, numeric 0/1,
  factor, or character with exactly two non-missing classes.

- target_auc:

  One target AUC or a numeric vector of target AUC values in `[0, 1]`.

- n_sim:

  Number of simulations. Default 10000.

- positive:

  Optional positive class label when `y` is a factor or character
  vector. By default, the second factor level is treated as the positive
  class.

- seed:

  Optional random seed.

- cor_method:

  Correlation method passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). One of
  `"pearson"`, `"spearman"`, or `"kendall"`. Default `"pearson"`.

- keep_x:

  Logical. If `TRUE`, also return the simulated `X` values as a matrix
  with one column per simulation. Default `FALSE`.

## Value

A list with components:

- `results` — A data.frame with one row per simulation per requested
  target and columns `target_auc`, `sim`, `auc`, and `correlation`.

- `target_auc` — The requested target AUC value(s).

- `mean_shift` — A data.frame with columns `target_auc` and
  `mean_shift`. For `target_auc = 1`, `mean_shift` is reported as `Inf`
  because perfect separation is handled as a special case rather than a
  finite normal shift.

- `n_sim` — Number of simulations performed per target AUC.

- `x_matrix` — Optional matrix of simulated predictors if
  `keep_x = TRUE` and `length(target_auc) == 1`.

- `x_matrix_list` — Optional named list of simulated predictor matrices
  if `keep_x = TRUE` and multiple target AUC values are supplied.

## Examples

``` r
set.seed(1)
y <- rbinom(1000, size = 1, prob = 0.3)
sim <- simulate_auc_correlation(y, target_auc = c(0.8, 0.9, 1), n_sim = 100)
head(sim$results)
#>   target_auc sim       auc correlation
#> 1        0.8   1 0.8107607   0.4961646
#> 2        0.8   2 0.8088182   0.4879859
#> 3        0.8   3 0.7929144   0.4710591
#> 4        0.8   4 0.8149577   0.4988884
#> 5        0.8   5 0.8116776   0.4981350
#> 6        0.8   6 0.7895729   0.4689970
aggregate(auc ~ target_auc, data = sim$results, mean)
#>   target_auc       auc
#> 1        0.8 0.7999534
#> 2        0.9 0.8990580
#> 3        1.0 1.0000000
```
