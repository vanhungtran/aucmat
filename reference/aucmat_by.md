# Stratified biomarker screening by group

Runs
[`aucmat()`](https://vanhungtran.github.io/aucmat/reference/aucmat.md)
within each level of a grouping variable and returns per-group results
plus a heterogeneity summary (I-squared, Cochran's Q). Useful for
assessing whether biomarker discrimination varies by batch, site, or
clinical stratum.

## Usage

``` r
aucmat_by(
  X,
  y,
  group,
  positive = NULL,
  ci = c("delong", "bootstrap", "none"),
  conf_level = 0.95,
  adjust = "BH",
  boot_n = 2000,
  min_group_size = 10L
)
```

## Arguments

- X:

  Numeric matrix.

- y:

  Binary outcome.

- group:

  Factor or vector coercible to factor. Screening is performed within
  each level.

- positive:

  Optional positive class label.

- ci, conf_level, adjust, boot_n:

  Passed to
  [`aucmat()`](https://vanhungtran.github.io/aucmat/reference/aucmat.md).

- min_group_size:

  Minimum number of observations per group (both classes required).
  Default 10.

## Value

A list of class `aucmat_by` with `results` (per-group data.frame),
`heterogeneity` (I-squared, Q, p per biomarker), `groups`, `settings`.

## Examples

``` r
# \donttest{
set.seed(1)
sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
  target_aucs = c(0.8, 0.7, 0.6), correlation = 0.3,
  structure = "exchangeable")
grp <- sample(c("SiteA", "SiteB"), 300, replace = TRUE)
res <- aucmat_by(as.matrix(sim$data[, 1:3]), sim$data$truth, grp,
  ci = "none")
print(res)
#> <aucmat_by>  2 groups: SiteA, SiteB
#> 
#> Heterogeneity summary:
#>  biomarker i_squared cochran_q q_p_value n_groups
#>         X1       NaN       NaN       NaN        2
#>         X2       NaN       NaN       NaN        2
#>         X3       NaN       NaN       NaN        2
# }
```
