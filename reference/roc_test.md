# Single-biomarker ROC test against chance or another biomarker

A focused wrapper around DeLong's test for one or two biomarkers on
common subjects. Returns a clean summary suitable for reporting.

## Usage

``` r
roc_test(
  x,
  y,
  X = NULL,
  x2 = NULL,
  positive = NULL,
  alternative = c("two.sided", "greater", "less"),
  method = c("delong", "bootstrap"),
  boot_n = 2000,
  conf_level = 0.95
)
```

## Arguments

- x:

  Numeric biomarker vector (or matrix column name if X supplied).

- y:

  Binary outcome.

- X:

  Optional numeric matrix. If supplied, `x` and `x2` are column names or
  indices.

- x2:

  Optional second biomarker for paired comparison.

- positive:

  Positive class label.

- alternative:

  `"two.sided"` (default), `"greater"`, or `"less"`.

- method:

  `"delong"` (default) or `"bootstrap"`.

- boot_n:

  Bootstrap replicates when `method = "bootstrap"`.

- conf_level:

  Confidence level.

## Value

A list of class `aucmat_roc_test` with `auc`, `se`, `ci`, `p_value`,
`n`, `method`.

## Examples

``` r
set.seed(1)
y <- rbinom(100, 1, 0.3)
x <- rnorm(100) + y * 1.5
roc_test(x, y)
#> <aucmat_roc_test>  delong  |  two.sided
#>   AUC = 0.8755, SE = 0.0346
#>   95% CI: [0.8076, 0.9433]  p = 0.0000
#>   n = 100 (32 + / 68 -)
```
