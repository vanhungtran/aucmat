# Screen a biomarker matrix against a binary outcome

Computes direction-preserving raw AUC, discrimination strength, and
inferential statistics for every biomarker in a numeric matrix against a
binary outcome. The package will never silently reverse a biomarker to
force its reported AUC above 0.5.

## Usage

``` r
aucmat(
  X,
  y,
  positive = NULL,
  ci = c("delong", "bootstrap", "none"),
  conf_level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  adjust = c("BH", "holm", "bonferroni", "none"),
  na_action = c("featurewise", "complete"),
  boot_n = 2000,
  seed = NULL,
  retain_data = FALSE,
  feature_metadata = NULL
)
```

## Arguments

- X:

  Numeric matrix or data.frame. Samples in rows, biomarkers in columns.
  Column names must be unique and non-empty.

- y:

  Binary outcome vector. Logical (`TRUE` = positive), numeric `0/1` (`1`
  = positive), factor, or character with exactly two observed levels.

- positive:

  Positive class label. Required when `y` is a character or factor with
  ambiguous ordering. Ignored for logical and numeric `0/1` outcomes.

- ci:

  Confidence interval method: `"delong"` (default), `"bootstrap"`, or
  `"none"`.

- conf_level:

  Confidence level in (0, 1). Default 0.95.

- alternative:

  Alternative hypothesis: `"two.sided"` (default), `"greater"`, or
  `"less"`.

- adjust:

  Multiplicity adjustment: `"BH"` (default), `"holm"`, `"bonferroni"`,
  or `"none"`.

- na_action:

  `"featurewise"` (default) uses all subjects with observed values for
  each biomarker independently. `"complete"` uses only subjects with
  complete observations across all biomarkers.

- boot_n:

  Number of bootstrap replicates when `ci = "bootstrap"`.

- seed:

  Optional integer seed for bootstrap reproducibility.

- retain_data:

  If `TRUE`, retain `X` and `y` in the result object for interactive
  plotting. Default `FALSE`.

- feature_metadata:

  Optional data.frame with metadata about biomarkers (rownames must
  match column names of `X`).

## Value

An object of class `aucmat_screen`, a list with components: `results`,
`sample_summary`, `settings`, `call`, and optionally `feature_metadata`,
`X`, `y`.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(100 * 20), nrow = 100)
colnames(X) <- paste0("bm", 1:20)
y <- rep(c(0, 1), each = 50)
fit <- aucmat(X, y)
print(fit)
#> <aucmat_screen>  20 biomarkers
#>   Outcome: 50 positive / 50 negative  (positive = pos)
#>   CI: delong  |  adjust: BH  |  na_action: featurewise
#> 
#> Top biomarkers by discrimination strength:
#>  rank biomarker auc_raw auc_strength   effect_direction    p_value   q_value
#>     1      bm20  0.6152       0.6152 higher_in_positive 0.04462509 0.4082926
#>     2       bm9  0.6144       0.6144 higher_in_positive 0.04287747 0.4082926
#>     3      bm10  0.6076       0.6076 higher_in_positive 0.06124389 0.4082926
#>     4       bm7  0.4300       0.5700  lower_in_positive 0.23079810 0.9733953
#>     5       bm2  0.5676       0.5676 higher_in_positive 0.24851964 0.9733953
#>     6       bm6  0.4468       0.5532  lower_in_positive 0.36163538 0.9733953
#>     7      bm13  0.4548       0.5452  lower_in_positive 0.44014349 0.9733953
#>     8       bm4  0.5448       0.5448 higher_in_positive 0.44492575 0.9733953
#>     9      bm15  0.4580       0.5420  lower_in_positive 0.47621631 0.9733953
#>    10       bm1  0.5408       0.5408 higher_in_positive 0.48669763 0.9733953
#>  status
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
#>      ok
head(as.data.frame(fit))
#>   biomarker auc_raw auc_strength   effect_direction n_used n_pos n_neg status
#> 1      bm20  0.6152       0.6152 higher_in_positive    100    50    50     ok
#> 2       bm9  0.6144       0.6144 higher_in_positive    100    50    50     ok
#> 3      bm10  0.6076       0.6076 higher_in_positive    100    50    50     ok
#> 4       bm7  0.4300       0.5700  lower_in_positive    100    50    50     ok
#> 5       bm2  0.5676       0.5676 higher_in_positive    100    50    50     ok
#> 6       bm6  0.4468       0.5532  lower_in_positive    100    50    50     ok
#>       prauc n_total n_missing missing_fraction  std_error  conf_low conf_high
#> 1 0.5859941     100         0                0 0.05736562 0.5027654 0.7276346
#> 2 0.5972794     100         0                0 0.05649657 0.5036687 0.7251313
#> 3 0.6349720     100         0                0 0.05748691 0.4949277 0.7202723
#> 4 0.4501436     100         0                0 0.05841582 0.3155071 0.5444929
#> 5 0.5370241     100         0                0 0.05858127 0.4527828 0.6824172
#> 6 0.4419427     100         0                0 0.05831713 0.3325005 0.5610995
#>      p_value   q_value rank warning
#> 1 0.04462509 0.4082926    1    <NA>
#> 2 0.04287747 0.4082926    2    <NA>
#> 3 0.06124389 0.4082926    3    <NA>
#> 4 0.23079810 0.9733953    4    <NA>
#> 5 0.24851964 0.9733953    5    <NA>
#> 6 0.36163538 0.9733953    6    <NA>
```
