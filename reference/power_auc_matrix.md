# Power and sample size for AUC comparisons

Computes the required sample size or achievable power for comparing AUCs
of correlated biomarkers via DeLong's test. Uses the binormal model to
translate target AUCs into mean separations and the DeLong variance
structure for correlated ROC curves.

## Usage

``` r
power_auc_matrix(
  target_aucs,
  null_auc = 0.5,
  correlation = 0,
  prevalence = 0.5,
  power = NULL,
  n = NULL,
  alpha = 0.05,
  alternative = c("two.sided", "greater", "less"),
  adjust = c("none", "BH", "bonferroni")
)
```

## Arguments

- target_aucs:

  Numeric vector of AUC values under the alternative. The null value
  (chance) is `null_auc`.

- null_auc:

  Null AUC value. Default 0.5.

- correlation:

  Between-biomarker correlation (scalar for exchangeable). Default 0.

- prevalence:

  Expected prevalence in (0, 1).

- power:

  Desired power. If `NULL` (default), compute power for given `n`. If
  supplied, compute required `n`.

- n:

  Sample size. Required when `power = NULL`. Ignored otherwise.

- alpha:

  Significance level. Default 0.05.

- alternative:

  `"two.sided"` (default), `"greater"`, or `"less"`.

- adjust:

  Multiplicity adjustment: `"none"` (default), `"BH"`, `"bonferroni"`.

## Value

A list of class `aucmat_power` with components `n`, `power`,
`target_aucs`, `effect_sizes`, `alpha`, `adjusted_alpha`.

## Examples

``` r
# Power for n=200, detecting AUC=0.70 vs chance
power_auc_matrix(target_aucs = c(0.70, 0.65), n = 200,
  prevalence = 0.3, correlation = 0.2)
#> <aucmat_power>
#>   Target AUCs: 0.700, 0.650
#>   N: 200  |  Prevalence: 0.30  |  Correlation: 0.20
#>   Alpha: 0.050  |  Adjustment: none  |  Alternative: two.sided
#> 
#>   Power: 1.000, 1.000

# Sample size for 80% power
# \donttest{
power_auc_matrix(target_aucs = 0.70, power = 0.80,
  prevalence = 0.3)
#> <aucmat_power>
#>   Target AUC: 0.700 (null = 0.500)
#>   N: 7  |  Prevalence: 0.30  |  Correlation: 0.00
#>   Alpha: 0.050  |  Adjustment: none  |  Alternative: two.sided
#> 
#>   Power: 0.800
# }
```
