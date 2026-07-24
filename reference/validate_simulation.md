# Validate a `simulate_auc_matrix()` specification by repeated simulation

Repeats a
[`simulate_auc_matrix()`](https://vanhungtran.github.io/aucmat/reference/simulate_auc_matrix.md)
call `times` times with independent seeds and reports how closely the
achieved AUCs and correlations track their targets, with Monte
Carlo-aware uncertainty. A single simulated draw is never sufficient
evidence that a simulator is calibrated; this function is the check.

## Usage

``` r
validate_simulation(
  ...,
  times = 100,
  tolerance = list(auc = 0.02, correlation = 0.05),
  seed = NULL
)
```

## Arguments

- ...:

  Arguments forwarded to
  [`simulate_auc_matrix()`](https://vanhungtran.github.io/aucmat/reference/simulate_auc_matrix.md)
  on every replicate (e.g. `n`, `prevalence`, `target_aucs`,
  `correlation`, `structure`, `feasibility`). The replicate stream is
  controlled by `validate_simulation()`'s own `seed` argument, below.

- times:

  Number of independent replicates. Default 100.

- tolerance:

  Named list with `auc` and `correlation` tolerances used for the
  target-interval hit rate. Default
  `list(auc = 0.02, correlation = 0.05)`.

- seed:

  Optional integer seed controlling the replicate stream. Each replicate
  uses seed `seed + replicate_index - 1` when supplied.

## Value

An object of class `aucmat_simulation_validation`, a list with
components `auc_bias`, `auc_rmse`, `auc_mc_se`, `auc_hit_rate`,
`corr_bias`, `corr_rmse`, `corr_mc_se`, `corr_hit_rate`, `n_requested`,
`n_ok`, `n_failed`, `settings`.

## Examples

``` r
# \donttest{
val <- validate_simulation(
  n = 200, prevalence = 0.3,
  target_aucs = c(0.8, 0.7), correlation = 0.3, structure = "exchangeable",
  times = 50, seed = 1
)
print(val)
#> <aucmat_simulation_validation>  50/50 successful replicates
#> 
#> AUC calibration:
#>  biomarker    bias   rmse  mc_se hit_rate
#>         X1  0.0007 0.0333 0.0048     0.50
#>         X2 -0.0134 0.0416 0.0056     0.38
#> 
#> Correlation calibration (upper triangle):
#>   Mean bias:     -0.0081
#>   Mean RMSE:     0.065
#>   Mean hit rate: 0.6
# }
```
