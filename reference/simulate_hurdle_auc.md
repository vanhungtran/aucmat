# Simulate zero-inflated biomarker data with controlled Hurdle-AUC

Generates biomarkers where a large fraction of values are exactly zero
(e.g., scRNA-seq, microbiome, detection-limited assays). Each biomarker
is generated via a two-stage hurdle model:

## Usage

``` r
simulate_hurdle_auc(
  n,
  prevalence,
  target_hurdle_aucs,
  zero_rate_neg,
  zero_rate_pos,
  nonzero_target_aucs = NULL,
  cv_between = 0.5,
  seed = NULL
)
```

## Arguments

- n:

  Number of observations.

- prevalence:

  Proportion of positive class.

- target_hurdle_aucs:

  Target combined Hurdle-AUC values in (0, 1).

- zero_rate_neg:

  Vector of zero rates in the negative class.

- zero_rate_pos:

  Vector of zero rates in the positive class.

- nonzero_target_aucs:

  Optional vector of target AUCs on nonzero values only. If NULL,
  derived from `target_hurdle_aucs`.

- cv_between:

  Coefficient of variation for nonzero values. Default 0.5.

- seed:

  Optional seed.

## Value

An object of class `aucmat_simulation` with `data`, `target_*`,
`achieved_*`, and hurdle-specific metadata.

## Details

1.  **Zero hurdle**: P(X = 0 \| y) = logistic(beta0 + beta1 \* y)

2.  **Magnitude**: X \| X \> 0 follows shifted log-normal with
    binormal-calibrated mean separation.

Values are clipped to \[0, Inf).

## Examples

``` r
# \donttest{
set.seed(1)
sim <- simulate_hurdle_auc(
  n = 500, prevalence = 0.3,
  target_hurdle_aucs = c(0.85, 0.72, 0.55),
  zero_rate_neg = c(0.55, 0.30, 0.80),
  zero_rate_pos = c(0.25, 0.10, 0.70)
)
# Compare standard vs hurdle AUC
fit_std <- aucmat(as.matrix(sim$data[, 1:3]), sim$data$truth, ci = "none")
fit_hur <- hurdle_auc(as.matrix(sim$data[, 1:3]), sim$data$truth)
data.frame(
  biomarker = paste0("X", 1:3),
  standard_AUC = round(fit_std$results$auc_strength, 3),
  hurdle_AUC   = round(fit_hur$results$hurdle_auc, 3)
)
#>   biomarker standard_AUC hurdle_AUC
#> 1        X1        0.912      0.007
#> 2        X2        0.861      0.466
#> 3        X3        0.609      0.247
# }
```
