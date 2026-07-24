# Hurdle-AUC for zero-inflated biomarkers

Computes a two-stage Hurdle-AUC for data where a large fraction of
values are exactly zero (e.g., scRNA-seq, microbiome, mass-spec with
detection limits). Stage 1 models zero vs non-zero via logistic
regression; Stage 2 computes standard AUC on the non-zero values only.

## Usage

``` r
hurdle_auc(X, y, positive = NULL, zero_threshold = 0)
```

## Arguments

- X:

  Numeric matrix (n x p). Zeros are treated as the hurdle.

- y:

  Binary outcome.

- positive:

  Optional positive class label.

- zero_threshold:

  Values at or below this are treated as zero. Default 0.

## Value

A list of class `aucmat_hurdle` with components: `results` (data.frame
with per-biomarker hurdle diagnostics), `stage1_models` (list of
logistic regression fits), `settings`.

## Examples

``` r
# \donttest{
set.seed(1)
sim <- simulate_hurdle_auc(n = 300, prevalence = 0.3,
  target_hurdle_aucs = c(0.85, 0.72, 0.55),
  zero_rate_neg = c(0.55, 0.30, 0.80),
  zero_rate_pos = c(0.25, 0.10, 0.70))
X <- as.matrix(sim$data[, 1:3])
y <- sim$data$truth
res <- hurdle_auc(X, y)
print(res)
#> <aucmat_hurdle>  3 biomarkers
#>   Zero threshold: 0 
#> 
#>  biomarker zero_rate_total  zero_auc nonzero_auc hurdle_auc n_nonzero
#>         X2       0.2400000 0.5841270   0.9982160  0.4003704       228
#>         X3       0.7766667 0.5468254   1.0000000  0.2122751        67
#>         X1       0.4766667 0.6341270   0.9963038  0.1390476       157
# }
```
