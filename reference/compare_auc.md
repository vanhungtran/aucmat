# Compare AUCs of selected biomarkers (paired, common-subject)

Performs DeLong-based or bootstrap paired comparisons of AUC between
biomarkers measured on the same subjects. Supports superiority,
non-inferiority, and equivalence (TOST) hypotheses.

## Usage

``` r
compare_auc(
  fit = NULL,
  X,
  y,
  reference = NULL,
  biomarkers = NULL,
  top_n = NULL,
  max_pairs = 100,
  alternative = c("two.sided", "greater", "less"),
  hypothesis = c("superiority", "noninferiority", "equivalence"),
  margin = NULL,
  conf_level = 0.95,
  method = c("delong", "bootstrap"),
  boot_n = 2000,
  seed = NULL,
  adjust = c("none", "BH", "holm", "bonferroni")
)
```

## Arguments

- fit:

  An `aucmat_screen` object (optional; stores outcome encoding).

- X:

  Numeric matrix (same as passed to
  [`aucmat()`](https://vanhungtran.github.io/aucmat/reference/aucmat.md)).

- y:

  Binary outcome.

- reference:

  Single biomarker name to compare all others against.

- biomarkers:

  Character vector of biomarker names.

- top_n:

  Take the top `top_n` biomarkers by `auc_strength` from `fit`.
  Post-selection inference is flagged as exploratory.

- max_pairs:

  Safety limit. Default 100.

- alternative:

  `"two.sided"` (default), `"greater"`, or `"less"`. `greater` means
  AUC_a \> AUC_b; `less` means AUC_a \< AUC_b.

- hypothesis:

  `"superiority"` (default), `"noninferiority"`, or `"equivalence"`.

- margin:

  Positive numeric margin for non-inferiority / equivalence.

- conf_level:

  Confidence level in (0, 1). Default 0.95.

- method:

  `"delong"` (default) or `"bootstrap"`.

- boot_n:

  Bootstrap replicates. Default 2000.

- seed:

  Optional integer seed.

- adjust:

  Multiplicity adjustment: `"BH"`, `"holm"`, `"bonferroni"`, or `"none"`
  (default).

## Value

A data.frame of class `aucmat_compare`.

## Details

The function refuses unbounded all-pairs comparisons.

## Examples

``` r
# \donttest{
set.seed(42)
sim <- simulate_auc_matrix(n=100, prevalence=0.3,
  target_aucs=c(0.85,0.75,0.65), correlation=0.3, structure="exchangeable")
X <- as.matrix(sim$data[,1:3]); y <- sim$data$truth
fit <- aucmat(X, y, ci="none")
compare_auc(fit, X, y, reference="X1")
#> <aucmat_compare>  2 pairwise comparisons
#>  biomarker_a biomarker_b     auc_a     auc_b   auc_diff  std_error    conf_low
#>           X1          X2 0.8338095 0.7847619 0.04904762 0.05359846 -0.05600343
#>           X1          X3 0.8338095 0.6166667 0.21714286 0.06799382  0.08387742
#>  conf_high     p_value q_value n_common n_pos n_neg  hypothesis margin
#>  0.1540987 0.360142363      NA      100    30    70 superiority      0
#>  0.3504083 0.001405264      NA      100    30    70 superiority      0
# }
```
