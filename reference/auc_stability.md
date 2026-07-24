# Bootstrap rank-stability analysis for biomarker screening

Resamples positive and negative subjects separately with replacement,
recomputes the complete feature screen for each replicate, and
aggregates rank distributions to assess the stability of leading
biomarkers.

## Usage

``` r
auc_stability(
  X,
  y,
  positive = NULL,
  times = 1000,
  top_k = c(10, 25, 50),
  seed = NULL,
  max_pairs_for_coselection = 20
)
```

## Arguments

- X:

  Numeric matrix (samples in rows, biomarkers in columns).

- y:

  Binary outcome vector.

- positive:

  Optional positive class label.

- times:

  Number of bootstrap replicates. Default 1000.

- top_k:

  Integer vector of set sizes for top-k probability reporting. Default
  `c(10, 25, 50)`.

- seed:

  Integer seed for reproducibility.

- max_pairs_for_coselection:

  Maximum number of top biomarkers for which pairwise co-selection
  frequencies are reported. Default 20.

## Value

An object of class `aucmat_stability`, a list with components:
`rank_summary`, `top_k_probs`, `coselection`, `settings`.

## Examples

``` r
# \donttest{
set.seed(42)
X <- matrix(rnorm(200*10), 200, 10, dimnames=list(NULL, paste0("bm",1:10)))
y <- rep(c(0,1), each=100)
stab <- auc_stability(X, y, times=50, seed=1)
print(stab)
#> <aucmat_stability>  50/50 successful replicates
#> Top biomarkers by median rank:
#>  biomarker rank_median rank_q25 rank_q75 auc_mean     auc_sd top1_freq
#>        bm9           1     1.00        2 0.613900 0.03687265      0.52
#>       bm10           3     2.00        4 0.585694 0.04279510      0.22
#>        bm7           4     3.00        5 0.577786 0.03752179      0.10
#>        bm1           5     3.00        7 0.560710 0.03636033      0.10
#>        bm3           6     4.00        8 0.548664 0.03205187      0.00
#>        bm6           6     4.25        8 0.548626 0.03121019      0.02
#>        bm2           7     6.00        9 0.534052 0.02595008      0.00
#>        bm4           7     5.00        9 0.545850 0.03650601      0.04
#>        bm5           7     5.25        9 0.538390 0.02971093      0.00
#>        bm8           8     6.00        9 0.535682 0.02734570      0.00
# }
```
