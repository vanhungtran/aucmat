# Screen biomarkers against a multiclass outcome

Computes one-versus-rest AUC for each biomarker-class pair, plus
Hand-Till multiclass AUC per biomarker. Hand-Till is the average AUC of
all pairwise class comparisons weighted by class prevalence.

## Usage

``` r
aucmat_multiclass(
  X,
  y,
  ci = c("none", "bootstrap"),
  conf_level = 0.95,
  boot_n = 2000,
  seed = NULL,
  adjust = c("none", "BH", "holm", "bonferroni")
)
```

## Arguments

- X:

  Numeric matrix (n x p).

- y:

  Outcome with 3+ classes. Factor, character, or integer.

- ci:

  `"none"` (default) or `"bootstrap"`. DeLong is not defined for
  multiclass.

- conf_level:

  Confidence level.

- boot_n:

  Bootstrap replicates when `ci = "bootstrap"`.

- seed:

  Optional seed.

- adjust:

  Multiplicity adjustment across all biomarker-class pairs.

## Value

A list of class `aucmat_multiclass` with components `ovr` (one-vs-rest),
`hand_till` (one value per biomarker), `classes`, `settings`.

## Examples

``` r
# \donttest{
set.seed(1)
# 3-class outcome with 2 biomarkers
n <- 300
y <- sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
X <- cbind(
  X1 = rnorm(n) + ifelse(y == "A", 1.5, ifelse(y == "B", 0.5, -1)),
  X2 = rnorm(n) + ifelse(y == "C", 2.0, 0)
)
colnames(X) <- c("bm1", "bm2")
res <- aucmat_multiclass(X, y, ci = "none")
print(res)
#> <aucmat_multiclass>  3 classes (A, B, C)
#>   2 biomarkers
#> 
#> Hand-Till multiclass AUC:
#>  biomarker auc_hand_till
#>        bm1     0.8564228
#>        bm2     0.2809321
#> 
#> Top one-vs-rest pairs:
#>  biomarker class    auc_ovr n_class
#>        bm2     C 0.87543703      49
#>        bm1     A 0.84334821     160
#>        bm2     B 0.41495347      91
#>        bm2     A 0.36607143     160
#>        bm1     B 0.36358378      91
#>        bm1     C 0.08561672      49
# }
```
