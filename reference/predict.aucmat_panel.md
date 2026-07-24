# Predict panel scores for new subjects

Applies a fitted panel's coefficients to new data, using the
centering/scaling stored from the original (training) fit.

## Usage

``` r
# S3 method for class 'aucmat_panel'
predict(object, newdata, ...)
```

## Arguments

- object:

  An `aucmat_panel` object from
  [`fit_auc_panel()`](https://vanhungtran.github.io/aucmat/reference/fit_auc_panel.md).

- newdata:

  Numeric matrix or data.frame containing at least the columns used to
  fit `object`.

- ...:

  Ignored.

## Value

A numeric vector of panel scores, one per row of `newdata`.

## Examples

``` r
# \donttest{
set.seed(1)
sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
  target_aucs = c(0.85, 0.75, 0.65, 0.60),
  correlation = 0.3, structure = "exchangeable")
X <- as.matrix(sim$data[, 1:4]); y <- sim$data$truth
panel <- fit_auc_panel(X[1:200, ], y[1:200], method = "ridge", n_folds = 0)
predict(panel, X[201:300, ])
#>   [1] -0.60426456 -1.87740605 -1.39787268  1.22442479  2.71114341  2.38053378
#>   [7]  3.74585471 -0.11773248 -3.46307358 -2.00021955 -0.17967577  0.94750887
#>  [13]  0.81581884  0.14071727  0.31019425 -1.25698687 -1.04132918 -1.60707686
#>  [19] -2.06911115  0.26349292  0.92721392  1.28804134  0.09957508 -0.08573197
#>  [25]  0.46865265  1.18154275  0.88970867  0.35698883  0.67505816 -1.43316713
#>  [31] -2.70860639 -1.37808841  2.68700911 -0.68560446 -2.65626096 -1.63490571
#>  [37] -1.09128772  3.75840074  1.04756562 -0.73729305  1.28332030 -1.47652118
#>  [43]  1.87659723  1.92243451  0.22364147 -0.89436562  3.58783771  0.90837668
#>  [49]  0.15547158 -0.71600932  0.56782378  2.91050321 -1.90076025  1.20604018
#>  [55]  1.70710660 -0.54957100 -0.09661353  0.71259946  1.71239902  2.53567697
#>  [61] -2.04880276 -3.03085241 -1.85518483 -0.69129303 -2.10944040  1.84174881
#>  [67]  1.16427075 -1.89474836  0.62729704 -2.04860654  3.45307342 -2.63676270
#>  [73] -0.57217526 -2.36761808  0.01856380  2.79571325  0.47642707 -1.65243619
#>  [79]  1.43926728  1.10460744 -0.99782617  2.51492647 -0.58538846 -0.64017403
#>  [85] -0.83338684  0.82663839  0.42650361  1.42654147 -0.53745679 -2.14630260
#>  [91]  0.17995935 -3.64781856 -0.49079664  0.94186731  1.33321407 -0.20511955
#>  [97]  2.93587614  3.87789879  1.69951730 -0.59227416
# }
```
