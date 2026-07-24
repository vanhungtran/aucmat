# Simulate correlated biomarkers with two-phase Copula-AUC generation

Generates biomarkers in two phases: (1) a class-conditional Gaussian
copula captures the full correlation structure, then (2) iterative
perturbation fine-tunes each biomarker's AUC against the outcome while
minimally disturbing correlations.

## Usage

``` r
simulate_auc_copula(
  n,
  prevalence,
  target_aucs,
  corr_matrix,
  n_iterations = 5,
  step_decay = 0.5,
  convergence_tol = 0.02,
  verify = TRUE,
  seed = NULL
)
```

## Arguments

- n:

  Number of observations.

- prevalence:

  Proportion of positive class, in (0, 1).

- target_aucs:

  Numeric vector of target AUC values, each in (0, 1).

- corr_matrix:

  A p x p target correlation matrix.

- n_iterations:

  Number of perturbation iterations. Default 5.

- step_decay:

  Rate at which perturbation step size decreases. Default 0.5.

- convergence_tol:

  Stop perturbing features where \|achieved - target\| \< tol. Default
  0.02.

- verify:

  Logical. If TRUE (default), verify achieved AUCs and correlations.

- seed:

  Optional integer seed for reproducibility.

## Value

An object of class `aucmat_simulation` with components `data`,
`target_aucs`, `achieved_aucs`, `requested_correlation`,
`achieved_correlation`, `phase1_aucs`, `n`, `prevalence`.

## Details

**Comparison to other aucmat simulators:**

- [`generate_data_analytical()`](https://vanhungtran.github.io/aucmat/reference/generate_data_analytical.md):

  Sequential binormal decomposition. Good AUC control, weak correlation
  control.

- [`simulate_auc_matrix()`](https://vanhungtran.github.io/aucmat/reference/simulate_auc_matrix.md):

  Class-conditional MVN. Good correlation control via parametrized
  structures, normal marginals only.

- [`generate_data_probit()`](https://vanhungtran.github.io/aucmat/reference/generate_data_probit.md):

  Latent probit model. Simultaneously controls AUC and correlation via
  (p+1)x(p+1) matrix.

- `simulate_auc_copula()`:

  **This function.** Gaussian copula + iterative AUC perturbation.
  Strong correlation preservation without the AUC-correlation tradeoff.

## Examples

``` r
# \donttest{
set.seed(42)
sim <- simulate_auc_copula(
  n = 500, prevalence = 0.3,
  target_aucs = c(0.85, 0.75, 0.65),
  corr_matrix = matrix(c(1.0, 0.4, 0.2, 0.4, 1.0, 0.3, 0.2, 0.3, 1.0), 3, 3)
)
data.frame(
  target   = sim$target_aucs,
  phase1   = round(sim$phase1_aucs, 3),
  achieved = round(sim$achieved_aucs, 3)
)
#>    target phase1 achieved
#> X1   0.85  0.853    0.853
#> X2   0.75  0.752    0.752
#> X3   0.65  0.651    0.651
# }
```
