# Generate correlated biomarkers with specified AUCs via latent probit model

Models the binary outcome as arising from a latent continuous variable
Y\* ~ N(0,1) with threshold tau = qnorm(1 - pi). Biomarkers (X_1, ...,
X_p, Y\*) are jointly multivariate normal with a full \$(p+1) \times
(p+1)\$ correlation matrix \$\Sigma\$.

## Usage

``` r
generate_data_probit(
  n,
  target_aucs,
  corr_matrix,
  prevalence,
  n_cal = 1e+05,
  verify = TRUE
)
```

## Arguments

- n:

  Number of observations.

- target_aucs:

  Numeric vector of target AUC values, each in (0, 1).

- corr_matrix:

  A \$p \times p\$ target correlation matrix for the biomarkers. Must be
  symmetric and positive definite.

- prevalence:

  Proportion of positive class, in (0, 1).

- n_cal:

  Calibration sample size for numerical root-finding of latent
  correlations. Larger values give more precise AUC targeting at the
  cost of calibration time. Default 100000.

- verify:

  Logical. If `TRUE` (default), empirically verify the achieved AUCs and
  correlations.

## Value

A list with components:

- data:

  A data.frame with columns `X1, ..., Xp` and `truth` (binary outcome,
  0/1).

- achieved_aucs:

  Empirical AUCs of each biomarker against the outcome.

- achieved_correlations:

  Empirical correlation matrix of the biomarkers.

- latent_rhos:

  The calibrated latent correlations used for generation.

- Sigma:

  The full \$(p+1) \times (p+1)\$ correlation matrix used.

## Details

This approach simultaneously controls **both** the between-biomarker
correlations **and** each biomarker's AUC against the outcome, subject
only to the positive-definiteness constraint on \$\Sigma\$.

**Theoretical basis**:

Under the latent probit model, a biomarker X_k and the binary outcome Y
are linked through rho_k = Cor(X_k, Y\*), the latent correlation. For a
given prevalence pi, the resulting AUC is a strictly increasing function
of rho_k.

\$\$AUC(rho_k, pi) = f(rho_k, pi)\$\$

The function f is not available in closed form (it involves the
bivariate normal CDF), so we calibrate rho_k numerically via bisection
with a large calibration sample.

**Relationship to the binormal model**:

The standard binormal model used in
[`generate_data_analytical()`](https://vanhungtran.github.io/aucmat/reference/generate_data_analytical.md)
assumes X_k \| Y=0 ~ N(mu_0, sigma^2) and X_k \| Y=1 ~ N(mu_1, sigma^2)
(equal variance). This links AUC and Cor(X_k, Y) through a single
parameter delta = (mu_1 - mu_0) / sigma.

The latent probit model **relaxes this constraint** by allowing the
correlation between biomarkers to be specified independently of the
biomarker-outcome AUC, through the full (p+1) x (p+1) matrix.

**Positive definiteness**:

The matrix \$\Sigma\$ must be positive definite. If the user-supplied
`corr_matrix` is incompatible with the calibrated latent correlations
\$\rho_k\$, the function uses the nearest positive-definite matrix (via
Higham's algorithm) and issues a warning.

## Examples

``` r
# \donttest{
# Three biomarkers with AUCs 0.9, 0.8, 0.7, moderately correlated
set.seed(42)
sim <- generate_data_probit(
  n = 500,
  target_aucs = c(0.9, 0.8, 0.7),
  corr_matrix = matrix(c(1.0, 0.3, 0.1, 0.3, 1.0, 0.2, 0.1, 0.2, 1.0), 3, 3),
  prevalence = 0.3
)

# Compare targets to achieved values
data.frame(
  Target  = sim$target_aucs,
  Achieved = round(sim$achieved_aucs, 4)
)
#>   Target Achieved
#> 1    0.9   0.9061
#> 2    0.8   0.7775
#> 3    0.7   0.7200
print(round(sim$achieved_correlations, 3))
#>       X1    X2    X3
#> X1 1.000 0.317 0.081
#> X2 0.317 1.000 0.162
#> X3 0.081 0.162 1.000
# }
```
