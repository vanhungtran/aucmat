# Methods: Generating Biomarker Data with Specified AUC

## 1. `generate_data_analytical()` — Correlated Multivariate Biomarkers

### 1.1 Problem

Given a binary outcome Y \in \\0, 1\\ with prevalence \pi = P(Y = 1),
generate p continuous biomarkers X_1, \ldots, X_p such that:

- Each X_k has a specified empirical AUC \theta_k against Y
- The biomarkers have a specified p \times p correlation matrix
  \mathbf{R}

### 1.2 Theoretical Foundation

Under the **binormal model**, if X_k \mid Y = 0 \sim
\mathcal{N}(\mu\_{k0}, \sigma^2) and X_k \mid Y = 1 \sim
\mathcal{N}(\mu\_{k1}, \sigma^2) (equal variance), the AUC is:

\text{AUC}\_k = \Phi\\\left(\frac{\mu\_{k1} -
\mu\_{k0}}{\sigma\sqrt{2}}\right)

where \Phi is the standard normal CDF. Equivalently, to achieve AUC
\theta_k, the required **standardized mean separation** is:

\delta_k = \frac{\mu\_{k1} - \mu\_{k0}}{\sigma} = \sqrt{2} \cdot
\Phi^{-1}(\theta_k)

### 1.3 Algorithm

The algorithm proceeds **sequentially from highest to lowest AUC**,
decomposing each new biomarker into a **correlated component**
(inherited from previously generated biomarkers) and a **residual
component** (carrying the remaining signal).

#### Step 1: Generate the outcome

Y_i \sim \text{Bernoulli}(\pi), \quad i = 1, \ldots, n

Let n_1 = \sum Y_i (positives) and n_0 = \sum (1 - Y_i) (negatives).

#### Step 2: Generate the first biomarker X_1 (highest AUC)

1.  Compute the required total separation: \delta_1 = \sqrt{2} \cdot
    \Phi^{-1}(\theta_1)

2.  Generate X_1 \sim \mathcal{N}(0, 1) and shift each class to achieve
    the target means: \mu_1^{(0)} = \delta_1 \cdot
    \left(-\frac{n_1}{n}\right), \quad \mu_1^{(1)} = \delta_1 \cdot
    \left(\frac{n_0}{n}\right)

3.  Center each class at its target mean, then standardize: X\_{1,i}
    \gets X\_{1,i} - \bar{X}\_1^{(Y_i)} + \mu_1^{(Y_i)} X_1 \gets
    \frac{X_1 - \bar{X}\_1}{\text{sd}(X_1)}

This ensures \mathbb{E}\[X_1 \mid Y = 1\] - \mathbb{E}\[X_1 \mid Y = 0\]
= \delta_1 and \text{Var}(X_1) = 1, satisfying the binormal AUC
equation.

#### Step 3: Sequentially generate X_k for k = 2, \ldots, p

For each subsequent biomarker, we decompose it into two orthogonal
parts:

X_k = X_k^{\text{(corr)}} + X_k^{\text{(res)}}

##### 3a. Correlated component

From the k-1 already-generated biomarkers \mathbf{X}\_{1:k-1}, we
construct the part of X_k that achieves the target correlations:

\boldsymbol{\beta} = \mathbf{R}\_{1:k-1, 1:k-1}^{-1} \cdot \mathbf{r}\_k

where \mathbf{r}\_k = (\rho\_{k,1}, \ldots, \rho\_{k,k-1})^\top is the
vector of target correlations between X_k and the previous biomarkers.

X_k^{\text{(corr)}} = \tilde{\mathbf{X}}\_{1:k-1} \boldsymbol{\beta}

where \tilde{\mathbf{X}} is the column-standardized previous data. This
component carries the correlation structure but has **its own
separation** between the two outcome classes:

\Delta^{\text{(corr)}} = \bar{X}\_k^{\text{(corr)}} \mid\_{Y=1} -
\bar{X}\_k^{\text{(corr)}} \mid\_{Y=0}

##### 3b. Residual component

The residual variance after accounting for the correlated part is:

\sigma^2\_{\text{res}} = 1 - \boldsymbol{\beta}^\top \mathbf{R}\_{1:k-1,
1:k-1} \boldsymbol{\beta}

(clamped to \geq 0 to handle floating-point errors).

The **remaining separation** the residual must provide to hit the target
AUC \theta_k:

\delta_k = \sqrt{2} \cdot \Phi^{-1}(\theta_k) \delta\_{\text{res}} =
\delta_k - \Delta^{\text{(corr)}}

We generate X_k^{\text{(res)}} \sim \mathcal{N}(0, \sigma\_{\text{res}})
and shift the classes:

\mu\_{\text{res}}^{(0)} = \delta\_{\text{res}} \cdot
\left(-\frac{n_1}{n}\right), \quad \mu\_{\text{res}}^{(1)} =
\delta\_{\text{res}} \cdot \left(\frac{n_0}{n}\right)

Each observation’s residual is centered at the class-appropriate mean.

##### 3c. Combine and standardize

X_k = X_k^{\text{(corr)}} + X_k^{\text{(res)}} X_k \gets \frac{X_k -
\bar{X}\_k}{\text{sd}(X_k)}

#### Step 4: Reorder and verify

Columns are returned to the user’s original AUC ordering. The achieved
AUCs and correlation matrix are empirically verified using
[`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html) and
[`cor()`](https://rdrr.io/r/stats/cor.html).

### 1.4 Key Properties

- **Sequential from highest AUC**: Ensures the strongest signal is
  established first; weaker signals inherit residual from stronger ones
- **Exact class means**: By explicitly setting per-class empirical
  means, the achieved AUC closely matches the target
- **Standardization at each step**: Each biomarker has unit variance, so
  the covariance matrix equals the correlation matrix
- **Floating-point guard**: Residual variance is clamped to \geq 0 to
  handle numerical edge cases where the correlation structure nearly
  fully determines the separation

------------------------------------------------------------------------

## 2. `generate_auc_vector()` — Single Vector with Exact Empirical AUC

### 2.1 Problem

Given a binary outcome y with n_1 positives and n_0 negatives, construct
a continuous score vector x such that the **empirical Mann-Whitney AUC**
exactly equals a specified value \theta.

### 2.2 Theoretical Foundation

The empirical AUC is a **rank-based statistic**. On a finite sample with
n_1 positives and n_0 negatives, the achievable AUC values form a
discrete grid with step size:

\Delta\_{\text{AUC}} = \frac{1}{n_1 \cdot n_0}

An AUC of \theta corresponds to a **target U-statistic** of:

U\_{\text{target}} = \theta \cdot n_1 \cdot n_0

The U-statistic is related to the sum of ranks of positive observations:

U = \sum\_{i \in \text{pos}} R_i - \frac{n_1(n_1 + 1)}{2}

where R_i is the rank of observation i among all n = n_1 + n_0
observations. Therefore:

\sum\_{i \in \text{pos}} R_i = U\_{\text{target}} + \frac{n_1(n_1 +
1)}{2}

### 2.3 Algorithm

1.  **Validate attainability**: Check that \theta \cdot n_1 n_0 is
    achievable (integer). If `unattainable = "nearest"`, round to the
    nearest attainable value and warn. If `"error"`, stop.

2.  **Determine rank subset**: Find a subset of n_1 ranks from \\1,
    \ldots, n\\ whose sum equals the target: S\_{\text{target}} =
    U\_{\text{target}} + \frac{n_1(n_1 + 1)}{2}

    Start with the minimal ranks \\1, 2, \ldots, n_1\\ (sum =
    n_1(n_1+1)/2, AUC = 0) and greedily increment from the highest rank
    downward until the target sum is reached.

3.  **Assign scores from the normal quantile function**:
    \text{scores}\_i = \Phi^{-1}\\\left(\frac{R_i - 0.5}{n}\right) where
    R_i is the rank assigned to observation i.

4.  **Shuffle within classes** (if `shuffle_within_class = TRUE`):
    Randomly permute scores within each outcome class, preserving the
    rank ordering between classes (and therefore the AUC).

### 2.4 Key Properties

- **Exact empirical AUC**: On the finite sample, the achieved AUC equals
  the target (to within the grid resolution)
- **Rank-determined**: The AUC depends only on the rank ordering, not on
  the numeric spacing of scores
- **Normal quantile scores**: Produces normally-distributed scores that
  are interpretable as continuous biomarkers
- **NA propagation**: Missing values in y are preserved as NA in x

------------------------------------------------------------------------

## 3. `generate_auc_cor_vector()` — AUC + Pearson Correlation

### 3.1 Problem

Construct a score vector that simultaneously achieves a target empirical
AUC \theta and a target Pearson correlation \rho with the binary
outcome.

### 3.2 Theoretical Foundation

The **key insight** is that:

1.  The empirical AUC depends **only on the rank ordering** of scores —
    monotone transformations preserve AUC
2.  The Pearson correlation depends on the **numeric spacing** of scores
    — monotone transformations change correlation

Therefore the method proceeds in two stages: **rank construction**
(fixes AUC) followed by **monotone transform tuning** (adjusts
correlation without changing AUC).

### 3.3 Algorithm

#### Stage 1: AUC-achieving base vector

Call `generate_auc_vector(y, target_auc)` to obtain x\_{\text{base}}
with the exact target AUC.

#### Stage 2: Box-Cox transform tuning

Apply a strictly increasing Box-Cox transformation to x\_{\text{base}}:

x\_{\text{transformed}} = \begin{cases}
\dfrac{(x\_{\text{shifted}})^\lambda - 1}{\lambda}, & \lambda \neq 0
\\\[8pt\] \ln(x\_{\text{shifted}}), & \lambda = 0 \end{cases}

where x\_{\text{shifted}} = x\_{\text{base}} - \min(x\_{\text{base}}) +
1 ensures strictly positive values.

The parameter \lambda is chosen by **grid search followed by local
optimization**:

1.  Evaluate \text{cor}(y, x\_\lambda) on a grid of \lambda \in
    \[\lambda\_{\min}, \lambda\_{\max}\] (default: \[-8, 8\], 401
    points)
2.  Find the grid point closest to the target correlation
3.  Use [`optimize()`](https://rdrr.io/r/stats/optimize.html) to refine
    \lambda by minimizing (\text{cor}(y, x\_\lambda) -
    \rho\_{\text{target}})^2

#### Handling infeasible targets

Not every (\theta, \rho) pair is achievable. The monotone transform can
only move the correlation within a bounded range \[\rho\_{\min},
\rho\_{\max}\]. If the target falls outside this range with
`cor_unattainable = "nearest"`, the closest achievable correlation is
returned with a warning.

### 3.4 Key Properties

- **AUC is exactly preserved**: The Box-Cox transform is strictly
  increasing, so ranks are unchanged
- **Correlation is approximately matched**: Tuned to within numerical
  tolerance
- **Achievable range**: Depends on n_1, n_0, \theta, and the specific
  rank configuration

------------------------------------------------------------------------

## 4. `simulate_auc_correlation()` — Monte Carlo AUC-Correlation Sampling

### 4.1 Purpose

For a **fixed** binary outcome y, simulate many independent score
vectors under the **normal homoscedastic binormal model** and record the
empirical AUC and correlation for each. This characterizes the joint
sampling distribution of (AUC, correlation) for a given prevalence and
target AUC.

### 4.2 Algorithm

For each simulation replicate and each target AUC \theta:

1.  If \theta \< 1 (imperfect separation): X \mid Y = 0 \sim
    \mathcal{N}(0, 1) X \mid Y = 1 \sim \mathcal{N}(\delta, 1) where
    \delta = \sqrt{2} \cdot \Phi^{-1}(\theta)

2.  If \theta = 1 (perfect separation):

    - Generate negatives: X \mid Y = 0 \sim \mathcal{N}(0, 1)
    - Generate positives: X \mid Y = 1 \sim \mathcal{N}(0, 1)
    - Shift all positive values above the maximum negative value

3.  Record the empirical AUC and \text{cor}(y, X)

### 4.3 Key Properties

- **Sampling distribution characterization**: Reveals the variability of
  empirical AUC and correlation for a given sample size and prevalence
- **Perfect separation handling**: When \theta = 1, uses a deterministic
  gap rather than infinite mean shift
- **Multiple targets**: Accepts a vector of target AUCs, returning
  per-target results

------------------------------------------------------------------------

## Reference

- Hanley JA, McNeil BJ (1982). The meaning and use of the area under a
  receiver operating characteristic (ROC) curve. *Radiology*, 143(1),
  29–36.
- DeLong ER, DeLong DM, Clarke-Pearson DL (1988). Comparing the areas
  under two or more correlated receiver operating characteristic curves.
  *Biometrics*, 44(3), 837–845.
