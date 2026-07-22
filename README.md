# aucmat

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/vanhungtran/aucmat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vanhungtran/aucmat/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Scalable and statistically principled ROC analysis for high-dimensional biomarker matrices in R.**

`aucmat` is a matrix-first R package for univariate screening and ranking of
molecular biomarkers (genes, proteins, metabolites) measured on the same
subjects. It reports **direction-preserving AUCs** with DeLong and
stratified-bootstrap inference, multiplicity adjustment, paired comparisons
(superiority, non-inferiority, equivalence), omnibus global testing, bootstrap
rank-stability analysis, and a multivariate simulation engine with validation.

## Installation

```r
remotes::install_github("vanhungtran/aucmat")
```

## Function Reference

| Category | Function | Description |
|----------|----------|-------------|
| **Screening** | `aucmat()` | Screen all biomarkers against a binary outcome |
| **Inference** | `compare_auc()` | Paired AUC comparisons (superiority, NI, equivalence) |
| | `compare_auc_global()` | Omnibus Wald test: H0 all AUCs equal |
| | `auc_stability()` | Bootstrap rank-stability analysis |
| **Simulation** | `simulate_auc_matrix()` | Class-conditional MVN with named correlation structures |
| | `validate_simulation()` | Repeated-draw calibration check |
| | `generate_data_probit()` | Latent probit: independent AUC + correlation |
| | `generate_data_analytical()` | Sequential binormal: fast, AUC-correlation linked |
| | `generate_auc_vector()` | Single score with exact empirical AUC |
| | `generate_auc_cor_vector()` | Single score with exact AUC + approximate Pearson r |
| | `simulate_auc_correlation()` | Monte Carlo sampling distribution of (AUC, r) |
| **Visualization** | `plot_auc_rank()` | Ordered discrimination strengths with CIs |
| | `plot_auc_volcano()` | Effect magnitude vs statistical evidence |
| | `plot_auc_forest()` | Selected AUCs with confidence intervals |
| | `plot_auc_stability()` | Bootstrap rank distributions |
| | `plot_roc_top()` | ROC curves for selected biomarkers |
| **S3 Methods** | `print()`, `summary()` | Display and summarise results |
| | `as.data.frame()` | Convert screening results to data.frame |
| | `plot()` | Default plot (rank plot) |
| | `subset()` | Filter screening results |

## Quick Start

```r
library(aucmat)

# 1. Generate correlated biomarkers with controlled AUCs
sim <- generate_data_probit(
  n = 500, prevalence = 0.3,
  target_aucs = c(0.90, 0.80, 0.70),
  corr_matrix  = matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), 3, 3)
)

# 2. Screen all biomarkers
fit <- aucmat(sim$data[, 1:3], sim$data$truth)
print(fit)

# 3. Visualize
plot_auc_rank(fit)       # ordered discrimination strengths
plot_auc_volcano(fit)    # effect vs evidence
plot_auc_forest(fit)     # multiple AUCs with CIs
plot_roc_top(fit, X = sim$data[, 1:3], y = sim$data$truth)

# 4. Compare AUCs — three comparison modes
compare_auc(fit, sim$data[, 1:3], sim$data$truth, top_n = 3)          # all pairs
compare_auc(fit, sim$data[, 1:3], sim$data$truth, reference = "X1")   # vs reference
compare_auc(fit, sim$data[, 1:3], sim$data$truth,
  biomarkers = c("X1", "X2"))                                         # named set

# 5. Global test — are all AUCs equal?
compare_auc_global(fit, sim$data[, 1:3], sim$data$truth)

# 6. Assess rank stability
stab <- auc_stability(sim$data[, 1:3], sim$data$truth, times = 500, seed = 42)
plot_auc_stability(stab)
```

---

## 1. Screening: `aucmat()`

`aucmat()` computes direction-preserving AUCs for every biomarker column
against a binary outcome. It supports DeLong or bootstrap inference,
multiplicity adjustment, and feature-wise or complete-case missing-data
handling.

```r
fit <- aucmat(X, y)                                       # default: DeLong, BH
fit <- aucmat(X, y, ci = "bootstrap", boot_n = 2000)      # bootstrap CIs
fit <- aucmat(X, y, ci = "none")                          # no inference (fastest)
fit <- aucmat(X, y, adjust = "bonferroni")                # stricter adjustment
fit <- aucmat(X, y, na_action = "complete")               # complete-case only
fit <- aucmat(X, y, retain_data = TRUE)                   # keep data for plotting
```

### Result Columns

| Column | Description |
|--------|-------------|
| `biomarker` | Column name |
| `auc_raw` | AUC on observed direction (may be < 0.5) |
| `auc_strength` | 0.5 + \|auc_raw − 0.5\| — always ≥ 0.5 |
| `effect_direction` | `higher_in_positive` or `lower_in_positive` |
| `std_error` | DeLong or bootstrap standard error |
| `conf_low`, `conf_high` | Confidence interval |
| `p_value`, `q_value` | Raw and multiplicity-adjusted p-values |
| `rank` | Ordered by `auc_strength` descending |
| `n_used`, `n_pos`, `n_neg` | Sample counts per biomarker |
| `n_missing`, `missing_fraction` | Missing data diagnostics |
| `status` | `ok`, `constant`, `insufficient_positive`, etc. |
| `warning` | `small_class_counts`, `high_missingness`, etc. |

### S3 Methods

```r
print(fit)               # compact top-N display
summary(fit)             # class counts, statuses, missingness, multiplicity
as.data.frame(fit)       # raw results table
plot(fit)                # alias for plot_auc_rank(fit)
subset(fit, q_value < 0.05)   # filter by significance
subset(fit, auc_strength > 0.8)  # filter by discrimination
```

---

## 2. Visualization

Every plot returns a `ggplot2` object and labels only a limited subset of
biomarkers (via `ggrepel`) to keep wide-matrix views readable.

### 2.1 Rank Plot — `plot_auc_rank()`

Ordered discrimination strengths with optional DeLong CI error bars.
Colour indicates effect direction (blue = higher in positives, red = lower).

```r
plot_auc_rank(fit, n_label = 20, show_ci = TRUE)
```

![Rank plot](man/figures/README-rank.png)

### 2.2 Volcano Plot — `plot_auc_volcano()`

Effect magnitude (\|AUC − 0.5\|) against statistical evidence (−log₁₀ q).
Biomarkers with q < cutoff are highlighted in red.

```r
plot_auc_volcano(fit, n_label = 20, q_cutoff = 0.05)
```

![Volcano plot](man/figures/README-volcano.png)

### 2.3 Forest Plot — `plot_auc_forest()`

Point estimates with confidence intervals for selected biomarkers.

```r
plot_auc_forest(fit, n = 8)                          # top 8 by auc_strength
plot_auc_forest(fit, biomarkers = c("X1", "X2"))     # specific set
```

![Forest plot](man/figures/README-forest.png)

### 2.4 ROC Curves — `plot_roc_top()`

Empirical ROC curves for deliberately selected biomarkers. Supports optional
bootstrap CI ribbons.

```r
plot_roc_top(fit, X = X, y = y)                            # top 6 by default
plot_roc_top(fit, X = X, y = y, biomarkers = c("X1", "X2"))
plot_roc_top(fit, X = X, y = y, add_ci = TRUE, boot_n = 500)  # with CI ribbons
```

![ROC curves](man/figures/README-roc.png)

### 2.5 Stability Plot — `plot_auc_stability()`

Bootstrap rank distributions with median (point) and IQR (error bar).

```r
plot_auc_stability(stab, n_label = 20)
```

![Stability plot](man/figures/README-stability.png)

---

## 3. Paired Comparisons: `compare_auc()`

`compare_auc()` tests AUC differences between biomarkers measured on the
**same subjects**. It supports three selection modes, three hypothesis types,
and DeLong or bootstrap inference.

### Comparison Modes

| Mode | Argument | Description |
|------|----------|-------------|
| Reference | `reference = "X1"` | Every other biomarker vs one reference |
| Top-N | `top_n = 5` | All pairs among the top N (flagged exploratory) |
| Named set | `biomarkers = c("a","b","c")` | All pairs within a specific list |

### Hypothesis Types

| Hypothesis | `hypothesis` | H0 | Use case |
|------------|-------------|-----|----------|
| Superiority | `"superiority"` | Δ = 0 | Is one biomarker better? |
| Non-inferiority | `"noninferiority"` | Δ ≤ −margin | Is one at least as good? |
| Equivalence (TOST) | `"equivalence"` | Δ ≤ −m or Δ ≥ m | Are they practically the same? |

With `alternative = "two.sided"` (default), `"greater"`, or `"less"`.

```r
# Two-sided superiority (default)
compare_auc(fit, X, y, reference = "X1")

# Directional: H0: AUC_a ≤ AUC_b (biomarker a is no better than b)
compare_auc(fit, X, y, reference = "X1", alternative = "greater")

# Non-inferiority: H0: AUC_a is worse than AUC_b by at least 0.05
compare_auc(fit, X, y, reference = "X1",
  hypothesis = "noninferiority", margin = 0.05)

# Equivalence (TOST): are AUCs within ±0.15?
compare_auc(fit, X, y, reference = "X1",
  hypothesis = "equivalence", margin = 0.15)
```

**Output columns**: `biomarker_a`, `biomarker_b`, `auc_a`, `auc_b`,
`auc_diff`, `std_error`, `conf_low`, `conf_high`, `p_value`, `q_value`,
`n_common`, `n_pos`, `n_neg`, `hypothesis`, `margin`.

When `top_n` selects and tests on the same data, the result carries
`selection_status = "same_data"` and a warning that inference is exploratory.

### Multiplicity Adjustment

```r
compare_auc(fit, X, y, top_n = 5, adjust = "BH")
compare_auc(fit, X, y, reference = "X1", adjust = "holm")
compare_auc(fit, X, y, biomarkers = c("X1","X2","X3"), adjust = "bonferroni")
```

---

## 4. Global Test: `compare_auc_global()`

Tests H₀: AUC₁ = AUC₂ = … = AUCₚ on a common subject set using the
joint DeLong covariance matrix and a Wald statistic.

```r
global <- compare_auc_global(fit, X, y)              # all biomarkers
global <- compare_auc_global(fit, X, y,
  biomarkers = c("X1", "X2", "X3"))                   # selected set
print(global)
```

Output: χ² statistic, degrees of freedom, p-value, per-biomarker AUC
estimates, covariance matrix, contrast specification, and sample counts.
The default safety limit is 100 biomarkers; raise it with `max_biomarkers`.

With two biomarkers, the global Wald statistic equals the squared
paired-DeLong z statistic (within numerical tolerance).

---

## 5. Rank Stability: `auc_stability()`

Quantifies how much bootstrap-determined ranks change under sampling
variability. Resamples positive and negative subjects separately with
replacement.

```r
stab <- auc_stability(X, y, times = 1000, top_k = c(10, 25, 50), seed = 42)
print(stab)              # median rank, Q25/Q75, top-1 frequency, mean/SD AUC
head(stab$rank_summary)  # full rank-distribution summary
head(stab$top_k_probs)   # top-k selection probabilities
stab$coselection          # pairwise co-selection frequencies
```

---

## 6. Simulation

`aucmat` provides seven simulation functions. The core distinction is whether
AUC and between-biomarker correlation are **independent** (free parameters) or
**linked** (determined by the binormal constraint).

### 6.1 Simulator Comparison

| Function | AUC + Correlation | Approach | Speed |
|----------|-------------------|----------|-------|
| `simulate_auc_matrix()` | **Independent** | Class-conditional MVN, closed-form mean shift, named correlation structures | Fast |
| `generate_data_probit()` | **Independent** | Latent MVN, numerical ρ calibration | Slower |
| `generate_data_analytical()` | **Linked** | Sequential binormal decomposition | Fast |
| `generate_auc_vector()` | **AUC only** | Rank construction, exact finite-sample AUC | Instant |
| `generate_auc_cor_vector()` | **AUC exact, r approx** | Rank + Box-Cox tuning | Fast |
| `simulate_auc_correlation()` | Sampling distribution | Monte Carlo over binormal draws | Moderate |
| `validate_simulation()` | Calibration check | Repeated independent draws + diagnostics | Slow |

### 6.2 `simulate_auc_matrix()` — Recommended General-Purpose Simulator

Draws biomarkers class-conditionally: X\|Y=0 ~ MVN(μ₀, R), X\|Y=1 ~ MVN(μ₁, R).
The outcome Y is fixed **before** any biomarker is drawn — supplied outcomes
keep their exact values and row order; generated outcomes get exact class
counts by default.

```r
# Exchangeable: all off-diagonal entries equal
sim_ex <- simulate_auc_matrix(
  n = 500, prevalence = 0.3,
  target_aucs = c(0.9, 0.8, 0.7, 0.65),
  correlation = 0.4, structure = "exchangeable"
)

# AR(1): correlation decays with distance |i−j|
sim_ar1 <- simulate_auc_matrix(
  n = 500, prevalence = 0.3,
  target_aucs = c(0.9, 0.8, 0.7, 0.65),
  correlation = 0.6, structure = "ar1"
)

# Block: within-block and between-block correlations
sim_blk <- simulate_auc_matrix(
  n = 500, prevalence = 0.3,
  target_aucs = c(0.9, 0.8, 0.7, 0.65),
  structure = "block", block_sizes = c(2, 2),
  rho_within = 0.6, rho_between = 0.1
)

# User-supplied matrix
sim_usr <- simulate_auc_matrix(
  n = 500, prevalence = 0.3,
  target_aucs = c(0.9, 0.8, 0.7),
  correlation = matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), 3, 3),
  structure = "user"
)

# Supplied outcome (retains exact values and row order)
y_supplied <- rbinom(500, 1, 0.3)
sim_y <- simulate_auc_matrix(
  y = y_supplied,
  target_aucs = c(0.9, 0.8, 0.7), correlation = 0.3, structure = "exchangeable"
)
```

**Correlation structures:**

| `structure` | Parameters | Off-diagonal entries |
|-------------|-----------|---------------------|
| `"user"` | `correlation` = full matrix | As supplied |
| `"exchangeable"` | `correlation` = single value | Constant ρ everywhere |
| `"ar1"` | `correlation` = single value | ρ^{\|i−j\|} |
| `"block"` | `block_sizes`, `rho_within`, `rho_between` | Block-diagonal |

**Feasibility diagnostics:** When the requested correlation matrix is not
positive definite, `feasibility = "error"` (default) stops with a structured
condition. `feasibility = "nearest"` projects to the nearest valid matrix and
reports the adjustment magnitude in `$feasibility`.

### 6.3 `validate_simulation()` — Calibration Check

Repeats a `simulate_auc_matrix()` specification across many independent draws
and reports bias, RMSE, Monte Carlo standard error, and target-interval hit
rates for both AUCs and pairwise correlations. A single draw is never evidence
that a simulator is calibrated.

```r
val <- validate_simulation(
  n = 300, prevalence = 0.3,
  target_aucs = c(0.85, 0.75, 0.65),
  correlation = 0.3, structure = "exchangeable",
  times = 100, seed = 1
)
print(val)
# $auc_bias, $auc_rmse, $auc_mc_se, $auc_hit_rate per biomarker
# $correlation$bias, $correlation$rmse, $correlation$mc_se, $correlation$hit_rate
```

A bias that stays large relative to the Monte Carlo SE across replicates
signals a calibration problem, not sampling noise.

### 6.4 Other Simulators

```r
# Latent probit: independent AUC + correlation (numerical calibration)
generate_data_probit(n = 500, target_aucs = c(0.9, 0.8, 0.7),
  corr_matrix = diag(3) + 0.3 - diag(0.3, 3), prevalence = 0.3)

# Binormal sequential: fast, AUC-correlation linked
generate_data_analytical(n = 500, target_aucs = c(0.85, 0.75, 0.65),
  corr_matrix = diag(3), prevalence = 0.3)

# Single score with exact empirical AUC
generate_auc_vector(y, target_auc = 0.80)

# Single score with exact AUC + approximate Pearson r
generate_auc_cor_vector(y, target_auc = 0.80, target_cor = 0.45)

# Monte Carlo sampling distribution of (AUC, r) under binormal model
simulate_auc_correlation(y, target_auc = c(0.7, 0.8, 0.9), n_sim = 500)
```

### 6.5 How Simulation Works

Under the **binormal model**, biomarker values are normally distributed in
each class, shifted apart:

```
Negative class:  X | Y=0  ~  N(μ₀, σ²)
Positive class:  X | Y=1  ~  N(μ₁, σ²)
```

The AUC is determined by the standardized mean separation δ = (μ₁ − μ₀)/σ:

$$\text{AUC} = \Phi\left(\frac{\delta}{\sqrt{2}}\right)$$

To target a specific AUC, invert: δ = √2 · Φ⁻¹(AUC).

For multivariate simulation, biomarkers share a correlation matrix R.
`simulate_auc_matrix()` and `generate_data_analytical()` use this equal-variance
binormal relationship (fast, closed-form). `generate_data_probit()` uses a
latent probit model with numerical calibration for independent AUC-correlation
control.

---

## 7. Missing Data

Two strategies via `na_action`:

| `na_action` | Behaviour |
|-------------|-----------|
| `"featurewise"` (default) | Each biomarker uses all subjects with observed values — maximises data |
| `"complete"` | Only subjects with complete observations across all biomarkers — ensures comparable populations |

```r
fit <- aucmat(X, y, na_action = "featurewise")
```

Missing-outcome observations are always removed. Imputation should be done
**before** calling `aucmat()` — the package no longer performs internal
imputation.

---

## 8. Reproducibility

All functions with randomness accept a `seed` argument. With an explicit seed,
the global `.Random.seed` is restored on exit (present or absent), even on error:

```r
fit1 <- aucmat(X, y, ci = "bootstrap", boot_n = 500, seed = 123)
fit2 <- aucmat(X, y, ci = "bootstrap", boot_n = 500, seed = 123)
identical(fit1$results$auc_raw, fit2$results$auc_raw)  # TRUE
```

---

## 9. Documentation & Help

```r
# Package-level help
?aucmat

# Main functions
?aucmat
?compare_auc
?compare_auc_global
?auc_stability

# Simulation
?simulate_auc_matrix
?validate_simulation
?generate_data_probit
?generate_data_analytical
?generate_auc_vector
?generate_auc_cor_vector
?simulate_auc_correlation

# Visualization
?plot_auc_rank
?plot_auc_volcano
?plot_auc_forest
?plot_auc_stability
?plot_roc_top
```

## Vignettes

- [Introduction to aucmat](https://vanhungtran.github.io/aucmat/articles/introduction-to-aucmat.html) — screening, visualization, comparisons, stability, missing data
- [Simulating Biomarker Data](https://vanhungtran.github.io/aucmat/articles/simulating-biomarker-data.html) — mathematical foundations and comparison of all simulation engines

## Online Documentation

Full pkgdown site: <https://vanhungtran.github.io/aucmat>

## Citation

```bibtex
@software{tran_aucmat,
  author  = {Lucas VHH Tran},
  title   = {aucmat: Scalable and Statistically Principled ROC Analysis
             for High-Dimensional Biomarker Matrices},
  year    = {2026},
  version = {0.2.0},
  url     = {https://github.com/vanhungtran/aucmat}
}
```

## License

MIT © Lucas VHH Tran
