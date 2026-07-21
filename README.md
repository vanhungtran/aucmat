# aucmat

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/vanhungtran/aucmat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vanhungtran/aucmat/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Scalable and statistically principled ROC analysis for high-dimensional biomarker matrices in R.**

## What is aucmat?

`aucmat` is a matrix-first R package for univariate screening and ranking of
molecular biomarkers (genes, proteins, metabolites) measured on the same
subjects. It reports **direction-preserving AUCs** with proper multiplicity
adjustment and bootstrap rank-stability assessment.

## Installation

```r
remotes::install_github("vanhungtran/aucmat")
```

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

# 4. Compare AUCs - three approaches
compare_auc(fit, sim$data[, 1:3], sim$data$truth, top_n = 3)          # all pairs
compare_auc(fit, sim$data[, 1:3], sim$data$truth, reference = "X1")   # vs reference

# 5. Assess rank stability
stab <- auc_stability(sim$data[, 1:3], sim$data$truth, times = 500, seed = 42)
plot_auc_stability(stab)
```

## How Simulation Works

AUC measures: if you pick one positive and one negative subject at random,
how often does the positive have a higher biomarker value?

Under the binormal model, biomarker values are normally distributed in each
class with a mean shift. For target AUC = 0.8, the positive class mean is
shifted by about 1.19 standard deviations: `delta = sqrt(2) * qnorm(AUC)`.

| Function | AUC + Correlation | Speed |
|----------|-------------------|-------|
| `generate_data_probit()` | **Independent** (both free) | Slower (calibration) |
| `generate_data_analytical()` | **Linked** (binormal constraint) | Fast |

## Plots

### Rank Plot
Ordered discrimination strengths with CI error bars:

![Rank plot](man/figures/README-rank.png)

### Volcano Plot
Effect magnitude vs statistical evidence (-log10 q-value):

![Volcano plot](man/figures/README-volcano.png)

### Forest Plot
Multiple AUCs with DeLong confidence intervals:

![Forest plot](man/figures/README-forest.png)

### ROC Curves
ROC curves for selected biomarkers overlaid:

![ROC curves](man/figures/README-roc.png)

### Stability Plot
Bootstrap rank distributions with median and IQR:

![Stability plot](man/figures/README-stability.png)

## All Functions

### Simulation

```r
# Latent probit: independent AUC + correlation control
generate_data_probit(n = 500, target_aucs = c(0.9, 0.8, 0.7),
  corr_matrix = diag(3), prevalence = 0.3)

# Binormal sequential: fast, AUC-correlation linked
generate_data_analytical(n = 500, target_aucs = c(0.85, 0.75, 0.65),
  corr_matrix = diag(3), prevalence = 0.3)

# Single score with exact empirical AUC
generate_auc_vector(y, target_auc = 0.80)

# Single score with exact AUC + approximate Pearson r
generate_auc_cor_vector(y, target_auc = 0.80, target_cor = 0.45)

# Monte Carlo sampling distribution of (AUC, r)
simulate_auc_correlation(y, target_auc = c(0.7, 0.8, 0.9), n_sim = 500)
```

### Screening

```r
fit <- aucmat(X, y)                              # default: DeLong, BH
fit <- aucmat(X, y, ci = "bootstrap", boot_n = 2000)  # bootstrap CIs
fit <- aucmat(X, y, adjust = "bonferroni")       # stricter adjustment
fit <- aucmat(X, y, na_action = "complete")       # complete-case only
```

### S3 Methods

```r
print(fit)          # compact top-N display
summary(fit)        # class counts, missingness, multiplicity
as.data.frame(fit)  # raw results table
plot(fit)           # same as plot_auc_rank(fit)
subset(fit, q_value < 0.05)  # filter by q-value
```

### Visualization

```r
plot_auc_rank(fit)       # rank plot with CI error bars
plot_auc_volcano(fit)    # effect size vs -log10(q)
plot_auc_forest(fit)     # forest plot of multiple AUCs with CIs
plot_auc_stability(stab) # bootstrap rank distributions
plot_roc_top(fit, X, y)  # ROC curves for selected biomarkers
```

### Comparisons - Three Approaches

```r
# 1. All pairwise among top N
compare_auc(fit, X, y, top_n = 3)

# 2. All vs a reference biomarker
compare_auc(fit, X, y, reference = "IL13")

# 3. A specific named set
compare_auc(fit, X, y, biomarkers = c("IL13", "CCL17", "CCL22"))

# With multiplicity adjustment
compare_auc(fit, X, y, top_n = 5, adjust = "BH")
```

### Stability

```r
stab <- auc_stability(X, y, times = 1000, top_k = c(10, 25, 50), seed = 42)
print(stab)          # median ranks, top-1 frequencies
stab$top_k_probs     # top-k selection probabilities
stab$coselection     # pairwise co-selection frequencies
```

## Vignettes

- [Introduction to aucmat](https://vanhungtran.github.io/aucmat/articles/introduction-to-aucmat.html)
- [Simulating Biomarker Data](https://vanhungtran.github.io/aucmat/articles/simulating-biomarker-data.html)

## Citation

```bibtex
@software{tran_aucmat,
  author  = {Lucas VHH Tran},
  title   = {aucmat: Scalable and Statistically Principled ROC Analysis
             for High-Dimensional Biomarker Matrices},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/vanhungtran/aucmat}
}
```

## License

MIT (c) Lucas VHH Tran
