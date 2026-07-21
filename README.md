# aucmat

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/vanhungtran/aucmat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vanhungtran/aucmat/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Scalable and statistically principled ROC analysis for high-dimensional biomarker matrices in R.**

## What is aucmat?

`aucmat` is a matrix-first R package for univariate screening and ranking of
molecular biomarkers (genes, proteins, metabolites) measured on the same
subjects.  It reports **direction-preserving AUCs** — never silently
reversing a biomarker to force its AUC above 0.5 — with proper multiplicity
adjustment and bootstrap rank-stability assessment.

## Installation

```r
remotes::install_github("vanhungtran/aucmat")
```

## In 30 Seconds

```r
library(aucmat)

# Generate correlated biomarkers with controlled AUCs
sim <- generate_data_probit(
  n = 500, prevalence = 0.3,
  target_aucs = c(0.90, 0.80, 0.70),
  corr_matrix  = matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), 3, 3)
)

# Screen all biomarkers
fit <- aucmat(sim$data[, 1:3], sim$data$truth)
print(fit)

# Plot discrimination strengths
plot_auc_rank(fit)

# Compare top biomarkers
compare_auc(fit, sim$data[, 1:3], sim$data$truth, top_n = 3)

# Assess rank stability
stab <- auc_stability(sim$data[, 1:3], sim$data$truth, times = 500, seed = 42)
plot_auc_stability(stab)
```

## Features

| Category | Functions | Description |
|----------|-----------|-------------|
| **Matrix Screening** | `aucmat()` | Screen thousands of biomarkers with DeLong or bootstrap CIs |
| **Direction-Preserving** | `auc_raw`, `auc_strength` | Raw AUC on observed scale; discrimination strength always ≥ 0.5 |
| **Multiplicity** | BH, Holm, Bonferroni | Adjust p-values across all screened biomarkers |
| **Simulation** | `generate_data_probit()`, `generate_data_analytical()` | Latent probit (AUC + correlation independent) or binormal |
| **Comparisons** | `compare_auc()` | Paired DeLong comparisons on common subjects with safety limits |
| **Stability** | `auc_stability()` | Bootstrap rank distributions and top-k selection probabilities |
| **Visualization** | `plot_auc_rank()`, `plot_auc_volcano()`, `plot_auc_forest()`, `plot_auc_stability()`, `plot_roc_top()` | Publication-quality ggplot2 graphics |

## Two Simulation Engines

| | `generate_data_probit()` | `generate_data_analytical()` |
|---|---|---|
| **Approach** | Latent multivariate normal | Sequential binormal decomposition |
| **AUC control** | Target AUC via calibrated ρ | Target AUC via mean separation |
| **Correlation control** | Directly specified in Σ | Approximately achieved |
| **AUC–correlation link** | **Decoupled** (both free) | Linked (binormal constraint) |
| **Speed** | Slower (calibration) | Fast (closed form) |

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

MIT © Lucas VHH Tran
