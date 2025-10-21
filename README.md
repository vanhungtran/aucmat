aucmat: A Toolbox for Biomarkers Data Analysis
================

- [🧬 aucmat](#dna-aucmat)
  - [A Toolbox for Biomarkers Data Analysis](#a-toolbox-for-biomarkers-data-analysis)
- [📖 Overview](#open_book-overview)
- [🚀 Installation](#rocket-installation)
- [✨ Features](#features)
- [📊 Main Functions](#main-functions)
  - [`tableroc()`](#tableroc)
  - [`plot_roc_with_combos()`](#plot_roc_with_combos)
  - [`generate_data_analytical()`](#generate_data_analytical)
  - [`install_and_load()`](#install_and_load)
- [💡 Usage Examples](#usage-examples)
  - [Quick AUC Analysis](#quick-auc-analysis)
  - [ROC Curves with Combinations](#roc-curves-with-combinations)
  - [Synthetic Data Generation](#synthetic-data-generation)
  - [Missing Data Handling](#missing-data-handling)
- [🔧 Advanced Usage](#advanced-usage)
  - [Multiple Imputation Methods](#multiple-imputation-methods)
  - [Customizing ROC Plots](#customizing-roc-plots)
- [📚 Documentation](#documentation)
- [📄 License](#license)
- [📖 Citation](#citation)
- [🆘 Support](#support)
- [🤝 Contributing](#contributing)

<!-- README.md is auto-generated from README.Rmd -->

<div align="center">

# 🧬 aucmat

### A Toolbox for Biomarkers Data Analysis

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<br>

> *"Comprehensive ROC analysis and biomarker evaluation — all in R."*

</div>

------------------------------------------------------------------------

## 📖 Overview

**`aucmat`** is a comprehensive R package designed for biomarker data analysis with a focus on ROC curve analysis, AUC calculation, and advanced visualization. The package provides robust tools for handling missing data, generating synthetic datasets with specified properties, and creating publication-ready ROC plots with confidence intervals.

The package is particularly valuable for bioinformaticians, statisticians, and data scientists working with binary classification problems in biomedical research, offering both simple interfaces for quick analysis and advanced features for comprehensive biomarker evaluation.

## 🚀 Installation

The code of *aucmat* is freely available at <https://github.com/vanhungtran/aucmat>.

The following commands can be used to install this R package, and an R version >= 4.0.0 is required.

```r
# Install from GitHub
library(devtools)
install_github("vanhungtran/aucmat")

# Load the package
library(aucmat)
```

## ✨ Features

- **📈 ROC Curve Analysis**: Generate ROC curves with confidence intervals for single predictors or combinations
- **🎯 AUC Calculation**: Compute AUC values with DeLong or bootstrap confidence intervals
- **🧬 Data Simulation**: Generate synthetic datasets with specified AUC characteristics and correlation structures
- **🔧 Missing Data Handling**: Multiple imputation strategies including mean/median, KNN, MICE, and missForest
- **📊 Comprehensive Visualization**: Professional-quality ROC plots with customizable options
- **⚡ Package Management**: Intelligent package installation from CRAN, Bioconductor, and GitHub

## 📊 Main Functions

### `tableroc()`

**Primary Function**: Compute ROC AUC and confidence intervals for each numeric column in a dataset against a binary outcome.

This is the core function for quick biomarker screening and analysis.

**Key Features:**
- Supports multiple missing data imputation methods
- Provides both DeLong and bootstrap confidence intervals
- Handles various binary outcome formats (factor, logical, numeric 0/1)
- Optional ranking by AUC performance
- Comprehensive diagnostic information

```r
# Basic usage
library(aucmat)

# Create example data with missing values
set.seed(123)
X <- data.frame(
  biomarker1 = c(rnorm(95, mean = 1), rep(NA, 5)),
  biomarker2 = rnorm(100, mean = 0.5),
  biomarker3 = c(rep(NA, 3), rnorm(97, mean = -0.2))
)
y <- factor(rbinom(100, 1, 0.4), labels = c("Control", "Case"))

# Compute AUC table with median imputation
auc_results <- tableroc(
  X = X,
  y = y,
  na_impute = "median",
  rank = TRUE,
  ci_method = "delong"
)

print(auc_results)
```

### `plot_roc_with_combos()`

**Advanced ROC Analysis**: Generate ROC curves for single predictors or combinations with confidence intervals.

**Key Features:**
- Supports combinations of multiple predictors using logistic regression
- Bootstrap confidence intervals with ribbon visualization
- Customizable plot aesthetics and color palettes
- Option to skip CI computation for perfect classifiers (AUC = 1)
- Professional ggplot2-based visualizations

```r
# Generate ROC curves for predictor combinations
set.seed(42)
data <- data.frame(
  M1 = rnorm(200, mean = 0, sd = 1),
  M2 = rnorm(200, mean = 0.5, sd = 1.2),
  M3 = rnorm(200, mean = -0.3, sd = 0.8),
  M4 = rnorm(200, mean = 0.8, sd = 1.5),
  Outcome = sample(c(0, 1), 200, replace = TRUE, prob = c(0.6, 0.4))
)

# Generate ROC curves for combinations of 1-3 predictors
roc_results <- plot_roc_with_combos(
  data = data,
  outcome = "Outcome",
  predictors = c("M1", "M2", "M3", "M4"),
  combo_sizes = 1:3,
  add_ci = TRUE,
  boot_n = 1000,
  reorder_by_auc = TRUE,
  title = "ROC Analysis: Single and Combined Biomarkers"
)

# Display the plot
print(roc_results$plot)

# Access AUC values
print(roc_results$auc)
```

### `generate_data_analytical()`

**Synthetic Data Generation**: Create datasets with precise AUC values and correlation structures for simulation studies.

**Key Features:**
- Analytically generates data to match target AUC values exactly
- Supports custom correlation matrices between predictors
- Useful for method validation and simulation studies
- Returns achieved vs. target metrics for verification

```r
# Generate synthetic data with specific properties
target_aucs <- c(0.85, 0.75, 0.65, 0.55)
corr_matrix <- matrix(c(
  1.0, 0.4, 0.2, 0.1,
  0.4, 1.0, 0.3, 0.2,
  0.2, 0.3, 1.0, 0.4,
  0.1, 0.2, 0.4, 1.0
), nrow = 4)

sim_data <- generate_data_analytical(
  n = 1000,
  prevalence = 0.3,
  target_aucs = target_aucs,
  corr_matrix = corr_matrix
)

# Verify achieved results
cat("Target AUCs:", target_aucs, "\n")
cat("Achieved AUCs:", round(sim_data$achieved_aucs, 3), "\n")
```

### `install_and_load()`

**Package Management**: Intelligently install and load R packages from multiple sources.

**Key Features:**
- Automatic source detection (CRAN, Bioconductor, GitHub)
- Interactive GitHub repository search for unknown packages
- Graceful error handling with informative messages
- Support for direct GitHub installation (`username/repository`)

```r
# Install packages from various sources
packages_needed <- c(
  "ggplot2",              # CRAN
  "limma",                # Bioconductor  
  "tidyverse/dplyr",      # Direct from GitHub
  "pROC"                  # CRAN
)

install_and_load(packages_needed)
```

## 💡 Usage Examples

### Quick AUC Analysis

```r
library(aucmat)

# Load your biomarker data
# X: matrix/data.frame with biomarkers as columns
# y: binary outcome vector

# Quick screening of all biomarkers
results <- tableroc(X, y, rank = TRUE)
head(results)
```

### ROC Curves with Combinations

```r
# Compare single biomarkers vs combinations
roc_analysis <- plot_roc_with_combos(
  data = your_data,
  outcome = "disease_status", 
  predictors = c("biomarker1", "biomarker2", "biomarker3"),
  combo_sizes = 1:3,
  add_ci = TRUE
)

# View results
roc_analysis$plot
```

### Synthetic Data Generation

```r
# Create validation dataset
validation_data <- generate_data_analytical(
  n = 500,
  prevalence = 0.25,
  target_aucs = c(0.9, 0.8, 0.7),
  corr_matrix = diag(3)  # Independent predictors
)
```

### Missing Data Handling

```r
# Advanced missing data imputation
results_knn <- tableroc(X, y, na_impute = "knn", knn_k = 5)
results_mice <- tableroc(X, y, na_impute = "mice", mice_maxit = 10)
results_forest <- tableroc(X, y, na_impute = "missForest")

# Compare imputation methods
rbind(
  transform(results_knn, method = "KNN"),
  transform(results_mice, method = "MICE"),
  transform(results_forest, method = "missForest")
)
```

## 🔧 Advanced Usage

### Multiple Imputation Methods

The `tableroc()` function supports various imputation strategies:

```r
# Available imputation methods:
imputation_methods <- c(
  "none",              # No imputation
  "mean",              # Column means
  "median",            # Column medians  
  "constant",          # Constant value
  "median_by_class",   # Class-specific medians
  "halfmin",           # Half minimum positive value
  "quantile",          # Specified quantile
  "knn",               # K-nearest neighbors
  "mice",              # Multiple imputation (requires mice package)
  "missForest"         # Random forest (requires missForest package)
)

# Example with KNN imputation
results <- tableroc(
  X = X, 
  y = y,
  na_impute = "knn",
  knn_k = 10,
  knn_scale = TRUE,
  knn_weighted = TRUE
)
```

### Customizing ROC Plots

```r
# Customized ROC visualization
custom_roc <- plot_roc_with_combos(
  data = data,
  outcome = "outcome",
  predictors = biomarkers,
  combo_sizes = 1:2,
  add_ci = TRUE,
  boot_n = 2000,
  ci_points = 201,
  palette = c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00"),
  title = "Biomarker Performance Comparison",
  ribbon_alpha = 0.15,
  line_size = 1.2
)

# Further customize with ggplot2
library(ggplot2)
custom_roc$plot + 
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
```

## 📚 Documentation

Full documentation is available through R's help system:

```r
# View function documentation
?tableroc
?plot_roc_with_combos
?generate_data_analytical
?install_and_load

# Browse package vignettes
browseVignettes("aucmat")
```

## 📄 License

This package is released under the MIT License. See the `LICENSE` file for details.

## 📖 Citation

If you use `aucmat` in your research, please cite it as:

```
Tran, L. V. H. H. (2025). aucmat: A Toolbox for Biomarkers Data Analysis. 
R package version 0.0.0.9000. https://github.com/vanhungtran/aucmat
```

## 🆘 Support

For bug reports, feature requests, or questions, please open an issue on the [GitHub repository](https://github.com/vanhungtran/aucmat/issues).

## 🤝 Contributing

Contributions to `aucmat` are welcome! Please feel free to:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Commit your changes (`git commit -m 'Add amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

---

**Acknowledgments**: This package was developed with support from CK-CARE (Christine Kühne - Center for Allergy Research and Education).