aucmat: Advanced AUC Analysis and Visualization for R
================

- [aucmat: Advanced AUC Analysis and Visualization for R (comming
  soon!!!)](#aucmat-advanced-auc-analysis-and-visualization-for-r-comming-soon)
- [🧬 aucmat](#dna-aucmat)
  - [Work with AUC](#work-with-auc)
- [📖 Overview](#open_book-overview)
- [🚀 Installation](#rocket-installation)
- [Features](#features)
- [Main Functions](#main-functions)
  - [`plot_roc_with_combos()`](#plot_roc_with_combos)
  - [`generate_data_analytical()`](#generate_data_analytical)
  - [`tableroc()`](#tableroc)
- [Advanced Usage](#advanced-usage)
  - [Customizing ROC Plots](#customizing-roc-plots)
  - [Benchmarking Clustering
    Algorithms](#benchmarking-clustering-algorithms)
- [Package Structure](#package-structure)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)
- [Support](#support)
- [Acknowledgments](#acknowledgments)
- [🤝 Contributing](#handshake-contributing)

<!-- README.md is auto-generated from README.Rmd -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- README.md is auto-generated from README.Rmd -->

# aucmat: Advanced AUC Analysis and Visualization for R (comming soon!!!)

Lucas TRAN 10/09/2025

The `aucmat` package provides comprehensive tools for ROC curve
analysis, AUC calculation, and visualization in R. It’s designed for
bioinformaticians, statisticians, and data scientists working with
binary classification problems.

*aucmat*

<div align="center">

# 🧬 aucmat

### Work with AUC

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/yourusername/atcddd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourusername/atcddd/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/atcddd)](https://cran.r-project.org/package=atcddd)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<br>

> *“Work with AUC — all in R.”*

</div>

------------------------------------------------------------------------

## 📖 Overview

**`aucmat`** is an R package that simplifies working with

## 🚀 Installation

The code of *aucmat* is freely available at
<https://github.com/vanhungtran/aucmat>.

The following commands can be used to install this R package, and an R
version \>= 4.2.3 is required.

    library(devtools)
    install_github("vanhungtran/aucmat")

Load the package:

``` r
library(aucmat)
```

## Features

- **ROC Curve Analysis**: Generate ROC curves with confidence intervals
  for single predictors or combinations
- **AUC Calculation**: Compute AUC values with various statistical
  methods
- **Data Simulation**: Generate synthetic datasets with specified AUC
  characteristics
- **Missing Data Handling**: Multiple imputation strategies for handling
  missing values
- **Comprehensive Visualization**: Professional-quality ROC plots with
  customizable options

## Main Functions

### `plot_roc_with_combos()`

Generate ROC curves for single predictors or combinations with
confidence intervals.

``` r
library(aucmat)

# Generate example data
set.seed(42)
data <- data.frame(
  M1 = rnorm(200, mean = 0, sd = 1),
  M2 = rnorm(200, mean = 5, sd = 2),
  M3 = rnorm(200, mean = 10, sd = 3),
  M4 = rnorm(200, mean = 15, sd = 4),
  Outcome = sample(c(0, 1), 200, replace = TRUE)
)

# Generate ROC curves for all combinations of 1-4 predictors
roc_results <- plot_roc_with_combos(
  data = data,
  outcome = "Outcome",
  predictors = c("M1", "M2", "M3", "M4"),
  combo_sizes = 1:4,
  add_ci = TRUE,
  boot_n = 1000
)

# Display the plot
print(roc_results$plot)
```

### `generate_data_analytical()`

Generate synthetic datasets with specified AUC values and correlation
structure.

``` r
# Generate data with specific AUC targets
sim_data <- generate_data_analytical(
  n = 1000,
  prevalence = 0.3,
  target_aucs = c(0.9, 0.85, 0.8, 0.75),
  corr_matrix = matrix(c(
    1.0, 0.6, 0.3, 0.1,
    0.6, 1.0, 0.5, 0.2,
    0.3, 0.5, 1.0, 0.4,
    0.1, 0.2, 0.4, 1.0
  ), nrow = 4)
)

# Check achieved AUCs
print(sim_data$achieved_aucs)
```

### `tableroc()`

Compute AUC and confidence intervals for each numeric column in a data
frame.

``` r
# Create example data with missing values
X <- data.frame(
  M1 = c(rnorm(95), rep(NA, 5)),
  M2 = rnorm(100),
  M3 = c(rep(NA, 3), rnorm(97))
)
y <- factor(rbinom(100, 1, 0.5), labels = c("neg", "pos"))

# Compute AUC table with median imputation
auc_table <- tableroc(
  X = X,
  y = y,
  na_impute = "median",
  rank = TRUE
)

print(auc_table)
```

## Advanced Usage

### Customizing ROC Plots

``` r
# Customized ROC plot
roc_results <- plot_roc_with_combos(
  data = data,
  outcome = "Outcome",
  predictors = c("M1", "M2", "M3", "M4"),
  combo_sizes = 1:4,
  add_ci = TRUE,
  boot_n = 2000,
  ci_points = 101,
  palette = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"),
  title = "Custom ROC Analysis",
  legacy_axes = FALSE
)

print(roc_results$plot)
```

### Benchmarking Clustering Algorithms

``` r
# Benchmark different clustering algorithms
cluster_results <- benchmark_clusterings(
  data = iris,
  features = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
  ks = 2:4,
  algorithms = c("kmeans", "hclust_complete", "hclust_average"),
  distances = c("euclidean", "manhattan"),
  status_col = "Species"
)

print(cluster_results$plot)
```

## Package Structure

The package includes the following main components:

- **ROC Analysis**: Functions for generating and analyzing ROC curves
- **Data Simulation**: Tools for creating synthetic datasets with known
  properties
- **Statistical Utilities**: Helper functions for imputation,
  normalization, and statistical testing
- **Visualization**: Custom plotting functions for professional-quality
  graphics

## Documentation

Full documentation is available through R’s help system:

``` r
# View function documentation
?plot_roc_with_combos
?generate_data_analytical
?tableroc
?benchmark_clusterings
```

## Contributing

Contributions to `aucmat` are welcome! Please feel free to:

1.  Fork the repository
2.  Create a feature branch
3.  Make your changes
4.  Submit a pull request

## License

This package is released under the MIT License. See the `LICENSE` file
for details.

## Citation

If you use `aucmat` in your research, please cite it as:

    Tran, L. V. H. H. (2025). aucmat: Advanced AUC Analysis and Visualization for R. 
    R package version 0.1.0. https://github.com/vanhungtran/aucmat

## Support

For bug reports, feature requests, or questions, please open an issue on
the [GitHub repository](https://github.com/vanhungtran/aucmat/issues).

## Acknowledgments

This package was developed with support from CK-CARE (Christine Kühne -
Center for Allergy Research and Education).


    To render this README.Rmd to README.md, use:

    ```r
    rmarkdown::render("README.Rmd", output_format = "github_document")

This README.Rmd provides: 1. A comprehensive introduction to the package
2. Installation instructions 3. Examples of main functions 4. Advanced
usage examples 5. Package structure overview 6. Documentation and
contribution guidelines 7. License and citation information

The document uses R code chunks with `eval=FALSE` to show the code
without executing it, making it suitable for a GitHub README that
displays code examples without requiring execution.

------------------------------------------------------------------------

## 🤝 Contributing

We welcome contributions! Please see our
[CONTRIBUTING.md](CONTRIBUTING.md) guide for how to:

- Report bugs
- Suggest features
- Submit pull requests
- Improve documentation

------------------------------------------------------------------------
