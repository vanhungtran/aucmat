# aucmat

**aucmat** is an R package for computing AUC (Area Under the Curve) tables for multiple predictors using the `pROC` package, with robust missing-data handling and optional ranking.

## Features

- **Compute ROC AUC and confidence intervals**  
  Calculate ROC AUC and its confidence intervals (using DeLong or bootstrap methods) for each numeric predictor column.

- **Flexible missing-data strategies:**  
  Handle missing values with a variety of imputation approaches:
  - `none` (row-drop or per-column omission)
  - `mean`, `median`, `constant` value
  - `median_by_class` (stratified by outcome)
  - `halfmin`, `quantile`-based imputation
  - `knn` (base R implementation)
  - `mice` (if installed)
  - `missForest` (if installed)

- **Optional ranking and tidy output**  
  Rank predictors by AUC and display results in a tidy, easily readable format.

## Installation

You can install **aucmat** from GitHub using `devtools`:

```r
devtools::install_github("vanhungtran/aucmat")
```

## Usage

```r
library(aucmat)

# Example: Compute AUC table for predictors in a data frame
result <- aucmat(
  data = your_data,
  outcome_col = "Outcome",
  predictors = c("Predictor1", "Predictor2", "Predictor3"),
  na_strategy = "median_by_class",
  rank = TRUE
)

print(result)
```

## Arguments

- `data`: Data frame containing outcome and predictor columns.
- `outcome_col`: Name of the outcome (binary) column.
- `predictors`: Vector of numeric predictor column names.
- `na_strategy`: Missing data imputation strategy (see above).
- `rank`: Logical; if `TRUE`, rank predictors by AUC.

## Missing Data Strategies

| Strategy          | Description                                   |
|-------------------|-----------------------------------------------|
| none              | Drop rows or omit per-column                  |
| mean              | Impute missing with column mean               |
| median            | Impute missing with column median             |
| constant          | Impute missing with a constant value          |
| median_by_class   | Impute with median stratified by outcome      |
| halfmin           | Impute with half of the minimum value         |
| quantile          | Impute using a specified quantile value       |
| knn               | k-Nearest Neighbors imputation                |
| mice              | Multiple Imputation by Chained Equations      |
| missForest        | Random forest-based imputation                |

## License

This package is licensed under the MIT License.

## Authors

- [vanhungtran](https://github.com/vanhungtran)

## Contributing

Contributions and suggestions are welcome! Please submit issues or pull requests via [GitHub](https://github.com/vanhungtran/aucmat).

## References

- [pROC package documentation](https://cran.r-project.org/web/packages/pROC/index.html)
- [mice package documentation](https://cran.r-project.org/web/packages/mice/index.html)
- [missForest package documentation](https://cran.r-project.org/web/packages/missForest/index.html)
