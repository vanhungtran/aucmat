# aucmat 0.2.0

## New Features

### Inference Extensions
- `compare_auc()` now supports superiority, non-inferiority, and equivalence (TOST) hypotheses
- Configurable `alternative` (two.sided, greater, less) and `conf_level`
- DeLong and bootstrap methods for paired comparisons
- Post-selection warning when `top_n` is used from the same data
- New `compare_auc_global()`: omnibus Wald test for equality of multiple correlated AUCs via joint DeLong covariance
- Structured condition classes for errors and RNG safety (`with_seed`)

### Simulation Engine
- New `simulate_auc_matrix()`: class-conditional multivariate normal biomarker simulator with closed-form binormal mean-shift calibration to target AUCs
- Parametrized correlation structures: `user` (full matrix), `exchangeable`, `ar1`, `block`
- Feasibility diagnostics with `error`/`nearest` handling and nearest positive-definite projection for infeasible correlation targets
- New `validate_simulation()`: repeated-draw calibration harness reporting bias, RMSE, Monte Carlo standard error, and target-interval hit rates for both AUCs and pairwise correlations

### Shared Engines
- Joint DeLong covariance engine for correlated full-AUC inference
- Unified stratified resampling engine with minimum-success rules and failure accounting
- Shared nearest positive-definite correlation projection (`R/pd-projection.R`)

# aucmat 0.1.0

## New Features

### Matrix-First Architecture
- New primary function `aucmat()` for screening biomarker matrices against binary outcomes
- Direction-preserving AUC: `auc_raw` (observed direction), `auc_strength` (≥ 0.5), `effect_direction`
- DeLong and stratified-bootstrap confidence intervals
- BH, Holm, and Bonferroni multiplicity adjustment
- Feature-wise and complete-case missing data handling
- S3 class `aucmat_screen` with `print()`, `summary()`, `as.data.frame()`, `plot()`, `subset()` methods

### Latent Probit Simulator
- New `generate_data_probit()`: jointly controls AUCs and between-biomarker correlations via latent multivariate normal model
- Numerical calibration of latent correlations for target AUCs
- Automatic positive-definiteness repair

### Paired Comparisons & Stability
- New `compare_auc()`: bounded paired DeLong comparisons on common subjects with safety limits
- New `auc_stability()`: stratified bootstrap rank distributions and top-k selection probabilities

### Visualization
- `plot_auc_rank()`: ordered discrimination strengths with CI error bars
- `plot_auc_volcano()`: effect magnitude vs. statistical evidence
- `plot_auc_forest()`: selected AUCs with confidence intervals
- `plot_auc_stability()`: bootstrap rank distributions
- `plot_roc_top()`: ROC curves for selected biomarkers

## Deprecated

- `tableroc()` is deprecated in favor of `aucmat()`
- `plot_roc_with_combos()` is deprecated in favor of `plot_roc_top()`

## Removed

- `install_and_load()` utility function
- Internal imputation functions (impute your data before calling `aucmat()`)

## Internal

- Shared binary outcome normalization via `.normalize_binary_y()`
- Matrix AUC engine using Mann-Whitney rank formulation
- Component-based architecture with separated concerns (engine, inference, multiplicity, stability)
