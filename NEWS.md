# aucmat 0.2.0

## New Features

### Multivariable Panels: Combining Biomarkers
- Two new `fit_auc_panel()` combination methods: `"su_liu"`, a closed-form
  linear combination that maximizes AUC under the multivariate-normal
  (binormal) model (Su & Liu, 1993), with shrinkage toward the diagonal
  covariance for stability when biomarkers are collinear or p approaches n;
  and `"unweighted"`, a fitting-free equal-weight sum of direction-aligned
  standardized biomarkers, useful as an assumption-light baseline
- New `predict.aucmat_panel()`: apply a fitted panel's coefficients (and its
  stored centering/scaling) to new subjects
- New `plot_roc_panel()`: ROC curve for a fitted panel's honest
  cross-validated score overlaid with its top individual biomarkers -- the
  non-deprecated, honestly-scored successor to `plot_roc_with_combos()`
- New `compare_panel_auc()`: DeLong or bootstrap paired test of whether the
  combined panel actually outperforms each of its individual biomarkers,
  reusing the `compare_auc()` engine

## Bug Fixes

### `fit_auc_panel()` correctness
- **Nested cross-validation.** The reported `auc_cv` previously reused a
  penalty (`lambda.min`) selected once on the *full* dataset inside every
  outer CV fold, and centering/scaling were likewise computed once on the
  full data. Both are now re-derived from each fold's training portion
  only, so `auc_cv` no longer carries the resulting optimistic bias.
- **Elastic-net default.** `method = "elasticnet"` defaulted to `alpha = 0`,
  which is ridge regression -- elastic net never actually differed from
  ridge unless a user manually supplied `alpha`. The default is now 0.5.
- **RNG safety.** `fit_auc_panel(..., seed = ...)` previously left the
  caller's global `.Random.seed` mutated on exit, and the final (full-data)
  model's own internal penalty-selection randomness was not covered by
  `seed` at all, so results were not fully reproducible. Both the final fit
  and the honest CV loop now run inside a single `with_seed()` block.

## Changed

- `fit_auc_panel()`'s `fitted_model` component is now a list of
  `coefficients`/`means`/`sds` (used by `predict.aucmat_panel()`) instead of
  the raw `glmnet`/`glm` fit object, for all methods uniformly.
- `plot_roc_top()` and `plot_roc_panel()` now default to `add_ci = TRUE`
  (bootstrap confidence ribbons shown), matching `plot_roc_smooth()`'s
  existing default. Pass `add_ci = FALSE` for a faster, ribbon-free plot.

# aucmat 0.1.0

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
