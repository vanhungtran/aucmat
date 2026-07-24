# Package index

## Matrix Screening

Screen a biomarker matrix against a binary outcome

- [`aucmat()`](https://vanhungtran.github.io/aucmat/reference/aucmat.md)
  : Screen a biomarker matrix against a binary outcome

## Data Simulation

Generate synthetic biomarker data with controlled AUC and correlation

- [`simulate_auc_matrix()`](https://vanhungtran.github.io/aucmat/reference/simulate_auc_matrix.md)
  : Simulate a biomarker matrix with target AUCs and a parametrized
  correlation structure

- [`simulate_auc_copula()`](https://vanhungtran.github.io/aucmat/reference/simulate_auc_copula.md)
  : Simulate correlated biomarkers with two-phase Copula-AUC generation

- [`simulate_hurdle_auc()`](https://vanhungtran.github.io/aucmat/reference/simulate_hurdle_auc.md)
  : Simulate zero-inflated biomarker data with controlled Hurdle-AUC

- [`validate_simulation()`](https://vanhungtran.github.io/aucmat/reference/validate_simulation.md)
  :

  Validate a
  [`simulate_auc_matrix()`](https://vanhungtran.github.io/aucmat/reference/simulate_auc_matrix.md)
  specification by repeated simulation

- [`generate_data_probit()`](https://vanhungtran.github.io/aucmat/reference/generate_data_probit.md)
  : Generate correlated biomarkers with specified AUCs via latent probit
  model

- [`generate_data_analytical()`](https://vanhungtran.github.io/aucmat/reference/generate_data_analytical.md)
  : Generate Synthetic Data with Specific AUC and Correlation

- [`generate_auc_vector()`](https://vanhungtran.github.io/aucmat/reference/generate_auc_vector.md)
  : Generate a continuous score vector for a target AUC

- [`generate_auc_cor_vector()`](https://vanhungtran.github.io/aucmat/reference/generate_auc_cor_vector.md)
  : Generate a continuous score vector for target AUC and correlation

- [`simulate_auc_correlation()`](https://vanhungtran.github.io/aucmat/reference/simulate_auc_correlation.md)
  : Simulate repeated X vectors for a fixed y and record AUC and
  correlation

## Comparisons & Stability

Compare biomarker AUCs and assess ranking stability

- [`compare_auc()`](https://vanhungtran.github.io/aucmat/reference/compare_auc.md)
  : Compare AUCs of selected biomarkers (paired, common-subject)
- [`compare_auc_global()`](https://vanhungtran.github.io/aucmat/reference/compare_auc_global.md)
  : Global Wald test for equality of multiple correlated AUCs
- [`auc_stability()`](https://vanhungtran.github.io/aucmat/reference/auc_stability.md)
  : Bootstrap rank-stability analysis for biomarker screening
- [`roc_test()`](https://vanhungtran.github.io/aucmat/reference/roc_test.md)
  : Single-biomarker ROC test against chance or another biomarker

## Visualization

Publication-quality graphics – every ROC/panel-ROC plot shows bootstrap
confidence ribbons by default

- [`plot_auc_rank()`](https://vanhungtran.github.io/aucmat/reference/plot_auc_rank.md)
  : Rank plot: ordered discrimination strengths
- [`plot_auc_volcano()`](https://vanhungtran.github.io/aucmat/reference/plot_auc_volcano.md)
  : Volcano plot: effect magnitude against statistical evidence
- [`plot_auc_forest()`](https://vanhungtran.github.io/aucmat/reference/plot_auc_forest.md)
  : Forest plot: selected AUCs with confidence intervals
- [`plot_auc_stability()`](https://vanhungtran.github.io/aucmat/reference/plot_auc_stability.md)
  : Stability plot: rank distributions and top-k probabilities
- [`plot_roc_top()`](https://vanhungtran.github.io/aucmat/reference/plot_roc_top.md)
  : ROC curves for a small set of biomarkers
- [`plot_roc_smooth()`](https://vanhungtran.github.io/aucmat/reference/plot_roc_smooth.md)
  : Plot smoothed ROC curves for selected biomarkers
- [`plot_auc_pr()`](https://vanhungtran.github.io/aucmat/reference/plot_auc_pr.md)
  : Plot Precision-Recall curves for selected biomarkers
- [`plot_correlation_heatmap()`](https://vanhungtran.github.io/aucmat/reference/plot_correlation_heatmap.md)
  : Plot achieved vs requested correlation heatmap
- [`plot(`*`<aucmat_compare>`*`)`](https://vanhungtran.github.io/aucmat/reference/plot.aucmat_compare.md)
  : Forest plot of paired AUC differences
- [`plot(`*`<aucmat_simulation>`*`)`](https://vanhungtran.github.io/aucmat/reference/plot.aucmat_simulation.md)
  : Plot simulated biomarker data
- [`theme_aucmat()`](https://vanhungtran.github.io/aucmat/reference/theme_aucmat.md)
  : aucmat ggplot2 theme
- [`scale_colour_aucmat_direction()`](https://vanhungtran.github.io/aucmat/reference/scale_colour_aucmat_direction.md)
  : aucmat direction colour scale (discrete)
- [`scale_fill_aucmat_direction()`](https://vanhungtran.github.io/aucmat/reference/scale_fill_aucmat_direction.md)
  : aucmat fill colour scale (discrete)
- [`scale_colour_aucmat_significance()`](https://vanhungtran.github.io/aucmat/reference/scale_colour_aucmat_significance.md)
  : aucmat significance scale

## Multivariable Panels

Combine biomarkers into a single score, with honest nested
cross-validation

- [`fit_auc_panel()`](https://vanhungtran.github.io/aucmat/reference/fit_auc_panel.md)
  : Fit and validate a multivariable biomarker panel
- [`predict(`*`<aucmat_panel>`*`)`](https://vanhungtran.github.io/aucmat/reference/predict.aucmat_panel.md)
  : Predict panel scores for new subjects
- [`plot_roc_panel()`](https://vanhungtran.github.io/aucmat/reference/plot_roc_panel.md)
  : ROC curve for a fitted biomarker panel vs. its individual components
- [`compare_panel_auc()`](https://vanhungtran.github.io/aucmat/reference/compare_panel_auc.md)
  : Test whether a combined panel improves on its individual biomarkers

## Cross-Validation & Study Planning

- [`cv_aucmat()`](https://vanhungtran.github.io/aucmat/reference/cv_aucmat.md)
  : Cross-validated AUC screening
- [`power_auc_matrix()`](https://vanhungtran.github.io/aucmat/reference/power_auc_matrix.md)
  : Power and sample size for AUC comparisons
- [`aucmat_by()`](https://vanhungtran.github.io/aucmat/reference/aucmat_by.md)
  : Stratified biomarker screening by group
- [`aucmat_multiclass()`](https://vanhungtran.github.io/aucmat/reference/aucmat_multiclass.md)
  : Screen biomarkers against a multiclass outcome

## Hurdle-AUC

Two-stage model for zero-inflated biomarker data

- [`hurdle_auc()`](https://vanhungtran.github.io/aucmat/reference/hurdle_auc.md)
  : Hurdle-AUC for zero-inflated biomarkers
- [`simulate_hurdle_auc()`](https://vanhungtran.github.io/aucmat/reference/simulate_hurdle_auc.md)
  : Simulate zero-inflated biomarker data with controlled Hurdle-AUC
- [`plot_hurdle_diagnostics()`](https://vanhungtran.github.io/aucmat/reference/plot_hurdle_diagnostics.md)
  : Plot Hurdle-AUC diagnostics

## Example Data

Real biomarker data shipped with the package

- [`beatmg_baseline`](https://vanhungtran.github.io/aucmat/reference/beatmg_baseline.md)
  : BeatMG Olink proteomics: baseline (pre-treatment)
- [`beatmg_ontreatment`](https://vanhungtran.github.io/aucmat/reference/beatmg_ontreatment.md)
  : BeatMG Olink proteomics: on-treatment (timepoint 3)

## Deprecated

Functions retained for backward compatibility

- [`tableroc()`](https://vanhungtran.github.io/aucmat/reference/tableroc.md)
  : (Deprecated) AUC table for each numeric predictor column
- [`plot_roc_with_combos()`](https://vanhungtran.github.io/aucmat/reference/plot_roc_with_combos.md)
  : (Deprecated) Plot ROC curves for predictor combinations
