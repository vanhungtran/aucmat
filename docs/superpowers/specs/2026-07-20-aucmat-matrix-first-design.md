# Matrix-first `aucmat` design

Date: 2026-07-20
Status: Approved concept; implementation not started

## Decision summary

`aucmat` will be a matrix-first R package for scalable and statistically principled ROC analysis of high-dimensional molecular biomarker data. Its primary task is univariate screening and ranking of genes, proteins, metabolites, and similar features measured on the same subjects.

The first release and JSS manuscript will not attempt to cover the full biomarker-evaluation lifecycle. Survival outcomes, joint models, Bayesian ROC models, diagnostic meta-analysis, clinical-utility analysis, sample-size calculations, and automatic multivariable panel construction are explicitly deferred.

The working manuscript title is:

> **aucmat: Scalable and statistically principled ROC analysis for high-dimensional biomarker matrices in R**

The package's distinguishing contribution will be the combination of:

1. a matrix-native computational interface;
2. direction-preserving AUC reporting;
3. simultaneous inference and multiplicity control;
4. selected paired comparisons for correlated biomarkers;
5. bootstrap rank-stability analysis; and
6. reproducible validation using correlated biomarker simulations with controlled AUCs.

## Target users and primary use case

The primary user is a biomedical researcher or statistician with:

- an `n x p` numeric matrix `X`;
- samples in rows;
- biomarkers in columns; and
- a binary phenotype or clinical outcome `y` of length `n`.

The expected scale is hundreds to low thousands of samples and thousands to tens of thousands of biomarkers. The package must still behave correctly for small matrices.

The primary workflow is:

1. screen every biomarker against the binary outcome;
2. quantify uncertainty and test against `AUC = 0.5`;
3. adjust for multiple testing;
4. inspect missingness and quality flags;
5. assess the stability of leading biomarkers; and
6. perform paired comparisons for a deliberately selected, small set of biomarkers.

## Goals

Version 1 will provide:

- fast AUC estimation over a biomarker matrix;
- correct handling of ties and paired subjects;
- raw, direction-preserving AUC values;
- a separate direction-independent discrimination measure for ranking;
- DeLong and stratified-bootstrap uncertainty estimates;
- two-sided tests of `AUC = 0.5`;
- BH, Holm, and Bonferroni multiplicity adjustment;
- explicit feature-wise missingness and analysis counts;
- paired AUC comparisons on a common subject set;
- stratified-bootstrap rank stability and top-k selection probabilities;
- focused publication-quality visualizations;
- S3 result classes with standard methods; and
- numerical, statistical, and performance validation.

## Non-goals

Version 1 will not provide:

- survival or competing-risk outcomes;
- time-dependent ROC curves or joint models;
- Bayesian AUC or ROC inference;
- diagnostic meta-analysis;
- NRI, IDI, calibration, or decision-curve analysis;
- sample-size or power calculations;
- automatic enumeration of biomarker combinations;
- internal multivariable model fitting or panel selection;
- internal MICE, random-forest, KNN, or class-specific imputation; or
- package installation and dependency-management helpers.

Multivariable panel construction may become a later companion package. It must not be added to `aucmat` version 1 without a separate design because valid panel assessment requires nested feature selection and model validation.

## Statistical definitions

For each biomarker, `auc_raw` is the probability that a randomly selected positive subject has a larger observed biomarker value than a randomly selected negative subject, with half credit for ties.

The package will never silently reverse a biomarker to force its reported AUC above 0.5. It will report:

- `auc_raw`, on the observed scale and direction;
- `auc_strength = 0.5 + abs(auc_raw - 0.5)`;
- `effect_direction`, equal to `"higher_in_positive"`, `"lower_in_positive"`, or `"none"`; and
- two-sided inference based on departure from `AUC = 0.5`.

The default rank is descending `auc_strength`. Statistical evidence remains visible through `p_value` and `q_value`; it does not replace the effect-size ranking. Tied strengths receive tied ranks.

## Public API

### Matrix screening

The primary function is:

```r
fit <- aucmat(
  X,
  y,
  positive = "disease",
  ci = "delong",
  conf_level = 0.95,
  adjust = "BH",
  na_action = "featurewise"
)
```

`X` must be a numeric matrix or data frame with unique, nonempty column names. `y` must contain exactly two observed outcome levels. For logical outcomes, `TRUE` is positive. For numeric `0/1` outcomes, `1` is positive. Character and factor outcomes require an explicit `positive` value. Missing outcomes are removed once, before all feature analyses, and their count is reported. Infinite biomarker values cause an error naming the affected columns; `NA` and `NaN` follow the missing-data policy.

The supported `ci` values are `"delong"`, `"bootstrap"`, and `"none"`. The supported adjustment methods are `"BH"`, `"holm"`, `"bonferroni"`, and `"none"`.

### Result object

`aucmat()` returns an object of class `aucmat_screen`. It is a list containing:

- `results`: one row per biomarker;
- `sample_summary`: global outcome and exclusion counts;
- `settings`: resolved analysis options;
- `call`: the matched call; and
- `feature_metadata`: optional metadata supplied by the user.

The object does not retain the full biomarker matrix by default. An opt-in `retain_data = TRUE` option may retain it for interactive plotting. The default protects memory use and reduces accidental retention of sensitive molecular data.

The `results` table contains at least:

- `biomarker`;
- `n_total`, `n_used`, `n_positive`, and `n_negative`;
- `n_missing` and `missing_fraction`;
- `auc_raw` and `auc_strength`;
- `effect_direction`;
- `std_error`, `conf_low`, and `conf_high`;
- `p_value` and `q_value`;
- `rank`; and
- `status` and `warning`.

Standard methods are:

```r
print(fit)
summary(fit)
as.data.frame(fit)
plot(fit)
subset(fit, q_value < 0.05)
```

`print()` gives a compact analysis summary and leading biomarkers. `summary()` adds class counts, missingness, warning counts, and multiplicity summaries. `as.data.frame()` returns the feature-level results without hidden transformations.

### Selected paired comparisons

The comparison API is deliberately bounded:

```r
compare_auc(fit, X, y, reference = "IL13")
compare_auc(fit, X, y, biomarkers = c("IL13", "CCL17", "CCL22"))
compare_auc(fit, X, y, top_n = 10)
```

All comparisons are paired because columns are measured on the same subjects. Each pair is recomputed on the subjects with observed values for both biomarkers, ensuring a common analysis population. The output reports both common-sample AUCs, their difference, standard error, confidence interval, raw p-value, adjusted p-value, and common sample counts.

The function refuses an unbounded all-pairs comparison when the requested set would create more than a documented safety limit. Users must explicitly select biomarkers or increase the limit. This prevents accidental quadratic computation and uninterpretable multiplicity.

### Rank stability

The stability API is:

```r
stability <- auc_stability(
  X,
  y,
  positive = "disease",
  times = 1000,
  top_k = c(10, 25, 50),
  seed = 20260720
)
```

Each replicate resamples positive and negative subjects separately with replacement, recomputes the complete feature screen, and ranks by `auc_strength`. The result has class `aucmat_stability` and reports:

- bootstrap AUC summaries;
- median rank and rank intervals;
- probability of appearing in each requested top-k set;
- frequency of being the top-ranked biomarker; and
- pairwise co-selection frequencies for the leading biomarkers.

Random-number streams must be reproducible with or without supported parallel execution. The output records the seed, number of successful replicates, failed replicates, and resolved execution settings.

## Missing-data policy

The default is `na_action = "featurewise"`: each biomarker uses all subjects with an observed outcome and observed value for that biomarker. The output makes the resulting sample sizes explicit.

`na_action = "complete"` uses the common set of subjects observed across every selected biomarker. This is appropriate when directly comparable populations are required but may discard substantial data in wide matrices.

The package accepts externally imputed matrices but does not perform imputation internally in version 1. In particular, outcome-dependent or class-specific imputation is excluded because it can manufacture discrimination. Documentation must warn that imputation, filtering, and normalization should be defined before outcome analysis and included inside resampling whenever they are learned from data.

## Component architecture

The implementation will be divided into small internal components:

1. **Input normalization** validates `X`, `y`, feature names, the positive class, and missingness policy.
2. **Matrix AUC engine** computes feature-wise rank statistics, tie corrections, counts, and raw AUCs.
3. **Inference engine** computes DeLong or bootstrap uncertainty and tests.
4. **Multiplicity layer** adjusts valid p-values and preserves missing results.
5. **Comparison engine** recomputes bounded paired comparisons on common samples.
6. **Stability engine** performs stratified resampling, aggregates rank distributions, and manages reproducible parallel execution.
7. **S3 layer** owns user-facing classes and standard methods.
8. **Visualization layer** converts result objects into `ggplot2` objects without recomputing statistics.
9. **Simulation and validation layer** generates controlled correlated biomarker matrices and supports tests and manuscript replication.

The matrix AUC engine must not depend on plotting or S3 presentation code. Inference must consume documented engine outputs rather than reach into user-facing objects. This separation permits optimization without changing the public API.

## Data flow

For `aucmat(X, y)`:

1. validate dimensions, names, types, outcome levels, and finite values;
2. remove observations with missing outcomes and record exclusions;
3. establish feature-wise or global complete-case masks;
4. calculate counts, rank statistics, ties, and raw AUCs;
5. calculate uncertainty and two-sided p-values where estimable;
6. calculate direction, strength, multiplicity adjustment, and rank;
7. attach feature warnings and global settings; and
8. return an `aucmat_screen` object.

No statistical result is calculated inside a plotting method.

## Errors, warnings, and feature status

Global input failures stop immediately with actionable messages. These include mismatched dimensions, a nonbinary outcome, missing or duplicated feature names, nonnumeric biomarker columns, and infinite values.

Feature-specific failures do not abort a large matrix analysis. A feature instead receives a status and `NA` inferential results. Defined statuses include:

- `ok`;
- `constant`;
- `insufficient_positive`;
- `insufficient_negative`;
- `all_missing`; and
- `inference_failed`.

Fewer than two usable observations in either class makes AUC inference unavailable. Small class counts and high missing fractions generate recorded feature warnings rather than thousands of console messages. `summary()` consolidates these warnings.

## Visualization scope

Version 1 contains five focused plotting functions:

- `plot_auc_rank()` shows ordered discrimination strengths and optional intervals;
- `plot_auc_volcano()` shows effect magnitude against adjusted statistical evidence;
- `plot_auc_forest()` shows selected raw AUCs and confidence intervals;
- `plot_auc_stability()` shows rank distributions and top-k probabilities; and
- `plot_roc_top()` shows ROC curves for a deliberately small biomarker set.

Every function returns a `ggplot2` object. Matrix-wide plots must avoid unreadable feature labels by labeling only an explicit or automatically limited subset and reporting the labeling rule.

When the screening object does not retain data, `plot_roc_top()` requires `X` and `y` explicitly. It verifies biomarker names and outcome compatibility before plotting.

## Validation strategy

### Numerical correctness

- Compare small-matrix AUCs, DeLong intervals, and paired tests against `pROC`.
- Test analytically constructed examples with exact empirical AUCs.
- Test complete separation, reversed effects, ties, constants, duplicated values, and imbalance.
- Verify feature-wise and complete-case counts by hand-constructed fixtures.

### Statistical behavior

- Evaluate bias and interval coverage across sample sizes, prevalences, AUCs, ties, and distribution shapes.
- Verify type I error under null biomarkers.
- Verify false-discovery-rate behavior under mixtures of null and non-null correlated biomarkers.
- Evaluate rank and top-k stability under highly correlated signals.
- Treat simulation as software validation rather than as a claim of a new statistical estimator.

### Software behavior

- Unit-test every public method and internal statistical boundary.
- Snapshot compact printing and warnings without snapshotting volatile numerical noise.
- Verify deterministic serial and parallel bootstrap results for fixed seeds.
- Test that plotting functions return valid `ggplot2` objects.
- Run `R CMD check --as-cran` with zero errors, warnings, and avoidable notes before release.

### Performance

The reference benchmark is 10,000 biomarkers by 500 samples on an ordinary laptop. The base DeLong screen should complete without excessive memory use and should materially outperform a straightforward R loop calling `pROC` once per biomarker. Benchmarks will also cover 1,000 and 50,000 biomarkers, with runtime, peak memory, hardware, R version, and package versions recorded.

Bootstrap stability is explicitly optional and parallelizable. Documentation will distinguish base-screen runtime from resampling runtime.

## Real-data application

The principal application will use an atopic-dermatitis molecular matrix, preferably a public or distributable protein panel. The workflow will:

1. document preprocessing without outcome-dependent transformation;
2. screen all biomarkers;
3. adjust for multiplicity;
4. inspect missingness and quality flags;
5. estimate rank stability;
6. compare a small set of correlated Th2 biomarkers; and
7. show why the largest observed AUC is not necessarily the most reproducible biomarker.

If the current cohort cannot legally be distributed, the manuscript will use a public primary dataset and reserve the private cohort for a separate applied publication.

## JSS manuscript contribution

The paper will be organized around high-dimensional biomarker screening rather than around a catalog of unrelated ROC methods:

1. motivation and failure modes of one-biomarker-at-a-time workflows;
2. direction-preserving AUC definitions and simultaneous inference;
3. matrix computation and result-object design;
4. paired comparison and rank stability;
5. numerical and statistical validation;
6. runtime and memory benchmarks;
7. the atopic-dermatitis application;
8. comparison with existing R workflows; and
9. limitations and planned extensions.

The submission will include a CRAN-ready package, a JSS-style manuscript, and a single commented replication script that reproduces all reported results in a practical runtime.

## Migration from the current package

The current `tableroc()` behavior will be mapped onto the new engine or deprecated with a documented transition path. The existing outcome-dependent imputation options will not be part of the new primary workflow. Existing AUC-vector and correlated-biomarker simulation functions will be retained only after their contracts and numerical properties are tested against the new definitions.

`plot_roc_with_combos()` will not define the version 1 architecture because it fits and evaluates multivariable models on the same data. It may remain temporarily as a deprecated legacy function, but it will not be used in the JSS manuscript or recommended workflow.

## Delivery phases and acceptance criteria

### Phase 1: API and statistical contract

- `aucmat()` input rules and result schema are implemented and documented.
- Raw AUC, strength, direction, missingness, and feature statuses are unambiguous.
- S3 printing and conversion methods work on small fixtures.

### Phase 2: Inference and multiplicity

- DeLong and stratified-bootstrap inference match reference calculations.
- Adjusted p-values preserve missing and failed results.
- Tie, constant-feature, and imbalance tests pass.

### Phase 3: comparisons and stability

- Bounded paired comparisons use common subjects.
- Bootstrap rank summaries and top-k probabilities are reproducible.
- Serial and supported parallel runs agree for fixed seeds.

### Phase 4: visualization

- The five approved plot types return customizable `ggplot2` objects.
- Wide-matrix labeling remains readable and deterministic.

### Phase 5: validation and benchmarking

- Numerical reference tests, statistical simulations, and scale benchmarks are reproducible.
- The reference 10,000-by-500 workflow runs within documented resource bounds.

### Phase 6: release and manuscript

- Package checks pass on supported platforms.
- Public documentation and vignettes describe the focused workflow.
- The standalone replication script reproduces manuscript results.
- The package is submitted to CRAN before or alongside the JSS submission.

## Success criteria

The design succeeds when a researcher can provide a molecular biomarker matrix and binary outcome, obtain statistically transparent matrix-wide AUC results, identify unstable rankings, and reproduce selected paired comparisons without writing feature-wise loops or silently introducing outcome leakage. The JSS contribution must remain understandable as one coherent software system rather than seven loosely connected modules.
