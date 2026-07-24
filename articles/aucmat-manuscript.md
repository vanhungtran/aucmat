# aucmat: Matrix-First ROC Analysis

## Introduction

Receiver operating characteristic (ROC) analysis is the dominant
framework for evaluating continuous biomarkers against binary outcomes
in diagnostic medicine, prognostic modeling, and molecular screening
(Hanley and McNeil 1982; DeLong et al. 1988; Pepe 2004). The area under
the ROC curve (AUC) summarizes the probability that a randomly chosen
positive subject has a higher biomarker value than a randomly chosen
negative subject. Its interpretability – 0.5 for chance, 1.0 for perfect
separation – has made it the default metric in omics-scale biomarker
discovery (Robin et al. 2011).

Modern molecular profiling technologies routinely measure thousands of
features on a few hundred subjects. A proteomics panel may contain 5,000
proteins; a metabolomics screen 2,000 metabolites; a transcriptomics
experiment 20,000 genes. The core statistical task – “which biomarkers
discriminate cases from controls?” – requires computing and testing an
AUC for each feature, adjusting for multiplicity, comparing top
candidates, and assessing the stability of the rankings under sampling
variability. Despite decades of methodological development, no
widely-used R package integrates these steps into a single, reproducible
workflow.

The existing R ecosystem provides excellent point solutions. **pROC**
(Robin et al. 2011) is the gold standard for single-biomarker ROC
inference: it computes DeLong confidence intervals using the Sun & Xu
(2014) O(N log N) algorithm (Sun and Xu 2014), tests paired or unpaired
AUC differences via `roc.test()`, handles partial AUC with bootstrap
inference, provides binormal and density-smoothed ROC curves, computes
covariance between paired AUCs via
[`cov()`](https://rdrr.io/r/stats/cor.html), and supports multiclass AUC
via Hand & Till (Hand and Till 2001). However, pROC operates one
biomarker at a time – screening a matrix of 10,000 features requires an
explicit per-column loop. **ROCR** (Sing et al. 2005) offers flexible
visualization with 28 performance measures and cross-validation
averaging but provides no statistical inference (no confidence
intervals, no hypothesis tests for AUC differences). **precrec** (Saito
and Rehmsmeier 2017) adds precision-recall curves and fast C++-based AUC
computation, with confidence bands for multiple test sets and functions
for partial AUC (`part()`, `pauc()`) and AUC confidence intervals
(`auc_ci()`), but provides no DeLong inference, no multiplicity
adjustment, and no paired comparison tests. **colAUC** from caTools
(Tuszynski and Dietze 2024) computes column-wise AUCs efficiently via
the Wilcoxon rank-sum test and supports multiclass AUC via column means
– but provides no inference whatsoever (no standard errors, confidence
intervals, or p-values). **dtComb** (Yerlitas et al. 2025) provides 142
methods for combining two diagnostic markers into a composite score,
with 30+ cutoff optimization methods, resampling, and caret integration
– but is limited to exactly two biomarkers and does not perform
matrix-wide screening. **cancerclass** (Budczies et al. 2024) builds
validated classifiers from high-dimensional molecular data using
nearest-centroid prediction with multiple random validation, but does
not offer DeLong-based screening or bootstrap stability analysis.

aucmat fills this gap. Its design philosophy is *matrix-first*: the
primary function `aucmat(X, y)` accepts a numeric matrix \mathbf{X} \in
\mathbb{R}^{n \times p} and a binary outcome vector \mathbf{y} \in
\\0,1\\^n, and returns a complete screening table with point estimates,
standard errors, confidence intervals, raw and multiplicity-adjusted
p-values, and diagnostic metadata for every biomarker in a single call.
Downstream functions for comparison, stability, cross-validation, panel
construction, and power analysis consume the same S3 object, preserving
outcome encoding and inferential settings throughout.

## Methods

### Direction-preserving AUC

For biomarker X and binary outcome Y \in \\0,1\\, the AUC equals the
Mann-Whitney U probability: \text{AUC} = P(X_i \> X_j \mid Y_i = 1, Y_j
= 0) (Hanley and McNeil 1982). Most ROC software reports
\max(\text{AUC}, 1 - \text{AUC}), silently reversing the biomarker to
force the reported value above 0.5. aucmat deliberately preserves the
observed direction, reporting three values per biomarker:
\text{AUC}\_{\text{raw}} (the Mann-Whitney estimate on the original
numeric scale, which may be below 0.5), \text{AUC}\_{\text{strength}} =
0.5 + \|\text{AUC}\_{\text{raw}} - 0.5\| (for ranking), and
`effect_direction` (`higher_in_positive` or `lower_in_positive`). This
distinction matters biologically: a tumor suppressor gene with
\text{AUC}\_{\text{raw}} = 0.15 discriminates just as strongly as an
oncogene with \text{AUC}\_{\text{raw}} = 0.85, but with the opposite
biological interpretation.

### DeLong variance and joint covariance

The DeLong variance estimator (DeLong et al. 1988; Sun and Xu 2014) uses
placement values V\_{10}(X_i) = n_0^{-1} \sum\_{j:Y_j=0}
\[\mathbf{1}(X_j \< X_i) + 0.5 \cdot \mathbf{1}(X_j = X_i)\] for
positive subjects and V\_{01}(X_j) = n_1^{-1} \sum\_{i:Y_i=1}
\[\mathbf{1}(X_i \> X_j) + 0.5 \cdot \mathbf{1}(X_i = X_j)\] for
negative subjects. Let S\_{10} = (n_1 - 1)^{-1} \sum\_{i \in
\mathcal{I}\_1} (V\_{10}(X_i) - \bar{V}\_{10})^2 and S\_{01} = (n_0 -
1)^{-1} \sum\_{j \in \mathcal{I}\_0} (V\_{01}(X_j) - \bar{V}\_{01})^2 be
the sample variances of the placement values. The DeLong variance
estimator is \widehat{\text{Var}}(\widehat{\text{AUC}}) = S\_{10}/n_1 +
S\_{01}/n_0. For p biomarkers on common subjects, the p \times p joint
covariance matrix \widehat{\mathbf{\Sigma}} =
n_1^{-1}\widehat{\text{Cov}}(\mathbf{V}\_{10}) +
n_0^{-1}\widehat{\text{Cov}}(\mathbf{V}\_{01}) supports an omnibus Wald
test of H_0: \text{AUC}\_1 = \cdots = \text{AUC}\_p. Let \mathbf{C} be a
(p-1) \times p contrast matrix (default: comparing biomarkers 2,\ldots,p
against biomarker 1). The Wald statistic W =
(\mathbf{C}\hat{\boldsymbol{\theta}})^\top
(\mathbf{C}\widehat{\mathbf{\Sigma}}\mathbf{C}^\top)^{-1}
(\mathbf{C}\hat{\boldsymbol{\theta}}) \sim \chi^2\_{r} uses the
effective rank r of the contrast covariance matrix, determined via
eigenvalue tolerance. A generalized inverse is used on the stable
subspace only; rank-zero or numerically unstable cases return
`status = "non_testable_singular"`. For p = 2, W equals the squared
paired-DeLong z-statistic within numerical tolerance, serving as an
internal consistency check.

### Multiplicity adjustment and stability

Screening p biomarkers produces p hypothesis tests. aucmat applies
Benjamini-Hochberg (Benjamini and Hochberg 1995), Holm, or Bonferroni
correction via
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), reporting
adjusted p-values as q-values. For equivalence (TOST) (Schuirmann 1987),
the reported p\_{\text{equiv}} = \max(p\_{\text{lower}},
p\_{\text{upper}}) and multiplicity adjustment applies to this combined
value.

Bootstrap rank-stability analysis
([`auc_stability()`](https://vanhungtran.github.io/aucmat/reference/auc_stability.md))
resamples the data B times with replacement (stratified by outcome),
recomputes the full screening for each replicate, and reports median/IQR
rank distributions, top-k selection probabilities, and pairwise
co-selection frequencies – quantifying how reproducible the ranking is
under sampling variability (Efron and Tibshirani 1994).

### Hurdle-AUC for zero-inflated data

Standard AUC assumes continuous biomarkers following a shifted normal
distribution per class. This assumption breaks on zero-inflated data
(scRNA-seq, microbiome, detection-limited assays) where 60%+ of values
are exactly zero. Under the standard model, zeros would need to arise
from the tail of a normal distribution below zero – biologically
implausible. The Hurdle-AUC uses a two-stage model:

**Stage 1 – Zero hurdle:** P(X_j = 0 \mid Y) = \text{logistic}(\beta_0 +
\beta_1 Y) models the binary expression decision via logistic
regression.

**Stage 2 – Magnitude:** Standard AUC on the non-zero values captures
magnitude discrimination among expressed observations.

The composite score S_i for observation i is: S_i = \hat{p}\_i if
X\_{ij} = 0, and S_i = 0.5 \cdot \hat{p}\_i + 0.5 \cdot \tilde{x}\_{ij}
if X\_{ij} \> 0, where \hat{p}\_i is the predicted probability of being
non-zero from Stage 1, and \tilde{x}\_{ij} is the normalized biomarker
value.

Simulation confirms recovery from standard AUC \approx 0.53 to
Hurdle-AUC \approx 0.90 on zero-inflated biomarkers by acknowledging
that “expression on/off” and “expression level when on” are separate
biological processes.

### Multivariable panel construction

Screening ranks biomarkers individually; combining several into one
score is a distinct problem, addressed by
[`fit_auc_panel()`](https://vanhungtran.github.io/aucmat/reference/fit_auc_panel.md).
Four of its six methods are ridge, lasso, elastic-net, and unpenalized
logistic regression on centred and scaled predictors, the penalized
variants fit via glmnet (Friedman et al. 2010). A fifth, `su_liu`, is
the closed-form linear combination that maximizes AUC under the
common-covariance multivariate normal model (Su and Liu 1993): for
class-conditional means \boldsymbol{\mu}\_1, \boldsymbol{\mu}\_0 and
pooled within-class covariance \boldsymbol{\Sigma}, the AUC-maximizing
direction is \mathbf{w} = \boldsymbol{\Sigma}^{-1}
(\boldsymbol{\mu}\_1 - \boldsymbol{\mu}\_0) – the same direction as
Fisher’s linear discriminant, and a natural counterpart to the binormal
model already used by aucmat’s simulation engine. Because
\boldsymbol{\Sigma} is frequently ill-conditioned when biomarkers are
collinear or p approaches n, the pooled covariance is shrunk toward its
diagonal, \boldsymbol{\Sigma}\_\lambda =
(1-\lambda)\boldsymbol{\Sigma} +
\lambda\\\mathrm{diag}(\boldsymbol{\Sigma}) for a user-set \lambda \in
\[0,1\], with a reciprocal-condition-number guard that fails
informatively rather than returning an unstable solution. A sixth,
`unweighted`, sums direction-aligned, standardized biomarkers with equal
weight and fits no parameters at all – an assumption-light baseline for
judging whether the fitted methods earn their added complexity.

All six methods share one **nested** cross-validation procedure for the
reported out-of-fold AUC: centering, scaling, and (for the penalized
methods) penalty selection are re-derived from each outer fold’s
training portion alone, including for the final full-data model’s own
internal penalty search, so nothing derived from a held-out fold ever
contributes to the model evaluated on it.
[`compare_panel_auc()`](https://vanhungtran.github.io/aucmat/reference/compare_panel_auc.md)
reuses the DeLong or bootstrap paired-comparison engine from
[`compare_auc()`](https://vanhungtran.github.io/aucmat/reference/compare_auc.md)
to test the panel’s honest cross-validated score against each individual
biomarker used to build it, on the same subjects, and
[`plot_roc_panel()`](https://vanhungtran.github.io/aucmat/reference/plot_roc_panel.md)
plots the corresponding ROC curves together.

### Multivariate simulation

Under the equal-variance binormal model, X_k \mid Y \sim
\mathcal{N}(\mu\_{kY}, \sigma^2) with \delta_k = (\mu\_{1k} -
\mu\_{0k})/\sigma = \sqrt{2} \cdot \Phi^{-1}(\text{AUC}\_k). For p
biomarkers with correlation matrix \mathbf{R}, the class-conditional
distribution is \mathbf{X} \mid Y \sim
\mathcal{N}\_p(\boldsymbol{\mu}\_Y, \mathbf{R}), where
\boldsymbol{\mu}\_0 and \boldsymbol{\mu}\_1 are centred so the
n-weighted grand mean is zero.
[`simulate_auc_matrix()`](https://vanhungtran.github.io/aucmat/reference/simulate_auc_matrix.md)
generates data under this model with four parametrized correlation
structures (user, exchangeable, AR(1), block) and feasibility
diagnostics with optional nearest positive-definite projection.

[`simulate_auc_copula()`](https://vanhungtran.github.io/aucmat/reference/simulate_auc_copula.md)
uses a two-phase approach: Phase 1 applies a Gaussian copula
(rank-transform to uniform, probit to normal scores) that preserves the
target correlation structure; Phase 2 applies iterative mean
perturbation with decaying step size to fine-tune each biomarker’s AUC.
This separates correlation control (copula) from discriminative signal
(perturbation), avoiding the AUC-correlation tradeoff inherent in
sequential binormal methods.

[`simulate_hurdle_auc()`](https://vanhungtran.github.io/aucmat/reference/simulate_hurdle_auc.md)
implements the two-stage hurdle model: Bernoulli zeros with
class-specific rates, and log-normal magnitudes for expressed values
with binormal-calibrated mean separation, clipped to \[0, \infty).

[`validate_simulation()`](https://vanhungtran.github.io/aucmat/reference/validate_simulation.md)
repeats a specification across independent draws and reports bias, RMSE,
Monte Carlo standard error, and target-interval hit rates for both AUCs
and pairwise correlations – a single draw is never evidence that a
simulator is calibrated.

## Implementation

aucmat is implemented in pure R (~6,900 lines across 37 source files)
with minimized dependencies. The primary function
[`aucmat()`](https://vanhungtran.github.io/aucmat/reference/aucmat.md)
performs six steps in sequence: (1) input validation (matrix dimensions,
column names, outcome levels, infinite values); (2) outcome
normalization via `.normalize_binary_y()` (accepting logical, numeric
0/1, factor, or character outcomes); (3) matrix-wide AUC computation via
the Mann-Whitney rank formulation; (4) DeLong or bootstrap inference via
`.apply_inference()` with configurable alternative hypothesis
(`two.sided`, `greater`, `less`); (5) multiplicity adjustment via
`.adjust_pvalues()`; and (6) result assembly with ranking by
discrimination strength and feature-level warnings for small class
counts or high missingness.

All primary functions return S3 objects with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), and (where
appropriate)
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`subset()`](https://rdrr.io/r/base/subset.html) methods. The package
stores outcome encoding, positive-class designation, and sample masks in
the returned object, eliminating the need for downstream functions to
re-infer factor direction or re-derive class counts.

Structured error classes (inheriting from `error`) allow programmatic
handling of common failures: `aucmat_invalid_outcome`,
`aucmat_invalid_correlation`, `aucmat_infeasible_targets`,
`aucmat_insufficient_sample`, and `aucmat_resampling_failure`. An
RNG-safety mechanism (`with_seed()`) saves the global `.Random.seed`,
executes the computation, and restores the previous state on exit,
including when `.Random.seed` was previously absent and when an error
occurs.

The package imports pROC (ROC curve objects and smoothed ROC), mvtnorm
(multivariate normal generation), ggplot2 (visualization), and rlang
(tidy evaluation). Optional Suggests include GGally (pairs plots for
simulated data), glmnet (penalized panel models), matrixStats (fast
column-wise ranks), knitr and rmarkdown (vignette building), and ggrepel
(non-overlapping plot labels). The complete source, documentation, and
four vignettes are available at <https://vanhungtran.github.io/aucmat>.

## Usage

### Basic screening

``` r

library(aucmat)

sim <- simulate_auc_matrix(
  n = 500, prevalence = 0.3,
  target_aucs = c(0.9, 0.8, 0.7, 0.65),
  correlation = 0.3, structure = "exchangeable"
)
X <- as.matrix(sim$data[, 1:4])
y <- sim$data$truth

fit <- aucmat(X, y, ci = "delong", adjust = "BH")
print(fit)
summary(fit)
plot_auc_rank(fit, n_label = 4)
```

### Paired comparisons with hypothesis testing

``` r

# Two-sided superiority (default: H0: delta = 0)
compare_auc(fit, X, y, reference = "X1")

# Non-inferiority: H0: AUC_a <= AUC_b - margin
compare_auc(fit, X, y, reference = "X1",
  hypothesis = "noninferiority", margin = 0.05)

# Equivalence (TOST): H0: |AUC_a - AUC_b| >= margin
compare_auc(fit, X, y, reference = "X1",
  hypothesis = "equivalence", margin = 0.15, adjust = "BH")

# Omnibus Wald test: H0: AUC_1 = AUC_2 = ... = AUC_p
compare_auc_global(fit, X, y, biomarkers = paste0("X", 1:4))

# Single-biomarker DeLong test against chance
roc_test(X[, 1], y)
```

### Hurdle-AUC for zero-inflated data

``` r

# scRNA-seq: 55% zeros in controls, 25% in cases
sim <- simulate_hurdle_auc(
  n = 400, prevalence = 0.3,
  target_hurdle_aucs = c(0.85, 0.72),
  zero_rate_neg = c(0.55, 0.80),
  zero_rate_pos = c(0.25, 0.70)
)
fit_hur <- hurdle_auc(as.matrix(sim$data[, 1:2]), sim$data$truth)
print(fit_hur)
# Standard AUC ~0.53 recovers to Hurdle AUC ~0.90
```

### Cross-validation, panel construction, power analysis

``` r

# Cross-validated AUC (honest out-of-fold estimates)
cv <- cv_aucmat(X, y, n_folds = 5, seed = 42)

# Combine top biomarkers into a panel -- six methods, all honestly
# nested-cross-validated
top4 <- X[, head(cv$results$biomarker, 4)]
panel    <- fit_auc_panel(top4, y, method = "ridge",      n_folds = 5)
panel_su <- fit_auc_panel(top4, y, method = "su_liu",     n_folds = 5)
panel_uw <- fit_auc_panel(top4, y, method = "unweighted", n_folds = 5)

# Does the panel actually beat its individual components?
compare_panel_auc(panel, top4, y)
plot_roc_panel(panel, top4, y)

# Sample size for 80% power to detect AUC = 0.70
power_auc_matrix(
  target_aucs = 0.70, power = 0.80,
  prevalence = 0.3, alpha = 0.05
)
```

## Real-Data Application

The examples above use simulated data with known ground truth. aucmat
also ships real biomarker data: plasma Olink proteomics from BeatMG
(NeuroNEXT NN102), a randomized, double-blind, placebo-controlled trial
of rituximab – an anti-CD20 monoclonal antibody that depletes B cells –
for myasthenia gravis. `data(beatmg_baseline)` and
`data(beatmg_ontreatment)` provide 730 Olink Inflammation-panel
proteins, quality-controlled and reshaped into subject-by-protein
matrices with randomized arm as the outcome, at the pre-treatment
baseline visit (n = 48) and the latest well-powered on-treatment
collection (n = 44).

``` r

data(beatmg_baseline)
data(beatmg_ontreatment)
protein_cols <- setdiff(names(beatmg_ontreatment), c("SampleID", "arm"))

fit_baseline <- aucmat(as.matrix(beatmg_baseline[, protein_cols]),
  beatmg_baseline$arm, positive = "Rituximab", ci = "delong", adjust = "BH")
fit_rtx <- aucmat(as.matrix(beatmg_ontreatment[, protein_cols]),
  beatmg_ontreatment$arm, positive = "Rituximab", ci = "delong", adjust = "BH")
```

Randomization implies no true baseline association between the proteome
and treatment arm; screening the pre-treatment matrix is a real-data
negative control. Zero of 730 proteins reach q \< 0.05 at baseline.
Screening on-treatment finds eight, all with AUC below 0.5 –
systematically lower under rituximab:

| Protein            | AUC   | 95% CI          | q       |
|--------------------|-------|-----------------|---------|
| FCRL2              | 0.035 | \[-0.01, 0.08\] | 1.1e-79 |
| CD22               | 0.054 | \[-0.01, 0.12\] | 9.2e-41 |
| TNFRSF13C (BAFF-R) | 0.066 | \[-0.01, 0.14\] | 4.0e-26 |
| IL17F              | 0.072 | \[0.00, 0.14\]  | 5.9e-30 |
| TNFRSF13B (TACI)   | 0.147 | \[0.04, 0.26\]  | 5.7e-08 |
| SPINK2             | 0.159 | \[0.04, 0.28\]  | 1.1e-06 |
| CD79B              | 0.199 | \[0.07, 0.33\]  | 6.7e-04 |
| CD72               | 0.219 | \[0.08, 0.36\]  | 5.6e-03 |

Proteins significant at q \< 0.05, rituximab vs. placebo on-treatment
(Timepoint 3, n = 44). {.table}

Six of the eight (FCRL2, CD22, TNFRSF13C/BAFF-R, TNFRSF13B/TACI, CD79B,
CD72) are canonical B-cell surface or B-cell-survival-signaling
proteins. Finding exactly these markers suppressed in a blinded,
randomized comparison is the expected pharmacodynamic signature of
anti-CD20 B-cell depletion, not a post hoc association fit to noise –
and the null baseline screen rules out a simple technical explanation,
since a batch or plate artifact would appear at every timepoint rather
than emerge only on-treatment. Bootstrap rank-stability analysis
([`auc_stability()`](https://vanhungtran.github.io/aucmat/reference/auc_stability.md),
1,000 replicates) shows FCRL2, CD22, IL17F, and TNFRSF13C are in the top
8 in 98.2-99.8% of resamples, while CD79B and CD72 are markedly less
stable (30.8% and 15.8%) and warrant more caution. This demonstrates
rituximab’s expected pharmacodynamic effect on circulating B-cell
markers; it is not a clinical efficacy result – the shipped datasets
carry randomized arm only, not the trial’s clinical endpoints. Full
validation detail and figures are in the `real-data-beatmg` vignette.

## Benchmarking

**Screening throughput.** On a desktop workstation (AMD Ryzen 3.5 GHz,
32 GB RAM, R 4.5.1, Windows 11), aucmat screens 100, 1,000, and 10,000
biomarkers (n = 500, DeLong inference) in 0.54, 5.02, and 59.76 seconds
respectively. The equivalent per-column pROC loop
([`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html) then
[`pROC::ci.auc()`](https://rdrr.io/pkg/pROC/man/ci.auc.html) with
`method = "delong"`) on the same data requires 0.68, 6.42, and an
estimated 64 seconds. The speed advantage comes from avoiding repeated
outcome encoding and from column-wise rank computation. At small p (\<
100), the difference is negligible; for wide matrices (p \> 1{,}000),
the per-column overhead of
[`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html) and
[`pROC::ci.auc()`](https://rdrr.io/pkg/pROC/man/ci.auc.html) becomes
meaningful.

**Simulation calibration.** Across 27 simulation conditions (n = 100,
300, 500; prevalence averaged over 0.2, 0.3, 0.5; target AUC = 0.65,
0.75, 0.85; 200 replicates each), DeLong 95% confidence interval
coverage ranged from 0.898 (at n = 100, AUC = 0.65) to 0.957 (at n =
500, AUC = 0.65), converging to the nominal 0.95 as sample size
increased. The coverage is slightly below nominal for small n and
extreme AUC values, consistent with the known finite-sample behavior of
the DeLong variance estimator.

| AUC  | n=100 | n=300 | n=500 |
|------|-------|-------|-------|
| 0.65 | 0.898 | 0.937 | 0.957 |
| 0.75 | 0.915 | 0.943 | 0.933 |
| 0.85 | 0.883 | 0.923 | 0.905 |

Coverage of DeLong 95% CI, averaged across 27 conditions (200 reps
each). {.table}

## Comparison with existing R packages

The table below compares aucmat (v0.2.0) with five established R
packages for ROC analysis and biomarker evaluation. Features are
verified against each package’s CRAN documentation and published
literature as of July 2026.

| Feature | aucmat | pROC | precrec | ROCR | colAUC | dtComb |
|:---|:---|:---|:---|:---|:---|:---|
| Matrix-first screening (single call) | \+ | \- | \- | \- | \+ | \- |
| Direction-preserving AUC (\<0.5 reported) | \+ | \- | \- | \- | \- | \- |
| DeLong confidence intervals | \+ | \+ | \- | \- | \- | \- |
| Stratified bootstrap CI | \+ | \+ | \- | \- | \- | \+ |
| Multiplicity adjustment (BH/Holm/Bonf) | \+ | \- | \- | \- | \- | \- |
| One-sided alternative tests | \+ | \+ | \- | \- | \- | \- |
| Paired biomarker comparisons | \+ | \+ | \- | \- | \- | \+ |
| Non-inferiority / equivalence (TOST) | \+ | \- | \- | \- | \- | \- |
| Omnibus Wald test (3+ correlated AUCs) | \+ | \- | \- | \- | \- | \- |
| Bootstrap rank-stability analysis | \+ | \- | \- | \- | \- | \- |
| Partial AUC | \+ | \+ | \+ | \- | \- | \- |
| Multivariate simulation engine | \+ | \- | \- | \- | \- | \- |
| Simulation calibration diagnostics | \+ | \- | \- | \- | \- | \- |
| Cross-validated AUC screening | \+ | \- | \- | \- | \- | \+ |
| Multivariable panel combination (6 methods) | \+ | \- | \- | \- | \- | \+ |
| Honest nested CV for panel scores | \+ | \- | \- | \- | \- | \- |
| Panel vs. component significance test | \+ | \- | \- | \- | \- | \- |
| Hurdle-AUC (zero-inflated data) | \+ | \- | \- | \- | \- | \- |
| Precision-Recall curves | \+ | \- | \+ | \- | \- | \- |
| Power / sample size computation | \+ | \+ | \- | \- | \- | \- |
| Multiclass AUC (Hand-Till) | \+ | \+ | \- | \- | \+ | \- |
| ROC smoothing (binormal/density) | \+ | \+ | \- | \- | \- | \- |
| Group-stratified screening | \+ | \- | \- | \- | \- | \- |

Feature comparison of R packages for ROC analysis and biomarker
screening. \texttt{+} = supported; \texttt{-} = not supported or
requires manual implementation. Verified against CRAN documentation as
of July 2026. pROC supports unpaired comparisons via
\texttt{roc.test(paired=FALSE)}; aucmat requires common subjects for all
comparisons. {.table}

**What aucmat uniquely provides.** Six capabilities are, to our
knowledge, not available in any other single R package: (i) matrix-first
screening with direction-preserving AUC and multiplicity adjustment in
one call; (ii) non-inferiority and equivalence (TOST) tests for paired
AUC comparisons; (iii) an omnibus Wald test for three or more correlated
AUCs via the joint DeLong covariance matrix; (iv) bootstrap
rank-stability analysis with top-k selection probabilities and pairwise
co-selection frequencies; (v) a two-stage Hurdle-AUC model for
zero-inflated data that recovers discriminative signal where standard
AUC fails; and (vi) multivariable panel combination with honestly nested
cross-validation and a dedicated significance test of whether the
combined panel outperforms its individual components.

**When to use what.** For a single biomarker with smoothed ROC curves,
pROC remains the gold standard, offering binormal and density smoothing,
partial AUC with bootstrap inference, and Venkatraman’s permutation test
for whole-curve comparison (Venkatraman and Begg 2000). For speed on
extremely wide matrices (\>100k columns), colAUC’s Wilcoxon-based
column-wise computation is recommended, though it provides no inference.
For precision-recall analysis of highly imbalanced data, precrec’s C++
backend offers fast and accurate PR curves with confidence bands for
multiple test sets. For an exhaustive search over composite scores for
exactly two biomarkers, dtComb’s library of 142 methods covers virtually
all published linear, non-linear, mathematical, and machine-learning
combination approaches; aucmat’s own
[`fit_auc_panel()`](https://vanhungtran.github.io/aucmat/reference/fit_auc_panel.md)
targets a smaller set of statistically transparent combination methods
(penalized regression, a closed-form AUC-maximizing linear combination,
and an assumption-light unweighted baseline) for panels of any size,
always honestly nested-cross-validated, with a built-in test of whether
combining helped. For end-to-end screening of hundreds to tens of
thousands of biomarkers with full statistical inference, multiplicity
adjustment, comparison testing, stability assessment, validation, and
study planning in a single coherent workflow, aucmat is the recommended
tool.

## Discussion

aucmat provides a unified, statistically principled framework for
biomarker screening against binary outcomes in R. Its matrix-first
design reduces the code burden for high-dimensional ROC analysis from
dozens of bespoke scripting steps to a handful of function calls with
consistent S3 interfaces. Cross-validated screening and multivariable
panel models address same-data optimism – the latter via a fully nested
procedure that re-derives every outcome-dependent step, including
preprocessing, from each fold’s training portion alone – while power
analysis tools support prospective study design.

**Limitations.** The current version (0.2.0) is restricted to binary
outcomes for the primary screening workflow (multiclass is available
experimentally via
[`aucmat_multiclass()`](https://vanhungtran.github.io/aucmat/reference/aucmat_multiclass.md)).
The DeLong variance estimator requires common subjects for paired
comparisons; unpaired comparisons between independent cohorts are not
yet implemented (pROC provides this via `roc.test(paired=FALSE)`).
Partial AUC computation is implemented as an internal engine
(`R/partial-auc.R`) but not yet exposed through the primary
[`aucmat()`](https://vanhungtran.github.io/aucmat/reference/aucmat.md)
interface. The package targets in-memory dense matrices; sparse and
on-disk backends are deferred. The DeLong variance estimator shows
slightly liberal coverage at small sample sizes (coverage \approx 0.90
at n = 100 for AUC = 0.65), consistent with its known finite-sample
behavior; bootstrap intervals are recommended when n \< 100 per class.

**Future work.** Version 0.2.0 closed the honest-nested-cross-validation
gap for panel models identified in the previous release and added two
combination methods (a closed-form AUC-maximizing linear combination and
an unweighted baseline) alongside the existing penalized-regression
panels. The development roadmap targets three further extensions:
version 0.3 will extend the framework to multiclass outcomes with
one-versus-rest, one-versus-one, and Hand-Till estimators; version 0.4
will incorporate cumulative/dynamic time-dependent AUCs for censored
survival outcomes (Heagerty et al. 2000; Blanche et al. 2013); and
version 0.5 will introduce sparse-matrix and on-disk backends for
genomics-scale data exceeding available RAM.

## Acknowledgements

We thank Xavier Robin and colleagues for pROC, which provided both the
statistical foundation and the ROC curve objects that aucmat builds
upon. Hadley Wickham’s ggplot2 and the mvtnorm authors (Genz, Bretz, et
al.) are gratefully acknowledged for the visualization and simulation
backends. The Sun & Xu (2014) fast DeLong algorithm implemented in pROC
since version 1.9 informs aucmat’s DeLong variance computation.

Thank KUHNE foundation for financial support.

## References

Benjamini, Yoav, and Yosef Hochberg. 1995. “Controlling the False
Discovery Rate: A Practical and Powerful Approach to Multiple Testing.”
*Journal of the Royal Statistical Society: Series B* 57 (1): 289–300.

Blanche, Paul, Jean-François Dartigues, and Hélène Jacqmin-Gadda. 2013.
“Estimating and Comparing Time-Dependent Areas Under Receiver Operating
Characteristic Curves for Censored Event Times with Competing Risks.”
*Statistics in Medicine* 32 (30): 5381–97.

Budczies, Jan, Daniel Kosztyla, Christian von Törne, et al. 2024.
*cancerclass: Development and Validation of Diagnostic Tests from
High-Dimensional Molecular Data*.

DeLong, Elizabeth R, David M DeLong, and Daniel L Clarke-Pearson. 1988.
“Comparing the Areas Under Two or More Correlated Receiver Operating
Characteristic Curves: A Nonparametric Approach.” *Biometrics* 44 (3):
837–45.

Efron, Bradley, and Robert J Tibshirani. 1994. *An Introduction to the
Bootstrap*. CRC Press.

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.

Hand, David J, and Robert J Till. 2001. “A Simple Generalisation of the
Area Under the ROC Curve for Multiple Class Classification Problems.”
*Machine Learning* 45 (2): 171–86.

Hanley, James A, and Barbara J McNeil. 1982. “The Meaning and Use of the
Area Under a Receiver Operating Characteristic (ROC) Curve.” *Radiology*
143 (1): 29–36.

Heagerty, Patrick J, Thomas Lumley, and Margaret S Pepe. 2000.
“Time-Dependent ROC Curves for Censored Survival Data and a Diagnostic
Marker.” *Biometrics* 56 (2): 337–44.

Pepe, Margaret S. 2004. *The Statistical Evaluation of Medical Tests for
Classification and Prediction*. Oxford University Press.

Robin, Xavier, Natacha Turck, Alexandre Hainard, et al. 2011. “pROC: An
Open-Source Package for R and S+ to Analyze and Compare ROC Curves.”
*BMC Bioinformatics* 12: 77. <https://doi.org/10.1186/1471-2105-12-77>.

Saito, Takaya, and Marc Rehmsmeier. 2017. “precrec: Fast and Accurate
Precision-Recall and ROC Curve Calculations in R.” *Bioinformatics* 33
(1): 145–47. <https://doi.org/10.1093/bioinformatics/btw570>.

Schuirmann, Donald J. 1987. “A Comparison of the Two One-Sided Tests
Procedure and the Power Approach for Assessing the Equivalence of
Average Bioavailability.” *Journal of Pharmacokinetics and
Biopharmaceutics* 15 (6): 657–80.

Sing, Tobias, Oliver Sander, Niko Beerenwinkel, and Thomas Lengauer.
2005. “ROCR: Visualizing Classifier Performance in R.” *Bioinformatics*
21 (20): 3940–41. <https://doi.org/10.1093/bioinformatics/bti623>.

Su, Jack Q, and Jun S Liu. 1993. “Linear Combinations of Multiple
Diagnostic Markers.” *Journal of the American Statistical Association*
88 (424): 1350–55. <https://doi.org/10.1080/01621459.1993.10476417>.

Sun, Xu, and Weichao Xu. 2014. “Fast Implementation of DeLong’s
Algorithm for Comparing the Areas Under Correlated Receiver Operating
Characteristic Curves.” *IEEE Signal Processing Letters* 21 (11):
1389–93. <https://doi.org/10.1109/LSP.2014.2337313>.

Tuszynski, Jarek, and Michael Dietze. 2024. *caTools: Tools: Moving
Window Statistics, GIF, Base64, ROC AUC, Etc.*
<https://doi.org/10.32614/CRAN.package.caTools>.

Venkatraman, E S, and Colin B Begg. 2000. “A Permutation Test to Compare
Receiver Operating Characteristic Curves.” *Biometrics* 56 (4): 1134–38.

Yerlitas, Serra Ilayda, Serra Bersan Gengec, Necla Kochan, Gozde Erturk
Zararsiz, Selcuk Korkmaz, and Gokmen Zararsiz. 2025. “dtComb: A
Comprehensive R Library and Web Tool for Combining Diagnostic Tests.”
*The R Journal* 17 (1): 1–20. <https://doi.org/10.32614/RJ-2025-036>.
