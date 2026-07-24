# BeatMG Olink proteomics: on-treatment (timepoint 3)

Companion to
[beatmg_baseline](https://vanhungtran.github.io/aucmat/reference/beatmg_baseline.md),
measured at the latest well-powered on-treatment collection (Olink
`Timepoint == 3`; `Timepoint == 4` has only 5 samples after QC and is
not shipped). Screening this against
[beatmg_baseline](https://vanhungtran.github.io/aucmat/reference/beatmg_baseline.md)
shows rituximab's pharmacodynamic effect on circulating B-cell markers
(e.g. `FCRL2`, `CD22`, `TNFRSF13C`, `TNFRSF13B`, `CD79B`, `CD72` are all
lower under rituximab, consistent with B-cell depletion) emerging where
none was present at baseline.

## Usage

``` r
beatmg_ontreatment
```

## Format

A data frame with 44 rows (subjects) and 732 columns; same structure as
[beatmg_baseline](https://vanhungtran.github.io/aucmat/reference/beatmg_baseline.md).

## Source

See
[beatmg_baseline](https://vanhungtran.github.io/aucmat/reference/beatmg_baseline.md).
Build script: `data-raw/beatmg.R`.

## Details

This reflects the trial's randomized biomarker panel only – it does not
include the trial's clinical efficacy endpoints (MG-ADL, QMG,
prednisone-sparing outcome), which are not part of this data release.

## See also

[beatmg_baseline](https://vanhungtran.github.io/aucmat/reference/beatmg_baseline.md),
[`aucmat()`](https://vanhungtran.github.io/aucmat/reference/aucmat.md),
[`auc_stability()`](https://vanhungtran.github.io/aucmat/reference/auc_stability.md).
