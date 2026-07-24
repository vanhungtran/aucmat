# BeatMG Olink proteomics: baseline (pre-treatment)

Plasma Olink NPX (Normalized Protein eXpression) values for 730 proteins
(Inflammation + Inflammation_II panels), measured at the pre-treatment
baseline visit of BeatMG (NeuroNEXT NN102), a randomized, double-blind,
placebo-controlled trial of rituximab for myasthenia gravis. Rituximab
is an anti-CD20 monoclonal antibody that depletes B cells.

## Usage

``` r
beatmg_baseline
```

## Format

A data frame with 48 rows (subjects) and 732 columns: `SampleID`
(numeric subject identifier), `arm` (factor, `"Placebo"` or
`"Rituximab"`, the randomized arm), and 730 further columns of NPX
values, one per protein, named by Olink `Assay` (gene symbol).

## Source

Olink Explore NPX export and trial demographics for BeatMG (NeuroNEXT
NN102). Quality-controlled: Olink `CONTROL` samples removed,
sample-level `QC_Warning == "PASS"` required, assay-level
`Assay_Warning == "EXCLUDED"` removed (assays flagged `WARN` are
retained). Build script: `data-raw/beatmg.R`.

## Details

At baseline, randomization implies no true association is expected
between the biomarker panel and treatment arm; this dataset is useful as
a real-data negative control / Type-I-error check alongside
[beatmg_ontreatment](https://vanhungtran.github.io/aucmat/reference/beatmg_ontreatment.md).

## See also

[beatmg_ontreatment](https://vanhungtran.github.io/aucmat/reference/beatmg_ontreatment.md)
for the on-treatment comparison,
[`aucmat()`](https://vanhungtran.github.io/aucmat/reference/aucmat.md)
for screening.
