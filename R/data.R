# ==============================================================================
# Documentation for shipped datasets. Objects live in data/*.rda, built by
# data-raw/beatmg.R from the raw BeatMG trial exports.
# ==============================================================================

#' BeatMG Olink proteomics: baseline (pre-treatment)
#'
#' Plasma Olink NPX (Normalized Protein eXpression) values for 730 proteins
#' (Inflammation + Inflammation_II panels), measured at the pre-treatment
#' baseline visit of BeatMG (NeuroNEXT NN102), a randomized, double-blind,
#' placebo-controlled trial of rituximab for myasthenia gravis. Rituximab is
#' an anti-CD20 monoclonal antibody that depletes B cells.
#'
#' At baseline, randomization implies no true association is expected
#' between the biomarker panel and treatment arm; this dataset is useful as
#' a real-data negative control / Type-I-error check alongside
#' [beatmg_ontreatment].
#'
#' @format A data frame with 48 rows (subjects) and 732 columns: `SampleID`
#'   (numeric subject identifier), `arm` (factor, `"Placebo"` or
#'   `"Rituximab"`, the randomized arm), and 730 further columns of NPX
#'   values, one per protein, named by Olink `Assay` (gene symbol).
#'
#' @source Olink Explore NPX export and trial demographics for BeatMG
#'   (NeuroNEXT NN102). Quality-controlled: Olink `CONTROL` samples removed,
#'   sample-level `QC_Warning == "PASS"` required, assay-level
#'   `Assay_Warning == "EXCLUDED"` removed (assays flagged `WARN` are
#'   retained). Build script: `data-raw/beatmg.R`.
#'
#' @seealso [beatmg_ontreatment] for the on-treatment comparison,
#'   [aucmat()] for screening.
"beatmg_baseline"

#' BeatMG Olink proteomics: on-treatment (timepoint 3)
#'
#' Companion to [beatmg_baseline], measured at the latest well-powered
#' on-treatment collection (Olink `Timepoint == 3`; `Timepoint == 4` has
#' only 5 samples after QC and is not shipped). Screening this against
#' [beatmg_baseline] shows rituximab's pharmacodynamic effect on circulating
#' B-cell markers (e.g. `FCRL2`, `CD22`, `TNFRSF13C`, `TNFRSF13B`, `CD79B`,
#' `CD72` are all lower under rituximab, consistent with B-cell depletion)
#' emerging where none was present at baseline.
#'
#' This reflects the trial's randomized biomarker panel only -- it does not
#' include the trial's clinical efficacy endpoints (MG-ADL, QMG,
#' prednisone-sparing outcome), which are not part of this data release.
#'
#' @format A data frame with 44 rows (subjects) and 732 columns; same
#'   structure as [beatmg_baseline].
#'
#' @source See [beatmg_baseline]. Build script: `data-raw/beatmg.R`.
#'
#' @seealso [beatmg_baseline], [aucmat()], [auc_stability()].
"beatmg_ontreatment"
