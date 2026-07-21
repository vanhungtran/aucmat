# ==============================================================================
# aucmat() — Matrix-first biomarker screening with direction-preserving AUC
# ==============================================================================

#' Screen a biomarker matrix against a binary outcome
#'
#' Computes direction-preserving raw AUC, discrimination strength, and
#' inferential statistics for every biomarker in a numeric matrix against
#' a binary outcome.  The package will never silently reverse a biomarker
#' to force its reported AUC above 0.5.
#'
#' @param X Numeric matrix or data.frame. Samples in rows, biomarkers in
#'   columns.  Column names must be unique and non-empty.
#' @param y Binary outcome vector.  Logical (`TRUE` = positive), numeric
#'   `0/1` (`1` = positive), factor, or character with exactly two
#'   observed levels.
#' @param positive Positive class label.  Required when `y` is a character
#'   or factor with ambiguous ordering.  Ignored for logical and numeric
#'   `0/1` outcomes.
#' @param ci Confidence interval method: `"delong"` (default),
#'   `"bootstrap"`, or `"none"`.
#' @param conf_level Confidence level in (0, 1).  Default 0.95.
#' @param adjust Multiplicity adjustment: `"BH"` (default), `"holm"`,
#'   `"bonferroni"`, or `"none"`.
#' @param na_action `"featurewise"` (default) uses all subjects with
#'   observed values for each biomarker independently.  `"complete"` uses
#'   only subjects with complete observations across all biomarkers.
#' @param boot_n Number of bootstrap replicates when `ci = "bootstrap"`.
#' @param seed Optional integer seed for bootstrap reproducibility.
#' @param retain_data If `TRUE`, retain `X` and `y` in the result object
#'   for interactive plotting.  Default `FALSE`.
#' @param feature_metadata Optional data.frame with metadata about
#'   biomarkers (rownames must match column names of `X`).
#'
#' @return An object of class `aucmat_screen`, a list with components:
#'   `results`, `sample_summary`, `settings`, `call`, and optionally
#'   `feature_metadata`, `X`, `y`.
#'
#' @export
#' @importFrom stats p.adjust pnorm qnorm sd var complete.cases
#' @importFrom rlang .data
#'
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(100 * 20), nrow = 100)
#' colnames(X) <- paste0("bm", 1:20)
#' y <- rep(c(0, 1), each = 50)
#' fit <- aucmat(X, y)
#' print(fit)
#' head(as.data.frame(fit))
aucmat <- function(X,
                   y,
                   positive       = NULL,
                   ci             = c("delong", "bootstrap", "none"),
                   conf_level     = 0.95,
                   adjust         = c("BH", "holm", "bonferroni", "none"),
                   na_action      = c("featurewise", "complete"),
                   boot_n         = 2000,
                   seed           = NULL,
                   retain_data    = FALSE,
                   feature_metadata = NULL) {

  ci        <- match.arg(ci)
  adjust    <- match.arg(adjust)
  na_action <- match.arg(na_action)

  # ---- 1. Validate and normalise inputs ----
  if (!is.matrix(X) && !is.data.frame(X))
    stop("X must be a numeric matrix or data.frame.")
  X <- as.matrix(X)
  if (!is.numeric(X))
    stop("All columns of X must be numeric.")
  if (nrow(X) == 0L)
    stop("X has 0 rows.")
  cn <- colnames(X)
  if (is.null(cn) || anyDuplicated(cn) || any(cn == ""))
    stop("Column names of X must be unique and non-empty.")
  if (any(is.infinite(X)))
    stop("X contains infinite values. Columns: ",
         paste(cn[apply(X, 2, function(v) any(is.infinite(v)))], collapse = ", "))

  if (length(y) != nrow(X))
    stop("Length of y (", length(y), ") must match nrow(X) (", nrow(X), ").")

  # Normalise outcome
  y_norm <- .normalize_binary_y(y, positive)

  # Remove observations with missing outcomes
  keep_y <- !is.na(y_norm)
  if (!all(keep_y)) {
    X <- X[keep_y, , drop = FALSE]
    y_norm <- y_norm[keep_y]
  }
  if (nrow(X) == 0L)
    stop("All rows removed because y was entirely NA.")

  pos <- y_norm == levels(y_norm)[2L]
  neg <- !pos
  n_total <- nrow(X)

  if (sum(pos) < 1L || sum(neg) < 1L)
    stop("y must contain at least one positive and one negative observation.")

  # ---- 2. Missing-data masks ----
  if (na_action == "complete") {
    complete_mask <- stats::complete.cases(X)
    if (!all(complete_mask)) {
      X <- X[complete_mask, , drop = FALSE]
      pos <- pos[complete_mask]
      neg <- neg[complete_mask]
    }
  }

  n_missing <- colSums(is.na(X))
  missing_fraction <- n_missing / n_total

  # ---- 3. Matrix AUC computation ----
  results <- .compute_matrix_auc(X, pos, neg)

  results$n_total          <- n_total
  results$n_missing        <- n_missing
  results$missing_fraction <- missing_fraction

  # ---- 4. Inference ----
  if (!is.null(seed)) set.seed(seed)

  results <- .apply_inference(results, X, pos, neg,
                               ci = ci, conf_level = conf_level,
                               boot_n = boot_n)

  # ---- 5. Multiplicity adjustment ----
  results$q_value <- .adjust_pvalues(results$p_value, method = adjust)

  # ---- 6. Rank by auc_strength descending ----
  ord <- order(results$auc_strength, decreasing = TRUE, na.last = TRUE)
  results <- results[ord, , drop = FALSE]
  results$rank <- seq_len(nrow(results))
  rownames(results) <- NULL

  # ---- 7. Feature warnings ----
  results$warning <- NA_character_
  lo_pos <- results$n_pos < 10L & results$status == "ok"
  lo_neg <- results$n_neg < 10L & results$status == "ok"
  results$warning[lo_pos & lo_neg] <- "small_class_counts"
  results$warning[lo_pos & !lo_neg] <- "small_positive"
  results$warning[!lo_pos & lo_neg] <- "small_negative"
  hi_miss <- results$missing_fraction > 0.5 & results$status == "ok"
  results$warning[hi_miss] <- ifelse(
    is.na(results$warning[hi_miss]),
    "high_missingness",
    paste(results$warning[hi_miss], "high_missingness", sep = ";")
  )

  # ---- 8. Assemble result object ----
  sample_summary <- list(
    n_total     = n_total,
    n_used      = sum(keep_y),
    n_excluded  = sum(!keep_y),
    n_positive  = sum(pos),
    n_negative  = sum(neg),
    positive_class = levels(y_norm)[2L],
    negative_class = levels(y_norm)[1L]
  )

  settings <- list(
    ci         = ci,
    conf_level = conf_level,
    adjust     = adjust,
    na_action  = na_action,
    boot_n     = boot_n
  )

  out <- list(
    results          = results,
    sample_summary   = sample_summary,
    settings         = settings,
    call             = match.call(),
    feature_metadata = feature_metadata
  )

  if (isTRUE(retain_data)) {
    out$X <- X
    out$y <- y_norm
  }

  class(out) <- "aucmat_screen"
  out
}
