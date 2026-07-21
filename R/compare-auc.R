# ==============================================================================
# compare_auc() — Bounded paired AUC comparisons on common subjects
# ==============================================================================

#' Compare AUCs of selected biomarkers (paired, common-subject)
#'
#' Performs DeLong-based paired comparisons of AUC between biomarkers.
#' Each pair is recomputed on the subjects with observed values for both
#' biomarkers, ensuring a common analysis population.
#'
#' The function refuses unbounded all-pairs comparisons.  Users must
#' explicitly select biomarkers or increase `max_pairs`.
#'
#' @param fit An `aucmat_screen` object.
#' @param X Numeric matrix (same as passed to [aucmat()]).
#' @param y Binary outcome (same as passed to [aucmat()]).
#' @param reference Single biomarker name to compare all others against.
#' @param biomarkers Character vector of biomarker names.  If length > 2,
#'   all pairwise comparisons among the set are computed.
#' @param top_n Take the top `top_n` biomarkers by `auc_strength` from
#'   `fit` and compare all pairs.  Requires `fit` to have a `$results`
#'   table.
#' @param max_pairs Safety limit on the number of pairwise comparisons.
#'   Default 100.  An error is raised if more pairs would be computed.
#' @param adjust Multiplicity adjustment for the comparison p-values:
#'   `"BH"`, `"holm"`, `"bonferroni"`, or `"none"` (default).
#'
#' @return A data.frame of class `aucmat_compare` with columns:
#'   `biomarker_a`, `biomarker_b`, `auc_a`, `auc_b`, `auc_diff`,
#'   `std_error`, `conf_low`, `conf_high`, `p_value`, `q_value`,
#'   `n_common`, `n_pos`, `n_neg`.
#'
#' @export
#' @importFrom stats p.adjust pnorm
compare_auc <- function(fit, X, y,
                         reference  = NULL,
                         biomarkers = NULL,
                         top_n      = NULL,
                         max_pairs  = 100,
                         adjust     = c("none", "BH", "holm", "bonferroni")) {

  adjust <- match.arg(adjust)

  # Build the list of biomarker pairs to compare
  if (!is.null(reference)) {
    if (length(reference) != 1L)
      stop("reference must be a single biomarker name.")
    if (is.null(biomarkers))
      biomarkers <- setdiff(colnames(X), reference)
    pairs <- data.frame(a = reference, b = biomarkers, stringsAsFactors = FALSE)
  } else if (!is.null(top_n)) {
    if (!inherits(fit, "aucmat_screen"))
      stop("fit must be an aucmat_screen object when using top_n.")
    top <- head(fit$results$biomarker, top_n)
    top <- top[!is.na(top)]
    if (length(top) < 2L) stop("top_n must select at least 2 biomarkers.")
    pairs <- t(utils::combn(top, 2))
    pairs <- data.frame(a = pairs[, 1], b = pairs[, 2], stringsAsFactors = FALSE)
  } else if (!is.null(biomarkers)) {
    if (length(biomarkers) < 2L)
      stop("biomarkers must contain at least 2 names.")
    pairs <- t(utils::combn(biomarkers, 2))
    pairs <- data.frame(a = pairs[, 1], b = pairs[, 2], stringsAsFactors = FALSE)
  } else {
    stop("One of reference, biomarkers, or top_n must be provided.")
  }

  npairs <- nrow(pairs)
  if (npairs > max_pairs) {
    stop(sprintf(
      "%d pairwise comparisons requested but max_pairs = %d.  ",
      npairs, max_pairs),
      "Select fewer biomarkers or increase max_pairs.")
  }

  # Normalise y
  y_norm <- .normalize_binary_y(y, positive = NULL)
  pos <- y_norm == levels(y_norm)[2L]
  neg <- !pos

  X <- as.matrix(X)

  out <- lapply(seq_len(npairs), function(k) {
    a <- pairs$a[k]; b <- pairs$b[k]
    if (!a %in% colnames(X)) stop("Biomarker not found: ", a)
    if (!b %in% colnames(X)) stop("Biomarker not found: ", b)

    # common subjects: observed for both biomarkers
    use <- !is.na(X[, a]) & !is.na(X[, b]) & !is.na(y_norm)
    xa <- X[use, a]; xb <- X[use, b]
    pa <- pos[use]; na <- neg[use]
    np <- sum(pa); nn <- sum(na)

    if (np < 2L || nn < 2L) {
      return(data.frame(
        biomarker_a = a, biomarker_b = b,
        auc_a = NA_real_, auc_b = NA_real_, auc_diff = NA_real_,
        std_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
        p_value = NA_real_, q_value = NA_real_,
        n_common = np + nn, n_pos = np, n_neg = nn,
        stringsAsFactors = FALSE
      ))
    }

    auc_a <- .compute_single_auc(xa, pa, na)$auc_raw
    auc_b <- .compute_single_auc(xb, pa, na)$auc_raw
    auc_diff <- auc_a - auc_b

    # DeLong variance of the difference requires joint placement values.
    # We compute the covariance term using the DeLong structure.
    se_diff <- .delong_diff_se(xa, xb, pa, na)

    z <- stats::qnorm(0.975)  # 95 % CI for the difference
    data.frame(
      biomarker_a = a, biomarker_b = b,
      auc_a = auc_a, auc_b = auc_b, auc_diff = auc_diff,
      std_error = se_diff,
      conf_low  = auc_diff - z * se_diff,
      conf_high = auc_diff + z * se_diff,
      p_value   = if (is.na(se_diff) || se_diff <= 0) NA_real_
                  else 2 * stats::pnorm(abs(auc_diff) / se_diff,
                                        lower.tail = FALSE),
      q_value   = NA_real_,
      n_common  = np + nn, n_pos = np, n_neg = nn,
      stringsAsFactors = FALSE
    )
  })

  res <- do.call(rbind, out)

  # multiplicity adjustment on valid p-values
  valid <- !is.na(res$p_value)
  if (any(valid) && adjust != "none") {
    res$q_value[valid] <- stats::p.adjust(res$p_value[valid], method = adjust)
  }

  class(res) <- c("aucmat_compare", class(res))
  res
}

#' DeLong standard error for the difference of two paired AUCs
#'
#' @param xa,xb Numeric vectors (same length, no NAs) for two biomarkers.
#' @param pos,neg Logical vectors.
#' @return Numeric scalar standard error of AUC_a - AUC_b.
#' @noRd
#' @keywords internal
.delong_diff_se <- function(xa, xb, pos, neg) {
  np <- sum(pos); nn <- sum(neg)
  if (np < 2L || nn < 2L) return(NA_real_)

  # placement values for biomarker a
  V10_a <- vapply(which(pos), function(i) {
    mean((xa[neg] < xa[i]) + 0.5 * (xa[neg] == xa[i]))
  }, numeric(1L))
  V01_a <- vapply(which(neg), function(j) {
    mean((xa[pos] > xa[j]) + 0.5 * (xa[pos] == xa[j]))
  }, numeric(1L))

  # placement values for biomarker b
  V10_b <- vapply(which(pos), function(i) {
    mean((xb[neg] < xb[i]) + 0.5 * (xb[neg] == xb[i]))
  }, numeric(1L))
  V01_b <- vapply(which(neg), function(j) {
    mean((xb[pos] > xb[j]) + 0.5 * (xb[pos] == xb[j]))
  }, numeric(1L))

  # difference in placement values
  D10 <- V10_a - V10_b
  D01 <- V01_a - V01_b

  S10 <- stats::var(D10)
  S01 <- stats::var(D01)

  se <- sqrt(S10 / np + S01 / nn)
  if (!is.finite(se) || se <= 0) NA_real_ else se
}
