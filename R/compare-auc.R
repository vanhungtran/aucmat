# ==============================================================================
# compare_auc() -- Bounded paired AUC comparisons on common subjects
# ==============================================================================

#' Compare AUCs of selected biomarkers (paired, common-subject)
#'
#' Performs DeLong-based or bootstrap paired comparisons of AUC between
#' biomarkers measured on the same subjects.  Supports superiority,
#' non-inferiority, and equivalence (TOST) hypotheses.
#'
#' The function refuses unbounded all-pairs comparisons.
#'
#' @param fit An `aucmat_screen` object (optional; stores outcome encoding).
#' @param X Numeric matrix (same as passed to [aucmat()]).
#' @param y Binary outcome.
#' @param reference Single biomarker name to compare all others against.
#' @param biomarkers Character vector of biomarker names.
#' @param top_n Take the top `top_n` biomarkers by `auc_strength` from `fit`.
#'   Post-selection inference is flagged as exploratory.
#' @param max_pairs Safety limit.  Default 100.
#' @param alternative `"two.sided"` (default), `"greater"`, or `"less"`.
#'   `greater` means AUC_a > AUC_b; `less` means AUC_a < AUC_b.
#' @param hypothesis `"superiority"` (default), `"noninferiority"`, or
#'   `"equivalence"`.
#' @param margin Positive numeric margin for non-inferiority / equivalence.
#' @param conf_level Confidence level in (0, 1).  Default 0.95.
#' @param method `"delong"` (default) or `"bootstrap"`.
#' @param boot_n Bootstrap replicates.  Default 2000.
#' @param seed Optional integer seed.
#' @param adjust Multiplicity adjustment: `"BH"`, `"holm"`, `"bonferroni"`,
#'   or `"none"` (default).
#'
#' @return A data.frame of class `aucmat_compare`.
#' @export
#' @importFrom stats p.adjust pnorm qnorm
compare_auc <- function(fit = NULL, X, y,
                         reference   = NULL,
                         biomarkers  = NULL,
                         top_n       = NULL,
                         max_pairs   = 100,
                         alternative = c("two.sided", "greater", "less"),
                         hypothesis  = c("superiority", "noninferiority",
                                         "equivalence"),
                         margin      = NULL,
                         conf_level  = 0.95,
                         method      = c("delong", "bootstrap"),
                         boot_n      = 2000,
                         seed        = NULL,
                         adjust      = c("none", "BH", "holm", "bonferroni")) {

  alternative <- match.arg(alternative)
  hypothesis  <- match.arg(hypothesis)
  method      <- match.arg(method)
  adjust      <- match.arg(adjust)

  # Validate hypothesis/margin consistency
  if (hypothesis %in% c("noninferiority", "equivalence")) {
    if (is.null(margin) || !is.numeric(margin) || length(margin) != 1L ||
        is.na(margin) || margin <= 0) {
      stop("A positive numeric margin is required for non-inferiority or equivalence.")
    }
  }
  if (hypothesis == "equivalence" && alternative != "two.sided") {
    stop("Equivalence testing requires alternative = 'two.sided'.")
  }

  selection_status <- "prespecified"

  # Build the list of biomarker pairs
  if (!is.null(reference)) {
    if (length(reference) != 1L)
      stop("reference must be a single biomarker name.")
    if (is.null(biomarkers))
      biomarkers <- setdiff(colnames(X), reference)
    pairs <- data.frame(a = reference, b = biomarkers, stringsAsFactors = FALSE)
  } else if (!is.null(top_n)) {
    if (!is.null(fit) && inherits(fit, "aucmat_screen")) {
      top <- head(fit$results$biomarker, top_n)
      top <- top[!is.na(top)]
      if (length(top) < 2L) stop("top_n must select at least 2 biomarkers.")
      pairs <- t(utils::combn(top, 2))
      pairs <- data.frame(a = pairs[, 1], b = pairs[, 2],
                          stringsAsFactors = FALSE)
      selection_status <- "same_data"
    } else {
      stop("fit must be an aucmat_screen object when using top_n.")
    }
  } else if (!is.null(biomarkers)) {
    if (length(biomarkers) < 2L)
      stop("biomarkers must contain at least 2 names.")
    pairs <- t(utils::combn(biomarkers, 2))
    pairs <- data.frame(a = pairs[, 1], b = pairs[, 2],
                        stringsAsFactors = FALSE)
  } else {
    stop("One of reference, biomarkers, or top_n must be provided.")
  }

  npairs <- nrow(pairs)
  if (npairs > max_pairs) {
    stop(sprintf("%d pairwise comparisons requested but max_pairs = %d.  ",
                 npairs, max_pairs),
         "Select fewer biomarkers or increase max_pairs.")
  }

  # Normalise y
  y_norm <- .normalize_binary_y(y)
  pos <- y_norm == levels(y_norm)[2L]
  neg <- !pos

  X <- as.matrix(X)

  with_seed(seed, {
    out <- lapply(seq_len(npairs), function(k) {
      a <- pairs$a[k]; b <- pairs$b[k]
      if (!a %in% colnames(X)) stop("Biomarker not found: ", a)
      if (!b %in% colnames(X)) stop("Biomarker not found: ", b)

      use <- !is.na(X[, a]) & !is.na(X[, b]) & !is.na(y_norm)
      xa <- X[use, a]; xb <- X[use, b]
      pa <- pos[use]; na <- neg[use]
      np <- sum(pa); nn <- sum(na)

      if (np < 2L || nn < 2L) {
        return(.empty_compare_row(a, b, np, nn))
      }

      auc_a <- .compute_single_auc(xa, pa, na)$auc_raw
      auc_b <- .compute_single_auc(xb, pa, na)$auc_raw
      auc_diff <- auc_a - auc_b

      if (method == "delong") {
        se_diff <- .delong_diff_se(xa, xb, pa, na)
      } else {
        bs <- bootstrap_auc_distribution(xa, pa, na, boot_n = boot_n)
        se_diff <- bs$std_error
      }

      if (is.na(se_diff) || se_diff <= 0) {
        return(.empty_compare_row(a, b, np, nn))
      }

      z <- stats::qnorm((1 + conf_level) / 2)
      z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)  # for one-sided NI

      # Build result based on hypothesis type
      if (hypothesis == "superiority") {
        res <- .compare_superiority(auc_diff, se_diff, alternative,
                                     conf_level, a, b, auc_a, auc_b,
                                     np, nn)
      } else if (hypothesis == "noninferiority") {
        res <- .compare_noninferiority(auc_diff, se_diff, margin,
                                        conf_level, a, b, auc_a, auc_b,
                                        np, nn)
      } else {
        res <- .compare_equivalence(auc_diff, se_diff, margin,
                                     conf_level, a, b, auc_a, auc_b,
                                     np, nn)
      }
      res
    })
  })

  res <- do.call(rbind, out)

  # Multiplicity adjustment
  valid <- !is.na(res$p_value)
  if (any(valid) && adjust != "none") {
    res$q_value[valid] <- stats::p.adjust(res$p_value[valid], method = adjust)
  }

  attr(res, "selection_status") <- selection_status
  class(res) <- c("aucmat_compare", class(res))
  res
}

# ---- Internal comparison helpers ----

.empty_compare_row <- function(a, b, np, nn) {
  data.frame(
    biomarker_a = a, biomarker_b = b,
    auc_a = NA_real_, auc_b = NA_real_, auc_diff = NA_real_,
    std_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
    p_value = NA_real_, q_value = NA_real_,
    n_common = np + nn, n_pos = np, n_neg = nn,
    hypothesis = NA_character_, margin = NA_real_,
    stringsAsFactors = FALSE
  )
}

.compare_superiority <- function(delta, se, alternative, conf_level,
                                  a, b, auc_a, auc_b, np, nn) {
  z <- stats::qnorm((1 + conf_level) / 2)
  z_stat <- delta / se

  if (alternative == "two.sided") {
    p_val <- 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
    ci <- c(delta - z * se, delta + z * se)
  } else if (alternative == "greater") {
    p_val <- stats::pnorm(z_stat, lower.tail = FALSE)
    ci <- c(delta - z * se, Inf)
  } else {  # less
    p_val <- stats::pnorm(z_stat, lower.tail = TRUE)
    ci <- c(-Inf, delta + z * se)
  }

  data.frame(
    biomarker_a = a, biomarker_b = b,
    auc_a = auc_a, auc_b = auc_b, auc_diff = delta,
    std_error = se, conf_low = ci[1], conf_high = ci[2],
    p_value = p_val, q_value = NA_real_,
    n_common = np + nn, n_pos = np, n_neg = nn,
    hypothesis = "superiority", margin = 0,
    stringsAsFactors = FALSE
  )
}

.compare_noninferiority <- function(delta, se, margin, conf_level,
                                     a, b, auc_a, auc_b, np, nn) {
  # H0: delta <= -margin  vs H1: delta > -margin
  z_alpha <- stats::qnorm(1 - (1 - conf_level))
  z_stat <- (delta + margin) / se
  p_val <- stats::pnorm(z_stat, lower.tail = FALSE)
  ci_low <- delta - z_alpha * se  # one-sided lower bound

  data.frame(
    biomarker_a = a, biomarker_b = b,
    auc_a = auc_a, auc_b = auc_b, auc_diff = delta,
    std_error = se, conf_low = ci_low, conf_high = Inf,
    p_value = p_val, q_value = NA_real_,
    n_common = np + nn, n_pos = np, n_neg = nn,
    hypothesis = "noninferiority", margin = margin,
    stringsAsFactors = FALSE
  )
}

.compare_equivalence <- function(delta, se, margin, conf_level,
                                  a, b, auc_a, auc_b, np, nn) {
  # TOST: H0_lower: delta <= -margin, H0_upper: delta >= margin
  z_alpha <- stats::qnorm(1 - (1 - conf_level))
  z_lower <- (delta + margin) / se
  z_upper <- (delta - margin) / se
  p_lower <- stats::pnorm(z_lower, lower.tail = FALSE)
  p_upper <- stats::pnorm(z_upper, lower.tail = TRUE)
  p_equiv <- max(p_lower, p_upper)

  # 1 - 2*alpha CI
  z_2a <- stats::qnorm(1 - (1 - conf_level))
  ci <- c(delta - z_2a * se, delta + z_2a * se)

  data.frame(
    biomarker_a = a, biomarker_b = b,
    auc_a = auc_a, auc_b = auc_b, auc_diff = delta,
    std_error = se, conf_low = ci[1], conf_high = ci[2],
    p_value = p_equiv, q_value = NA_real_,
    n_common = np + nn, n_pos = np, n_neg = nn,
    hypothesis = "equivalence", margin = margin,
    stringsAsFactors = FALSE
  )
}

#' @export
print.aucmat_compare <- function(x, n = 10L, ...) {
  ss <- attr(x, "selection_status")
  if (!is.null(ss) && ss == "same_data") {
    cat("NOTE: top_n selection from same data -- inference is exploratory.\n")
  }
  cat("<aucmat_compare>  ", nrow(x), " pairwise comparisons\n", sep = "")
  n_show <- min(n, nrow(x))
  print.data.frame(x[seq_len(n_show), , drop = FALSE], row.names = FALSE)
  if (nrow(x) > n_show) cat("... ", nrow(x) - n_show, " more rows\n", sep = "")
  invisible(x)
}

# ---- DeLong SE for paired difference ----

.delong_diff_se <- function(xa, xb, pos, neg) {
  np <- sum(pos); nn <- sum(neg)
  if (np < 2L || nn < 2L) return(NA_real_)

  pos_idx <- which(pos); neg_idx <- which(neg)

  V10_a <- vapply(pos_idx, function(i) {
    mean((xa[neg] < xa[i]) + 0.5 * (xa[neg] == xa[i]))
  }, numeric(1L))
  V01_a <- vapply(neg_idx, function(j) {
    mean((xa[pos] > xa[j]) + 0.5 * (xa[pos] == xa[j]))
  }, numeric(1L))
  V10_b <- vapply(pos_idx, function(i) {
    mean((xb[neg] < xb[i]) + 0.5 * (xb[neg] == xb[i]))
  }, numeric(1L))
  V01_b <- vapply(neg_idx, function(j) {
    mean((xb[pos] > xb[j]) + 0.5 * (xb[pos] == xb[j]))
  }, numeric(1L))

  D10 <- V10_a - V10_b
  D01 <- V01_a - V01_b
  S10 <- stats::var(D10)
  S01 <- stats::var(D01)
  se <- sqrt(S10 / np + S01 / nn)
  if (!is.finite(se) || se <= 0) NA_real_ else se
}
