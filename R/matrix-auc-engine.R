# ==============================================================================
# Matrix AUC Engine
#
# Computes feature-wise Mann-Whitney rank-statistic AUC values for a numeric
# biomarker matrix against a binary outcome.  Designed to be called from
# aucmat(); not exported.
# ==============================================================================

#' Compute raw AUC, strength, and effect direction for one biomarker
#'
#' Uses the Mann-Whitney U rank formulation: half credit for ties.
#'
#' @param x Numeric vector of biomarker values.
#' @param pos Logical vector; TRUE for positive-class observations.
#' @param neg Logical vector; TRUE for negative-class observations.
#'
#' @return A list with elements `auc_raw`, `auc_strength`,
#'   `effect_direction` (numeric: 1, -1, or 0), `n_pos`, `n_neg`,
#'   `n_used`, `status`.
#' @noRd
#' @keywords internal
.compute_single_auc <- function(x, pos, neg) {
  use   <- !is.na(x)
  p_use <- pos & use
  n_use <- neg & use
  np    <- sum(p_use)
  nn    <- sum(n_use)
  n_used <- np + nn

  if (np < 2L || nn < 2L) {
    status <- if (np < 2L && nn < 2L) "insufficient_both"
              else if (np < 2L) "insufficient_positive"
              else "insufficient_negative"
    return(list(
      auc_raw          = NA_real_,
      auc_strength     = NA_real_,
      effect_direction = NA_real_,
      n_pos            = np,
      n_neg            = nn,
      n_used           = n_used,
      status           = status
    ))
  }

  x_use <- x[use]

  # Detect constant or near-constant features
  if (stats::sd(x_use) == 0 || length(unique(x_use)) <= 1L) {
    return(list(
      auc_raw          = NA_real_,
      auc_strength     = NA_real_,
      effect_direction = NA_real_,
      n_pos            = np,
      n_neg            = nn,
      n_used           = n_used,
      status           = "constant"
    ))
  }
  p_sub <- match(which(p_use), which(use))

  r  <- rank(x_use, ties.method = "average")
  U  <- sum(r[p_sub]) - np * (np + 1) / 2
  auc_raw <- as.numeric(U / (np * nn))

  auc_strength    <- 0.5 + abs(auc_raw - 0.5)
  effect_direction <- if (auc_raw > 0.5) 1 else if (auc_raw < 0.5) -1 else 0

  list(
    auc_raw          = auc_raw,
    auc_strength     = auc_strength,
    effect_direction = effect_direction,
    n_pos            = np,
    n_neg            = nn,
    n_used           = n_used,
    status           = "ok"
  )
}

#' Matrix-wide AUC computation
#'
#' Applies `.compute_single_auc()` to every column of a numeric matrix.
#'
#' @param X Numeric matrix or data.frame (biomarkers in columns).
#' @param pos Logical vector; TRUE for positive class.
#' @param neg Logical vector; TRUE for negative class.
#'
#' @return A data.frame with columns `biomarker`, `auc_raw`, `auc_strength`,
#'   `effect_direction`, `n_used`, `n_pos`, `n_neg`, `status`.
#' @noRd
#' @keywords internal
.compute_matrix_auc <- function(X, pos, neg) {
  X <- as.matrix(X)
  cn <- colnames(X)
  if (is.null(cn)) cn <- paste0("V", seq_len(ncol(X)))

  out_list <- lapply(seq_len(ncol(X)), function(j) {
    .compute_single_auc(X[, j], pos, neg)
  })

  # Build data.frame from list elements
  df <- data.frame(
    biomarker        = cn,
    auc_raw          = vapply(out_list, `[[`, numeric(1L), "auc_raw"),
    auc_strength     = vapply(out_list, `[[`, numeric(1L), "auc_strength"),
    effect_direction = vapply(out_list, `[[`, numeric(1L), "effect_direction"),
    n_used           = vapply(out_list, `[[`, integer(1L), "n_used"),
    n_pos            = vapply(out_list, `[[`, integer(1L), "n_pos"),
    n_neg            = vapply(out_list, `[[`, integer(1L), "n_neg"),
    status           = vapply(out_list, `[[`, character(1L), "status"),
    stringsAsFactors = FALSE
  )

  # Map numeric direction to labels
  dir_map <- c("1" = "higher_in_positive", "-1" = "lower_in_positive",
               "0" = "none")
  df$effect_direction <- dir_map[as.character(df$effect_direction)]
  df$effect_direction[is.na(df$effect_direction)] <- NA_character_

  df
}
