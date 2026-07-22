# ==============================================================================
# aucmat_multiclass() — Matrix screening for multiclass outcomes
#
# v0.4: one-versus-rest, one-versus-one, and Hand-Till multiclass AUC.
# ==============================================================================

#' Screen biomarkers against a multiclass outcome
#'
#' Computes one-versus-rest AUC for each biomarker-class pair, plus
#' Hand-Till multiclass AUC per biomarker.  Hand-Till is the average AUC
#' of all pairwise class comparisons weighted by class prevalence.
#'
#' @param X Numeric matrix (n x p).
#' @param y Outcome with 3+ classes.  Factor, character, or integer.
#' @param ci `"none"` (default) or `"bootstrap"`.  DeLong is not defined
#'   for multiclass.
#' @param conf_level Confidence level.
#' @param boot_n Bootstrap replicates when `ci = "bootstrap"`.
#' @param seed Optional seed.
#' @param adjust Multiplicity adjustment across all biomarker-class pairs.
#'
#' @return A list of class `aucmat_multiclass` with components `ovr`
#'   (one-vs-rest), `hand_till` (one value per biomarker), `classes`,
#'   `settings`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' # 3-class outcome with 2 biomarkers
#' n <- 300
#' y <- sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
#' X <- cbind(
#'   X1 = rnorm(n) + ifelse(y == "A", 1.5, ifelse(y == "B", 0.5, -1)),
#'   X2 = rnorm(n) + ifelse(y == "C", 2.0, 0)
#' )
#' colnames(X) <- c("bm1", "bm2")
#' res <- aucmat_multiclass(X, y, ci = "none")
#' print(res)
#' }
aucmat_multiclass <- function(X, y, ci = c("none", "bootstrap"),
                               conf_level = 0.95, boot_n = 2000,
                               seed = NULL,
                               adjust = c("none", "BH", "holm", "bonferroni")) {
  ci     <- match.arg(ci)
  adjust <- match.arg(adjust)

  X <- as.matrix(X)
  p <- ncol(X)
  cn <- colnames(X)
  if (is.null(cn)) cn <- paste0("X", seq_len(p))

  if (!is.factor(y)) y <- factor(y)
  classes <- levels(y)
  K <- length(classes)
  if (K < 3L) stop("Need at least 3 classes for multiclass AUC. Use aucmat() for binary outcomes.")
  if (any(table(y) < 3L)) stop("Each class must have at least 3 observations.")

  n_total <- length(y)

  # ---- One-vs-rest AUC for each biomarker-class pair ----
  ovr_rows <- expand.grid(biomarker = cn, class = classes,
    stringsAsFactors = FALSE)
  ovr_rows$auc_ovr <- NA_real_
  ovr_rows$n_class <- NA_integer_

  for (i in seq_len(nrow(ovr_rows))) {
    bm  <- ovr_rows$biomarker[i]
    cls <- ovr_rows$class[i]
    x   <- X[, bm]
    use <- !is.na(x) & !is.na(y)
    pos <- y[use] == cls
    neg <- !pos
    np <- sum(pos); nn <- sum(neg)
    ovr_rows$n_class[i] <- np
    if (np >= 2L && nn >= 2L) {
      ovr_rows$auc_ovr[i] <- .compute_single_auc(x[use], pos, neg)$auc_raw
    }
  }

  # Multiplicity across all OVR pairs
  if (adjust != "none") {
    all_p <- vapply(seq_len(nrow(ovr_rows)), function(i) {
      bm <- ovr_rows$biomarker[i]; cls <- ovr_rows$class[i]
      x <- X[, bm]; use <- !is.na(x) & !is.na(y)
      pos <- y[use] == cls; neg <- !pos
      if (sum(pos) < 2L || sum(neg) < 2L) return(NA_real_)
      se <- try(unname(.delong_variance(x[use], pos, neg)["std_error"]), silent = TRUE)
      if (inherits(se, "try-error") || is.na(se) || se <= 0) return(NA_real_)
      auc_r <- .compute_single_auc(x[use], pos, neg)$auc_raw
      2 * stats::pnorm(abs(auc_r - 0.5) / se, lower.tail = FALSE)
    }, numeric(1L))
    ovr_rows$p_value <- all_p
    ovr_rows$q_value <- .adjust_pvalues(all_p, method = adjust)
  }

  # ---- Hand-Till multiclass AUC per biomarker ----
  hand_till <- data.frame(
    biomarker = cn,
    auc_hand_till = NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(p)) {
    x <- X[, j]
    use <- !is.na(x) & !is.na(y)
    xu <- x[use]; yu <- y[use]

    # All pairwise class AUCs
    pairwise_aucs <- numeric(0)
    weights       <- numeric(0)
    for (a in seq_len(K - 1L)) {
      for (b in seq(a + 1L, K)) {
        idx_ab <- which(yu %in% classes[c(a, b)])
        if (length(idx_ab) < 4L) next
        x_ab <- xu[idx_ab]
        y_ab <- yu[idx_ab]
        pos_ab <- y_ab == classes[a]
        neg_ab <- !pos_ab
        np <- sum(pos_ab); nn <- sum(neg_ab)
        if (np < 2L || nn < 2L) next
        auc_ab <- .compute_single_auc(x_ab, pos_ab, neg_ab)$auc_raw
        pairwise_aucs <- c(pairwise_aucs, auc_ab)
        weights <- c(weights, np + nn)
      }
    }
    if (length(pairwise_aucs) > 0L) {
      hand_till$auc_hand_till[j] <- stats::weighted.mean(
        pairwise_aucs, weights, na.rm = TRUE)
    }
  }

  out <- list(
    ovr       = ovr_rows,
    hand_till = hand_till,
    classes   = classes,
    settings  = list(ci = ci, conf_level = conf_level, adjust = adjust,
      n_total = n_total, n_classes = K)
  )
  class(out) <- "aucmat_multiclass"
  out
}

#' @export
print.aucmat_multiclass <- function(x, n = 10L, ...) {
  cat("<aucmat_multiclass>  ", x$settings$n_classes, " classes (",
    paste(x$classes, collapse = ", "), ")\n", sep = "")
  cat("  ", nrow(x$hand_till), " biomarkers\n\n", sep = "")

  cat("Hand-Till multiclass AUC:\n")
  df <- x$hand_till[order(x$hand_till$auc_hand_till, decreasing = TRUE, na.last = TRUE), ]
  n_show <- min(n, nrow(df))
  top <- df[seq_len(n_show), ]
  print(top, row.names = FALSE)
  if (nrow(df) > n_show) cat("... ", nrow(df) - n_show, " more\n", sep = "")

  cat("\nTop one-vs-rest pairs:\n")
  ovr <- x$ovr[order(x$ovr$auc_ovr, decreasing = TRUE, na.last = TRUE), ]
  n_show2 <- min(n, nrow(ovr))
  print(ovr[seq_len(n_show2), ], row.names = FALSE)
  invisible(x)
}

#' @export
plot.aucmat_multiclass <- function(x, ...) {
  # Heatmap: biomarker x class of one-vs-rest AUC
  ovr <- x$ovr
  ovr <- ovr[!is.na(ovr$auc_ovr), ]

  ggplot2::ggplot(ovr, ggplot2::aes(x = .data$class, y = .data$biomarker,
    fill = .data$auc_ovr)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0.5, limits = c(0, 1), name = "AUC (OVR)") +
    ggplot2::geom_text(ggplot2::aes(
      label = formatC(.data$auc_ovr, 2, format = "f")), size = 3) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = "Class", y = NULL,
      title = "One-vs-rest AUC by biomarker and class")
}
