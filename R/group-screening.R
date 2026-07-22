# ==============================================================================
# aucmat_by() — Stratified biomarker screening by group
#
# Screens biomarkers within each level of a grouping variable (batch, site,
# stratum), returning per-group results and between-group heterogeneity.
# ==============================================================================

#' Stratified biomarker screening by group
#'
#' Runs [aucmat()] within each level of a grouping variable and returns
#' per-group results plus a heterogeneity summary (I-squared, Cochran's Q).
#' Useful for assessing whether biomarker discrimination varies by batch,
#' site, or clinical stratum.
#'
#' @param X Numeric matrix.
#' @param y Binary outcome.
#' @param group Factor or vector coercible to factor.  Screening is
#'   performed within each level.
#' @param positive Optional positive class label.
#' @param ci,conf_level,adjust,boot_n Passed to [aucmat()].
#' @param min_group_size Minimum number of observations per group
#'   (both classes required).  Default 10.
#'
#' @return A list of class `aucmat_by` with `results` (per-group data.frame),
#'   `heterogeneity` (I-squared, Q, p per biomarker), `groups`, `settings`.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_auc_matrix(n = 300, prevalence = 0.3,
#'   target_aucs = c(0.8, 0.7, 0.6), correlation = 0.3,
#'   structure = "exchangeable")
#' grp <- sample(c("SiteA", "SiteB"), 300, replace = TRUE)
#' res <- aucmat_by(as.matrix(sim$data[, 1:3]), sim$data$truth, grp,
#'   ci = "none")
#' print(res)
#' }
aucmat_by <- function(X, y, group, positive = NULL,
                       ci = c("delong", "bootstrap", "none"),
                       conf_level = 0.95, adjust = "BH",
                       boot_n = 2000, min_group_size = 10L) {
  ci <- match.arg(ci)

  X <- as.matrix(X)
  if (!is.factor(group)) group <- factor(group)
  grp_levels <- levels(group)
  n_groups <- length(grp_levels)
  p <- ncol(X)

  if (n_groups < 2L) stop("Need at least 2 groups.")

  y_norm <- .normalize_binary_y(y, positive)
  pos <- y_norm == levels(y_norm)[2L]
  neg <- !pos

  # Per-group screening
  group_fits <- vector("list", n_groups)
  names(group_fits) <- grp_levels

  for (g in seq_len(n_groups)) {
    idx <- which(group == grp_levels[g])
    if (length(idx) < min_group_size) {
      group_fits[[g]] <- NULL
      next
    }
    Xg <- X[idx, , drop = FALSE]
    yg <- y[idx]
    group_fits[[g]] <- aucmat(Xg, yg, positive = positive,
      ci = ci, conf_level = conf_level, adjust = adjust, boot_n = boot_n)
  }

  # Per-biomarker heterogeneity: Q = sum(w_i * (theta_i - theta_bar)^2)
  cn <- colnames(X)
  if (is.null(cn)) cn <- paste0("X", seq_len(p))

  heterogeneity <- data.frame(
    biomarker   = cn,
    i_squared   = NA_real_,
    cochran_q   = NA_real_,
    q_p_value   = NA_real_,
    n_groups    = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(p)) {
    aucs <- numeric(n_groups)
    ses  <- numeric(n_groups)
    n_valid <- 0L

    for (g in seq_len(n_groups)) {
      if (is.null(group_fits[[g]])) next
      res_g <- group_fits[[g]]$results
      row_j <- which(res_g$biomarker == cn[j])
      if (length(row_j) != 1L) next
      n_valid <- n_valid + 1L
      aucs[n_valid] <- res_g$auc_raw[row_j]
      ses[n_valid]  <- res_g$std_error[row_j]
    }

    if (n_valid < 2L) next

    aucs <- aucs[seq_len(n_valid)]
    ses  <- ses[seq_len(n_valid)]

    # Inverse-variance weights
    w <- 1 / ses^2
    w[!is.finite(w)] <- 0
    theta_bar <- sum(w * aucs) / sum(w)

    # Cochran's Q
    Q <- sum(w * (aucs - theta_bar)^2)
    df <- n_valid - 1L

    # I-squared
    I2 <- max(0, (Q - df) / Q * 100)

    heterogeneity$i_squared[j] <- I2
    heterogeneity$cochran_q[j] <- Q
    heterogeneity$q_p_value[j] <- stats::pchisq(Q, df, lower.tail = FALSE)
    heterogeneity$n_groups[j] <- n_valid
  }

  out <- list(
    results       = heterogeneity,
    heterogeneity = heterogeneity,
    groups        = grp_levels,
    group_fits    = group_fits,
    settings      = list(
      ci = ci, conf_level = conf_level, adjust = adjust,
      min_group_size = min_group_size, n_groups = n_groups
    )
  )
  class(out) <- "aucmat_by"
  out
}

#' @export
print.aucmat_by <- function(x, n = 10L, ...) {
  cat("<aucmat_by>  ", x$settings$n_groups, " groups: ",
    paste(x$groups, collapse = ", "), "\n", sep = "")
  cat("\nHeterogeneity summary:\n")
  n_show <- min(n, nrow(x$results))
  top <- x$results[order(x$results$i_squared, decreasing = TRUE), ]
  top <- top[seq_len(n_show), ]
  print(top[, c("biomarker", "i_squared", "cochran_q", "q_p_value", "n_groups")],
    row.names = FALSE)
  if (nrow(x$results) > n_show)
    cat("... ", nrow(x$results) - n_show, " more\n", sep = "")
  invisible(x)
}

#' @export
plot.aucmat_by <- function(x, n_label = 15L, ...) {
  df <- x$results
  df <- df[order(df$i_squared, decreasing = TRUE), ]
  df$label <- ""
  n_lab <- min(n_label, nrow(df))
  df$label[seq_len(n_lab)] <- df$biomarker[seq_len(n_lab)]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$cochran_q,
    y = .data$i_squared)) +
    ggplot2::geom_point(ggplot2::aes(
      colour = .data$q_p_value < 0.05), size = 2, alpha = 0.7) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "grey60", "TRUE" = "#B2182B"),
      labels = c("FALSE" = "NS", "TRUE" = "p < 0.05"), name = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Cochran's Q", y = "I-squared (%)",
      title = "Between-group heterogeneity in biomarker AUCs")
}
