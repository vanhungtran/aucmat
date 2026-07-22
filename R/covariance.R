# ==============================================================================
# Joint DeLong Covariance Engine
#
# Computes the full (p x p) covariance matrix of correlated AUC estimates
# for biomarkers measured on the same subjects.  Supports the global Wald
# test and correlated inference.
# ==============================================================================

#' Joint DeLong covariance matrix for multiple biomarkers
#'
#' Computes the DeLong placement-value covariance matrix for a set of
#' biomarkers on a common subject set.  All biomarkers must be observed
#' on the same subjects (common complete cases).
#'
#' @param X Numeric matrix (n × p), no missing values.
#' @param pos Logical vector; TRUE for positive-class.
#' @param neg Logical vector; TRUE for negative-class.
#'
#' @return A list with `covariance` (p × p matrix), `aucs` (p-vector),
#'   `n_pos`, `n_neg`.
#' @noRd
#' @keywords internal
joint_delong_covariance <- function(X, pos, neg) {
  X <- as.matrix(X)
  p <- ncol(X)
  n_pos <- sum(pos)
  n_neg <- sum(neg)

  if (n_pos < 2L || n_neg < 2L) {
    return(list(covariance = matrix(NA_real_, p, p),
                aucs = rep(NA_real_, p),
                n_pos = n_pos, n_neg = n_neg,
                status = "insufficient_sample"))
  }

  pos_idx <- which(pos)
  neg_idx <- which(neg)

  # Storage for placement values
  # V10: (n_pos × p) — placement of each positive subject per biomarker
  # V01: (n_neg × p) — placement of each negative subject per biomarker
  V10 <- matrix(NA_real_, n_pos, p)
  V01 <- matrix(NA_real_, n_neg, p)
  aucs <- numeric(p)

  for (j in seq_len(p)) {
    x <- X[, j]

    # Per-subject placement values
    v10 <- vapply(pos_idx, function(i) {
      mean((x[neg_idx] < x[i]) + 0.5 * (x[neg_idx] == x[i]))
    }, numeric(1L))
    v01 <- vapply(neg_idx, function(jj) {
      mean((x[pos_idx] > x[jj]) + 0.5 * (x[pos_idx] == x[jj]))
    }, numeric(1L))

    V10[, j] <- v10
    V01[, j] <- v01
    aucs[j] <- mean(v10)
  }

  # Covariance components
  # S10: p × p matrix from V10, scaled by 1/n_pos
  # S01: p × p matrix from V01, scaled by 1/n_neg
  if (n_pos > 1L) {
    S10 <- stats::cov(V10) / n_pos   # cov divides by (n_pos-1), then we divide by n_pos
  } else {
    S10 <- matrix(0, p, p)
  }
  if (n_neg > 1L) {
    S01 <- stats::cov(V01) / n_neg
  } else {
    S01 <- matrix(0, p, p)
  }

  covariance <- S10 + S01

  list(covariance = covariance, aucs = aucs,
       n_pos = n_pos, n_neg = n_neg,
       status = "ok")
}

#' Global Wald test for equality of multiple correlated AUCs
#'
#' Tests H0: AUC_1 = AUC_2 = ... = AUC_p on a common subject set
#' using the joint DeLong covariance matrix.
#'
#' @param fit An `aucmat_screen` object (optional, for outcome encoding).
#' @param X Numeric matrix.
#' @param y Binary outcome.
#' @param biomarkers Character vector of biomarker names to include.
#'   Default: all.
#' @param max_biomarkers Safety limit.  Default 100.
#'
#' @return A list of class `aucmat_global_test` with components:
#'   `statistic`, `df`, `p_value`, `aucs`, `covariance`, `contrasts`,
#'   `n_common`, `n_pos`, `n_neg`, `status`.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' sim <- simulate_auc_matrix(n=100, prevalence=0.3,
#'   target_aucs=c(0.85,0.75,0.65), correlation=0.3, structure="exchangeable")
#' X <- as.matrix(sim$data[,1:3]); y <- sim$data$truth
#' compare_auc_global(X=X, y=y)
#' }
#' @export
compare_auc_global <- function(fit = NULL, X, y,
                                biomarkers = NULL,
                                max_biomarkers = 100) {
  # Normalise outcome
  if (!is.null(fit) && inherits(fit, "aucmat_screen")) {
    y_norm <- .normalize_binary_y(y)
  } else {
    y_norm <- .normalize_binary_y(y)
  }

  X <- as.matrix(X)

  if (is.null(biomarkers)) {
    biomarkers <- colnames(X)
  }
  biomarkers <- biomarkers[!is.na(biomarkers)]

  p <- length(biomarkers)
  if (p < 2L) stop("At least 2 biomarkers are required.")
  if (p > max_biomarkers) {
    stop(sprintf("%d biomarkers requested but max_biomarkers = %d.  ",
                 p, max_biomarkers),
         "Select fewer biomarkers or increase max_biomarkers.")
  }

  # Verify biomarkers exist
  missing <- setdiff(biomarkers, colnames(X))
  if (length(missing) > 0L) {
    stop("Biomarkers not found: ", paste(missing, collapse = ", "))
  }

  # Common complete cases
  X_sub <- X[, biomarkers, drop = FALSE]
  keep <- complete.cases(X_sub) & !is.na(y_norm)
  X_clean <- X_sub[keep, , drop = FALSE]
  y_clean <- y_norm[keep]

  pos <- y_clean == levels(y_clean)[2L]
  neg <- !pos
  n_common <- nrow(X_clean)
  n_pos <- sum(pos)
  n_neg <- sum(neg)

  if (n_pos < 2L || n_neg < 2L) {
    out <- list(statistic = NA_real_, df = NA_integer_,
                p_value = NA_real_,
                aucs = rep(NA_real_, p),
                covariance = matrix(NA_real_, p, p),
                contrasts = NULL,
                n_common = n_common, n_pos = n_pos, n_neg = n_neg,
                status = "insufficient_sample")
    class(out) <- "aucmat_global_test"
    return(out)
  }

  # Joint covariance
  jc <- joint_delong_covariance(X_clean, pos, neg)
  aucs <- jc$aucs
  Sigma <- jc$covariance

  # Default contrasts: compare biomarkers 2:p with biomarker 1
  C <- cbind(-1, diag(p - 1L))
  rownames(C) <- biomarkers[-1L]
  colnames(C) <- biomarkers

  # Wald statistic: (C * aucs)' * inv(C * Sigma * C') * (C * aucs)
  contrast_est <- as.vector(C %*% aucs)
  contrast_cov <- C %*% Sigma %*% t(C)

  # Effective rank via eigenvalue tolerance
  eig <- eigen(contrast_cov, symmetric = TRUE)
  tol <- max(eig$values) * p * .Machine$double.eps
  rank <- sum(eig$values > tol)

  if (rank == 0L) {
    out <- list(statistic = NA_real_, df = NA_integer_,
                p_value = NA_real_,
                aucs = aucs, covariance = Sigma,
                contrasts = C,
                n_common = n_common, n_pos = n_pos, n_neg = n_neg,
                status = "non_testable_singular")
    class(out) <- "aucmat_global_test"
    return(out)
  }

  # Use generalized inverse on the stable subspace
  inv_eig <- eig$values
  inv_eig[inv_eig > tol] <- 1 / inv_eig[inv_eig > tol]
  inv_eig[inv_eig <= tol] <- 0
  cov_inv <- eig$vectors %*% diag(inv_eig) %*% t(eig$vectors)

  stat <- as.numeric(t(contrast_est) %*% cov_inv %*% contrast_est)
  p_value <- stats::pchisq(stat, df = rank, lower.tail = FALSE)

  out <- list(statistic = stat, df = rank, p_value = p_value,
              aucs = stats::setNames(aucs, biomarkers),
              covariance = Sigma, contrasts = C,
              n_common = n_common, n_pos = n_pos, n_neg = n_neg,
              status = "ok")
  class(out) <- "aucmat_global_test"
  out
}

#' @export
print.aucmat_global_test <- function(x, ...) {
  cat("<aucmat_global_test>\n")
  cat("  H0: all AUCs equal\n")
  cat("  n =", x$n_common, " (", x$n_pos, "+ / ", x$n_neg, "-)\n", sep = "")
  if (x$status == "ok") {
    cat(sprintf("  Wald chi2(%d) = %.3f, p = %.4f\n",
                x$df, x$statistic, x$p_value))
  } else {
    cat("  Status:", x$status, "\n")
  }
  cat("\nAUC estimates:\n")
  print(round(x$aucs, 4))
  invisible(x)
}
