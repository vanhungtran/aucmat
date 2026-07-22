# ==============================================================================
# auc_stability() — Stratified bootstrap rank-stability analysis
# ==============================================================================

#' Bootstrap rank-stability analysis for biomarker screening
#'
#' Resamples positive and negative subjects separately with replacement,
#' recomputes the complete feature screen for each replicate, and
#' aggregates rank distributions to assess the stability of leading
#' biomarkers.
#'
#' @param X Numeric matrix (samples in rows, biomarkers in columns).
#' @param y Binary outcome vector.
#' @param positive Optional positive class label.
#' @param times Number of bootstrap replicates.  Default 1000.
#' @param top_k Integer vector of set sizes for top-k probability
#'   reporting.  Default `c(10, 25, 50)`.
#' @param seed Integer seed for reproducibility.
#' @param max_pairs_for_coselection Maximum number of top biomarkers for
#'   which pairwise co-selection frequencies are reported.  Default 20.
#'
#' @return An object of class `aucmat_stability`, a list with components:
#'   `rank_summary`, `top_k_probs`, `coselection`, `settings`.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' X <- matrix(rnorm(200*10), 200, 10, dimnames=list(NULL, paste0("bm",1:10)))
#' y <- rep(c(0,1), each=100)
#' stab <- auc_stability(X, y, times=50, seed=1)
#' print(stab)
#' }
#' @export
auc_stability <- function(X, y, positive = NULL,
                           times   = 1000,
                           top_k   = c(10, 25, 50),
                           seed    = NULL,
                           max_pairs_for_coselection = 20) {

  if (!is.null(seed)) set.seed(seed)

  # Normalise and get class indices
  y_norm <- .normalize_binary_y(y, positive)
  levs <- levels(y_norm)
  pos_idx <- which(y_norm == levs[2L])
  neg_idx <- which(y_norm == levs[1L])
  np <- length(pos_idx)
  nn <- length(neg_idx)

  if (np < 3L || nn < 3L)
    stop("Need at least 3 positive and 3 negative observations for bootstrap stability.")

  X <- as.matrix(X)
  p <- ncol(X)
  cn <- colnames(X)

  # Storage for bootstrap ranks
  rank_matrix <- matrix(NA_integer_, nrow = times, ncol = p)
  colnames(rank_matrix) <- cn
  auc_matrix   <- matrix(NA_real_, nrow = times, ncol = p)
  colnames(auc_matrix) <- cn

  n_ok <- 0L

  for (b in seq_len(times)) {
    bp <- sample(pos_idx, np, replace = TRUE)
    bn <- sample(neg_idx, nn, replace = TRUE)
    b_idx <- c(bp, bn)
    Xb <- X[b_idx, , drop = FALSE]
    yb <- y_norm[b_idx]

    res_b <- try(.compute_matrix_auc(Xb, yb == levs[2L], yb == levs[1L]),
                 silent = TRUE)
    if (inherits(res_b, "try-error")) next

    n_ok <- n_ok + 1L
    auc_matrix[n_ok, ] <- res_b$auc_strength
    rank_matrix[n_ok, ] <- rank(-res_b$auc_strength, ties.method = "min")
  }

  if (n_ok == 0L) stop("All bootstrap replicates failed.")

  # Trim unused rows
  if (n_ok < times) {
    rank_matrix <- rank_matrix[seq_len(n_ok), , drop = FALSE]
    auc_matrix   <- auc_matrix[seq_len(n_ok), , drop = FALSE]
  }

  # ---- Aggregate rank summaries ----
  rank_median  <- apply(rank_matrix, 2, stats::median, na.rm = TRUE)
  rank_q25     <- apply(rank_matrix, 2, stats::quantile, 0.25, na.rm = TRUE)
  rank_q75     <- apply(rank_matrix, 2, stats::quantile, 0.75, na.rm = TRUE)
  auc_mean     <- colMeans(auc_matrix, na.rm = TRUE)
  auc_sd       <- apply(auc_matrix, 2, stats::sd, na.rm = TRUE)

  rank_summary <- data.frame(
    biomarker    = cn,
    rank_median  = rank_median,
    rank_q25     = rank_q25,
    rank_q75     = rank_q75,
    auc_mean     = auc_mean,
    auc_sd       = auc_sd,
    top1_freq    = colMeans(rank_matrix == 1L, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  rank_summary <- rank_summary[order(rank_summary$rank_median), ]
  rownames(rank_summary) <- NULL

  # ---- Top-k probabilities ----
  top_k <- sort(unique(as.integer(top_k)))
  top_k <- top_k[top_k >= 1L & top_k <= p]

  # Build top-k probability table
  tk_biomarker <- character(0L)
  tk_k         <- integer(0L)
  tk_prob      <- numeric(0L)
  for (k in top_k) {
    probs <- as.numeric(colMeans(rank_matrix <= k, na.rm = TRUE))
    names(probs) <- NULL
    tk_biomarker <- c(tk_biomarker, cn)
    tk_k         <- c(tk_k, rep(as.integer(k), length(cn)))
    tk_prob      <- c(tk_prob, probs)
  }
  top_k_probs <- data.frame(
    biomarker = tk_biomarker,
    k         = tk_k,
    prob      = tk_prob,
    stringsAsFactors = FALSE
  )
  top_k_probs <- top_k_probs[order(top_k_probs$k, -top_k_probs$prob), ]
  rownames(top_k_probs) <- NULL

  # ---- Pairwise co-selection ----
  n_top_cs <- min(max_pairs_for_coselection, p)
  top_biomarkers <- rank_summary$biomarker[seq_len(n_top_cs)]

  co_matrix <- matrix(0, nrow = n_top_cs, ncol = n_top_cs)
  rownames(co_matrix) <- top_biomarkers
  colnames(co_matrix) <- top_biomarkers

  for (k in top_k) {
    in_top <- rank_matrix <= k
    for (i in seq_len(n_top_cs)) {
      for (j in seq_len(n_top_cs)) {
        if (i >= j) next
        co_matrix[i, j] <- co_matrix[i, j] +
          mean(in_top[, top_biomarkers[i]] & in_top[, top_biomarkers[j]],
               na.rm = TRUE)
        co_matrix[j, i] <- co_matrix[i, j]
      }
    }
  }

  coselection <- list(
    biomarkers = top_biomarkers,
    matrix     = co_matrix,
    top_k      = top_k
  )

  settings <- list(
    times    = times,
    n_ok     = n_ok,
    n_failed = times - n_ok,
    top_k    = top_k,
    seed     = seed
  )

  out <- list(
    rank_summary = rank_summary,
    top_k_probs  = top_k_probs,
    coselection  = coselection,
    settings     = settings
  )
  class(out) <- "aucmat_stability"
  out
}
