
#' AUC table for each numeric predictor column
#'
#' Compute ROC AUC and confidence intervals for each numeric column of a
#' matrix or data.frame against a binary outcome using \pkg{pROC}.
#' Provides robust missing-data handling for X, including row-dropping,
#' per-column omission, and multiple imputation strategies.
#'
#' Optional imputation methods (via `na_impute`):
#' - "none": no imputation (row-drop if `na_rm = TRUE`, or per-column omission if `na_rm = FALSE`)
#' - "mean", "median", "constant"
#' - "median_by_class": within-outcome medians
#' - "halfmin": half of minimum positive value
#' - "quantile": column quantile (via `na_quantile`)
#' - "knn": base-R Euclidean KNN (with safety cap)
#' - "mice": single imputation (if installed)
#' - "missForest": random forest imputation (if installed)
#'
#' @param X data.frame or matrix of predictors (only numeric columns are used).
#' @param y Binary outcome (factor with 2 levels, logical, or numeric 0/1).
#' @param round.digit Integer. Rounding digits for AUC/CIs. Default 4.
#' @param rank Logical. If TRUE, sort by AUC desc and add rank. Default FALSE.
#' @param ci_level Confidence level for CI. Default 0.95.
#' @param ci_method "delong" or "bootstrap". Default "delong".
#' @param boot_n Bootstrap reps if method = "bootstrap". Default 2000.
#' @param positive Positive class label for factor y (second level used by default).
#' @param direction One of "auto", "<", ">". Passed to pROC::roc. Default "auto".
#' @param na_rm If TRUE and na_impute = "none", drop rows with any NA in numeric X (after removing NA in y).
#' @param na_impute One of "none","mean","median","constant","median_by_class","halfmin","quantile","knn","mice","missForest".
#' @param na_constant Constant for "constant" imputation. Default 0.
#' @param na_quantile Quantile in (0,1) for "quantile" imputation. Default 0.25.
#' @param knn_k Neighbors for KNN. Default 5.
#' @param knn_scale Standardize columns before distances. Default TRUE.
#' @param knn_weighted Inverse-distance weights in KNN. Default TRUE.
#' @param knn_max_n Safety cap for KNN (O(n^2) memory). Default 5000.
#' @param mice_m Number of imputations for mice (coerced to 1). Default 1.
#' @param mice_maxit Iterations for mice. Default 5.
#' @param mice_method Method for mice. Default "pmm".
#' @param mice_seed Optional seed for mice.
#' @param missforest_maxiter Max iterations for missForest. Default 10.
#' @param missforest_ntree Trees for missForest. Default 100.
#' @param track_na If TRUE, include NA diagnostics (n_na, n_imputed, na_handling). Default TRUE.
#'
#' @return data.frame with rows = numeric predictors; columns:
#'   biomarker, auc, ci_low, ci_high, ci_level, ci_method, n_pos, n_neg, direction,
#'   na_handling, n_na, n_imputed, and rank (if rank=TRUE).
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(M1 = c(rnorm(95), rep(NA, 5)), M2 = rnorm(100))
#' y <- factor(rbinom(100, 1, 0.5), labels = c("neg", "pos"))
#' tableroc(X, y, na_impute = "median", rank = TRUE)
#'
#' @export
#' @importFrom pROC roc auc ci.auc
#' @importFrom stats complete.cases median dist sd quantile
tableroc <- function(
    X, y,
    round.digit = 4,
    rank = FALSE,
    ci_level = 0.95,
    ci_method = c("delong", "bootstrap"),
    boot_n = 2000,
    positive = NULL,
    direction = c("auto", "<", ">"),
    na_rm = TRUE,
    na_impute = c("none", "mean", "median", "constant", "median_by_class",
                  "halfmin", "quantile", "knn", "mice", "missForest"),
    na_constant = 0,
    na_quantile = 0.25,
    knn_k = 5,
    knn_scale = TRUE,
    knn_weighted = TRUE,
    knn_max_n = 5000,
    mice_m = 1,
    mice_maxit = 5,
    mice_method = "pmm",
    mice_seed = NULL,
    missforest_maxiter = 10,
    missforest_ntree = 100,
    track_na = TRUE
) {
  if (is.matrix(X)) X <- as.data.frame(X)
  if (!is.data.frame(X)) stop("X must be a data.frame or matrix.")
  if (nrow(X) == 0L) stop("X has 0 rows.")
  if (length(y) != nrow(X)) stop("Length of y must match nrow(X).")
  if (!is.numeric(round.digit) || length(round.digit) != 1L || round.digit < 0) stop("round.digit must be non-negative integer.")
  if (!is.numeric(ci_level) || length(ci_level) != 1L || ci_level <= 0 || ci_level >= 1) stop("ci_level must be in (0,1).")
  if (!is.numeric(na_quantile) || length(na_quantile) != 1L || na_quantile <= 0 || na_quantile >= 1) stop("na_quantile must be in (0,1).")
  if (!is.numeric(knn_k) || knn_k < 1) stop("knn_k must be >= 1.")
  if (!is.numeric(knn_max_n) || knn_max_n < 1) stop("knn_max_n must be >= 1.")
  if (!is.numeric(mice_m) || mice_m < 1) stop("mice_m must be >= 1.")
  if (!is.numeric(mice_maxit) || mice_maxit < 0) stop("mice_maxit must be >= 0.")
  if (!is.numeric(missforest_maxiter) || missforest_maxiter < 1) stop("missforest_maxiter must be >= 1.")
  if (!is.numeric(missforest_ntree) || missforest_ntree < 1) stop("missforest_ntree must be >= 1.")

  ci_method <- match.arg(ci_method)
  direction <- match.arg(direction)
  na_impute <- match.arg(na_impute)

  # Numeric-only predictors
  num_cols <- vapply(X, is.numeric, logical(1))
  X_num <- X[, num_cols, drop = FALSE]
  if (ncol(X_num) == 0L) stop("X must contain at least one numeric column.")

  # Drop rows with NA in y
  keep_y <- !is.na(y)
  X_num <- X_num[keep_y, , drop = FALSE]
  y <- y[keep_y]
  if (nrow(X_num) == 0L) stop("All rows removed because y was NA.")

  # Diagnostics
  pre_na_counts <- vapply(X_num, function(col) sum(is.na(col)), integer(1))
  na_handling_label <- switch(
    na_impute,
    none            = if (isTRUE(na_rm)) "drop_rows" else "per_column",
    mean            = "impute_mean",
    median          = "impute_median",
    constant        = "impute_constant",
    median_by_class = "impute_median_by_class",
    halfmin         = "impute_halfmin",
    quantile        = "impute_quantile",
    knn             = "impute_knn",
    mice            = "impute_mice",
    missForest      = "impute_missForest"
  )

  # NA handling
  if (na_impute == "none") {
    if (isTRUE(na_rm)) {
      keep <- stats::complete.cases(X_num)
      X_work <- X_num[keep, , drop = FALSE]
      y_work <- y[keep]
    } else {
      X_work <- X_num
      y_work <- y
    }
    imputed_counts <- integer(ncol(X_num))
  } else {
    imp_res <- impute_numeric_matrix(
      X_num,
      method = na_impute,
      constant = na_constant,
      y = y,
      quantile_p = na_quantile,
      knn_k = knn_k,
      knn_scale = knn_scale,
      knn_weighted = knn_weighted,
      knn_max_n = knn_max_n,
      mice_m = mice_m,
      mice_maxit = mice_maxit,
      mice_method = mice_method,
      mice_seed = mice_seed,
      missforest_maxiter = missforest_maxiter,
      missforest_ntree = missforest_ntree
    )
    X_work <- imp_res$data
    y_work <- y
    imputed_counts <- imp_res$n_imputed
  }
  if (nrow(X_work) == 0L) stop("No data left after NA handling.")

  # Normalize outcome
  y_fac <- normalize_binary_response(y_work, positive)
  levs <- levels(y_fac)
  total_pos <- sum(y_fac == levs[2L])
  total_neg <- sum(y_fac == levs[1L])
  if (total_pos == 0L || total_neg == 0L) stop("y must contain at least one positive and one negative sample.")

  # Compute ROC per column
  res <- vector("list", ncol(X_work))
  cn <- colnames(X_work)
  for (j in seq_along(res)) {
    pred <- X_work[[j]]
    unique_pred <- unique(pred[!is.na(pred)])
    if (length(unique_pred) <= 1L) {
      res[[j]] <- data.frame(
        biomarker   = cn[j],
        auc         = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        ci_level    = ci_level, ci_method = ci_method,
        n_pos       = 0L, n_neg = 0L,
        direction   = NA_character_,
        na_handling = na_handling_label,
        n_na        = if (track_na) pre_na_counts[j] else NA_integer_,
        n_imputed   = if (track_na) imputed_counts[j] else NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    roc_obj <- try(roc(
      response  = y_fac,
      predictor = pred,
      levels    = levels(y_fac),
      direction = if (direction == "auto") "auto" else direction,
      quiet     = TRUE,
      na.rm     = TRUE
    ), silent = TRUE)

    if (inherits(roc_obj, "try-error") || is.null(roc_obj)) {
      res[[j]] <- data.frame(
        biomarker   = cn[j],
        auc         = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        ci_level    = ci_level, ci_method = ci_method,
        n_pos       = 0L, n_neg = 0L,
        direction   = NA_character_,
        na_handling = na_handling_label,
        n_na        = if (track_na) pre_na_counts[j] else NA_integer_,
        n_imputed   = if (track_na) imputed_counts[j] else NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    used_cases    <- if (!is.null(roc_obj$cases))    sum(!is.na(roc_obj$cases))    else sum(!is.na(pred) & y_fac == levs[2L])
    used_controls <- if (!is.null(roc_obj$controls)) sum(!is.na(roc_obj$controls)) else sum(!is.na(pred) & y_fac == levs[1L])

    auc_val <- suppressWarnings(as.numeric(auc(roc_obj)))
    ci_obj <- try({
      if (ci_method == "bootstrap") {
        ci.auc(roc_obj, conf.level = ci_level, method = ci_method, boot.n = boot_n)
      } else {
        ci.auc(roc_obj, conf.level = ci_level, method = ci_method)
      }
    }, silent = TRUE)
    if (inherits(ci_obj, "try-error") || is.null(ci_obj) || length(ci_obj) < 3L) {
      ci_low <- NA_real_; ci_high <- NA_real_
    } else {
      ci_low <- unname(ci_obj[1L]); ci_high <- unname(ci_obj[3L])
    }

    res[[j]] <- data.frame(
      biomarker   = cn[j],
      auc         = auc_val,
      ci_low      = ci_low,
      ci_high     = ci_high,
      ci_level    = ci_level,
      ci_method   = ci_method,
      n_pos       = used_cases,
      n_neg       = used_controls,
      direction   = if (!is.null(roc_obj$direction)) as.character(roc_obj$direction) else NA_character_,
      na_handling = na_handling_label,
      n_na        = if (track_na) pre_na_counts[j] else NA_integer_,
      n_imputed   = if (track_na) imputed_counts[j] else NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, res)

  if (isTRUE(rank)) {
    ord <- order(out$auc, decreasing = TRUE, na.last = TRUE)
    out <- out[ord, , drop = FALSE]
    out$rank <- seq_len(nrow(out))
    rownames(out) <- NULL
  }

  cols_to_round <- intersect(c("auc", "ci_low", "ci_high"), names(out))
  if (length(cols_to_round)) {
    out[cols_to_round] <- lapply(out[cols_to_round], function(col) {
      ifelse(is.na(col), NA, round(col, digits = round.digit))
    })
  }

  out
}

# Internal helpers (not exported) ------------------------------------------------

#' @noRd
#' @keywords internal
normalize_binary_response <- function(y, positive) {
  if (is.logical(y)) {
    return(factor(y, levels = c(FALSE, TRUE), labels = c("neg", "pos")))
  }
  if (is.numeric(y)) {
    u <- unique(y); u <- u[!is.na(u)]
    if (length(u) > 0L && all(u %in% c(0, 1))) {
      return(factor(y, levels = c(0, 1), labels = c("neg", "pos")))
    }
  }
  if (!is.factor(y)) y <- factor(y)
  if (nlevels(y) != 2L) stop("y must have exactly 2 unique values/classes (excluding NA).")
  lev <- levels(y)
  if (!is.null(positive)) {
    if (!(positive %in% lev)) stop("`positive` must be one of levels(y).")
    neg <- setdiff(lev, positive)
    if (length(neg) != 1L) stop("Could not infer negative level from y.")
    factor(y, levels = c(neg[1L], positive))
  } else {
    y
  }
}

#' @noRd
#' @keywords internal
impute_numeric_matrix <- function(
    df, method, constant = 0, y = NULL, quantile_p = 0.25,
    knn_k = 5, knn_scale = TRUE, knn_weighted = TRUE, knn_max_n = 5000,
    mice_m = 1, mice_maxit = 5, mice_method = "pmm", mice_seed = NULL,
    missforest_maxiter = 10, missforest_ntree = 100
) {
  stopifnot(is.data.frame(df))
  if (!all(vapply(df, is.numeric, logical(1)))) stop("All columns of df must be numeric for imputation.")

  if (method %in% c("mean", "median", "constant", "halfmin", "quantile")) {
    if (method == "mean")    return(impute_mean(df))
    if (method == "median")  return(impute_median(df))
    if (method == "constant")return(impute_constant(df, constant))
    if (method == "halfmin") return(impute_halfmin(df))
    if (method == "quantile")return(impute_quantile(df, quantile_p, constant))
  }

  if (method == "median_by_class") {
    if (is.null(y)) stop("y is required for 'median_by_class' imputation.")
    return(impute_by_class(df, y = y, stat = "median", constant = constant))
  }

  if (method == "knn") {
    return(impute_knn(df, k = knn_k, scale = knn_scale, weighted = knn_weighted,
                      max_n = knn_max_n, constant = constant))
  }

  if (method == "mice") {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Package 'mice' is required for na_impute = 'mice'. Please install it or choose another method.")
    }
    if (mice_m != 1) {
      warning("Only single imputation supported for 'mice'. Using mice_m = 1.")
      mice_m <- 1
    }
    old_seed <- .Random.seed
    on.exit({
      if (!is.null(old_seed)) .Random.seed <<- old_seed
    }, add = TRUE)
    if (!is.null(mice_seed)) set.seed(mice_seed)

    mice_mice     <- getExportedValue("mice", "mice")
    mice_complete <- getExportedValue("mice", "complete")

    imp <- mice_mice(df, m = mice_m, maxit = mice_maxit, method = mice_method, printFlag = FALSE)
    comp <- mice_complete(imp, 1)
    n_imp <- vapply(df, function(v) sum(is.na(v)), integer(1))
    return(list(data = comp, n_imputed = n_imp))
  }

  if (method == "missForest") {
    if (!requireNamespace("missForest", quietly = TRUE)) {
      stop("Package 'missForest' is required for na_impute = 'missForest'. Please install it or choose another method.")
    }
    # Drop all-NA and all-constant columns
    is_all_na <- vapply(df, function(x) all(is.na(x)), logical(1))
    is_constant <- vapply(df, function(x) {
      x_non_na <- x[!is.na(x)]
      length(x_non_na) > 0 && length(unique(x_non_na)) == 1
    }, logical(1))
    dropme <- is_all_na | is_constant
    if (all(dropme)) {
      stop("No usable columns for missForest: all columns are all-NA or constant.")
    }
    if (any(dropme)) {
      df <- df[, !dropme, drop = FALSE]
    }

    missForest_fun <- getExportedValue("missForest", "missForest")
    res <- missForest_fun(df, maxiter = missforest_maxiter, ntree = missforest_ntree, verbose = FALSE)
    n_imp <- vapply(df, function(v) sum(is.na(v)), integer(1))
    return(list(data = as.data.frame(res$ximp), n_imputed = n_imp))
  }

  stop("Unknown imputation method: ", method)
}

#' @noRd
#' @keywords internal
impute_mean <- function(df) {
  df_out <- df; n_imp <- integer(ncol(df))
  for (j in seq_len(ncol(df))) {
    v <- df[[j]]; idx <- is.na(v); n_imp[j] <- sum(idx)
    if (any(idx)) { fill <- mean(v, na.rm = TRUE); if (is.na(fill)) fill <- 0; v[idx] <- fill }
    df_out[[j]] <- v
  }
  list(data = df_out, n_imputed = n_imp)
}

#' @noRd
#' @keywords internal
impute_median <- function(df) {
  df_out <- df; n_imp <- integer(ncol(df))
  for (j in seq_len(ncol(df))) {
    v <- df[[j]]; idx <- is.na(v); n_imp[j] <- sum(idx)
    if (any(idx)) { fill <- stats::median(v, na.rm = TRUE); if (is.na(fill)) fill <- 0; v[idx] <- fill }
    df_out[[j]] <- v
  }
  list(data = df_out, n_imputed = n_imp)
}

#' @noRd
#' @keywords internal
impute_constant <- function(df, constant = 0) {
  df_out <- df; n_imp <- integer(ncol(df))
  for (j in seq_len(ncol(df))) {
    v <- df[[j]]; idx <- is.na(v); n_imp[j] <- sum(idx)
    if (any(idx)) v[idx] <- constant
    df_out[[j]] <- v
  }
  list(data = df_out, n_imputed = n_imp)
}

#' @noRd
#' @keywords internal
impute_halfmin <- function(df) {
  df_out <- df; n_imp <- integer(ncol(df))
  for (j in seq_len(ncol(df))) {
    v <- df[[j]]; idx <- is.na(v); n_imp[j] <- sum(idx)
    if (any(idx)) {
      pos <- v[v > 0 & !is.na(v)]
      fill <- if (length(pos) > 0L) min(pos) / 2 else min(v, na.rm = TRUE) / 2
      if (!is.finite(fill)) fill <- 0
      v[idx] <- fill
    }
    df_out[[j]] <- v
  }
  list(data = df_out, n_imputed = n_imp)
}

#' @noRd
#' @keywords internal
impute_quantile <- function(df, p = 0.25, constant = 0) {
  df_out <- df; n_imp <- integer(ncol(df))
  for (j in seq_len(ncol(df))) {
    v <- df[[j]]; idx <- is.na(v); n_imp[j] <- sum(idx)
    if (any(idx)) {
      qv <- stats::quantile(v, probs = p, na.rm = TRUE, type = 7, names = FALSE)
      if (is.na(qv)) qv <- constant
      v[idx] <- qv
    }
    df_out[[j]] <- v
  }
  list(data = df_out, n_imputed = n_imp)
}

#' @noRd
#' @keywords internal
impute_by_class <- function(df, y, stat = c("median", "mean"), constant = 0) {
  stat <- match.arg(stat)
  if (!is.factor(y)) y <- factor(y)
  if (nlevels(y) != 2L) stop("'y' must be a 2-level factor for class-based imputation.")
  df_out <- df; n_imp <- integer(ncol(df)); levs <- levels(y)
  for (j in seq_len(ncol(df))) {
    v <- df[[j]]; idx <- is.na(v); n_imp[j] <- sum(idx)
    if (!any(idx)) { df_out[[j]] <- v; next }
    if (stat == "median") {
      fill0 <- stats::median(v[y == levs[1L]], na.rm = TRUE)
      fill1 <- stats::median(v[y == levs[2L]], na.rm = TRUE)
    } else {
      fill0 <- mean(v[y == levs[1L]], na.rm = TRUE)
      fill1 <- mean(v[y == levs[2L]], na.rm = TRUE)
    }
    glob <- if (stat == "median") stats::median(v, na.rm = TRUE) else mean(v, na.rm = TRUE)
    if (is.na(glob)) glob <- constant
    if (is.na(fill0)) fill0 <- glob
    if (is.na(fill1)) fill1 <- glob
    idx0 <- idx & (y == levs[1L]); idx1 <- idx & (y == levs[2L])
    v[idx0] <- fill0; v[idx1] <- fill1
    df_out[[j]] <- v
  }
  list(data = df_out, n_imputed = n_imp)
}

#' @noRd
#' @keywords internal
impute_knn <- function(df, k = 5, scale = TRUE, weighted = TRUE, max_n = 5000, constant = 0) {
  n <- nrow(df); p <- ncol(df)
  if (n > max_n) {
    stop("KNN imputation requires an n x n distance matrix. n = ", n,
         " exceeds knn_max_n = ", max_n, ". Use a different method or increase knn_max_n.")
  }
  mat <- as.matrix(df)
  if (scale) {
    means <- colMeans(mat, na.rm = TRUE)
    sds <- apply(mat, 2, stats::sd, na.rm = TRUE); sds[!is.finite(sds) | sds == 0] <- 1
    mat_scaled <- sweep(mat, 2, means, "-"); mat_scaled <- sweep(mat_scaled, 2, sds, "/")
  } else {
    mat_scaled <- mat
  }
  col_means_scaled <- colMeans(mat_scaled, na.rm = TRUE)
  mat_scaled_imp <- mat_scaled
  for (j in seq_len(p)) {
    m <- col_means_scaled[j]; if (is.na(m)) m <- 0
    na_j <- is.na(mat_scaled_imp[, j]); if (any(na_j)) mat_scaled_imp[na_j, j] <- m
  }
  D <- as.matrix(stats::dist(mat_scaled_imp, method = "euclidean"))
  df_out <- df; n_imp <- integer(p)
  for (j in seq_len(p)) {
    v <- df[[j]]; na_rows <- which(is.na(v)); n_imp[j] <- length(na_rows)
    if (length(na_rows) == 0L) { df_out[[j]] <- v; next }
    candidates <- which(!is.na(v))
    if (length(candidates) == 0L) { v[na_rows] <- constant; df_out[[j]] <- v; next }
    for (i in na_rows) {
      cand_d <- D[i, candidates]
      ord <- order(cand_d, decreasing = FALSE, na.last = NA)
      if (length(ord) == 0L) { v[i] <- constant; next }
      k_use <- min(k, length(ord)); nn_idx <- candidates[ord[seq_len(k_use)]]
      nn_vals <- v[nn_idx]
      if (weighted) {
        w <- 1 / (cand_d[ord[seq_len(k_use)]] + 1e-8)
        v[i] <- sum(w * nn_vals) / sum(w)
      } else {
        v[i] <- mean(nn_vals)
      }
    }
    df_out[[j]] <- v
  }
  list(data = as.data.frame(df_out), n_imputed = n_imp)
}
