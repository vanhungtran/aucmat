# ==============================================================================
# Deprecated Functions
#
# These functions are retained for backward compatibility but will be removed
# in a future version.  Use aucmat() instead of tableroc().
# ==============================================================================

#' @title (Deprecated) AUC table for each numeric predictor column
#' @description `tableroc()` is deprecated.  Use [aucmat()] for matrix
#'   screening.  Internal imputation is no longer supported; impute your
#'   data before passing it to `aucmat()`.
#' @export
#' @keywords internal
tableroc <- function(X, y, ...) {
  .Deprecated("aucmat", package = "aucmat",
    msg = paste(
      "tableroc() is deprecated. Use aucmat() for matrix screening.",
      "Internal imputation options are no longer supported;",
      "impute your data before calling aucmat()."
    ))
  # Forward compatible arguments to aucmat()
  dots <- list(...)
  args <- list(X = X, y = y)
  if (!is.null(dots$positive))   args$positive   <- dots$positive
  if (!is.null(dots$ci_method))  args$ci         <- dots$ci_method
  if (!is.null(dots$ci_level))   args$conf_level <- dots$ci_level
  if (!is.null(dots$boot_n))     args$boot_n     <- dots$boot_n
  if (!is.null(dots$na_rm)) {
    args$na_action <- if (isTRUE(dots$na_rm)) "complete" else "featurewise"
  }
  args$adjust <- "none"
  do.call(aucmat, args)
}

#' @title (Deprecated) Plot ROC curves for predictor combinations
#' @description `plot_roc_with_combos()` is deprecated and will be removed
#'   in a future version.  Use [plot_roc_top()] for focused ROC
#'   visualisation of selected biomarkers.
#' @export
#' @keywords internal
plot_roc_with_combos <- function(
    data, outcome, predictors,
    show_auc_labels = FALSE, label_fpr = 0.25, label_size = 3.5,
    label_nudge_x = 0.02, combo_sizes = NULL, positive = NULL,
    add_ci = TRUE, conf_level = 0.95, boot_n = 2000, ci_points = 201,
    legacy_axes = FALSE, ribbon_alpha = 0.20, line_size = 0.9,
    palette = NULL, reorder_by_auc = TRUE, max_curves = Inf,
    title = "ROC: Sensitivity vs 1 - Specificity with CI ribbons",
    seed = NULL, skip_ci_for_auc1 = TRUE, auc1_eps = 1e-12,
    auc_digits = 3, legend_position = "right", legend_ncol = 1,
    use_ggrepel = FALSE, show_auc_sd = TRUE, sd_digits = NULL,
    show_chance_line = FALSE) {

  .Deprecated("plot_roc_top",
    msg = paste(
      "plot_roc_with_combos() is deprecated.",
      "Use plot_roc_top() for focused ROC visualisation of selected biomarkers."
    ))

  stopifnot(is.data.frame(data), is.character(outcome), is.character(predictors))
  if (is.null(combo_sizes)) combo_sizes <- seq_along(predictors)
  if (!all(c(outcome, predictors) %in% names(data))) {
    stop("Some columns are missing in data.")
  }

  subsets <- unlist(lapply(combo_sizes, function(k) {
    if (k <= 0 || k > length(predictors)) return(list())
    utils::combn(predictors, k, simplify = FALSE)
  }), recursive = FALSE)

  build_roc_for_subset <- function(vars) {
    nm <- paste(vars, collapse = " + ")
    mf <- stats::model.frame(
      stats::as.formula(paste(outcome, "~", paste(vars, collapse = " + "))),
      data = data, na.action = stats::na.omit
    )
    y_val <- mf[[outcome]]
    mf_pred <- mf[, vars, drop = FALSE]
    if (is.numeric(y_val) || is.logical(y_val)) y_val <- factor(y_val)
    if (!is.factor(y_val) || length(levels(y_val)) != 2)
      stop("Outcome must be binary.")

    lv <- levels(y_val)
    pos <- if (is.null(positive)) lv[2] else positive
    if (!pos %in% lv) stop("Provided 'positive' not found in outcome levels.")
    neg <- setdiff(lv, pos)[1]

    y_glm <- stats::relevel(y_val, ref = pos)
    df_fit <- stats::model.frame(stats::as.formula(paste("y_glm ~ .")),
      data = cbind(y_glm = y_glm, mf_pred))
    fit <- suppressWarnings(stats::glm(y_glm ~ ., data = df_fit,
      family = stats::binomial()))
    prob_pos <- stats::predict(fit, newdata = mf_pred, type = "response")

    roc_obj <- pROC::roc(response = y_val, predictor = prob_pos,
      levels = c(neg, pos), direction = ">", quiet = TRUE)
    list(name = nm, roc = roc_obj)
  }

  roc_items <- lapply(subsets, build_roc_for_subset)
  names(roc_items) <- vapply(roc_items, `[[`, "", "name")
  roc_list <- stats::setNames(lapply(roc_items, `[[`, "roc"), names(roc_items))

  aucs <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)), numeric(1))
  ord <- if (isTRUE(reorder_by_auc)) order(aucs, decreasing = TRUE) else seq_along(aucs)
  if (is.finite(max_curves) && max_curves < length(ord)) ord <- ord[seq_len(max_curves)]
  roc_list <- roc_list[ord]
  aucs <- aucs[ord]

  is_perfect <- (abs(aucs - 1) <= auc1_eps)
  auc_sds <- rep(NA_real_, length(roc_list))
  auc_cis <- matrix(NA_real_, nrow = length(roc_list), ncol = 2)

  if (isTRUE(add_ci)) {
    if (!is.null(seed)) set.seed(seed)
    which_auc_ci <- seq_along(roc_list)
    if (isTRUE(skip_ci_for_auc1)) which_auc_ci <- which(!is_perfect)
    if (length(which_auc_ci) > 0) {
      for (i in which_auc_ci) {
        r <- roc_list[[i]]
        auc_ci_obj <- try(pROC::ci.auc(r, conf.level = conf_level,
          boot.n = boot_n), silent = TRUE)
        if (!inherits(auc_ci_obj, "try-error")) {
          auc_cis[i, ] <- c(auc_ci_obj[1], auc_ci_obj[3])
          z_val <- stats::qnorm((1 + conf_level) / 2)
          auc_sds[i] <- (auc_ci_obj[3] - auc_ci_obj[1]) / (2 * z_val)
        }
      }
    }
  }

  auc_fmt <- function(x) formatC(x, auc_digits, format = "f")
  if (is.null(sd_digits)) sd_digits <- auc_digits

  legend_labels <- vapply(seq_along(roc_list), function(i) {
    nm <- names(roc_list)[i]
    auc_str <- auc_fmt(aucs[i])
    if (isTRUE(add_ci) && isTRUE(show_auc_sd)) {
      if (is_perfect[i] && isTRUE(skip_ci_for_auc1)) {
        paste0(nm, " (AUC=", auc_str, ", no CI)")
      } else if (!is.na(auc_sds[i])) {
        sd_str <- formatC(auc_sds[i], digits = sd_digits, format = "f")
        paste0(nm, " (AUC=", auc_str, " +/- ", sd_str, ")")
      } else {
        paste0(nm, " (AUC=", auc_str, ")")
      }
    } else {
      paste0(nm, " (AUC=", auc_str, ")")
    }
  }, character(1))

  if (is.null(palette)) {
    palette <- if (requireNamespace("scales", quietly = TRUE)) {
      scales::hue_pal()(length(roc_list))
    } else grDevices::rainbow(length(roc_list))
  } else if (length(palette) < length(roc_list)) {
    palette <- rep(palette, length.out = length(roc_list))
  }

  ggroc_obj <- pROC::ggroc(roc_list, legacy.axes = legacy_axes, size = line_size)
  p <- ggroc_obj +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  if (isTRUE(show_chance_line)) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0,
      linetype = "dashed", alpha = 0.7, color = "grey50")
  }
  p <- p +
    ggplot2::labs(x = "1 - Specificity (FPR)", y = "Sensitivity", title = title) +
    ggplot2::scale_colour_manual(values = palette, labels = legend_labels,
      name = "Model (AUC)") +
    ggplot2::theme(legend.position = legend_position) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = legend_ncol,
      override.aes = list(fill = NA)))

  ci_df <- NULL
  if (isTRUE(add_ci)) {
    if (!is.null(seed)) set.seed(seed)
    spec_grid <- seq(1, 0, length.out = ci_points)
    which_ci <- seq_along(roc_list)
    if (isTRUE(skip_ci_for_auc1)) which_ci <- which(!is_perfect)
    if (length(which_ci)) {
      ci_list <- lapply(which_ci, function(i) {
        r <- roc_list[[i]]
        pROC::ci.se(r, specificities = spec_grid, conf.level = conf_level,
          boot.n = boot_n)
      })
      ci_df <- do.call(rbind, lapply(seq_along(ci_list), function(j) {
        i <- which_ci[j]
        ci <- ci_list[[j]]
        spec <- as.numeric(rownames(ci))
        x_coord <- if (isTRUE(legacy_axes)) (1 - spec) else spec
        o <- order(x_coord)
        data.frame(name = names(roc_list)[i], x = x_coord[o],
          lower = ci[o, 1], upper = ci[o, 3], stringsAsFactors = FALSE)
      }))
      if (!is.null(ci_df) && nrow(ci_df) > 0) {
        ci_df$name <- factor(ci_df$name, levels = names(roc_list))
      }
    }
  }

  if (!is.null(ci_df) && nrow(ci_df)) {
    p <- p +
      ggplot2::geom_ribbon(data = ci_df,
        ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper,
          fill = .data$name),
        alpha = ribbon_alpha, inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = palette, guide = "none")
  }

  if (isTRUE(show_auc_labels)) {
    lbls <- do.call(rbind, lapply(seq_along(roc_list), function(i) {
      r <- roc_list[[i]]
      spec_req <- 1 - label_fpr
      cp <- try(pROC::coords(r, x = spec_req, input = "specificity",
        ret = c("specificity", "sensitivity"), transpose = FALSE), silent = TRUE)
      if (inherits(cp, "try-error")) return(NULL)
      if (is.matrix(cp)) {
        if (nrow(cp) < 1) return(NULL)
        cp_row <- cp[1, , drop = TRUE]
      } else if (is.numeric(cp) && length(cp) >= 1) {
        cp_row <- cp
      } else {
        return(NULL)
      }
      if (any(is.na(cp_row))) return(NULL)
      spec_val <- unname(cp_row["specificity"])
      x_val <- if (isTRUE(legacy_axes)) (1 - spec_val) else spec_val
      data.frame(name = names(roc_list)[i], x = x_val,
        sens = unname(cp_row["sensitivity"]),
        label = paste0("AUC=", auc_fmt(aucs[i])), stringsAsFactors = FALSE)
    }))
    if (!is.null(lbls) && nrow(lbls) > 0) {
      lbls$name <- factor(lbls$name, levels = names(roc_list))
      if (isTRUE(use_ggrepel) && requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(data = lbls,
          ggplot2::aes(x = .data$x, y = .data$sens, label = .data$label,
            colour = .data$name),
          inherit.aes = FALSE, size = label_size, show.legend = FALSE)
      } else {
        p <- p + ggplot2::geom_text(data = lbls,
          ggplot2::aes(x = .data$x, y = .data$sens, label = .data$label,
            colour = .data$name),
          inherit.aes = FALSE, size = label_size,
          position = ggplot2::position_nudge(x = label_nudge_x),
          show.legend = FALSE)
      }
    }
  }

  list(plot = p, roc_list = roc_list, auc = aucs,
    auc_sd = stats::setNames(auc_sds, names(roc_list)),
    ci_df = ci_df, skipped_ci = names(roc_list)[is_perfect])
}
