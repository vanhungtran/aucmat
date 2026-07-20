#' Plot ROC curves for predictor combinations with optional CI ribbons and AUC labels
#'
#' Build ROC curves for combinations of predictors by fitting a binomial GLM
#' for each subset, predicting probabilities, and plotting ROC curves with
#' optional confidence ribbons and AUC labels.
#'
#' @param data data.frame containing outcome and predictor columns
#' @param outcome character; name of the outcome column (binary)
#' @param predictors character vector of predictor column names
#' @param show_auc_labels logical; add AUC text labels on the plot
#' @param label_fpr numeric; FPR position (0-1) where labels are placed
#' @param label_size numeric; label text size
#' @param label_nudge_x numeric; horizontal nudge for labels
#' @param combo_sizes integer vector of combination sizes to evaluate (NULL => all)
#' @param positive optional positive class label
#' @param add_ci logical; add CI ribbons
#' @param conf_level numeric; CI level
#' @param boot_n integer; bootstrap reps for CI
#' @param ci_points integer; grid points for CI ribbons
#' @param legacy_axes logical; pass to pROC::ggroc
#' @param ribbon_alpha numeric; alpha for CI ribbons
#' @param line_size numeric; ROC line width
#' @param palette character vector of colours
#' @param reorder_by_auc logical; reorder curves by AUC
#' @param max_curves maximum number of curves to show
#' @param title plot title
#' @param seed optional seed for CI reproducibility
#' @param skip_ci_for_auc1 logical; skip CI for (near-)perfect AUCs
#' @param auc1_eps numeric tolerance for detecting AUC == 1
#' @param legend_position legend position (default "right")
#' @param legend_ncol integer number of legend columns
#' @param use_ggrepel logical; if TRUE and `ggrepel` available, use geom_text_repel for labels
#' @param show_auc_sd logical; if TRUE and add_ci is TRUE, show AUC ± SD in legend (default TRUE)
#' @param sd_digits integer; number of digits for SD in legend (default: same as auc_digits)
#'
#' @return list with `plot` (ggplot), `roc_list`, `auc`, `ci_df`, `skipped_ci`, `auc_sd`
#' @importFrom pROC ggroc roc auc ci.se coords ci.auc
#' @importFrom ggplot2 ggplot geom_ribbon geom_text theme labs coord_equal coord_cartesian geom_abline scale_colour_manual scale_fill_manual guides guide_legend theme_minimal
#' @importFrom rlang .data
#' @export
#' @examples
#' \donttest{
#' # Generate example data
#' set.seed(123)
#' sim_data <- generate_data_analytical(
#'   n = 200,
#'   prevalence = 0.3,
#'   target_aucs = c(0.8, 0.7, 0.6),
#'   corr_matrix = diag(3)
#' )
#'
#' # Get predictor names
#' preds <- setdiff(names(sim_data$data), "truth")
#'
#' # Create ROC plot with legend showing AUC ± SD
#' res <- plot_roc_with_combos(
#'   data = sim_data$data,
#'   outcome = "truth",
#'   predictors = preds,
#'   combo_sizes = 1:2,  # Single predictors and pairs
#'   add_ci = TRUE,
#'   show_auc_labels = TRUE,
#'   auc_digits = 4,
#'   show_auc_sd = TRUE   # Show SD in legend (default)
#'   # sd_digits defaults to auc_digits (4 in this case)
#' )
#'
#' # Display the plot (legend shows "Model (AUC=0.xxxx ± 0.xxxx)")
#' print(res$plot)
#'
#' # Access AUC values and SDs
#' res$auc
#' res$auc_sd
#' }
#'
plot_roc_with_combos <- function(
  data,
  outcome,
  predictors,
  show_auc_labels = FALSE,
  label_fpr = 0.25,
  label_size = 3.5,
  label_nudge_x = 0.02,
  combo_sizes = NULL,
    positive = NULL,
    add_ci = TRUE,
    conf_level = 0.95,
    boot_n = 2000,
    ci_points = 201,
    legacy_axes = FALSE,              # FALSE => x = 1 - Specificity (FPR), left->right
    ribbon_alpha = 0.20,
    line_size = 0.9,
    palette = NULL,
    reorder_by_auc = TRUE,
    max_curves = Inf,
    title = "ROC: Sensitivity vs 1 - Specificity with CI ribbons",
    seed = NULL,
  skip_ci_for_auc1 = TRUE,          # <- skip CI when AUC ≈ 1
  auc1_eps = 1e-12,                 # <- tolerance for detecting perfect AUC
  auc_digits = 3,
  legend_position = "right",
  legend_ncol = 1,
  use_ggrepel = FALSE,
  show_auc_sd = TRUE,
  sd_digits = NULL,
  show_chance_line = FALSE
) {
  stopifnot(is.data.frame(data), is.character(outcome), is.character(predictors))
  # Default combo_sizes to all possible sizes if not provided
  if (is.null(combo_sizes)) combo_sizes <- seq_along(predictors)
  if (!all(c(outcome, predictors) %in% names(data))) {
    stop("Some columns are missing in data.")
  }

  # Build combinations of predictors
  subsets <- unlist(lapply(combo_sizes, function(k) {
    if (k <= 0 || k > length(predictors)) return(list())
    combn(predictors, k, simplify = FALSE)
  }), recursive = FALSE)

  # Helper: ROC from subset via binomial GLM -> predicted probabilities
  build_roc_for_subset <- function(vars) {
    nm <- paste(vars, collapse = " + ")
    mf <- stats::model.frame(
      stats::as.formula(paste(outcome, "~", paste(vars, collapse = " + "))),
      data = data, na.action = stats::na.omit
    )
    # Extract response and predictor data separately to avoid leaking the
    # outcome column into the predictors (which would produce perfect
    # predictions and incorrect AUCs).
    y <- mf[[outcome]]
    mf_pred <- mf[, vars, drop = FALSE]
    if (is.numeric(y) || is.logical(y)) y <- factor(y)
    if (!is.factor(y) || length(levels(y)) != 2) stop("Outcome must be binary.")

    lv <- levels(y)
    pos <- if (is.null(positive)) lv[2] else positive
    if (!pos %in% lv) stop("Provided 'positive' not found in outcome levels.")
    neg <- setdiff(lv, pos)[1]

    y_glm <- stats::relevel(y, ref = pos)  # make 'pos' the success
    # Build a fit data.frame including the response to avoid referencing
    # variables outside the data argument of glm. Predict on the predictor
    # frame to obtain predicted probabilities.
    df_fit <- stats::model.frame(stats::as.formula(paste("y_glm ~ .")),
                   data = cbind(y_glm = y_glm, mf_pred))
    fit <- suppressWarnings(stats::glm(y_glm ~ ., data = df_fit, family = stats::binomial()))
    prob_pos <- stats::predict(fit, newdata = mf_pred, type = "response")

    roc_obj <- pROC::roc(response = y, predictor = prob_pos,
                         levels = c(neg, pos), direction = ">", quiet = TRUE)
    list(name = nm, roc = roc_obj)
  }

  roc_items <- lapply(subsets, build_roc_for_subset)
  names(roc_items) <- vapply(roc_items, `[[`, "", "name")
  roc_list <- stats::setNames(lapply(roc_items, `[[`, "roc"), names(roc_items))

  # AUCs and ordering
  aucs <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)), numeric(1))
  ord <- if (isTRUE(reorder_by_auc)) order(aucs, decreasing = TRUE) else seq_along(aucs)
  if (is.finite(max_curves) && max_curves < length(ord)) ord <- ord[seq_len(max_curves)]
  roc_list <- roc_list[ord]
  aucs <- aucs[ord]

  # Compute AUC CIs and SDs for legend labels (if add_ci is TRUE)
  is_perfect <- (abs(aucs - 1) <= auc1_eps)
  auc_sds <- rep(NA_real_, length(roc_list))
  auc_cis <- matrix(NA_real_, nrow = length(roc_list), ncol = 2)

  if (isTRUE(add_ci)) {
    if (!is.null(seed)) set.seed(seed)
    # Compute AUC CIs for non-perfect curves
    which_auc_ci <- seq_along(roc_list)
    if (isTRUE(skip_ci_for_auc1)) which_auc_ci <- which(!is_perfect)

    if (length(which_auc_ci) > 0) {
      for (i in which_auc_ci) {
        r <- roc_list[[i]]
        auc_ci_obj <- try(pROC::ci.auc(r, conf.level = conf_level, boot.n = boot_n), silent = TRUE)
        if (!inherits(auc_ci_obj, "try-error")) {
          auc_cis[i, ] <- c(auc_ci_obj[1], auc_ci_obj[3])
          # Estimate SD from CI: SD ≈ (upper - lower) / (2 * z)
          # For 95% CI, z ≈ 1.96
          z_val <- stats::qnorm((1 + conf_level) / 2)
          auc_sds[i] <- (auc_ci_obj[3] - auc_ci_obj[1]) / (2 * z_val)
        }
      }
    }
  }

  # Format AUC values for legend labels using requested digits
  auc_fmt <- function(x) formatC(x, auc_digits, format = "f")

  # Determine SD digits (default to same as auc_digits for consistency)
  if (is.null(sd_digits)) {
    sd_digits <- auc_digits
  }

  # Build legend labels with AUC ± SD format
  legend_labels <- vapply(seq_along(roc_list), function(i) {
    nm <- names(roc_list)[i]
    auc_str <- auc_fmt(aucs[i])

    if (isTRUE(add_ci) && isTRUE(show_auc_sd)) {
      if (is_perfect[i] && isTRUE(skip_ci_for_auc1)) {
        # Perfect AUC, no CI computed
        paste0(nm, " (AUC=", auc_str, ", no CI)")
      } else if (!is.na(auc_sds[i])) {
        # AUC with SD - use same digits as AUC for consistency
        sd_str <- formatC(auc_sds[i], digits = sd_digits, format = "f")
        paste0(nm, " (AUC=", auc_str, " \u00b1 ", sd_str, ")")
      } else {
        # CI failed, just show AUC
        paste0(nm, " (AUC=", auc_str, ")")
      }
    } else {
      # No CI requested or show_auc_sd is FALSE
      paste0(nm, " (AUC=", auc_str, ")")
    }
  }, character(1))

  # Palette
  if (is.null(palette)) {
    palette <- if (requireNamespace("scales", quietly = TRUE)) {
      scales::hue_pal()(length(roc_list))
    } else grDevices::rainbow(length(roc_list))
  } else if (length(palette) < length(roc_list)) {
    palette <- rep(palette, length.out = length(roc_list))
  }

  # CI ribbons: vertical uncertainty in Sensitivity across specificities
  ci_df <- NULL
  if (isTRUE(add_ci)) {
    if (!is.null(seed)) set.seed(seed)
    spec_grid <- seq(1, 0, length.out = ci_points)  # decreasing so FPR increases left->right

    # Compute only for non-perfect AUCs if skipping
    which_ci <- seq_along(roc_list)
    if (isTRUE(skip_ci_for_auc1)) which_ci <- which(!is_perfect)

    if (length(which_ci)) {
      ci_list <- lapply(which_ci, function(i) {
        r <- roc_list[[i]]
        # ci.se bootstraps; skips pROC's warning by not calling for AUC==1 curves
        pROC::ci.se(r, specificities = spec_grid, conf.level = conf_level, boot.n = boot_n)
      })
        # Bind with names and use appropriate x coordinate
        # When legacy.axes = FALSE, ggroc uses specificity on x-axis (despite label)
        # When legacy.axes = TRUE, ggroc uses 1-specificity (FPR) on x-axis
      ci_df <- do.call(rbind, lapply(seq_along(ci_list), function(j) {
        i <- which_ci[j]
        ci <- ci_list[[j]]
        spec <- as.numeric(rownames(ci))
        x_coord <- if (isTRUE(legacy_axes)) (1 - spec) else spec
        o <- order(x_coord)
        data.frame(
          name  = names(roc_list)[i],
          x     = x_coord[o],
          lower = ci[o, 1],
          upper = ci[o, 3],
          stringsAsFactors = FALSE
        )
      }))
        # Ensure factor levels/order match the roc_list so colours/fills align
        if (!is.null(ci_df) && nrow(ci_df) > 0) {
          ci_df$name <- factor(ci_df$name, levels = names(roc_list))
        }
    }
  }

  # Start with ggroc to get proper ROC curves with correct aesthetics
  ggroc_obj <- pROC::ggroc(roc_list, legacy.axes = legacy_axes, size = line_size)

  # Build on top of ggroc plot
  p <- ggroc_obj +
    ggplot2::theme_minimal() +
    ggplot2::coord_equal() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  if (isTRUE(show_chance_line)) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         alpha = 0.7, color = "grey50")
  }
  p <- p + ggplot2::labs(x = "1 - Specificity (FPR)", y = "Sensitivity", title = title) +
    ggplot2::scale_colour_manual(
      values = palette,
      labels = legend_labels,
      name = "Model (AUC)"
    ) +
    ggplot2::theme(legend.position = legend_position) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = legend_ncol, override.aes = list(fill = NA)))

  # Add CI ribbons if available
  if (!is.null(ci_df) && nrow(ci_df)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = ci_df,
        ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper, fill = .data$name),
        alpha = ribbon_alpha,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(values = palette, guide = "none")
  }

  # Optionally add text labels showing AUC near the curves. We compute a
  # sensible label position by interpolating the ROC at a requested FPR
  # (1 - specificity) via pROC::coords().
  if (isTRUE(show_auc_labels)) {
    lbls <- do.call(rbind, lapply(seq_along(roc_list), function(i) {
      r <- roc_list[[i]]
      # target specificity corresponding to requested FPR
      spec_req <- 1 - label_fpr
      # coords will interpolate if needed; may return vector or matrix
      cp <- try(pROC::coords(r, x = spec_req, input = "specificity",
                            ret = c("specificity", "sensitivity"), transpose = FALSE), silent = TRUE)
      if (inherits(cp, "try-error")) return(NULL)
      # Handle possible shapes: numeric vector or 1-row matrix
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
      # Use same coordinate system as the plot
      x_val <- if (isTRUE(legacy_axes)) (1 - spec_val) else spec_val
      data.frame(
        name = names(roc_list)[i],
        x = x_val,
        sens = unname(cp_row["sensitivity"]),
        label = paste0("AUC=", auc_fmt(aucs[i])),
        stringsAsFactors = FALSE
      )
    }))
    if (!is.null(lbls) && nrow(lbls) > 0) {
      # match factor levels to ensure colours align
      lbls$name <- factor(lbls$name, levels = names(roc_list))
      if (isTRUE(use_ggrepel) && requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = lbls,
          ggplot2::aes(x = .data$x, y = .data$sens, label = .data$label, colour = .data$name),
          inherit.aes = FALSE,
          size = label_size,
          show.legend = FALSE
        )
      } else {
        p <- p + ggplot2::geom_text(
          data = lbls,
          ggplot2::aes(x = .data$x, y = .data$sens, label = .data$label, colour = .data$name),
          inherit.aes = FALSE,
          size = label_size,
          position = ggplot2::position_nudge(x = label_nudge_x),
          show.legend = FALSE
        )
      }
    }
  }

  list(
    plot = p,
    roc_list = roc_list,
    auc = aucs,
    auc_sd = stats::setNames(auc_sds, names(roc_list)),
    ci_df = ci_df,
    skipped_ci = names(roc_list)[is_perfect]
  )
}




