# Requires: pROC, ggplot2
plot_roc_with_combos <- function(
    data,
    outcome,
    predictors,
    combo_sizes = 1:length(predictors),
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
    auc1_eps = 1e-12                  # <- tolerance for detecting perfect AUC
) {
  stopifnot(is.data.frame(data), is.character(outcome), is.character(predictors))
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
    y <- mf[[outcome]]
    if (is.numeric(y) || is.logical(y)) y <- factor(y)
    if (!is.factor(y) || length(levels(y)) != 2) stop("Outcome must be binary.")

    lv <- levels(y)
    pos <- if (is.null(positive)) lv[2] else positive
    if (!pos %in% lv) stop("Provided 'positive' not found in outcome levels.")
    neg <- setdiff(lv, pos)[1]

    y_glm <- stats::relevel(y, ref = pos)  # make 'pos' the success
    fit <- suppressWarnings(stats::glm(y_glm ~ ., data = mf, family = stats::binomial()))
    prob_pos <- stats::predict(fit, type = "response")

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

  # Legend labels (tag AUC == 1 with "[no CI]" if skipping CIs)
  is_perfect <- (abs(aucs - 1) <= auc1_eps)
  legend_labels <- if (isTRUE(add_ci) && isTRUE(skip_ci_for_auc1)) {
    paste0(names(roc_list), " (AUC=", formatC(aucs, 3, format = "f"),
           ifelse(is_perfect, ", no CI", ""), ")")
  } else {
    paste0(names(roc_list), " (AUC=", formatC(aucs, 3, format = "f"), ")")
  }

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
      # Bind with names and transform x to FPR
      ci_df <- do.call(rbind, lapply(seq_along(ci_list), function(j) {
        i <- which_ci[j]
        ci <- ci_list[[j]]
        fpr <- 1 - as.numeric(rownames(ci))
        o <- order(fpr)
        data.frame(
          name  = names(roc_list)[i],
          fpr   = fpr[o],
          lower = ci[o, 1],
          upper = ci[o, 3],
          stringsAsFactors = FALSE
        )
      }))
    }
  }

  # Base plot: x = 1 - Specificity (FPR), y = Sensitivity
  p <- pROC::ggroc(roc_list, legacy.axes = legacy_axes, size = line_size) +
    ggplot2::theme_minimal() +
    ggplot2::coord_equal() +
    # Use coord_cartesian to avoid adding duplicate x/y scales (no "Scale for x..." warning)
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         alpha = 0.7, color = "grey50") +
    ggplot2::labs(x = "1 - Specificity (FPR)", y = "Sensitivity", title = title) +
    ggplot2::scale_colour_manual(values = palette, labels = legend_labels, name = NULL)

  if (!is.null(ci_df) && nrow(ci_df)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = ci_df,
        ggplot2::aes(x = fpr, ymin = lower, ymax = upper, fill = name),
        alpha = ribbon_alpha,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(values = palette, guide = "none")
  }

  list(plot = p, roc_list = roc_list, auc = aucs, ci_df = ci_df, skipped_ci = names(roc_list)[is_perfect])
}




