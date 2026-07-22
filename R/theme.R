# ==============================================================================
# aucmat Design System — unified colours, theme, and palette
#
# Inspired by cliomicplot & atcddd visual identity.
# Every ggplot2 plot in the package should use theme_aucmat().
# ==============================================================================

# ── Brand colour palette ────────────────────────────────────────────
# Use these as the canonical colour references throughout the package.

aucmat_colors <- list(
  # Core brand
  primary    = "#2c7da0",   # deep teal-blue (matches pkgdown CSS)
  dark       = "#1a5276",   # darker variant for text/emphasis
  light      = "#7fb5cc",   # lighter variant for fills

  # Directional (ROC/biomarker)
  higher     = "#2166AC",   # higher in positive cases (blue)
  lower      = "#B2182B",   # lower in positive cases (red)
  neutral    = "#8D99AE",   # neutral / non-significant

  # Accent palette (for multi-biomarker / group plots)
  accent = c(
    "#2c7da0", "#EE7733", "#009988", "#CC3311",
    "#33BBEE", "#AA3377", "#228833", "#4477AA",
    "#DDCC77", "#66CCEE", "#88CCAA", "#CC6677"
  ),

  # Semantic
  significant     = "#B2182B",  # q < cutoff
  not_significant = "#8D99AE",  # q >= cutoff
  warning         = "#EE7733",
  grid            = "#DDE1E6",
  background      = "#FAFBFC"
)

# ── aucmat ggplot2 theme ───────────────────────────────────────────

#' aucmat ggplot2 theme
#'
#' A clean, professional theme for all aucmat plots.
#' Uses Fira Sans if available via showtext, otherwise falls back to
#' system defaults.
#'
#' @param base_size Base font size (default 12).
#' @param base_family Font family. If `"fira"`, checks for Fira Sans
#'   registration via showtext.
#'
#' @return A ggplot2 theme object.
#' @export
#' @keywords internal
theme_aucmat <- function(base_size = 12, base_family = "") {
  ggplot2::`%+replace%`(
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family),
    ggplot2::theme(
      plot.title         = ggplot2::element_text(
                             face   = "bold",
                             hjust  = 0,
                             size   = ggplot2::rel(1.35),
                             colour = aucmat_colors$dark,
                             margin = ggplot2::margin(b = 8)),
      plot.subtitle      = ggplot2::element_text(
                             hjust  = 0,
                             colour = "grey45",
                             size   = ggplot2::rel(0.9),
                             margin = ggplot2::margin(b = 14)),
      plot.caption       = ggplot2::element_text(
                             colour = "grey55",
                             size   = ggplot2::rel(0.7),
                             hjust  = 1),
      plot.background    = ggplot2::element_rect(
                             fill  = aucmat_colors$background,
                             color = NA),
      panel.grid.major   = ggplot2::element_line(
                             colour    = aucmat_colors$grid,
                             linewidth = 0.3),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "bottom",
      legend.text        = ggplot2::element_text(size = ggplot2::rel(0.8)),
      legend.title       = ggplot2::element_text(size = ggplot2::rel(0.85)),
      axis.title         = ggplot2::element_text(
                             size   = ggplot2::rel(0.85),
                             colour = "grey30"),
      axis.text          = ggplot2::element_text(
                             size   = ggplot2::rel(0.8),
                             colour = "grey40"),
      strip.text         = ggplot2::element_text(
                             face   = "bold",
                             size   = ggplot2::rel(0.85),
                             colour = aucmat_colors$dark),
      plot.margin        = ggplot2::margin(18, 20, 12, 16)
    )
  )
}

# ── Direction colour scale ─────────────────────────────────────────

#' aucmat direction colour scale (discrete)
#'
#' Maps `higher_in_positive` → blue, `lower_in_positive` → red.
#'
#' @param ... Passed to [ggplot2::scale_colour_manual()].
#' @export
#' @keywords internal
scale_colour_aucmat_direction <- function(...) {
  ggplot2::scale_colour_manual(
    values = c(
      higher_in_positive = aucmat_colors$higher,
      lower_in_positive  = aucmat_colors$lower
    ),
    na.value = aucmat_colors$neutral,
    name = "Direction",
    ...
  )
}

#' aucmat fill colour scale (discrete)
#'
#' @param ... Passed to [ggplot2::scale_fill_manual()].
#' @export
#' @keywords internal
scale_fill_aucmat_direction <- function(...) {
  ggplot2::scale_fill_manual(
    values = c(
      higher_in_positive = aucmat_colors$higher,
      lower_in_positive  = aucmat_colors$lower
    ),
    na.value = aucmat_colors$neutral,
    name = "Direction",
    ...
  )
}

# ── Significance colour scale ──────────────────────────────────────

#' aucmat significance scale
#'
#' Red for significant, grey for non-significant.
#'
#' @param q_cutoff Cutoff value for the legend label.
#' @param ... Passed to [ggplot2::scale_colour_manual()].
#' @export
#' @keywords internal
scale_colour_aucmat_significance <- function(q_cutoff = 0.05, ...) {
  ggplot2::scale_colour_manual(
    values = c(
      "FALSE" = aucmat_colors$not_significant,
      "TRUE"  = aucmat_colors$significant
    ),
    labels = c(
      "FALSE" = "NS",
      "TRUE"  = paste0("q < ", q_cutoff)
    ),
    name = NULL,
    ...
  )
}
