# ==============================================================================
# Multiplicity Adjustment Layer
#
# Applies standard p-value adjustments (BH, Holm, Bonferroni) to the complete
# set of feature p-values while preserving NA entries for features where
# inference was not possible.
# ==============================================================================

#' Adjust p-values for multiplicity across a screening result
#'
#' Calls `stats::p.adjust()` on the valid (non-NA) subset of p-values and
#' writes the adjusted values back in their original positions.
#'
#' @param p_values Numeric vector of raw p-values, possibly containing NAs.
#' @param method One of `"BH"`, `"holm"`, `"bonferroni"`, or `"none"`.
#'
#' @return Numeric vector of adjusted p-values (same length as `p_values`).
#' @noRd
#' @keywords internal
.adjust_pvalues <- function(p_values, method = c("BH", "holm", "bonferroni", "none")) {
  method <- match.arg(method)
  if (method == "none") return(p_values)

  valid <- !is.na(p_values)
  if (!any(valid)) return(p_values)

  q <- rep(NA_real_, length(p_values))
  q[valid] <- stats::p.adjust(p_values[valid], method = method)
  q
}
