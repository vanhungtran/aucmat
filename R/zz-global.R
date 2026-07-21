## Register commonly used global variables for ggplot2/data.frame aesthetics
## This avoids R CMD check NOTE: no visible binding for global variable
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # plot_auc_rank / plot_auc_volcano / plot_auc_forest
    "rank", "auc_strength", "auc_raw", "effect_direction",
    "conf_low", "conf_high", "label", "effect", "neg_log10_q",
    "sig", "biomarker",
    # plot_auc_stability
    "rank_median", "rank_q25", "rank_q75",
    # plot_roc_top (inherited from pROC::ggroc)
    "fpr", "lower", "upper", "name",
    # deprecated plot_roc_with_combos
    "colour", "fill"
  ))
}
