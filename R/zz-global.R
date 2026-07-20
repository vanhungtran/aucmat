## Register commonly used global variables for ggplot/data.frame aesthetics
## This avoids R CMD check NOTE: no visible binding for global variable
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("fpr", "lower", "upper", "name"))
}
