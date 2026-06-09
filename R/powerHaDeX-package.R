#' powerHaDeX: Power analysis for HDX-MS experiments
#'
#' The R-package powerHaDeX provides tools for simulating and analyzing data
#' coming from HDX-MS experiments, along with the possibility of comparing the
#' power of tests verifying differences in deuteration levels.
#'
#' @import data.table
#' @importFrom stats AIC aggregate anova dbinom lm logLik pt qt rnorm sd weighted.mean coefficients
#' @importFrom nlme corCompSymm corAR1
#' @importFrom expm %^%
#' @importFrom Rcpp sourceCpp
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom ggplot2 ggplot aes geom_segment ylab xlab theme_minimal theme element_text ggtitle facet_wrap
#' @importFrom glmnet glmnet
#' @importFrom checkmate assert checkLogical checkChoice checkFALSE
#' @author Michal Burdukiewicz, Krystyna Grzesiak, Mateusz Staniak
#' @keywords internal
"_PACKAGE"

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "Test", "State_1", "State_2", "State", "Sequence",
        "Significant_difference", "Center", "z", "PF",
        "exp_mass", "Inten", "Experimental_state", "File",
        "Exposure", "avg_exp_mass", "Rep", "Mass", "Mz",
        "err_avg_mass", "err_deut_uptake", "deut_uptake",
        "Transformation", "Num_replicates", "Num_states",
        "Num_timepoints", "Time", "Intensity", "Charge"
    ))
}
