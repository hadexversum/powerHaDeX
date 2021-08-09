#' powerHaDeX
#'
#' @description The R-package powerHaDeX for simulating and analyzing data coming
#' from HDX-MS experiments along with the possibility of comparing the power of
#' the tests verifying differences in deuteration levels.
#'
#' @importFrom stats AIC aggregate anova dbinom lm
#' logLik pt rnorm sd weighted.mean coefficients
#' @importFrom nlme corCompSymm corAR1
#' @importFrom plyr .
#' @importFrom expm %^%
#' @importFrom Rcpp sourceCpp
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom ggplot2 ggplot aes geom_segment ylab xlab theme_minimal theme element_text ggtitle facet_wrap
#' @importFrom stats qt
#' @importFrom glmnet glmnet
#' @author Michal Burdukiewicz, Krystyna Grzesiak, Mateusz Staniak
#' @docType package
#' @name powerHaDeX
NULL
