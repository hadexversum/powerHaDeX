

#' Get truncated lines
#' @description This function prepares the function space for semiparametric
#' model using truncated lines.
#' @param x a vector. Variable which undergoes the truncating
#' @param kappa a vector of knots in which the spline model will break
#' @keywords internal


truncated_lines <- function(x, kappa){
    (x - kappa)*(x > kappa)
}


#' Semiparametric test for differences in deuteration
#' @description This function performs the semiparametric test for differences
#' in deuteration levels. Its input and output are compatible with
#' the function \code{\link[powerHaDeX]{calculate_hdx_power}}.
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @details This function uses \code{\link[powerHaDeX]{truncated_lines}}. The knots
#' considered in the testing procedure are chosen using ridge regression.
#' @importFrom data.table data.table
#' @importFrom glmnet glmnet
#' @export

semiparametric <- function(data, significance_level = 0.05) {

    States = unique(data$State)
    data[["id"]] <- paste0(data$Rep, data$Charge, data$Experimental_state)
    Test = aic = loglik = Test_statistic = p_value = NA

    knots <- unique(setdiff(data$Exposure, c(max(data$Exposure), min(data$Exposure))))
    X <- sapply(knots, function(k) {
        truncated_lines(data$Exposure, k)
    })
    colnames(X) <- c(paste0("knot_", as.character(knots)))

    cv_fit <- glmnet(X, data[["Mass"]], alpha = 0, lambda = 0.001)
    coefs <- coefficients(cv_fit)
    X_reduced <- cbind(intercept = 1, X)[, which(abs(coefs) >= 2*10^(-5))]


    model = lmerTest::lmer(Mass ~ Exposure*State + (1|id) + (1|Exposure) + X_reduced,
                           data = data,
                           REML = FALSE)
    model_reduced = lmerTest::lmer(Mass ~ Exposure + (1|id) + (1|Exposure) + X_reduced,
                                   data = data,
                                   REML = FALSE)
    result = anova(model, model_reduced)
    aic = AIC(model)
    loglik = as.numeric(logLik(model))
    Test_statistic = result$Chisq[2]
    p_value = result$`Pr(>Chisq)`[2]

    Test = "RIDGE_knots_random_intercept_id_exposure"

    data.table::data.table(Test = Test,
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = Test_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = NA,
                           Transformation = "identity",
                           AIC = aic,
                           logLik = loglik)
}
