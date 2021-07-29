#' Houde's test for deuteration curves
#' @description This function performs Damian Houde's confidence intervals test
#' for differences in deuteration levels. Its input and output are compatible with
#' the function \code{\link[powerHaDeX]{calculate_hdx_power}}.
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @returns This function returns a data table compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#' @references Houde, Damian, Steven A Berkowitz, and John R Engen (2011).
#' “The utility of hydrogen/deuterium exchange mass spectrometry in biopharmaceutical
#' comparabilitystudies”. In:Journal of pharmaceutical sciences100.6, pp. 2071–2086.
#' @seealso \code{\link[powerHaDeX]{calculate_hdx_power}} for estimation of power
#' of tests for differences in deuteration levels.
#' @importFrom data.table data.table
#' @export

houde <- function(data, significance_level = 0.05) {

    Sequence = State = Exposure = Rep = Experimental_state = Mass = err_avg_mass = avg_exp_mass = err_deut_uptake = deut_uptake = NULL

    States = unique(data$State)
    confidence_limit = 1 - significance_level

    alpha <- 1 - confidence_limit
    t_value <- qt(c(alpha/2, 1 - alpha/2), df = 2)[2]


    data_exp <- data[, .(avg_exp_mass = mean(Mass)),
                     by = list(Sequence, State, Exposure, Rep, Experimental_state)]

    data_exp <- data_exp[, .(deut_uptake = mean(avg_exp_mass), err_avg_mass = sd(avg_exp_mass)/sqrt(length(Rep))),
                         by = list(State, Sequence, Exposure, Experimental_state)]

    data_exp[, err_deut_uptake := sqrt(err_avg_mass^2 + err_avg_mass[Exposure == 0]^2),
             by = list(State, Sequence, Experimental_state)]


    calc_dat <- data_exp[, .(diff_deut_uptake = deut_uptake[Experimental_state == "A"] -
                                 deut_uptake[Experimental_state == "B"],
                             err_diff_deut_uptake = sqrt(err_deut_uptake[Experimental_state == "A"]^2 +
                                                             err_deut_uptake[Experimental_state == "B"]^2)),
                         by = list(Sequence, Exposure)]



    avg_difference <- mean(calc_dat$diff_deut_uptake)

    x_threshold <- t_value * mean(calc_dat[["err_diff_deut_uptake"]], na.rm = TRUE)/sqrt(length(calc_dat))

    data.table(Test = "Houde",
               State_1 = States[1],
               State_2 = States[2],
               Test_statistic = NA,
               P_value = NA,
               Significant_difference = abs(avg_difference) > x_threshold,
               Time = NA,
               Transformation = NA,
               AIC = NA,
               logLik = NA)
}



#' HDX-Analyzer model
#' @description This function performs the test based on the simplest linear models for deuteration
#' curves containing time, state of the protein and the interaction term. Its input
#' and output are compatible with the function \code{\link[powerHaDeX]{calculate_hdx_power}}.
#' @inheritParams houde
#' @returns This function returns a data table compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#' @seealso Liu, Sanmin et al. (2011). “HDX-analyzer: a novel package for statistical
#' analysis of protein structure dynamics”. In:BMC bioinformatics12.1, pp. 1–10.
#' @importFrom data.table data.table
#' @export

hdx_analyzer = function(data, significance_level = 0.05) {

    States = unique(data$State)

    Time = c("continuous", "categorical", "continuous")
    Transformation = c("identity", "identity", "log")

    aic = rep(NA, 3)
    loglik = rep(NA, 3)
    F_statistic = rep(NA, 3)
    p_value = rep(NA, 3)

    # continuous, identity
    model = lm(Mass ~ Exposure*State, data = data)
    model_reduced = lm(Mass ~ Exposure, data = data)
    result = anova(model, model_reduced)
    aic[1] = AIC(model)
    loglik[1] = logLik(model)
    F_statistic[1] = result$`F`[2]
    p_value[1] = result$`Pr(>F)`[2]

    # categorical, identity
    model = lm(Mass ~ factor(Exposure)*State, data = data)
    model_reduced = lm(Mass ~ factor(Exposure), data = data)
    result = anova(model, model_reduced)
    aic[2] = AIC(model)
    loglik[2] = logLik(model)
    F_statistic[2] = result$`F`[2]
    p_value[2] = result$`Pr(>F)`[2]

    # continuous, log
    model = lm(Mass ~ log(Exposure + 1)*State, data = data)
    model_reduced = lm(Mass ~ log(Exposure + 1), data = data)
    result = anova(model, model_reduced)
    aic[3] = AIC(model)
    loglik[3] = logLik(model)
    F_statistic[3] = result$`F`[2]
    p_value[3] = result$`Pr(>F)`[2]

    data.table::data.table(Test = "Deuteros lm",
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = F_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = Time,
                           Transformation = Transformation,
                           AIC = aic,
                           logLik = loglik)
}

#' MEMHDX model
#' @description This function performs the test based on a linear mixed effects
#' model used in MEMHDX tools. Its input and output are compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#' @inheritParams houde
#' @returns This function returns a data table compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#' @seealso Hourdel, Véronique et al. (July 2016). “MEMHDX: an interactive tool
#' to expedite thestatistical validation and visualization of large HDX-MS data sets”.
#' In:Bioinformatics32.22, pp. 3413–3419.issn: 1367-4803.
#' @import lmerTest
#' @export

memhdx_model = function(data, significance_level = 0.05) {

    States = unique(data$State)

    Time = c("continuous", "categorical", "continuous")
    Transformation = c("identity", "identity", "log")

    aic = rep(NA, 3)
    loglik = rep(NA, 3)
    Test_statistic = rep(NA, 3)
    p_value = rep(NA, 3)

    # continuous, identity
    model = lmerTest::lmer(Mass ~ Exposure*State + (1|Rep),
                           data = data,
                           REML = FALSE)
    model_reduced = lmerTest::lmer(Mass ~ Exposure + (1|Rep),
                                   data = data,
                                   REML = FALSE)
    result = anova(model, model_reduced)
    aic[1] = AIC(model)
    loglik[1] = logLik(model)
    Test_statistic[1] = result$Chisq[2]
    p_value[1] = result$`Pr(>Chisq)`[2]


    # categorical, identity
    model = lmerTest::lmer(Mass ~ factor(Exposure)*State + (1|Rep),
                           data = data,
                           REML = FALSE)
    model_reduced = lmerTest::lmer(Mass ~ factor(Exposure) + (1|Rep),
                                   data = data,
                                   REML = FALSE)
    result = anova(model, model_reduced)
    aic[2] = AIC(model)
    loglik[2] = logLik(model)
    Test_statistic[2] = result$Chisq[2]
    p_value[2] = result$`Pr(>Chisq)`[2]

    # continuous, log
    model = lmerTest::lmer(Mass ~ log(Exposure + 1)*State + (1|Rep),
                           data = data,
                           REML = FALSE)
    model_reduced = lmerTest::lmer(Mass ~ log(Exposure+1) + (1|Rep),
                                   data = data,
                                   REML = FALSE)
    result = anova(model, model_reduced)
    aic[3] = AIC(model)
    loglik[3] = logLik(model)
    Test_statistic[3] = result$Chisq[2]
    p_value[3] = result$`Pr(>Chisq)`[2]

    data.table(Test = "MEMHDX lmm",
               State_1 = States[1],
               State_2 = States[2],
               Test_statistic = Test_statistic,
               P_value = p_value,
               Significant_difference = (p_value <= significance_level),
               Time = Time,
               Transformation = Transformation,
               AIC = aic,
               logLik = loglik)
}

#' Get truncated lines
#' @description This function prepares the function space for semiparametric
#' model using truncated lines.
#' @param x a vector. Variable which undergoes the truncating
#' @param kappa a vector of knots in which the spline model will break
#' @keywords internal

truncated_lines <- function(x, kappa){
    (x - kappa)*(x > kappa)
}

#' Semiparametric test for differences in deuteration levels
#' @description This function performs the semiparametric test for differences
#' in deuteration levels. Its input and output are compatible with
#' the function \code{\link[powerHaDeX]{calculate_hdx_power}}.
#' @inheritParams houde
#' @details This function uses \code{\link[powerHaDeX]{truncated_lines}}. The knots
#' considered in the testing procedure are chosen using ridge regression.
#' @returns This function returns a data table compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
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






