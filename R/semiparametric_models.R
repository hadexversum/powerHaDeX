#' Truncated lines
#' @param x vector
#' @param kappa a number. Knot
#' @importFrom data.table data.table
#'
truncated_lines <- function(x, kappa){
    (x - kappa)*(x > kappa)
}


#' Semiparametric for deuteration curves
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @importFrom data.table data.table
#' @export

S1 <- function(data, significance_level = 0.05) {

    States = unique(data$State)
    Transformation = c("identity", "log")
    aic = loglik = Test_statistic = p_value = rep(NA, 2)

    knots <- unique(setdiff(data$Exposure, c(max(data$Exposure), min(data$Exposure))))

    X <- sapply(1:length(knots), function(knot) {
        truncated_lines(data$Exposure, knots[knot])
    })

    # identity
    model <- lm(Mass ~ Exposure + State + X,
                data = data)
    model_reduced <- lm(Mass ~ Exposure + X,
                        data = data)
    result = anova(model, model_reduced)
    aic[1] = AIC(model)
    loglik[1] = logLik(model)
    Test_statistic[1] = result$`F`[2]
    p_value[1] = result$`Pr(>F)`[2]

    # log
    model = lm(Mass ~ log(Exposure + 1) + State + log(X + 1),
               data = data)
    model_reduced = lm(Mass ~ log(Exposure + 1) + log(X + 1),
                       data = data)
    result = anova(model, model_reduced)
    aic[2] = AIC(model)
    loglik[2] = logLik(model)
    Test_statistic[2] = result$`F`[2]
    p_value[2] = result$`Pr(>F)`[2]

    data.table::data.table(Test = "S1",
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = Test_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = NA,
                           Transformation = Transformation,
                           AIC = aic,
                           logLik = loglik)
}

#' Semiparametric for deuteration curves
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @importFrom data.table data.table
#' @export

S2 <- function(data, significance_level = 0.05) {

    States = unique(data$State)
    Transformation = c("identity", "log")
    aic = loglik = Test_statistic = p_value = rep(NA, 2)

    knots <- unique(setdiff(data$Exposure, c(max(data$Exposure), min(data$Exposure))))

    X <- sapply(1:length(knots), function(knot) {
        truncated_lines(x, knots[knot])
    })
    # identity
    model <- lm(Mass ~ Exposure*State + X,
                data = data)
    model_reduced <- lm(Mass ~ Exposure + X,
                        data = data)
    result = anova(model, model_reduced)
    aic[1] = AIC(model)
    loglik[1] = logLik(model)
    Test_statistic[1] = result$`F`[2]
    p_value[1] = result$`Pr(>F)`[2]

    # log
    model = lm(Mass ~ log(Exposure + 1)*State + log(X + 1),
               data = data)
    model_reduced = lm(Mass ~ log(Exposure + 1) + log(X + 1),
                       data = data)
    result = anova(model, model_reduced)
    aic[2] = AIC(model)
    loglik[2] = logLik(model)
    Test_statistic[2] = result$`F`[2]
    p_value[2] = result$`Pr(>F)`[2]

    data.table::data.table(Test = "S2",
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = Test_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = NA,
                           Transformation = Transformation,
                           AIC = aic,
                           logLik = loglik)
}


#' Semiparametric for deuteration curves
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @importFrom data.table data.table
#' @export


S3 <- function(data, significance_level = 0.05) {

    States = unique(data$State)
    Transformation = c("identity", "log")
    aic = loglik = Test_statistic = p_value = rep(NA, 2)

    # identity
    model <- gam(Mass~s(Exposure) + Exposure*State,
                 data = data,
                 method="REML")
    model_reduced <- gam(Mass ~ s(Exposure) + Exposure,
                         data=data,
                         method="REML")

    result = anova.gam(model, model_reduced, test = "Chisq")
    aic[1] = AIC(model)
    Test_statistic[1] = result$Deviance[2]
    p_value[1] = result$`Pr(>Chi)`[2]
    loglik[1] = as.numeric(logLik(model))

    # log
    model <- gam(Mass~s(log(Exposure + 1)) + log(Exposure + 1)*State,
                 data = data,
                 method="REML")
    model_reduced <- gam(Mass ~ s(log(Exposure+ 1)) + log(Exposure + 1),
                         data=data,
                         method="REML")

    result = anova.gam(model, model_reduced, test = "Chisq")
    aic[2] = AIC(model)
    Test_statistic[2] = result$Deviance[2]
    p_value[2] = result$`Pr(>Chi)`[2]
    loglik[2] = as.numeric(logLik(model))

    data.table::data.table(Test = "S3",
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = Test_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = NA,
                           Transformation = Transformation,
                           AIC = aic,
                           logLik = loglik)
}


#' Semiparametric for deuteration curves
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @importFrom data.table data.table
#' @export


S4 <- function(data, significance_level = 0.05) {

    States = unique(data$State)
    Transformation = c("identity", "log")
    aic = loglik = Test_statistic = p_value = rep(NA, 2)

    #identity
    model_reduced <- gam(Mass ~ s(Exposure),
                         data=data,
                         method="REML")
    model <- gam(Mass~s(Exposure) + State,
                 data = data,
                 method="REML")

    result = anova.gam(model, model_reduced, test = "Chisq")
    aic[1] = AIC(model)
    Test_statistic[1] = result$Deviance[2]
    p_value[1] = result$`Pr(>Chi)`[2]
    loglik[1] = as.numeric(logLik(model))

    #log
    model_reduced <- gam(Mass ~ s(log(Exposure + 1)),
                         data=data,
                         method="REML")
    model <- gam(Mass~s(log(Exposure + 1)) + State,
                 data = data,
                 method="REML")

    result = anova.gam(model, model_reduced, test = "Chisq")
    aic[2] = AIC(model)
    Test_statistic[2] = result$Deviance[2]
    p_value[2] = result$`Pr(>Chi)`[2]
    loglik[2] = as.numeric(logLik(model))

    data.table::data.table(Test = "S4",
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = Test_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = NA,
                           Transformation = Transformation,
                           AIC = aic,
                           logLik = loglik)
}


#' Semiparametric for deuteration curves
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @importFrom data.table data.table
#' @export


S5 <- function(data, significance_level = 0.05) {

    States = unique(data$State)
    Transformation = c("identity", "log")
    aic = loglik = Test_statistic = p_value = rep(NA, 2)

    #identity
    model_reduced = gam(Mass ~ s(Exposure) + Exposure+ s(id, bs = "re"),
                        data=data,
                        method="REML")
    model = gam(Mass ~ s(Exposure) + Exposure*State + s(id, bs = "re"),
                data=data,
                method="REML")

    result = anova.gam(model, model_reduced, test = "Chisq")
    aic[1] = AIC(model)
    Test_statistic[1] = result$Deviance[2]
    p_value[1] = result$`Pr(>Chi)`[2]
    loglik[1] = as.numeric(logLik(model))

    #log
    model_reduced = gam(Mass ~ s(log(Exposure + 1)) + log(Exposure + 1)+ s(id, bs = "re"),
                        data=data,
                        method="REML")
    model = gam(Mass ~ s(log(Exposure + 1)) + log(Exposure + 1)*State + s(id, bs = "re"),
                data=data,
                method="REML")

    result = anova.gam(model, model_reduced, test = "Chisq")
    aic[2] = AIC(model)
    Test_statistic[2] = result$Deviance[2]
    p_value[2] = result$`Pr(>Chi)`[2]
    loglik[2] = as.numeric(logLik(model))

    data.table::data.table(Test = "S5",
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = Test_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = NA,
                           Transformation = Transformation,
                           AIC = aic,
                           logLik = loglik)
}





#' Semiparametric for deuteration curves
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @importFrom data.table data.table
#' @export


S6 <- function(data, significance_level = 0.05) {

    States = unique(data$State)
    Transformation = c("identity", "log")
    aic = loglik = Test_statistic = p_value = rep(NA, 2)

    #identity
    model_reduced = gam(Mass~s(Exposure) + Exposure*State +
                            s(id, bs="re") + s(id, Exposure, bs="re"),
                        data=example_data,
                        method="REML")
    model = gam(Mass~s(Exposure) + Exposure*State +
                    s(id, bs="re") + s(id, Exposure, bs="re"),
                data=example_data,
                method="REML")

    result = anova.gam(model, model_reduced, test = "Chisq")
    aic[1] = AIC(model)
    Test_statistic[1] = result$Deviance[2]
    p_value[1] = result$`Pr(>Chi)`[2]
    loglik[1] = as.numeric(logLik(model))

    #log
    model_reduced = gam(Mass~s(log(Exposure + 1)) + log(Exposure + 1) +
                            s(id, bs="re") + s(id, Exposure, bs="re"),
                        data=example_data,
                        method="REML")

    model = gam(Mass~s(log(Exposure + 1)) + log(Exposure + 1)*State +
                    s(id, bs="re") + s(id, Exposure, bs="re"),
                data=example_data,
                method="REML")

    result = anova.gam(model, model_reduced, test = "Chisq")
    aic[2] = AIC(model)
    Test_statistic[2] = result$Deviance[2]
    p_value[2] = result$`Pr(>Chi)`[2]
    loglik[2] = as.numeric(logLik(model))

    data.table::data.table(Test = "S6",
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = Test_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = NA,
                           Transformation = Transformation,
                           AIC = aic,
                           logLik = loglik)
}


