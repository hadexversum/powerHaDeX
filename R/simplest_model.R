#' Simplest linear models for deuteration curves containing time, state of the protein and the interaction term,.
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @importFrom data.table data.table
#' @export
deuteros = function(data, significance_level = 0.05) {

  States = unique(data$State)

  Time = c("continuous", "categorical", "continuous")
  Transformation = c("identity", "identity", "log")

  aic = rep(NA, 3)
  loglik = rep(NA, 3)
  F_statistic = rep(NA, 3)
  p_value = rep(NA, 3)

  # continuous, identity
  model = lm(Mass ~ Exposure*State, data = data)
  result = anova(model)
  aic[1] = AIC(model)
  loglik[1] = logLik(model)
  F_statistic[1] = result$`F value`[2]
  p_value[1] = result$`Pr(>F)`[2]

  # categorical, identity
  model = lm(Mass ~ factor(Exposure)*State, data = data)
  result = anova(model)
  aic[2] = AIC(model)
  loglik[2] = logLik(model)
  F_statistic[2] = result$`F value`[2]
  p_value[2] = result$`Pr(>F)`[2]

  # continuous, log
  model = lm(Mass ~ log(Exposure + 1)*State, data = data)
  result = anova(model)
  aic[3] = AIC(model)
  loglik[3] = logLik(model)
  F_statistic[3] = result$`F value`[2]
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

