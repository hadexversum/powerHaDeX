#' Simplest linear models for deuteration curves
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @importFrom data.table data.table
#' @export
deuteros = function(data, significance_level = 0.05) {

  Time = c("continuous", "categorical", "continuous")
  Transformation = c("identity", "identity", "log")

  aic = rep(NA, 3)
  loglik = rep(NA, 3)
  p_value = rep(NA, 3)

  # continuous, identity
  model = lm(Mass ~ Exposure*State, data = data)
  result = anova(model)
  aic[1] = AIC(model)
  loglik[1] = logLik(model)
  p_value[1] = result$`Pr(>F)`[2]

  # categorical, identity
  model = lm(Mass ~ factor(Exposure)*State, data = data)
  result = anova(model)
  aic[2] = AIC(model)
  loglik[2] = logLik(model)
  p_value[2] = result$`Pr(>F)`[2]

  # continuous, log
  model = lm(log(Mass + 1) ~ Exposure*State, data = data)
  result = anova(model)
  aic[3] = AIC(model)
  loglik[3] = logLik(model)
  p_value[3] = result$`Pr(>F)`[2]

  data.table::data.table(Test = "Deuteros", #TODO: name model
                         P_value = p_value,
                         Significant_difference = (p_value <= significance_level),
                         Time = Time,
                         Transformation = Transformation,
                         AIC = aic,
                         logLik = loglik)
}

