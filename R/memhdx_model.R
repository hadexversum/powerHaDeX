library(lme4)
library(lmerTest)

memhdx_model = function(data, significance_level = 0.05) {
  
  States = unique(data$State)
  
  Time = c("continuous", "categorical", "continuous")
  Transformation = c("identity", "identity", "log")
  
  aic = rep(NA, 3)
  loglik = rep(NA, 3)
  Test_statistic = rep(NA, 3)
  p_value = rep(NA, 3)
  
  # continuous, identity
  model = lmer(Mass ~ State + Exposure*State + Exposure + (1|Rep), data = data)
  result = anova(model)
  aic[1] = AIC(model)
  loglik[1] = logLik(model)
  Test_statistic[1] = result$`F value`[1]
  p_value[1] = result$`Pr(>F)`[1]
  
  
  # categorical, identity
  model = lmer(Mass ~ State + factor(Exposure)*State + factor(Exposure) + (1|Rep), data = data)
  result = anova(model)
  aic[2] = AIC(model)
  loglik[2] = logLik(model)
  Test_statistic[2] = result$`F value`[1]
  p_value[2] = result$`Pr(>F)`[1]
  
  # continuous, log
  model = lmer(Mass ~ State + log(Exposure+1)*State + log(Exposure+1) + (1|Rep), data = data)
  result = anova(model)
  aic[3] = AIC(model)
  loglik[3] = logLik(model)
  Test_statistic[3] = result$`F value`[1]
  p_value[3] = result$`Pr(>F)`[1]
  
  
  data.table(Test = "MEMHDX lmm",
             State_1 = States[1],
             State_2 = States[2],
             F_statistic = F_statistic,
             P_value = p_value,
             Significant_difference = (p_value <= significance_level),
             Time = Time,
             Transformation = Transformation,
             AIC = aic,
             logLik = loglik)
}