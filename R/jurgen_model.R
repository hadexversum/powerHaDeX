library(lme4)
library(lmerTest)

jurgen_model = function(data, alpha = 0.05) {
  
  States = unique(data$State)
  
  Time = c("continuous", "categorical", "continuous")
  Transformation = c("identity", "identity", "log")
  
  aic = rep(NA, 3)
  loglik = rep(NA, 3)
  Test_statistic = rep(NA, 3)
  p_value = rep(NA, 3)
  
  # continuous, identity
  model = lmer(Mass ~ State * Exposure , 
               data = data,
               REML = FALSE)
  model_reduced = lmer(Mass ~ Exposure , 
                       data = data,
                       REML = FALSE)
  result = anova(model, model_reduced)
  aic[1] = AIC(model)
  Test_statistic[1] = result$Chisq[2]
  p_value[1] = result$`Pr(>Chisq)`[2]
  
  # categorical, identity
  model = lmer(Mass ~ State * factor(Exposure) , 
               data = data,
               REML = FALSE)
  model_reduced = lmer(Mass ~ factor(Exposure) , 
                       data = data,
                       REML = FALSE)
  result = anova(model, model_reduced)
  aic[2] = AIC(model)
  loglik[2] = logLik(model)
  Test_statistic[2] = result$Chisq[2]
  p_value[2] = result$`Pr(>Chisq)`[2]
  
  # continuous, log
  model = lmer(Mass ~ State * log(Exposure+1) , 
               data = data,
               REML = FALSE)
  model_reduced = lmer(Mass ~ log(Exposure+1) , 
                       data = data,
                       REML = FALSE)
  result = anova(model, model_reduced)
  aic[3] = AIC(model)
  loglik[3] = logLik(model)
  Test_statistic[3] = result$Chisq[2]
  p_value[3] = result$`Pr(>Chisq)`[2]
  
  
  data.frame(Test = "Jurgen Claesen lmm", 
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