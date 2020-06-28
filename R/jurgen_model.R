library(lme4)
library(lmerTest)

jurgen_model = function(data, alpha = 0.05) {
  
  Time = c("continuous", "categorical", "continuous")
  Transformation = c("identity", "identity", "log")
  
  aic = rep(NA, 3)
  loglik = rep(NA, 3)
  p_value = rep(NA, 3)
  
  # continuous, identity
  model = lmer(Mass ~ State + Exposure*State + Exposure + (Rep| Exposure) + (1|State), data = data)
  result = anova(model)
  aic[1] = AIC(model)
  loglik[1] = logLik(model)
  p_value[1] = result$`Pr(>F)`[1]
  
  # categorical, identity
  model = lmer(Mass ~ State + factor(Exposure)*State + factor(Exposure) + (Rep| Exposure) + (1|State), data = data)
  result = anova(model)
  aic[2] = AIC(model)
  loglik[2] = logLik(model)
  p_value[2] = result$`Pr(>F)`[1]
  
  # continuous, log
  model = lmer(log(Mass+1) ~ State + Exposure*State + Exposure + (Rep| Exposure) + (1|State), data = data)
  result = anova(model)
  aic[3] = AIC(model)
  loglik[3] = logLik(model)
  p_value[3] = result$`Pr(>F)`[1]

  
  if (p_value <= alpha) {
    conclusion = "different"
  } else {
    conclusion = "equal"
  }
  
  data.frame(Test = "Jurgen Cl", #TODO: name model
             Covariate = compare,
             p_value = p_value,
             conclusion = conclusion)
}