library(lme4)
library(lmerTest)

jurgen_model = function(data, compare = "Deuteration", alpha = 0.05) {
  if (compare != "Mass" & compare != "Deuteration") stop("Either Mass or Deuteration must be chosen.")
  
  data = calculate_deuteration(data, compare)
  
  model = lmer(Mass ~ State + Exposure*State + Exposure + (Rep| Exposure) + (1|State), data = data)
  result = anova(model)
  p_value = result$`F value`[1]
  
  if (p_value <= alpha) {
    conclusion = "different"
  } else {
    conclusion = "equal"
  }
  
  data.frame(Test = "model Jurgena", #TODO: name model
             Covariate = compare,
             p_value = p_value,
             conclusion = conclusion)
}