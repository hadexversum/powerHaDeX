simplest_model = function(data, compare = "Deuteration", alpha = 0.05) {
  if (compare != "Mass" & compare != "Deuteration") stop("Either Mass or Deuteration must be chosen.")
  
  data = calculate_deuteration(data, compare)
  
  model = lm(Mass ~ Exposure*State, data = data)
  result = anova(model)
  
  p_value = result$`Pr(>F)`[2]
  
  if (p_value <= alpha) {
    conclusion = "different"
  } else {
    conclusion = "equal"
  }
  
  data.frame(Test = "model 1", #TODO: name model
             Covariate = compare,
             p_value = p_value,
             conclusion = conclusion)
}