hdx_anova = function(data, time_point, compare = "Deuteration", alpha = 0.05) {
  if (compare != "Mass" & compare != "Deuteration") stop("Either Mass or Deuteration must be chosen.")

  States = unique(data$State)
  if (length(States) < 2) stop("More than one state must be chosen.")

  data = calculate_deuteration(data, compare)
  data = data[data$Exposure == time_point, ]

  result = aov(data[["Mass"]] ~ data[["State"]])
  p_value = summary.aov(result)[[1]][5][1, 1]
  if (p_value <= alpha) {
    conclusion = "different"
  } else {
    conclusion = "equal"
  }

  data.frame(Test = "One-way ANOVA",
             Covariate = compare,
             p_value = p_value,
             conclusion = conclusion)
}
