calculate_deuteration = function(data, compare) {
  if (compare == "Deuteration") {
    columns = colnames(data)[!(colnames(data) %in% c("Mass", "Exposure"))]
    data = split(data, data[, columns])
    data = lapply(data, function(x) cbind(x[ !(colnames(x) == "Mass")], Mass = x$Mass - x$Mass[x$Exposure == 0]))
    data = do.call("rbind", c(data, make.row.names = FALSE))
  }
  data
}

student = function(data, time_point, compare = "Deuteration", alpha = 0.05, ...) {
  if (compare != "Mass" & compare != "Deuteration") stop("Either Mass or Deuteration must be chosen.")
  States = unique(data$State)
  if (length(States) < 2) stop("More than one state must be chosen.")

  data = data[data$Exposure == time_point, ]
  data = calculate_deuteration(data, compare)
  tests = t(combn(States, 2))
  n_tests = length(tests[, 1])

  test_results = lapply(1:n_tests, function(n) {
    state1 = tests[n, 1]
    state2 = tests[n, 2]

    state1_data = data[["Mass"]][data$State == state1]
    state2_data = data[["Mass"]][data$State == state2]

    if (length(state1_data) > 1 & length(state2_data) > 1) {
      result = t.test(x = state1_data, y = state2_data, conf.level = 1 - alpha, ...)

      p_value = result$p.value
      estimate = result$estimate
      if (p_value <= alpha) {
        conclusion = "difference"
      } else {
        conclusion = "no difference"
      }
    } else {
      p_value = NA
      conclusion = NA
      estimate = c(NA, NA)
      warning("Sample size must be greater than 1.")
    }

    data.frame(
      test = "Welch Two Sample t-test",
      State_1 = state1,
      State_2 = state2,
      Estimated = unname(diff(estimate)),
      P_value = p_value,
      Significance_level = alpha,
      Decision = conclusion
    )
  })
  do.call("rbind", c(test_results, make.row.names = FALSE))
}
