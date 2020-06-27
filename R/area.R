


simplest_model = function(data, compare = "Deuteration", alpha = 0.05) {
  if (compare != "Mass" & compare != "Deuteration") stop("Either Mass or Deuteration must be chosen.")
  States = unique(data$State)
  if (length(States) < 2) stop("More than one state must be chosen.")

  data = calculate_deuteration(data, compare)
  tests = t(combn(States, 2))
  n_tests = length(tests[, 1])
  
  test_results = lapply(1:n_tests, function(n) {
    state1 = tests[n, 1]
    state2 = tests[n, 2]
    state1_data = data[data$State == state1,]
    state2_data = data[data$State == state2,]
    
    states_exposure = setdiff(intersect(state1_data$Exposure, state2_data$Exposure), 0)

    state1_data = state1_data[state1_data$Exposure %in% states_exposure,]
    state2_data = state2_data[state2_data$Exposure %in% states_exposure,]
    
    
    t_n = tail(Exposure, n = 1)
    t_1 = Exposure[1]
    times = length(Exposure)
    
    y_a = aggregate(state1_data$Mass, list(state1_data$Exposure), mean)$x
    y_b = aggregate(state2_data$Mass, list(state2_data$Exposure), mean)$x
    
    A_aver = log(t_n/t_1) * (1/times) * sum(y_a - y_b)
    
    A_trapez_a = 0.5*(y_a[1]*log(states_exposure[2]/states_exposure[1]) + 
                      y_a[times]*log(states_exposure[times]/states_exposure[times - 1])+
                      sum(y_a[2:(times-1)]*log(states_exposure[3:times]/states_exposure[1:(times-2)])))
    
    A_trapez_b = 0.5*(y_b[1]*log(states_exposure[2]/states_exposure[1]) + 
                        y_b[times]*log(states_exposure[times]/states_exposure[times - 1])+
                        sum(y_b[2:(times-1)]*log(states_exposure[3:times]/states_exposure[1:(times-2)])))
    
    A_trapez = A_trapez_b - A_trapez_a
    
    
    if (length(state1_data) > 1 & length(state2_data) > 1) {
      p_value = #####?? nie mam pojecia 
        
      if (p_value <= alpha) {
        conclusion = "difference"
      } else {
        conclusion = "no difference"
      }
    } else {
      p_value = NA
      conclusion = NA
      warning("Sample size must be greater than 1.")
    }
    
    data.frame(
      test = "Test based on the Area Between Exchange Curves",
      State_1 = state1,
      State_2 = state2,
      P_value = p_value,
      Significance_level = alpha,
      Decision = conclusion
    )
  })
  do.call("rbind", c(test_results, make.row.names = FALSE))

}