library(dplyr)
library(gtools)


data = readRDS("hdx_processed.RDS", refhook = NULL)


test_data = data %>% 
  group_by(Protein, Sequence, State, Exposure) %>% 
  mutate(Mass = mean(na.omit(Mass))) %>% 
  select(-Rep) %>% 
  unique() %>% 
  group_by(Protein, Sequence, State) %>% 
  mutate(Deuteration = Mass - Mass[Exposure == 0])


test_data_time_1 = data %>% 
  group_by(Protein, Sequence, State, Exposure) %>% 
  mutate(Mass = mean(na.omit(Mass))) %>% 
  select(-Rep) %>% 
  unique() %>% 
  group_by(Protein, Sequence, State) %>% 
  mutate(Deuteration = Mass - Mass[Exposure == 0]) %>% 
  filter(Exposure == 1)



student = function(data, compare = "Mass", alpha) {
  if(compare != "Mass" & compare != "Deuteration") stop("Either Mass or Deuteration must be chosen.")
  
  States = unique(data$State)
  if(length(States) < 2) stop("More than one state must be chosen.")
  
  tests = combinations(length(States), 2, States)
  n_tests = length(tests[,1])
  
  test_results = sapply(1:n_tests, function(n) {
    state1 = tests[n, 1]
    state2 = tests[n, 2]
    
    state1_data = data[[compare]][data$State == state1]
    state2_data = data[[compare]][data$State == state2]
    
    if(length(state1_data) > 1 & length(state2_data) > 1) {
      result = t.test(x = state1_data, y = state2_data, conf.level = 1-alpha)
      
      p_value = result$p.value
      estimate = result$estimate
      if(p_value <= alpha) {conclusion = "different"}else {conclusion = "equal"}
      
    }else {
      p_value = NA
      conclusion = NA
      warning("Sample size must be greater than 1.")
    }
    
    c("Welch Two Sample t-test", state1, state2, estimate, 
      significance_level = alpha, p_value, conclusion)
  })
  output = data.frame(t(test_results))
  colnames(output) = c("Test", "State_1", "State_2", "Estimate_State_1", 
                       "Estimate_State_2", "Significance_level", "p_value", 
                       "Conclusion")
  output
}


student(test_data_time_1)

