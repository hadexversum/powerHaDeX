library(dplyr)
library(gtools)


data = readRDS("hdx_processed.RDS", refhook = NULL)


test_data_time_1 = data %>% 
  group_by(Protein, Sequence, State, Exposure) %>% 
  mutate(Mass = mean(na.omit(Mass))) %>% 
  select(-Rep) %>% 
  unique() %>% 
  group_by(Protein, Sequence, State) %>% 
  mutate(Deuteration = Mass - Mass[Exposure == 0]) %>% 
  filter(Exposure == 1)


hdx_anova = function(data, compare = "Mass" ,alpha = 0.05) {
  if(compare != "Mass" & compare != "Deuteration") stop("Either Mass or Deuteration must be chosen.")
  
  States = unique(data$State)
  if(length(States) < 2) stop("More than one state must be chosen.")
  
  result = aov(data[[compare]] ~ data[["State"]])
  p_value = summary.aov(result)[[1]][5][1,1]
  if(p_value <= alpha) {conclusion = "different"}else {conclusion = "equal"}
  
  data.frame(Test = "One-way ANOVA", 
             Covariate = compare, 
             p_value = p_value,
             conclusion = conclusion)
  
}

hdx_anova(test_data_time_1)

