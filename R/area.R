#' Test based on area under the deuteration curve
#' @param data data.table of deuteration curves
#' @param significance_level significance level that will be used for testing
#' @importFrom data.table data.table
#' @export
auc_test = function(data, significance_level = 0.05) {
  States = unique(data$State)
  if (length(States) < 2) stop("More than one state must be chosen.")

  state1_data = data[data$State == States[1], ]
  state2_data = data[data$State == States[2], ]

  states_exposure = setdiff(intersect(state1_data$Exposure, state2_data$Exposure), 0)

  state1_data = state1_data[state1_data$Exposure %in% states_exposure, ]
  state2_data = state2_data[state2_data$Exposure %in% states_exposure, ]

  t_n = max(states_exposure)
  t_1 = min(states_exposure)
  times = length(data$Exposure)
  state1_masses = (max(state1_data$Mass) - state1_data$Mass) / max(state1_data$Mass)
  state2_masses = (max(state2_data$Mass) - state2_data$Mass) / max(state2_data$Mass)

  y_a = aggregate(state1_masses,
                  list(state1_data$Exposure), mean)$x
  y_b = aggregate(state2_masses,
                  list(state2_data$Exposure), mean)$x
  s_a = aggregate(state1_masses,
                  list(state1_data$Exposure), sd)$x
  s_b = aggregate(state2_masses,
                  list(state2_data$Exposure), sd)$x
  S_a = sqrt(log(t_n / t_1) * mean(s_a^2))
  S_b = sqrt(log(t_n / t_1) * mean(s_b^2))
  n_rep_a = length(unique(state1_data$Rep))
  n_rep_b = length(unique(state2_data$Rep))
  S = sqrt(((n_rep_a - 1)*S_a^2 + (n_rep_b - 1)*S_b^2 ) / (n_rep_a + n_rep_b - 2))
  A_aver = log(t_n/t_1) * (1/times) * sum(y_a - y_b)

  Test_statistic = A_aver / (S * sqrt(1/n_rep_a + 1/n_rep_b))

  if (length(state1_data) > 1 & length(state2_data) > 1) {
    P_value = pt(abs(Test_statistic), df = (n_rep_a + n_rep_b - 2))
  } else {
    P_value = NA
    warning("Sample size must be greater than 1.")
  }

  data.table::data.table(
    test = "AUC test",
    State_1 = as.character(States[1]),
    State_2 = as.character(States[2]),
    Test_statistic = Test_statistic,
    P_value = P_value,
    Time = "continuous",
    Transformation = "log",
    Significant_difference = P_value < significance_level,
    AIC = NA,
    logLik = NA
  )

}
