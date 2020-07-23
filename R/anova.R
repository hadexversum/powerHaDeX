#' ANOVA based on gls 
#' @param data data.table of deuteration curves at given time point
#' @param significance_level significance level for testing
#' @import nlme
#' @export

hdx_anova = function(data, significance_level = 0.05) {
  
  structure = c("compound symmetry", "AR(1)")
  
  aic = rep(NA, 2)
  loglik = rep(NA, 2)
  Test_statistic = rep(NA, 2)
  p_value = rep(NA, 2)
  
  # compound symmetry structure, categorical time
  model = nlme::gls(Mass ~ factor(Exposure) + State + factor(Exposure) : State,
                    data = data, method = "ML",
                    correlation = corCompSymm(form = ~1|Rep))
  
  model_reduced = nlme::gls(Mass ~ factor(Exposure) ,
                     data = data, method = "ML",
                     correlation = corCompSymm(form = ~1|Rep))
  
  result = anova(model, model_reduced)
  aic[1] = AIC(model)
  Test_statistic[1] = result$L.Ratio[2]
  p_value[1] = result$`p-value`[2]
  loglik[1] = logLik(model)
  
  # AR(1) structure, categorical
  model = nlme::gls(Mass ~ factor(Exposure) + State + factor(Exposure) : State,
                    data = data, method = "ML",
                    correlation = corAR1(form = ~1|Rep))
  
  model_reduced = nlme::gls(Mass ~ factor(Exposure) ,
                            data = data, method = "ML",
                            correlation = corAR1(form = ~1|Rep))
  
  result = anova(model, model_reduced)
  aic[2] = AIC(model)
  Test_statistic[2] = result$L.Ratio[2]
  p_value[2] = result$`p-value`[2]
  loglik[2] = logLik(model)

  
  data.frame(Test = "ANOVA",
             Structure = structure,
             Test_statistic = Test_statistic,
             P_value = p_value,
             Significant_difference = (p_value <= significance_level),
             AIC = aic,
             logLik = loglik)
}