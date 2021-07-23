#' Houde's test for deuteration curves
#' @description This function performs Damian Houde's confidence intervals test
#' for differences in deuteration levels. Its input and output are compatible with
#' the function \code{\link[powerHaDeX]{calculate_hdx_power}}.
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#' @importFrom data.table data.table
#' @export

houde <- function(data, significance_level = 0.05) {

    Sequence = State = Exposure = Rep = Experimental_state = Mass = err_avg_mass = avg_exp_mass = err_deut_uptake = deut_uptake = NULL

    States = unique(data$State)
    confidence_limit = 1 - significance_level

    alpha <- 1 - confidence_limit
    t_value <- qt(c(alpha/2, 1 - alpha/2), df = 2)[2]


    data_exp <- data[, .(avg_exp_mass = mean(Mass)),
                     by = list(Sequence, State, Exposure, Rep, Experimental_state)]

    data_exp <- data_exp[, .(deut_uptake = mean(avg_exp_mass), err_avg_mass = sd(avg_exp_mass)/sqrt(length(Rep))),
                         by = list(State, Sequence, Exposure, Experimental_state)]

    data_exp[, err_deut_uptake := sqrt(err_avg_mass^2 + err_avg_mass[Exposure == 0]^2),
             by = list(State, Sequence, Experimental_state)]


    calc_dat <- data_exp[, .(diff_deut_uptake = deut_uptake[Experimental_state == "A"] -
                                 deut_uptake[Experimental_state == "B"],
                             err_diff_deut_uptake = sqrt(err_deut_uptake[Experimental_state == "A"]^2 +
                                                             err_deut_uptake[Experimental_state == "B"]^2)),
                         by = list(Sequence, Exposure)]



    avg_difference <- mean(calc_dat$diff_deut_uptake)

    x_threshold <- t_value * mean(calc_dat[["err_diff_deut_uptake"]], na.rm = TRUE)/sqrt(length(calc_dat))

    data.table(Test = "Houde",
               State_1 = States[1],
               State_2 = States[2],
               Test_statistic = NA,
               P_value = NA,
               Significant_difference = abs(avg_difference) > x_threshold,
               Time = NA,
               Transformation = NA,
               AIC = NA,
               logLik = NA)
}


