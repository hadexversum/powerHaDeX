#' Houde's test for deuteration curves
#'
#' @importFrom data.table data.table
#'
#' @description This function performs Damian Houde's confidence intervals test
#' for differences in deuteration levels. Its input and output are compatible with
#' the function \code{\link[powerHaDeX]{calculate_hdx_power}}.
#'
#' @param data data.table with deuteration curves
#' @param significance_level significance level for tests
#'
#' @returns This function returns a data table compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#'
#' @references Houde, Damian, Steven A Berkowitz, and John R Engen (2011).
#' “The utility of hydrogen/deuterium exchange mass spectrometry in biopharmaceutical
#' comparabilitystudies”. In:Journal of pharmaceutical sciences100.6, pp. 2071–2086.
#'
#' @seealso
#' Other tests:
#'
#' - \code{\link[powerHaDeX]{test_hdx_analyzer}}
#'
#' - \code{\link[powerHaDeX]{test_memhdx_model}}
#'
#' -\code{\link[powerHaDeX]{test_semiparametric}}
#'
#' Or \code{\link[powerHaDeX]{calculate_hdx_power}} for estimation of power
#' of tests for differences in deuteration levels.
#'
#' @examples
#' theo_spectra_pf_100 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 100,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#' theo_spectra_pf_200 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 200,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#'
#' theo_spectra_two_states <- rbind(theo_spectra_pf_100, theo_spectra_pf_200)
#'
#' deuteration_curves_paired_states <- get_noisy_deuteration_curves(theo_spectra_two_states,
#'                                                                  n_replicates = 1,
#'                                                                  n_experiments = 1,
#'                                                                  reference = 100)[[1]][[1]]
#' test_houde(deuteration_curves_paired_states)
#'
#' @export

test_houde <- function(data, significance_level = 0.05) {

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



#' HDX-Analyzer model
#'
#' @description This function performs the test based on the simplest linear models for deuteration
#' curves containing time, state of the protein and the interaction term. Its input
#' and output are compatible with the function \code{\link[powerHaDeX]{calculate_hdx_power}}.
#'
#' @inheritParams test_houde
#'
#' @returns This function returns a data table compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#'
#' @references  Liu, Sanmin et al. (2011). “HDX-analyzer: a novel package for statistical
#' analysis of protein structure dynamics”. In:BMC bioinformatics12.1, pp. 1–10.
#'
#' @seealso
#' Other tests:
#'
#' - \code{\link[powerHaDeX]{test_houde}}
#'
#' - \code{\link[powerHaDeX]{test_memhdx_model}}
#'
#' -\code{\link[powerHaDeX]{test_semiparametric}}
#'
#' Or \code{\link[powerHaDeX]{calculate_hdx_power}} for estimation of power
#' of tests for differences in deuteration levels.
#'
#' @examples
#' theo_spectra_pf_100 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 100,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#' theo_spectra_pf_200 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 200,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#'
#' theo_spectra_two_states <- rbind(theo_spectra_pf_100, theo_spectra_pf_200)
#'
#' deuteration_curves_paired_states <- get_noisy_deuteration_curves(theo_spectra_two_states,
#'                                                                  n_replicates = 1,
#'                                                                  n_experiments = 1,
#'                                                                  reference = 100)[[1]][[1]]
#' test_hdx_analyzer(deuteration_curves_paired_states)
#'
#' @export

test_hdx_analyzer = function(data, significance_level = 0.05) {

    States = unique(data$State)

    Time = c("continuous", "categorical", "continuous")
    Transformation = c("identity", "identity", "log")

    aic = rep(NA, 3)
    loglik = rep(NA, 3)
    F_statistic = rep(NA, 3)
    p_value = rep(NA, 3)

    # continuous, identity
    model = lm(Mass ~ Exposure*State, data = data)
    model_reduced = lm(Mass ~ Exposure, data = data)
    result = anova(model, model_reduced)
    aic[1] = AIC(model)
    loglik[1] = logLik(model)
    F_statistic[1] = result$`F`[2]
    p_value[1] = result$`Pr(>F)`[2]

    # categorical, identity
    model = lm(Mass ~ factor(Exposure)*State, data = data)
    model_reduced = lm(Mass ~ factor(Exposure), data = data)
    result = anova(model, model_reduced)
    aic[2] = AIC(model)
    loglik[2] = logLik(model)
    F_statistic[2] = result$`F`[2]
    p_value[2] = result$`Pr(>F)`[2]

    # continuous, log
    model = lm(Mass ~ log(Exposure + 1)*State, data = data)
    model_reduced = lm(Mass ~ log(Exposure + 1), data = data)
    result = anova(model, model_reduced)
    aic[3] = AIC(model)
    loglik[3] = logLik(model)
    F_statistic[3] = result$`F`[2]
    p_value[3] = result$`Pr(>F)`[2]

    data.table::data.table(Test = "Deuteros lm",
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = F_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = Time,
                           Transformation = Transformation,
                           AIC = aic,
                           logLik = loglik)
}

#' MEMHDX model
#'
#' @importFrom lmerTest lmer
#'
#' @description This function performs the test based on a linear mixed effects
#' model used in MEMHDX tools. Its input and output are compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#'
#' @inheritParams test_houde
#'
#' @returns This function returns a data table compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#'
#' @references  Hourdel, Véronique et al. (July 2016). “MEMHDX: an interactive tool
#' to expedite thestatistical validation and visualization of large HDX-MS data sets”.
#' In:Bioinformatics32.22, pp. 3413–3419.issn: 1367-4803.
#'
#' @seealso
#' Other tests:
#'
#' - \code{\link[powerHaDeX]{test_houde}}
#'
#' - \code{\link[powerHaDeX]{test_hdx_analyzer}}
#'
#' -\code{\link[powerHaDeX]{test_semiparametric}}
#'
#' Or \code{\link[powerHaDeX]{calculate_hdx_power}} for estimation of power
#' of tests for differences in deuteration levels.
#'
#' @examples
#' theo_spectra_pf_100 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 100,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#' theo_spectra_pf_200 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 200,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#'
#' theo_spectra_two_states <- rbind(theo_spectra_pf_100, theo_spectra_pf_200)
#'
#' deuteration_curves_paired_states <- get_noisy_deuteration_curves(theo_spectra_two_states,
#'                                                                  n_replicates = 1,
#'                                                                  n_experiments = 1,
#'                                                                  reference = 100)[[1]][[1]]
#' test_memhdx_model(deuteration_curves_paired_states)
#'
#' @export

test_memhdx_model = function(data, significance_level = 0.05) {

    States = unique(data$State)

    Time = c("continuous", "categorical", "continuous")
    Transformation = c("identity", "identity", "log")

    aic = rep(NA, 3)
    loglik = rep(NA, 3)
    Test_statistic = rep(NA, 3)
    p_value = rep(NA, 3)

    # continuous, identity
    model = lmer(Mass ~ Exposure*State + (1|Rep),
                           data = data,
                           REML = FALSE)
    model_reduced = lmer(Mass ~ Exposure + (1|Rep),
                                   data = data,
                                   REML = FALSE)
    result = anova(model, model_reduced)
    aic[1] = AIC(model)
    loglik[1] = logLik(model)
    Test_statistic[1] = result$Chisq[2]
    p_value[1] = result$`Pr(>Chisq)`[2]


    # categorical, identity
    model = lmer(Mass ~ factor(Exposure)*State + (1|Rep),
                           data = data,
                           REML = FALSE)
    model_reduced = lmer(Mass ~ factor(Exposure) + (1|Rep),
                                   data = data,
                                   REML = FALSE)
    result = anova(model, model_reduced)
    aic[2] = AIC(model)
    loglik[2] = logLik(model)
    Test_statistic[2] = result$Chisq[2]
    p_value[2] = result$`Pr(>Chisq)`[2]

    # continuous, log
    model = lmer(Mass ~ log(Exposure + 1)*State + (1|Rep),
                           data = data,
                           REML = FALSE)
    model_reduced = lmer(Mass ~ log(Exposure+1) + (1|Rep),
                                   data = data,
                                   REML = FALSE)
    result = anova(model, model_reduced)
    aic[3] = AIC(model)
    loglik[3] = logLik(model)
    Test_statistic[3] = result$Chisq[2]
    p_value[3] = result$`Pr(>Chisq)`[2]

    data.table(Test = "MEMHDX lmm",
               State_1 = States[1],
               State_2 = States[2],
               Test_statistic = Test_statistic,
               P_value = p_value,
               Significant_difference = (p_value <= significance_level),
               Time = Time,
               Transformation = Transformation,
               AIC = aic,
               logLik = loglik)
}

#' Get truncated lines
#' @description This function prepares the function space for semiparametric
#' model using truncated lines.
#' @param x a vector. Variable which undergoes the truncating
#' @param knots a vector of knots in which the spline model will break
#' @keywords internal

truncated_lines <- function(x, knots){
    sapply(knots, function(kappa) {
        (x - kappa)*(x > kappa)
    })
}

#' Semiparametric test for differences in deuteration levels
#'
#' @importFrom glmnet glmnet
#'
#' @description This function performs the semiparametric test for differences
#' in deuteration levels. Its input and output are compatible with
#' the function \code{\link[powerHaDeX]{calculate_hdx_power}}.
#'
#' @inheritParams test_houde
#'
#' @details This function uses \code{\link[powerHaDeX]{truncated_lines}}. The knots
#' considered in the testing procedure are chosen using ridge regression.
#'
#' @returns This function returns a data table compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#'
#' @seealso
#' Other tests:
#'
#' - \code{\link[powerHaDeX]{test_houde}}
#'
#' - \code{\link[powerHaDeX]{test_hdx_analyzer}}
#'
#' -\code{\link[powerHaDeX]{test_memhdx_model}}
#'
#' Or \code{\link[powerHaDeX]{calculate_hdx_power}} for estimation of power
#' of tests for differences in deuteration levels.
#'
#' @examples
#' theo_spectra_pf_100 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 100,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#' theo_spectra_pf_200 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 200,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#'
#' theo_spectra_two_states <- rbind(theo_spectra_pf_100, theo_spectra_pf_200)
#'
#' deuteration_curves_paired_states <- get_noisy_deuteration_curves(theo_spectra_two_states,
#'                                                                  n_replicates = 1,
#'                                                                  n_experiments = 1,
#'                                                                  reference = 100)[[1]][[1]]
#' test_semiparametric(deuteration_curves_paired_states)
#'
#' @export

test_semiparametric <- function(data, significance_level = 0.05) {

    States = unique(data$State)
    data[["id"]] <- paste0(data$Rep, data$Charge, data$Experimental_state)
    Test = aic = loglik = Test_statistic = p_value = NA

    knots <- unique(setdiff(data$Exposure, c(max(data$Exposure), min(data$Exposure))))
    X <- truncated_lines(data$Exposure, knots)

    colnames(X) <- c(paste0("knot_", as.character(knots)))

    cv_fit <- glmnet(X, data[["Mass"]], alpha = 0, lambda = 0.001)
    coefs <- coefficients(cv_fit)
    X_reduced <- cbind(intercept = 1, X)[, which(as.logical(abs(coefs) >= 2*10^(-5)))]


    model = lmer(Mass ~ Exposure*State + (1|id) + (1|Exposure) + X_reduced,
                           data = data,
                           REML = FALSE)
    model_reduced = lmer(Mass ~ Exposure + (1|id) + (1|Exposure) + X_reduced,
                                   data = data,
                                   REML = FALSE)
    result = anova(model, model_reduced)
    aic = AIC(model)
    loglik = as.numeric(logLik(model))
    Test_statistic = result$Chisq[2]
    p_value = result$`Pr(>Chisq)`[2]

    Test = "RIDGE_knots_random_intercept_id_exposure"

    data.table::data.table(Test = Test,
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = Test_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = NA,
                           Transformation = "identity",
                           AIC = aic,
                           logLik = loglik)
}




#' Test based on area under the deuteration curve for differences in deuteration levels
#'
#' @inheritParams test_houde
#'
#' @returns This function returns a data table compatible with the function
#' \code{\link[powerHaDeX]{calculate_hdx_power}}.
#'
#' @seealso Mazur, Sharlyn J and Daniel P Weber (2017). “The area between exchange
#' curves as a measure of conformational differences in hydrogen-deuterium exchange
#' mass spectrometry studies”. In:Journal of the American Society for Mass Spectrometry
#' 28.5, pp. 978–981.
#'

test_auc_test = function(data, significance_level = 0.05) {
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
        Test = "AUC test",
        State_1 = as.character(States[1]),
        State_2 = as.character(States[2]),
        Test_statistic = Test_statistic,
        P_value = P_value,
        Significant_difference = P_value < significance_level,
        Time = "continuous",
        Transformation = "log",
        AIC = NA,
        logLik = NA
    )

}

