#' Calculate power of statistical tests for HDX experiments
#'
#' @importFrom data.table rbindlist uniqueN
#'
#' @description This function estimates power of statistical tests for HDX experiments.
#'
#' @param deuteration_curves list returned by the \code{\link[powerHaDeX]{get_noisy_deuteration_curves}}
#' @param tests lists of tests to perform. Each test function should take two parameters - data
#' (data_table containing replicated curves) and \code{significance_level}, and have particular
#' output - data frame of variables: \code{Test} (name of a test which should be displayed in
#' the final result), \code{State_1}, \code{State_2} (biological states of interest), \code{Test_statistic},
#' \code{P_value}, \code{Significant_difference} (the same as \code{p_value <= significance_level}), \code{Time}
#' (character, "continuous" or "categorical"), \code{Transformation} (character, transformation that is used
#' for exposure), \code{AIC}, \code{logLik}. For example see \code{\link[powerHaDeX]{test_houde}}.
#' @param significance_level significance level that will be used for testing. See \code{tests}
#' @param summarized logical. Indicates whether the power should be calculated. Default \code{TRUE}.
#'
#' @return list of data.tables with test result, optionally summarized with power.
#'
#' @seealso \code{\link[powerHaDeX]{test_houde}}, \code{\link[powerHaDeX]{test_semiparametric}},
#' \code{\link[powerHaDeX]{test_hdx_analyzer}}, \code{\link[powerHaDeX]{test_memhdx_model}}
#'
#' @examples
#'
#' theo_spectra_pf_100 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 100,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#'
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
#'                                                                  n_replicates = 4,
#'                                                                  n_experiments = 2,
#'                                                                  compare_pairs = TRUE,
#'                                                                  reference = "all")
#' calculate_hdx_power(deuteration_curves_paired_states,
#'                     tests = list(test_houde),
#'                     summarized = TRUE)
#'
#' @export
#'
calculate_hdx_power <- function(deuteration_curves, tests, significance_level = 0.05,
                                summarized = TRUE) {

    Significant_difference <- Sequence <- Rep <- State <- Exposure <- Transformation <- Test <-
        Num_replicates <- Num_states <- Num_timepoints <- State_2 <- Time <- State_1 <- NULL

    test_results <- lapply( deuteration_curves, function(curve) {

        single_curve <- lapply(curve, function(replicate_curve) {

            if(uniqueN(replicate_curve[["State"]]) > 2) stop("The data table should contain at most 2 different states.")

            if(uniqueN(replicate_curve[["State"]]) == 1 & uniqueN(replicate_curve[["Experimental_state"]]) == 2) {

                replicate_curve[["State"]] <- paste0(replicate_curve[["State"]],
                                                     replicate_curve[["Experimental_state"]])
                info <- replicate_curve[, list(Sequence = as.character(unique(Sequence)),
                                               Num_replicates = uniqueN(Rep),
                                               Num_states = 2L,
                                               Num_timepoints = uniqueN(Exposure))]
                type_one_err = TRUE

            }else {

                type_one_err = FALSE
                info = replicate_curve[, list(Sequence = as.character(unique(Sequence)),
                                              Num_replicates = uniqueN(Rep),
                                              Num_states = uniqueN(State),
                                              Num_timepoints = uniqueN(Exposure))]
            }

            replicate_curve[["id"]] <- as.numeric(paste0(replicate_curve[["Rep"]],
                                                         replicate_curve[["Charge"]],
                                                         match(replicate_curve[["Experimental_state"]], LETTERS)))

            all_tests <- lapply(tests, function(test) {
                test_for_replicate <- tryCatch({suppressMessages(suppressWarnings(
                    test(replicate_curve,
                         significance_level = 0.05)
                ))}, error = function(e) {
                    print(e)
                    data.table::data.table()
                })
                if(type_one_err & nrow(test_for_replicate) != 0) {

                    test_for_replicate[["State_1"]] <- substr(test_for_replicate[["State_1"]], 1 ,
                                                              nchar(test_for_replicate[["State_1"]]) - 1)
                    test_for_replicate[["State_2"]] <- substr(test_for_replicate[["State_2"]], 1 ,
                                                              nchar(test_for_replicate[["State_2"]]) - 1)
                }
                cbind(info, test_for_replicate)
            })
            data.table::rbindlist(all_tests)
        })

        data.table::rbindlist(single_curve)
    })

    if (summarized) {

        test_results <- lapply(test_results, function(test_result) {

            test_result[, list(Power = mean(Significant_difference)),
                        by = list(Sequence, Num_replicates, Num_states, Num_timepoints,
                                  Test, Transformation, Time, State_1, State_2)]
        })
    }

    rbindlist(test_results)

}
