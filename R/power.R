#' Calculate power of statistical tests for HDX experiments
#' @description This function estimates power of statistical tests for HDX experiments.
#'
#' @param deuteration_curves list returned by the \code{\link[powerHDX]{get_noisy_deuteration_curves}}
#' @param tests lists of tests to perform. Each test function should take two parameters - data
#' (data_table containing replicated curves) and \code{significance_level}, and have particular
#' output - data frame of variables: \code{Test} (name of a test which should be displayed in
#' the final result), \code{State_1}, \code{State_2} (biological states of interest), \code{Test_statistic},
#' \code{P_value}, \code{Significant_difference} (the same as \code{p_value <= significance_level}), \code{Time}
#' (character, "continuous" or "categorical"), \code{Transformation} (character, transformation that is used
#' for exposure), \code{AIC}, \code{logLik}.
#'
#' @param significance_level significance level that will be used for testing. See \code{tests}
#' @param summarized logical. Indicates whether the power should be calculated. Default \code{TRUE}
#' @return list of data.tables with test result, optionally summarized with power.
#' @importFrom data.table rbindlist uniqueN
#' @export
#'
calculate_hdx_power = function(deuteration_curves, tests, significance_level = 0.05,
                               summarized = TRUE) {
    Significant_difference = Sequence = Rep = State = Exposure = NULL

    test_results = lapply(
        deuteration_curves, function(curve) {
            single_curve = lapply(curve, function(replicate_curve) {

                if(uniqueN(replicate_curve[["State"]]) == 1 & uniqueN(replicate_curve[["Experimental_state"]]) == 2) {

                    replicate_curve[["State"]] = paste0(replicate_curve[["State"]],
                                                        replicate_curve[["Experimental_state"]])
                    info = replicate_curve[, list(Sequence = unique(Sequence),
                                                  Num_replicates = uniqueN(Rep),
                                                  Num_states = 2,
                                                  Num_timepoints = uniqueN(Exposure))]
                    type_one_err = TRUE

                }else {
                    type_one_err = FALSE
                    info = replicate_curve[, list(Sequence = unique(Sequence),
                                                  Num_replicates = uniqueN(Rep),
                                                  Num_states = uniqueN(State),
                                                  Num_timepoints = uniqueN(Exposure))]
                }
                replicate_curve[, id := as.numeric(paste0(Rep, Charge, match(Experimental_state, LETTERS)))]

                all_tests = lapply(tests, function(test) {
                    test_for_replicate = tryCatch({suppressMessages(suppressWarnings(
                        test(replicate_curve,
                             significance_level = 0.05)
                    ))}, error = function(e) {
                        print(e)
                        data.table::data.table()
                    })
                    if(type_one_err) {
                        test_for_replicate[, State_1 := substr(State_1, 1 , nchar(State_1) - 1)]
                        test_for_replicate[, State_2 := substr(State_2, 1 , nchar(State_2) - 1)]
                    }

                    cbind(info, test_for_replicate)
                })
                data.table::rbindlist(all_tests)
            })
            data.table::rbindlist(single_curve)
        }
    )

    if (summarized) {
        test_results = lapply(test_results, function(test_result) {
            grouping_columns = setdiff(colnames(test_result),
                                       c("Significant_difference", "P_value", "Estimated"))
            test_result[, plyr::.(Power = mean(Significant_difference)), by = grouping_columns]
        })
    }
    rbindlist(test_results)
}
