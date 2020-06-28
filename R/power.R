#' Calculate power of statistical tests for HDX experiments
#'
#' @param deuteration_curves list returned by the `get_noisy_deuteration_curves`
#' @param tests lists of tests to perform. Each test function should have
#' @param significance_level significance level that will be used for testing
#'
#' @return list of data.tables with test result, optionally summarized with power
#' @importFrom data.table rbindlist uniqueN
#' @export
#'
calculate_hdx_power = function(deuteration_curves, tests, significance_level = 0.05,
                               summarized = TRUE) {
    Significant_difference = Sequence = Rep = State = Exposure = NULL

    test_results = lapply(
        deuteration_curves, function(curve) {
            single_curve = lapply(curve, function(replicate_curve) {
                all_tests = lapply(tests, function(test) {
                    info = replicate_curve[, list(Sequence = unique(Sequence),
                                                  Num_replicates = uniqueN(Rep),
                                                  Num_states = uniqueN(State),
                                                  Num_timepoints = uniqueN(Exposure))]
                    cbind(info, test(replicate_curve, significance_level = 0.05))
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
            test_result[, .(Power = mean(Significant_difference)), by = grouping_columns]
        })
    }
    rbindlist(test_results)
}
