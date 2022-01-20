#' Apply tests for HaDeX data
#'
#' @description This function converts the data from HaDeX in order to make it
#' compatible with the input of test functions and perform the testing
#' procedures of provided tests.
#'
#' @param dat data.table. The data of hdx_data class from the HaDeX package.
#' @param states a character vector containing two states from provided `dat`
#' that should be tested. By default the first two states (if exist) from `dat`
#' are chosen.
#' @param tests a list of testing functions. In the `powerHaDeX` package the
#' following tests are implemented:
#'
#' - \code{\link[powerHaDeX]{test_houde}},
#'
#' - \code{\link[powerHaDeX]{test_hdx_analyzer}},
#'
#' - \code{\link[powerHaDeX]{test_memhdx_model}},
#'
#' -\code{\link[powerHaDeX]{test_semiparametric}}.
#'
#' @returns This function returns a data table of variables:
#'
#' - \code{Test} - name of test,
#'
#' - \code{State_1}, \code{State_2} - tested states from \code{states},
#'
#' - \code{Significant_difference} - \code{TRUE} or \code{FALSE}, indicating
#' whether the null hypothesis is rejected
#'
#' - \code{Sequence} - amino acid sequence that was tested
#'
#' @export


test_hadex_data <- function(dat, states = unique(dat[["State"]])[1:2],
                            tests = list(test_houde)) {

    if(length(states[!is.na(states)]) != 2) {
        stop("Two states must be provided for pairwise testing.")
    }
    if(any(!(states %in% unique(dat[["State"]])))) {
        stop("The provided states can not be found in the given data.")
    }

    dat <- dat[State %in% states]

    dat <- convert_hadex_data(dat)

    sequences <- unique(dat[["Sequence"]])

    rbindlist(lapply(sequences, function(one_seq) {

        by_seq <- dat[Sequence == one_seq]

        rbindlist(all_tests <- lapply(tests, function(test) {

            test_for_replicate <- tryCatch({suppressMessages(suppressWarnings(
                test(by_seq,
                     significance_level = 0.05)
            ))}, error = function(e) {
                print(e)
                data.table::data.table()
            })[, list(Test, State_1, State_2, Significant_difference)]

        }))[, Sequence := one_seq]

    }))

}


#' Convert HaDeX to powerHaDeX data type.
#'
#' @description This function converts the data from HaDeX in order to make it
#' compatible with the input of test functions.
#'
#' @param dat data.table. The data of hdx_data class from the HaDeX package.
#'
#' @returns This function returns a data table..
#'
#' @seealso
#' The tests:
#'
#' - \code{\link[powerHaDeX]{test_houde}},
#'
#' - \code{\link[powerHaDeX]{test_hdx_analyzer}},
#'
#' - \code{\link[powerHaDeX]{test_memhdx_model}},
#'
#' -\code{\link[powerHaDeX]{test_semiparametric}}.
#'
#' And \code{\link[powerHaDeX]{test_hadex_data}}.
#'
#' @keywords internal

convert_hadex_data <- function(dat) {

    proton_mass <- 1.00727647

    dat[,  `:=`(exp_mass = Center * z - z * proton_mass,
                Center = NULL)]

    dat <- dat[,
               list(avg_exp_mass = weighted.mean(exp_mass, Inten,
                                                 na.rm = TRUE)),
               by = c("Sequence", "Start", "End", "MHP", "MaxUptake", "State",
                      "Exposure", "Protein", "File", "z")]

    dat[, Experimental_state := ifelse(State == unique(dat[["State"]])[1],
                                       "A", "B")]
    dat <- dat[, list(Sequence, File, State, Exposure,
                      avg_exp_mass, z, Experimental_state)]

    setnames(dat,
             c("Sequence", "File", "State", "Exposure", "avg_exp_mass",
               "z",  "Experimental_state"),
             c("Sequence", "Rep", "State", "Exposure", "Mass", "Charge",
               "Experimental_state"))

    dat

}


