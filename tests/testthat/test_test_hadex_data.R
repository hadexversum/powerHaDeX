

test_that("test_hadex_data works", {

    hadex_dat <- readRDS("dat.RDS")[Sequence == "INITSSASQEGTRLN"]

    test_res <- test_hadex_data(hadex_dat, tests = list(test_houde, test_hdx_analyzer))

    expect_equal(test_res, structure(list(Test = c("Houde",
                                                   "Deuteros lm identity continuous",
                                                   "Deuteros lm identity categorical",
                                                   "Deuteros lm log continuous"),
                                          State_1 = c("CD160", "CD160", "CD160", "CD160"),
                                          State_2 = c("CD160_HVEM", "CD160_HVEM", "CD160_HVEM", "CD160_HVEM"),
                                          Significant_difference = c(TRUE, FALSE, TRUE, FALSE),
                                          Sequence = c("INITSSASQEGTRLN", "INITSSASQEGTRLN", "INITSSASQEGTRLN", "INITSSASQEGTRLN")),
                                     row.names = c(NA, -4L),
                                     class = c("data.table", "data.frame")))
})


test_that("test_hadex_data returns proper error", {

    hadex_dat <- readRDS("dat.RDS")[Sequence == "INITSSASQEGTRLN" & State == "CD160"]

    expect_error(test_hadex_data(hadex_dat, tests = list(test_houde, test_hdx_analyzer)),
                             "Two states must be provided for pairwise testing.")

    expect_error(test_hadex_data(hadex_dat,
                                 states = c("state", "another_state")),
                 "The provided states can not be found in the given data.")

    hadex_dat <- readRDS("dat.RDS")[Sequence == "INITSSASQEGTRLN"]

    expect_equal(test_hadex_data(hadex_dat,
                                 tests = list(test = function(data, significance_level = 0.05) stop("Error"))),
                 structure(list(Sequence = "INITSSASQEGTRLN"),
                           row.names = c(NA, -1L),
                           class = c("data.table", "data.frame")))

})

