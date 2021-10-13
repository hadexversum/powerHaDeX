

test_that("calculate_hdx_power works", {
    curves <- readRDS("noisy_curves.RDS")

    #summarized = TRUE
    output <- calculate_hdx_power(curves, list(test_hdx_analyzer))
    expect_equal(output, readRDS("power.RDS"))

    #summarized = FALSE
    output <- calculate_hdx_power(curves, list(test_hdx_analyzer), summarized = FALSE)
    expect_equal(output, readRDS("power_summarized.RDS"))
})

test_that("calculate_hdx_power returns error for more than two states", {
    curves <- readRDS("get_noisy_deuteration_curves_jointly.RDS")

    expect_error(calculate_hdx_power(curves, list(test_hdx_analyzer)),
                 "The data table should contain at most 2 different states.")
})


test_that("calculate_hdx_power returns empty data.table for invalid test", {
    curves <- readRDS("noisy_curves.RDS")
    invalid_test <- function(data, significance_level = 0.05) stop("This function couses an error.")

    expect_equal(calculate_hdx_power(curves, list(invalid_test), summarized = FALSE),
                 readRDS("empty_power_result.RDS"))
    expect_output(calculate_hdx_power(curves, list(invalid_test), summarized = FALSE),
                  "This function couses an error.")

})

