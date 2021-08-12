

test_that("calculate_hdx_power works", {
    curves = readRDS("noisy_curves.RDS")

    #summarized = TRUE
    output <- calculate_hdx_power(curves, list(test_hdx_analyzer))
    expect_equal(output, readRDS("power.RDS"))

    #summarized = FALSE
    output <- calculate_hdx_power(curves, list(test_hdx_analyzer), summarized = FALSE)
    expect_equal(output, readRDS("power_summarized.RDS"))
})

