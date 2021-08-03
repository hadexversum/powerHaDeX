

test_that("number of columns and types of columns of the returned data.table are correct", {
    set.seed(10)
    times = c(5, 30, 60, 100, 500, 900)

    spec1 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10, charge = 1:3, times = times)
    spec2 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 20, charge = 1:3, times = times)
    spec3 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 30, charge = 1:3, times = times)

    spectra <- rbind(spec1, spec2, spec3)

    curves <- get_noisy_deuteration_curves(spectra,
                                           compare_pairs = TRUE,
                                           reference = "all",
                                           n_replicates = 4,
                                           n_experiments = 10)
    #summarized = TRUE
    output <- calculate_hdx_power(curves, list(hdx_analyzer))
    expect_true(typeof(output) == "list")
    expect_identical(unname(sapply(output, class)), c("character", "integer", "integer", "integer",
                                                      "character", "character", "character",
                                                      "factor", "factor", "numeric"))
    expect_equal(ncol(output), 10)

    #summarized = FALSE
    output <- calculate_hdx_power(curves, list(hdx_analyzer), summarized = FALSE)
    expect_true(typeof(output) == "list")
    expect_identical(unname(sapply(output, class)), c("character", "integer", "integer", "integer",
                                                      "character", "factor", "factor", "numeric",
                                                      "numeric", "logical", "character", "character",
                                                      "numeric", "numeric" ))
    expect_equal(ncol(output), 14)
})







