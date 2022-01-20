

test_that("simulate_theoretical_spectra works", {
    set.seed(10)
    times <- c(5, 30, 60, 100, 500, 900)

    #classic
    spec1 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10,
                                          charge = 1:3, times = times)
    expect_equal(spec1, readRDS("spectrum.RDS"))

    #time 0
    spec2 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10,
                                          charge = 1:3, times = 0)
    expect_equal(spec2, readRDS("spectrum_0.RDS"))

    #random charge
    spec3 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10,
                                          charge = NULL, times = times)
    expect_equal(spec3, readRDS("spectrum_rand_charge.RDS"))

    #rcpp algorithm
    spec4 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10,
                                          charge = NULL, times = times, use_markov = FALSE)
    expect_equal(spec4, readRDS("spectrum_rcpp.RDS"))

    #more than one PF
    spec5 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 2:7,
                                          charge = 1:3, times = times)
    expect_equal(spec5, readRDS("spectrum_pfs.RDS"))
})



test_that("Spectra for provided parameters are simulated", {
    set.seed(10)
    times <- c(5, 30, 60, 100, 500, 900)
    spec1 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10,
                                          charge = 1:3, times = times)

    expect_identical(unique(spec1[["Exposure"]]), c(0, times)) #time points
    expect_identical(unique(spec1[["Charge"]]), 1:3) #charges
    expect_identical(unique(spec1[["Sequence"]]), "PPAQHI") #sequence
    expect_identical(unique(spec1[["PH"]]), 7.5) #pH
    expect_identical(unique(spec1[["PF"]]), 10) # protection factor
})



test_that("Spectrum at time 0 is correct", {

    expect_output(simulate_theoretical_spectra("PPAQHI", protection_factor = 10,
                                               charge = 1:3, times = 0),
                  "There is no deuteration before given time point.",
                  fixed = TRUE)
    #one time point
    set.seed(10)
    #returned
    spec1 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10,
                                          charge = 1:3, times = 0)

    #expected
    intensity <- get_approx_isotopic_distribution(strsplit("PPAQHI", "")[[1]],
                                                  min_probability = 1e-4)[[2]]
    h2_o_mass <- 18.01056
    peptide_mass <- sum(AAmonoMass[strsplit("PPAQHI", "")[[1]]]) + h2_o_mass

    spec_expected <- data.table::data.table(Exposure = 0,
                                            PH = 7.5,
                                            Intensity = rep(intensity,
                                                            each = 3),
                                            Mz = rep(peptide_mass / 1:3 + 1.007276,
                                                     length(intensity)),
                                            Charge = 1:3,
                                            Sequence = "PPAQHI",
                                            PF = 10)
    expect_equal(spec1, spec_expected)

    #more than one time point
    #returned
    spec1 <- simulate_theoretical_spectra("SPADKTNVKAAWGKVGA",
                                          protection_factor = 100,
                                          charge = 1:3,
                                          times = c(10, 30, 50))
    spec1_time0 <- spec1[Exposure == 0, ]

    #expected
    intensity <- get_approx_isotopic_distribution(strsplit("SPADKTNVKAAWGKVGA",
                                                           "")[[1]],
                                                  min_probability = 1e-4)[[2]]
    h2_o_mass <- 18.01056
    peptide_mass <- sum(AAmonoMass[strsplit("SPADKTNVKAAWGKVGA", "")[[1]]]) + h2_o_mass

    spec_expected <- data.table::data.table(Exposure = 0,
                                            PH = 7.5,
                                            Intensity = rep(intensity,
                                                            each = 3),
                                            Mz = rep(peptide_mass / 1:3 + 1.007276,
                                                     length(intensity)),
                                            Charge = 1:3,
                                            Sequence = "SPADKTNVKAAWGKVGA",
                                            PF = 100)
    expect_equal(spec1_time0, spec_expected)
})


test_that("Sequence validation works.", {
    set.seed(10)
    times <- c(5, 30, 60, 100, 500, 900)

    expect_error(simulate_theoretical_spectra("ILOVETOTOIV",
                                              protection_factor = 10,
                                              charge = 1:3,
                                              times = times),
                 "The sequence is invalid. There is no interpretation for: O.")
})






