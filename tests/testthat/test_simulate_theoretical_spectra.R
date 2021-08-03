
test_that("number, names and types of columns of the returned data.table are correct", {
    set.seed(10)
    times = c(5, 30, 60, 100, 500, 900)
    spec1 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10, charge = 1:3, times = times)

    expect_true(typeof(spec1) == "list")
    expect_identical(unname(sapply(spec1, class)), c("numeric", "numeric", "numeric",
                                                     "numeric", "integer", "character",
                                                     "numeric"))

    expect_identical(names(spec1), c("Exposure", "PH", "Intensity", "Mz", "Charge",
                                     "Sequence", "PF"))
    expect_equal(ncol(spec1), 7)

    #time 0
    spec1 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10, charge = 1:3, times = 0)
    expect_true(typeof(spec1) == "list")
    expect_identical(unname(sapply(spec1, class)), c("numeric", "numeric", "numeric",
                                                     "numeric", "integer", "character",
                                                     "numeric"))
    expect_identical(names(spec1), c("Exposure", "PH", "Intensity", "Mz", "Charge",
                                     "Sequence", "PF"))
    expect_equal(ncol(spec1), 7)
})




test_that("Spectra for provided parameters are simulated", {
    set.seed(10)
    times = c(5, 30, 60, 100, 500, 900)
    spec1 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10, charge = 1:3, times = times)

    expect_identical(unique(spec1[["Exposure"]]), c(0, times)) #time points
    expect_identical(unique(spec1[["Charge"]]), 1:3) #charges
    expect_identical(unique(spec1[["Sequence"]]), "PPAQHI") #sequence
    expect_identical(unique(spec1[["PH"]]), 7.5) #pH
    expect_identical(unique(spec1[["PF"]]), 10) # protection factor
})



test_that("Spectrum at time 0 is correct", {
    #one time point
    #returned
    spec1 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10, charge = 1:3, times = 0)

    #expected
    intensity <- get_approx_isotopic_distribution(strsplit("PPAQHI", "")[[1]], min_probability = 1e-4)[[2]]
    h2_o_mass <- 18.01056
    peptide_mass <- sum(AAmonoMass[strsplit("PPAQHI", "")[[1]]]) + h2_o_mass

    spec_expected <- data.table::data.table(Exposure = 0,
                                            PH = 7.5,
                                            Intensity = rep(intensity, each = 3),
                                            Mz = rep(peptide_mass / 1:3 + 1.007276, length(intensity)),
                                            Charge = 1:3,
                                            Sequence = "PPAQHI",
                                            PF = 10)
    expect_true(all.equal(spec1, spec_expected))


    #more than one time point
    #returned
    spec1 <- simulate_theoretical_spectra("SPADKTNVKAAWGKVGA",
                                          protection_factor = 100,
                                          charge = 1:3,
                                          times = c(10, 30, 50))
    spec1_time0 <- spec1[Exposure == 0, ]

    #expected
    intensity <- get_approx_isotopic_distribution(strsplit("SPADKTNVKAAWGKVGA", "")[[1]],
                                                  min_probability = 1e-4)[[2]]
    h2_o_mass <- 18.01056
    peptide_mass <- sum(AAmonoMass[strsplit("SPADKTNVKAAWGKVGA", "")[[1]]]) + h2_o_mass

    spec_expected <- data.table::data.table(Exposure = 0,
                                            PH = 7.5,
                                            Intensity = rep(intensity, each = 3),
                                            Mz = rep(peptide_mass / 1:3 + 1.007276, length(intensity)),
                                            Charge = 1:3,
                                            Sequence = "SPADKTNVKAAWGKVGA",
                                            PF = 10)
    expect_true(all.equal(spec1_time0, spec_expected))
})





