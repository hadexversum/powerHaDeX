
test_that("get_noisy_deuteration_curves works", {
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
    expect_equal(curves, readRDS("noisy_curves.RDS"))
})


test_that("get_undeuterated_mass works", {
    spectra = readRDS("spectra_3_states.RDS")
    expect_equal(get_undeuterated_mass(spectra), 661.35474)
})


test_that("get_spectra_list works", {
    spectra = readRDS("spectra_3_states.RDS")

    expect_equal(get_spectra_list(spectra), readRDS("spectra_list.RDS"))
    expect_equal(get_spectra_list(spectra, compare_pairs = TRUE, reference = 'all'),
                 readRDS("spectra_list1.RDS"))
})

test_that("get_paired_spectra works", {
    spectra = readRDS("spectra_3_states.RDS")
    expect_equal(get_paired_spectra(spectra, reference = 'all'),
                 readRDS("spectra_list1.RDS"))
})

test_that("add_noise_to_spectra works", {
    set.seed(10)
    undeuterated_mass = 661.35474
    spectra = readRDS("spectra_list.RDS")
    spectra = add_noise_to_spectra(spectra, n_replicates = 4, n_experiments = 10,
                                   undeuterated_mass, mass_deviations = 50,
                                   intensity_deviations = NULL)
    expect_equal(spectra, readRDS("spectra_noise_added.RDS"))
})

test_that("make_experimental_design works", {
    set.seed(10)
    undeuterated_mass = 661.35474
    spectra = readRDS("spectra_list1.RDS")
    spectra = make_experimental_design(spectra, n_replicates = 4)

    expect_equal(spectra, readRDS("spectra_exp.RDS"))
})


test_that("make_noisy_spectra works", {
    set.seed(10)
    spectra_exp_design = readRDS("spectra_exp.RDS")

    spectra = make_noisy_spectra(spectra_exp_design, n_experiments = 10,
                                 undeuterated_mass = 661.35474,
                                 mass_deviations = 50, intensity_deviations = NULL)

    expect_equal(spectra, readRDS("make_noisy_spectra.RDS"))
})




