
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


test_that("add_noise_to_one_spectrum works", {
    set.seed(10)
    spectra_exp_design = readRDS("spectra_exp.RDS")[[6]]

    spectra = add_noise_to_one_spectrum(spectra_exp_design, undeuterated_mass = 661.35474,
                                        mass_deviations = 50, intensity_deviations = NULL)

    expect_equal(spectra, readRDS("add_noise_to_one_spectrum.RDS"))
})


test_that("get_deuteration_curves_from_spectra works", {
    set.seed(10)
    spectra = readRDS("spectra_noise_added.RDS")
    curves = get_deuteration_curves_from_spectra(spectra)
    expect_equal(curves, readRDS("get_deuteration_curves_from_spectra.RDS"))
})


test_that("get_deuteration_curve_single_spectrum works", {
    set.seed(10)
    curves = readRDS("spectra_noise_added.RDS")[[1]][[1]]
    curves = get_deuteration_curve_single_spectrum(curves)
    expect_equal(curves, readRDS("get_deuteration_curve_single_spectrum.RDS"))
})


test_that("add_noise_to_curves works", {
    set.seed(10)
    curves = readRDS("get_deuteration_curves_from_spectra.RDS")
    curves = add_noise_to_curves(curves)
    expect_equal(curves, readRDS("add_noise_to_curves.RDS"))
})


test_that("add_noise_to_single_curve works", {
    set.seed(10)
    curves = readRDS("get_deuteration_curves_from_spectra.RDS")[[1]][[1]]
    curves = add_noise_to_single_curve(curves)
    expect_equal(curves, readRDS("add_noise_to_single_curve.RDS"))
})


test_that("fix_columns_names_types works", {
    set.seed(10)
    curves = readRDS("add_noise_to_curves.RDS")
    curves = fix_columns_names_types(curves)
    expect_equal(curves, readRDS("fix_columns_names_types.RDS"))
})


test_that("get_relative_mass works", {
    set.seed(10)
    mass = rnorm(20)
    time = 0:19
    rel_mass = mass - mass[time == 0]
    expect_equal(rel_mass, c(0, -0.20299871301089, -1.39007672086434, -0.617913886725544,
                             0.275798955625682, 0.371048129758341, -1.22682234637131, -0.382422188412689,
                             -1.64541885264492, -0.275224565065818, 1.0830333321453, 0.737035337085511,
                             -0.256979726960545, 0.968698532471564, 0.722643957441997, 0.0706010955539952,
                             -0.973690027094204, -0.213896555609065, 0.906775091152254, 0.464232353894785
    ))
})

