
test_that(("calculate_peptide_mass works"), {
    sequence <- c("P", "P", "A", "Q", "H", "I")
    pep_mass <- calculate_peptide_mass(sequence)
    expect_equal(pep_mass, 661.35474)
})


test_that("get_approx_isotopic_distribution works", {
    sequence <- c("P", "P", "A", "Q", "H", "I")
    min_probability <- 1e-3
    iso_dist <- get_approx_isotopic_distribution(sequence, min_probability)
    expect_equal(iso_dist, list(mass = 661.35474,
                                isotopic_distribution = c(0.682463359996824,
                                                          0.252250381517293,
                                                          0.0551007838121806,
                                                          0.00888872718735956,
                                                          0.0011563382635169),
                                max_ND = 4,
                                n_exchangeable = 6L))
})

