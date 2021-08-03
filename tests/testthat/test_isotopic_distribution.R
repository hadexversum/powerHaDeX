
test_that("get_approx_isotopic_distribution returns proper output", {
    sequence = c("P", "P", "A", "Q", "H", "I")
    min_probability = 1e-3
    iso_dist <- get_approx_isotopic_distribution(sequence, min_probability)

    expect_true(length(iso_dist) == 4)
    expect_equal(names(iso_dist), c("mass", "isotopic_distribution", "max_ND", "n_exchangeable"))
    expect_true(class(iso_dist) == "list")
    expect_equal(unname(sapply(iso_dist, class)), c("numeric", "numeric", "numeric", "integer"))
    expect_true(length(iso_dist[[1]]) == 1)
    expect_true(length(iso_dist[[3]]) == 1)
    expect_true(length(iso_dist[[4]]) == 1)
})

test_that("get_approx_isotopic_distribution works", {
    sequence = c("S", "P", "A", "D", "K", "T", "N", "V", "K", "A", "A", "W", "G", "K", "V", "G", "A")
    min_probability = 1e-3
    iso_dist <- get_approx_isotopic_distribution(sequence, min_probability)

    expect_true(iso_dist[[1]] == 1698.90527)
    expect_equal(iso_dist[[2]], c(0.381926443964, 0.352219070740, 0.177810926702,
                                  0.064141646884, 0.018364492284, 0.004412363138))
    expect_true(iso_dist[[3]] == 5)
    expect_true(iso_dist[[4]] == 17)
})

