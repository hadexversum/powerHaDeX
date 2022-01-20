
test_that("experimental file is created", {
    set.seed(44)

    peptides <- data.frame(sequence = c("FPTTKTY", "LVRKDLQN"),
                           protection_factor = c(10, 100))
    res <- create_experimental_file(peptides, charge = 1:3)

    expected <- readRDS("experimental_file.RDS")
    expect_equal(res, expected)
})


test_that("experimental returns error", {
    set.seed(44)

    peptides <- data.frame()
    expect_error(create_experimental_file(peptides, charge = 1:3),
                 "at least one sequence")
})
