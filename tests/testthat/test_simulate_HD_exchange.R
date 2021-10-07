
#get_exchange_probabilities
test_that("get_exchange_probabilities works", {
    HD_rate <- c(0, 0, 58.6518885939, 117.0259030127, 142.6053257511, 0.4411805721)
    DH_rate <- c(0, 0, 293.955039644, 586.517412896, 670.603937240, 2.130998827)
    time_step <- 0.0001
    protection_factor <- 100
    probs <- get_exchange_probabilities(HD_rate, DH_rate, time_step, protection_factor)
    expect_equal(probs, list(HD = c(0, 0, 5.86501686055119e-05, 0.00011701905574879,
                                    0.000142595158095005, 4.41180474819447e-07),
                             DH = c(0, 0, 0.000293911839094441,
                                    0.000586345445180547, 0.000670379132674093, 2.13099655643756e-06
                             ),
                             HH = c(1, 1, 0.999941349831394, 0.999882980944251, 0.999857404841905,
                                    0.999999558819525),
                             DD = c(1, 1, 0.999706088160906, 0.999413654554819,
                                    0.999329620867326, 0.999997869003444)))
})

test_that("get_exchange_probabilities returns error", {
    HD_rate <- c(0, 0, 58.6518885939, 117.0259030127, 142.6053257511, 0.4411805721)
    DH_rate <- c(0, 0, 293.955039644, 586.517412896, 670.603937240, 2.130998827)
    time_step <- 0.0001
    protection_factor <- 0
    expect_error(get_exchange_probabilities(HD_rate, DH_rate, time_step, protection_factor))
})


#get_HD_matrices
test_that("get_HD_matrices works", {
    set.seed(17)
    sequence <- c("P", "P", "A", "Q", "H", "I")
    transition_probs = list(HD = c(0, 0, 0.0000586502, 0.0001170191, 0.0001425952, 0.0000004412),
                            DH = c(0, 0, 0.000293912, 0.000586345, 0.000670379, 0.000002131),
                            HH = c(1, 1, 0.9999413498, 0.9998829809, 0.9998574048, 0.9999995588),
                            DD = c(1, 1, 0.9997060882, 0.9994136546, 0.9993296209, 0.9999978690))
    experiment_times <- seq(0, 30, 0.0001)
    times_to_record <- c(10, 20, 30)
    n_molecules <- 100
    matrices <- get_HD_matrices(sequence, transition_probs, experiment_times,
                               times_to_record, n_molecules = 100)
    expect_equal(matrices, readRDS("matrices.RDS"))

    #close times of measurements
    experiment_times <- seq(0, 30, 0.01)
    times_to_record <- c(10, 10.0001)
    n_molecules <- 100
    matrices <- get_HD_matrices(sequence, transition_probs, experiment_times,
                               times_to_record, n_molecules = 100)
    expect_equal(matrices, readRDS("matrices_close_times.RDS"))
})


#get_HD_matrices_using_markov
test_that("get_HD_matrices_using_markov returns proper output", {
    set.seed(17)
    sequence <- c("P", "P", "A", "Q", "H", "I")
    transition_probs <- list(HD = c(0, 0, 0.0000586502, 0.0001170191, 0.0001425952, 0.0000004412),
                            DH = c(0, 0, 0.000293912, 0.000586345, 0.000670379, 0.000002131),
                            HH = c(1, 1, 0.9999413498, 0.9998829809, 0.9998574048, 0.9999995588),
                            DD = c(1, 1, 0.9997060882, 0.9994136546, 0.9993296209, 0.9999978690))
    experiment_times <- seq(0, 30, 0.0001)
    times_to_record <- c(10, 20, 30)
    n_molecules <- 100
    steps_between_time_points <- ceiling(times_to_record/0.0001)

    matrices <- get_HD_matrices_using_markov(sequence, transition_probs, steps_between_time_points,
                                            n_molecules)
    expect_equal(matrices, readRDS("matrices_markov.RDS"))
})


