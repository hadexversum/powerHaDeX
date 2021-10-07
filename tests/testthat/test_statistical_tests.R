

test_that("test_houde works", {
    curves <- readRDS("noisy_curves.RDS")[[1]][[1]]

    expect_equal(test_houde(curves), structure(list(Test = "Houde",
                                                    State_1 = structure(1L, .Label = c("10", "20"), class = "factor"),
                                                    State_2 = structure(2L, .Label = c("10", "20"), class = "factor"),
                                                    Test_statistic = NA,
                                                    P_value = NA,
                                                    Significant_difference = TRUE,
                                                    Time = NA,
                                                    Transformation = NA,
                                                    AIC = NA,
                                                    logLik = NA),
                                               row.names = c(NA, -1L),
                                               class = c("data.table", "data.frame")))
})


test_that("test_hdx_analyzer works", {
    curves <- readRDS("noisy_curves.RDS")[[1]][[1]]

    expect_equal(test_hdx_analyzer(curves), structure(list(Test = c("Deuteros lm", "Deuteros lm", "Deuteros lm"),
                                                           State_1 = structure(c(1L, 1L, 1L), class = "factor", .Label = c("10", "20")),
                                                           State_2 = structure(c(2L, 2L, 2L), class = "factor", .Label = c("10", "20")),
                                                           Test_statistic = c(1.90596540277125, 30.4891372765257, 5.00315365023482),
                                                           P_value = c(0.151957845006165, 3.63467693806125e-26, 0.00777797073256792),
                                                           Significant_difference = c(FALSE, TRUE, TRUE),
                                                           Time = c("continuous", "categorical", "continuous"),
                                                           Transformation = c("identity", "identity", "log"),
                                                           AIC = c(168.992007835247, -404.05519750975, 8.47535423232984),
                                                           logLik = c(-79.4960039176237, 217.027598754875, 0.762322883835081)),
                                                      row.names = c(NA, -3L), class = c("data.table", "data.frame")))
})


test_that("test_memhdx_model works", {
    set.seed(10)
    curves <- readRDS("noisy_curves.RDS")[[1]][[1]]

    expect_equal(test_memhdx_model(curves), structure(list(Test = c("MEMHDX lmm", "MEMHDX lmm", "MEMHDX lmm"),
                                                           State_1 = structure(c(1L, 1L, 1L), class = "factor", .Label = c("10", "20")),
                                                           State_2 = structure(c(2L, 2L, 2L), class = "factor", .Label = c("10", "20")),
                                                           Test_statistic = c(3.86021412427414, 148.417835464572, 9.94982009055384),
                                                           P_value = c(0.145132659439139, 8.72405171449471e-29, 0.00690914042131481),
                                                           Significant_difference = c(FALSE, TRUE, TRUE),
                                                           Time = c("continuous", "categorical", "continuous"),
                                                           Transformation = c("identity", "identity", "log"),
                                                           AIC = c(170.992007835247, -404.387996166469, 10.4753542323301),
                                                           logLik = c(-79.4960039176237, 218.193998083235, 0.762322883834969)),
                                                      row.names = c(NA, -3L), class = c("data.table", "data.frame")))
})


test_that("test_semiparametric works", {
    set.seed(10)
    curves <- readRDS("noisy_curves.RDS")[[1]][[1]]
    expect_equal(test_semiparametric(curves), structure(list(Test = "RIDGE_knots_random_intercept_id_exposure",
                                                             State_1 = structure(1L, .Label = c("10", "20"), class = "factor"),
                                                             State_2 = structure(2L, .Label = c("10", "20"), class = "factor"),
                                                             Test_statistic = 18.8790884443026,
                                                             P_value = 7.9516642216805e-05,
                                                             Significant_difference = TRUE,
                                                             Time = NA,
                                                             Transformation = "identity",
                                                             AIC = -367.11483443455,
                                                             logLik = 195.557417217275),
                                                        row.names = c(NA, -1L), class = c("data.table", "data.frame")))
})

test_that("truncated_lines works", {
    set.seed(10)
    x <- rnorm(30)
    knots <- seq(min(x), max(x), 1)
    expect_equal(truncated_lines(x, knots), structure(c(2.20403300911136, 2.00103429610047, 0.813956288247019,
                                                        1.58611912238581, 2.47983196473704, 2.5750811388697, 0.977210662740044,
                                                        1.82161082069867, 0.55861415646644, 1.92880844404554, 3.28706634125666,
                                                        2.94106834619687, 1.94705328215081, 3.17273154158292, 2.92667696655335,
                                                        2.27463410466535, 1.23034298201715, 1.99013645350229, 3.11080810026361,
                                                        2.66826536300614, 1.58897620144932, 0, 1.51042090029441, 0.0662256462593565,
                                                        0.920088816638627, 1.81162528301483, 1.49773140778161, 1.31312801145184,
                                                        2.08352583194471, 1.93150630806707, 1.20403300911136, 1.00103429610047,
                                                        0, 0.586119122385812, 1.47983196473704, 1.5750811388697, 0, 0.821610820698668,
                                                        0, 0.928808444045539, 2.28706634125666, 1.94106834619687, 0.947053282150812,
                                                        2.17273154158292, 1.92667696655335, 1.27463410466535, 0.230342982017153,
                                                        0.990136453502291, 2.11080810026361, 1.66826536300614, 0.588976201449324,
                                                        0, 0.510420900294414, 0, 0, 0.811625283014828, 0.497731407781613,
                                                        0.313128011451841, 1.08352583194471, 0.931506308067068, 0.204033009111357,
                                                        0.00103429610046663, 0, 0, 0.479831964737038, 0.575081138869698,
                                                        0, 0, 0, 0, 1.28706634125666, 0.941068346196868, 0, 1.17273154158292,
                                                        0.926676966553354, 0.274634104665352, 0, 0, 1.11080810026361,
                                                        0.668265363006141, 0, 0, 0, 0, 0, 0, 0, 0, 0.0835258319447139,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.287066341256661, 0, 0, 0.17273154158292,
                                                        0, 0, 0, 0, 0.11080810026361, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ), .Dim = c(30L, 4L)))
})


