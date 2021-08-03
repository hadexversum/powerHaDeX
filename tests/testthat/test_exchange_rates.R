
#get_F_const
test_that("get_F_const returns proper output", {
    gas_constant = 1.9858775
    temp_kelvin = 291.15
    F_const <- get_F_const(temp_kelvin, gas_constant)

    expect_true(typeof(F_const) == "list")
    expect_true(length(F_const) == 3)
    expect_equal(unname(sapply(F_const, class)), c("numeric", "numeric", "numeric"))
})


test_that("get_F_const works", {
    gas_constant = 1.9858775

    temp_kelvin = 291.15
    F_const <- as.numeric(unlist(get_F_const(temp_kelvin, gas_constant)))
    expect_true(all.equal(F_const, c(0.8582290117, 0.8305682088, 0.81262475)))

    temp_kelvin = 300
    F_const <- get_F_const(temp_kelvin, gas_constant)
    expect_equal(as.numeric(unlist(F_const)), c(1.753153117, 1.977274469, 2.142389176))
})

test_that("get_F_const returns error", {
    gas_constant = 1.9858775
    temp_kelvin = 0
    expect_error(get_F_const(temp_kelvin, gas_constant))
})

#get_poly_const
test_that("get_poly_const returns error for invalid input", {
    expect_error(get_poly_const("blahblah", "HD"))
    expect_error(get_poly_const("poly", "exchange"))

    expect_error(get_poly_const("HDD"))
    expect_error(get_poly_const("oligo", "poly"))
})


test_that("get_poly_const returns proper output", {
    #exchange HD
    poly_const = get_poly_const("poly", "HD")
    oligo_const = get_poly_const("oligo", "HD")

    expect_true(typeof(poly_const) == "list")
    expect_true(length(poly_const) == 3)
    expect_equal(unname(sapply(poly_const, class)), c("numeric", "numeric", "numeric"))

    expect_true(typeof(oligo_const) == "list")
    expect_true(length(oligo_const) == 3)
    expect_equal(unname(sapply(oligo_const, class)), c("numeric", "numeric", "numeric"))

    expect_equal(names(poly_const), c("Ka", "Kb", "Kw"))
    expect_equal(names(oligo_const), c("Ka", "Kb", "Kw"))

    #exchange DH
    poly_const = get_poly_const("poly", "DH")
    oligo_const = get_poly_const("oligo", "DH")

    expect_true(typeof(poly_const) == "list")
    expect_true(length(poly_const) == 3)
    expect_equal(unname(sapply(poly_const, class)), c("numeric", "numeric", "numeric"))

    expect_true(typeof(oligo_const) == "list")
    expect_true(length(oligo_const) == 3)
    expect_equal(unname(sapply(oligo_const, class)), c("numeric", "numeric", "numeric"))

    expect_equal(names(poly_const), c("Ka", "Kb", "Kw"))
    expect_equal(names(oligo_const), c("Ka", "Kb", "Kw"))
})


test_that("get_poly_const works", {
    #HD
    poly_const = get_poly_const("poly", "HD")
    oligo_const = get_poly_const("oligo", "HD")

    poly_const = as.numeric(unlist(poly_const))
    oligo_const = as.numeric(unlist(oligo_const))

    expect_equal(poly_const, c(0.6947823, 187003076, 0.0005270463))
    expect_equal(poly_const, c(1.625791, 252454152, 0.0008353683))

    #DH
    poly_const = get_poly_const("poly", "DH")
    oligo_const = get_poly_const("oligo", "DH")

    poly_const = as.numeric(unlist(poly_const))
    oligo_const = as.numeric(unlist(oligo_const))

    expect_equal(poly_const, c(0.4186477, 123551707, 0.0004186477))
    expect_equal(poly_const, c(0.9796357, 166794804, 0.0006635567))
})


#get_pkc
test_that("get_pkc returns proper output", {
    gas_constant = 1.9858775
    temp_kelvin = 291.15
    #HD
    pkc_const_HD = get_pkc(temp_kelvin, gas_constant, exchange = "HD")
    expect_true(typeof(pkc_const_HD) == "list")
    expect_true(length(pkc_const_HD) == 3)
    expect_equal(unname(sapply(pkc_const_HD, class)), c("numeric", "numeric", "numeric"))
    expect_equal(names(pkc_const_HD), c("asp", "glu", "his"))
    #DH
    pkc_const_DH = get_pkc(temp_kelvin, gas_constant, exchange = "DH")
    expect_true(typeof(pkc_const_DH) == "list")
    expect_true(length(pkc_const_DH) == 3)
    expect_equal(unname(sapply(pkc_const_DH, class)), c("numeric", "numeric", "numeric"))
    expect_equal(names(pkc_const_DH), c("asp", "glu", "his"))
})


test_that("get_pkc works", {

    expect_equal(as.numeric(unlist(get_pkc(300, 1.9858775, "HD"))),
                 c(4.422311601, 4.867523464, 6.987337008))
    expect_equal(as.numeric(unlist(get_pkc(300, 1.9858775, "DH"))),
                 c(3.814619137, 4.267523464, 6.567337008))
})

test_that("get_pkc returns error in the case of invalid output", {
    expect_error(get_pkc(0, 1.9858775, "HD"))
    expect_error(get_pkc(300, 1.9858775, "blah"))
})


#get_exchange_constants
test_that("get_pkc returns proper output", {
    pH = 7
    pkc_consts = list(asp = 4.444469934,
                      glu = 4.891520938,
                      his = 7.153524502)
    k_consts = list(Ka = 0.6947823058,
                    Kb = 187003075.7,
                    Kw = 0.0005270462767)

    constants = get_exchange_constants(pH, pkc_consts, k_consts)
    expect_equal(dim(constants), c(28, 4))
    expect_equal(class(constants), c("matrix", "array"))
    expect_equal(typeof(constants), "double")
    expect_false(all(is.na(constants)))
})


test_that("get_exchange_constants works", {
    pH = 7
    pkc_consts = list(asp = 4.444469934,
                      glu = 4.891520938,
                      his = 7.153524502)
    k_consts = list(Ka = 0.6947823058,
                    Kb = 187003075.7,
                    Kw = 0.0005270462767)

    constants = matrix(c(0.0000000000,  0.0000000000,  0.000000000000,  0.00000000000,
                         -0.5900000000, -0.3200000000,  0.076712254282,  0.22000000000,
                         0.8988123106,  0.5790342220, -0.289554452432, -0.17398516009,
                         -0.9000000000, -0.1200000000,  0.690000000000,  0.60000000000,
                         -0.5800000000, -0.1300000000,  0.490000000000,  0.32000000000,
                         -0.5400000000, -0.4600000000,  0.620000000000,  0.55000000000,
                         -0.7400000000, -0.5800000000,  0.550000000000, 0.46000000000,
                         -0.2200000000,  0.2181760471,  0.267251569286, 0.17000000000,
                         -0.8966718129,  0.3075189990, -0.494750659376, -0.14179534052,
                         -0.6000000000, -0.2700000000,  0.240000000000,  0.39000000000,
                         -0.4700000000, -0.2700000000,  0.060000000000,  0.20000000000,
                         -0.2961582545, -0.2261567011,  0.605773780975,  0.65717231213,
                         -0.9100000000, -0.5900000000, -0.730000000000, -0.23000000000,
                         -0.5700000000, -0.1300000000, -0.576252727722, -0.21000000000,
                         -0.5600000000, -0.2900000000, -0.040000000000,  0.12000000000,
                         -0.6400000000, -0.2800000000, -0.008954842653,  0.11000000000,
                         -0.5200000000, -0.4300000000, -0.235859464059,  0.06313158663,
                         0.0000000000, -0.1947734720,  0.000000000000, -0.24000000000,
                         0.0000000000, -0.8544165343,  0.000000000000,  0.60000000000,
                         -0.4379922777, -0.3885189346,  0.370000000000,  0.29955028561,
                         -0.7900000000, -0.4680731257, -0.066257979840,  0.20000000000,
                         -0.4000000000, -0.4400000000, -0.410000000000, -0.11000000000,
                         -0.4100000000, -0.3700000000, -0.270000000000,  0.05000000000,
                         -0.7390222734, -0.3000000000, -0.701934483300, -0.14000000000,
                         0.0000000000, -1.3200000000, 0.000000000000,  1.62000000000,
                         0.9570460868,  0.0000000000, -1.800000000000,  0.00000000000,
                         0.5119392952,  0.0000000000, -0.577243550643,  0.00000000000,
                         0.0000000000,  0.2930000000,  0.000000000000, -0.19700000000),
                       byrow = TRUE, ncol = 4)

    expect_equal(get_exchange_constants(pH, pkc_consts, k_consts), constants)
})



#get_exchange_rates
test_that("get_exchange_rates returns proper output", {
    sequence = c("P", "P", "A", "Q", "H", "I")
    pH = 9
    temperature = 15
    if_corr = 0

    kcHD = get_exchange_rates(sequence, exchange = "HD", mol_type = "poly")
    expect_true(length(kcHD) == 6)
    expect_true(class(kcHD) == "numeric")
    expect_true(typeof(kcHD) == "double")

    kcDH = get_exchange_rates(sequence, exchange = "DH", mol_type = "poly")
    expect_true(length(kcDH) == 6)
    expect_true(class(kcDH) == "numeric")
    expect_true(typeof(kcDH) == "double")
})


test_that("get_exchange_rates returns error in the case of invalid output", {
    sequence = c("P", "P", "A", "Q", "H", "I")
    pH = 9
    temperature = -273.15
    if_corr = 0
    expect_error(get_exchange_rates(sequence, exchange = "HD", mol_type = "poly", temperature = temperature))
    expect_error(get_exchange_rates(sequence, exchange = "blah", mol_type = "poly"))
    expect_error(get_exchange_rates(sequence, exchange = "DH", mol_type = "blahblah"))
    expect_error(get_exchange_rates(sequence, if_corr = 17))
})


test_that("get_exchange_rates works", {
    sequence = c("P", "P", "A", "Q", "H", "I")

    get_exchange_rates(sequence)

    expect_equal(get_exchange_rates(sequence),
                 c(0, 0, 58.6518885939, 117.0259030127, 142.6053257511, 0.4411805721))

    expect_equal(get_exchange_rates(sequence, exchange = "DH"),
                 c(0, 0, 293.955039644, 586.517412896, 670.603937240, 2.130998827))

})

