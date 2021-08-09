
#get_F_const
test_that("get_F_const works", {
    F_const <- get_F_const(291.15, 1.9858775)
    expect_equal(F_const, list(Fta = 0.858229011688035,
                               Ftb = 0.830568208785807,
                               Ftw = 0.812624750043961))
    F_const <- get_F_const(300, 1.9858775)
    expect_equal(F_const, list(Fta = 1.75315311725701,
                               Ftb = 1.97727446855674,
                               Ftw = 2.14238917600377))
})

test_that("get_F_const returns error", {
    expect_error(get_F_const(0, 1.9858775))
})


#get_poly_const
test_that("get_poly_const returns error for invalid input", {
    expect_error(get_poly_const("blahblah", "HD"))
    expect_error(get_poly_const("poly", "exchange"))

    expect_error(get_poly_const("HDD"))
    expect_error(get_poly_const("oligo", "poly"))
})

test_that("get_poly_const works", {
    #exchange HD
    poly_const = get_poly_const("poly", "HD")
    expect_equal(poly_const, list(Ka = 0.694782305783893,
                                  Kb = 187003075.716994,
                                  Kw = 0.00052704627669473))
    oligo_const = get_poly_const("oligo", "HD")
    expect_equal(oligo_const, list(Ka = 1.62579059553431,
                                   Kb = 252454152.217942,
                                   Kw = 0.000835368348561147))

    #exchange DH
    poly_const = get_poly_const("poly", "DH")
    expect_equal(poly_const, list(Ka = 0.41864773858493,
                                  Kb = 123551706.883486,
                                  Kw = 0.00041864773858493))
    oligo_const = get_poly_const("oligo", "DH")
    expect_equal(oligo_const, list(Ka = 0.979635708288736,
                                   Kb = 166794804.292706,
                                   Kw = 0.000663556665657114))
})


#get_pkc
test_that("get_pkc returns proper output", {
    gas_constant = 1.9858775
    temp_kelvin = 291.15
    #HD
    pkc_const_HD = get_pkc(temp_kelvin, gas_constant, exchange = "HD")
    expect_equal(pkc_const_HD, list(asp = 4.44446993365174,
                                    glu = 4.89152093814484,
                                    his = 7.15352450238806))
    #DH
    pkc_const_DH = get_pkc(temp_kelvin, gas_constant, exchange = "DH")
    expect_equal(pkc_const_DH, list(asp = 3.83589113630567,
                                    glu = 4.29152093814484,
                                    his = 6.73352450238806))
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
    expect_equal(constants, structure(c(0, -0.59, 0.898812310601082, -0.9, -0.58, -0.54,
                                        -0.74, -0.22, -0.896671812875742, -0.6, -0.47, -0.296158254546754,
                                        -0.91, -0.57, -0.56, -0.64, -0.52, 0, 0, -0.437992277698594,
                                        -0.79, -0.4, -0.41, -0.739022273362575, 0, 0.957046086760691,
                                        0.511939295200356, 0, 0, -0.32, 0.57903422200311, -0.12, -0.13,
                                        -0.46, -0.58, 0.218176047120386, 0.307518999008732, -0.27, -0.27,
                                        -0.22615670112472, -0.59, -0.13, -0.29, -0.28, -0.43, -0.194773472023435,
                                        -0.854416534276379, -0.388518934646472, -0.468073125742265, -0.44,
                                        -0.37, -0.3, -1.32, 0, 0, 0.293, 0, 0.0767122542818456, -0.289554452431933,
                                        0.69, 0.49, 0.62, 0.55, 0.267251569286023, -0.494750659376026,
                                        0.24, 0.06, 0.605773780975473, -0.73, -0.576252727721677, -0.04,
                                        -0.00895484265292644, -0.235859464059171, 0, 0, 0.37, -0.0662579798400606,
                                        -0.41, -0.27, -0.701934483299758, 0, -1.8, -0.577243550643321,
                                        0, 0, 0.22, -0.173985160086455, 0.6, 0.32, 0.55, 0.46, 0.17,
                                        -0.141795340523512, 0.39, 0.2, 0.657172312134108, -0.23, -0.21,
                                        0.12, 0.11, 0.0631315866300978, -0.24, 0.6, 0.299550285605933,
                                        0.2, -0.11, 0.05, -0.14, 1.62, 0, 0, -0.197), .Dim = c(28L, 4L
                                        )))
})


#get_exchange_rates
test_that("get_exchange_rates returns proper output", {
    sequence = c("P", "P", "A", "Q", "H", "I")
    pH = 9
    temperature = 15
    if_corr = 0

    kcHD = get_exchange_rates(sequence, exchange = "HD", mol_type = "poly")
    expect_equal(kcHD, c(0, 0, 58.6518885939027, 117.025903012734, 142.605325751056,
                         0.441180572094274))

    kcDH = get_exchange_rates(sequence, exchange = "DH", mol_type = "poly")
    expect_equal(kcDH, c(0, 0, 293.955039643536, 586.51741289567, 670.603937240368,
                         2.13099882731867))
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
