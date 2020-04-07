#gas constant(cal/(K * mol)); the value in spreadsheet is 1.987

get_F_const = function(temp_kelvin, gas_constant) {
  Ea_A = 14000 #unit: cal/mol
  Ea_B = 17000
  Ea_W = 19000
  Q7 = (1 / temp_kelvin - 1 / 293) / gas_constant # 1/(dT * R)
  list(Fta = exp( - Q7 * Ea_A),
       Ftb = exp( - Q7 * Ea_B),
       Ftw = exp( - Q7 * Ea_W))
}

get_poly_const = function(mol_type) {
  if(!(mol_type %in% c("poly", "oligo"))) stop("Mol type must be poly / oligo")
  Ka_poly = (10^(1.4)) / 60
  Ka_oligo = Ka_poly * 2.34
  Kb_poly = (10^(9.87)) / 60
  Kb_oligo = Kb_poly * 1.35
  Kw_poly = (10^( - 1.6)) / 60
  Kw_oligo = Kw_poly * 1.585
  if(mol_type == "poly") {
    list(Ka = Ka_poly,
         Kb = Kb_poly,
         Kw = Kw_poly)
  } else {
    list(Ka = Ka_oligo,
         Kb = Kb_oligo,
         Kw = Kw_oligo)
  }
}


get_pkc = function(temp_kelvin, gas_constant, ) {
  Ea_Asp = 960 #unit: cal/mol
  Ea_Glu = 1083
  Ea_His = 7500
  pKc_Asp = -log10(10^(-3.87) * exp(-1 * Ea_Asp * ((1 / temp_kelvin - 1 / 278) / gas_constant)))
  pKc_Glu = -log10(10^(-4.33) * exp(-1 * Ea_Glu * ((1 / temp_kelvin - 1 / 278) / gas_constant)))
  pKc_His = -log10(10^(-7) * exp(-1 * Ea_His * ((1 / temp_kelvin - 1 / 278) / gas_constant)))
  list(asp = pKc_Asp,
       glu = pKc_Glu,
       his = pKc_His)
}

get_exchange_constants = function(ph, pkc_consts, k_consts) {
  constants = matrix(
    c(0, 0, 0, 0,
      -0.59, -0.32, 0.0767122542818456, 0.22,
      0, 0, 0, 0, #Asp; calculate later
      -0.90, -0.12, 0.69, 0.60,
      -0.58, -0.13, 0.49, 0.32,
      -0.54, -0.46, 0.62, 0.55,
      -0.74, -0.58, 0.55, 0.46,
      -0.22, 0.218176047120386, 0.267251569286023, 0.17,
      0, 0, 0, 0, #Glu; calculate later
      -0.6, -0.27, 0.24, 0.39,
      -0.47, -0.27, 0.06, 0.2,
      0, 0, 0, 0, #His; calculate later
      -0.91, -0.59, -0.73, -0.23,
      -0.57, -0.13, -0.576252727721677, -0.21,
      -0.56, -0.29, -0.040, 0.12,
      -0.64, -0.28, -0.00895484265292644, 0.11,
      -0.52, -0.43, -0.235859464059171, 0.0631315866300978,
      0, -0.194773472023435, 0, -0.24,
      0, -0.854416534276379, 0, 0.6,
      -0.437992277698594, -0.388518934646472, 0.37, 0.299550285605933,
      -0.79, -0.468073125742265, -0.0662579798400606, 0.20,
      -0.40, -0.44, -0.41, -0.11,
      -0.41, -0.37, -0.27, 0.050,
      -0.739022273362575, -0.30, -0.701934483299758, -0.14,
      0, -1.32,0,1.62,
      0, 0, -1.8, 0, #C -Term; col1 calculate later
      0, 0, 0, 0, # -NHMe; col1&3 calulate later
      0, 0.293, 0, -0.197
    ), ncol = 4, byrow = TRUE)

  constants[3,1] = log10(10^(-0.9 - ph)/(10^(-pkc_consts[["asp"]]) + 10^(- ph)) + 10^(0.9 -pkc_consts[["asp"]])/(10^(-pkc_consts[["asp"]]) + 10^(- ph)))
  constants[3,2] = log10(10^(-0.12 - ph)/(10^(-pkc_consts[["asp"]]) + 10^(- ph)) + 10^(0.58 -pkc_consts[["asp"]])/(10^(-pkc_consts[["asp"]]) + 10^(- ph)))
  constants[3,3] = log10(10^(0.69 - ph)/(10^(-pkc_consts[["asp"]]) + 10^(- ph)) + 10^(-0.3 -pkc_consts[["asp"]])/(10^(-pkc_consts[["asp"]]) + 10^(- ph)))
  constants[3,4] = log10(10^(0.6 - ph)/(10^(-pkc_consts[["asp"]]) + 10^(- ph)) + 10^(-0.18 -pkc_consts[["asp"]])/(10^(-pkc_consts[["asp"]]) + 10^(- ph)))

  constants[9,1] = log10(10^(-0.6 - ph)/(10^(-pkc_consts[["glu"]]) + 10^(- ph)) + 10^(-0.9 -pkc_consts[["glu"]])/(10^(-pkc_consts[["glu"]]) + 10^(- ph)))
  constants[9,2] = log10(10^(-0.27 - ph)/(10^(-pkc_consts[["glu"]]) + 10^(- ph)) + 10^(0.31 -pkc_consts[["glu"]])/(10^(-pkc_consts[["glu"]]) + 10^(- ph)))
  constants[9,3] = log10(10^(0.24 - ph)/(10^(-pkc_consts[["glu"]]) + 10^(- ph)) + 10^(-0.51 -pkc_consts[["glu"]])/(10^(-pkc_consts[["glu"]]) + 10^(- ph)))
  constants[9,4] = log10(10^(0.39 - ph)/(10^(-pkc_consts[["glu"]]) + 10^(- ph)) + 10^(-0.15 -pkc_consts[["glu"]])/(10^(-pkc_consts[["glu"]]) + 10^(- ph)))

  constants[12,1] = log10(10^(-0.8 - ph)/(10^(-pkc_consts[["his"]]) + 10^(- ph)) + 10^(0 -pkc_consts[["his"]])/(10^(-pkc_consts[["his"]]) + 10^(- ph)))
  constants[12,2] = log10(10^(-0.51 - ph)/(10^(-pkc_consts[["his"]]) + 10^(- ph)) + 10^(0 -pkc_consts[["his"]])/(10^(-pkc_consts[["his"]]) + 10^(- ph)))
  constants[12,3] = log10(10^(0.8 - ph)/(10^(-pkc_consts[["his"]]) + 10^(- ph)) + 10^(-0.1 -pkc_consts[["his"]])/(10^(-pkc_consts[["his"]]) + 10^(- ph)))
  constants[12,4] = log10(10^(0.83 - ph)/(10^(-pkc_consts[["his"]]) + 10^(- ph)) + 10^(0.14 -pkc_consts[["his"]])/(10^(-pkc_consts[["his"]]) + 10^(- ph)))

  constants[26,1] = log10(10^(0.05 - ph)/(10^(-pkc_consts[["glu"]]) + 10^(- ph)) + 10^(0.96 -pkc_consts[["glu"]])/(10^(-pkc_consts[["glu"]]) + 10^(- ph)));

  constants[27,1] = log10(135.5/(k_consts[["Ka"]] * 60))
  constants[27,3] = log10(2970000000/(k_consts[["Kb"]] * 60))
  constants
}

fbmme_dh = function(sequence, ph = 9, temperature = 15, mol_type = "poly",
                    gas_constant = 1.9858775) {
  D = 10 ^ (-ph) #[D + ]
  OD = 10 ^ (ph - 14.17) #[OD - ]
  temp_kelvin = temperature + 273.15  #convert to kelvin scale
  F_consts = get_F_const(temp_kelvin, gas_constant)
  poly_consts = get_poly_const(mol_type)
  pkc_consts = get_pkc(temp_kelvin, gas_constant)

  AAs = strsplit('ARDdNCsGEeQHILKMFPpSTWYVncma', "")[[1]]
  constants = get_exchange_constants(ph, pkc_consts, poly_consts)
  sequence = c("n", sequence, "c")
  N = length(sequence) #residue numbers of the input peptide/protein plus 2(N -  & C - term)
  kcDH = rep(0, N) # TODO: assert N > 2

  for(i in 1:N) {
    if(i == 1 || sequence[i] == 'P' || sequence[i] == 'p' || sequence[i] == 'a') {
      next()
    } else {
      if(i %in% c(1, 2, N) || sequence[i] %in% c("P", "a", "p")) {
        Fa = 0
        Fb = 0
      } else {
        j = which(AAs == sequence[i])
        k = which(AAs == sequence[i - 1])
        if(i - 2 <= 1 && i + 1 == N && sequence[i - 1] != "a" && sequence[i] != "m"){
          Fa = 10 ^ (constants[j, 1] + constants[k, 2] + constants[25, 2] + constants[26, 1])
          Fb = 10 ^ (constants[j, 3] + constants[k, 4] + constants[25, 4] + constants[26, 3])
        } else {
          if(i - 2 <= 1 && sequence[i - 1] != "a") {
            Fa = 10 ^ (constants[j, 1] + constants[k, 2] + constants[25, 2])
            Fb = 10 ^ (constants[j, 3] + constants[k, 4] + constants[25, 4])
          } else {
            if(i + 1 == N && sequence[i] != "m") {
              Fa = 10 ^ (constants[j, 1] + constants[k, 2] + constants[26, 1])
              Fb = 10 ^ (constants[j, 3] + constants[k, 4] + constants[26, 3])
            } else {
              Fa = 10^(constants[j, 1] + constants[k, 2])
              Fb = 10^(constants[j, 3] + constants[k, 4])
            }
          }
        }
      }
      kcDH[i] = Fa * poly_consts[["Ka"]] * D * F_consts[["Fta"]] +
        Fb * poly_consts[["Kb"]] * OD * F_consts[["Ftb"]] +
        Fb * poly_consts[["Kw"]] * F_consts[["Ftw"]]
    }
  }
  kcDH[2:(N - 1)]
}

