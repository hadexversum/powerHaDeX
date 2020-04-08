source("fbmme_hd.R")
source("fbmme_dh.R")
source("pepinfo.R")

library(phonTools)
library(signal)
library(BRAIN)

sequence = "KITEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELAKDPRIAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTRITK"
sequence = "KITEGKLVIWINGDKGYNGLAEV"
charge = 1
pf = 1
time_p = 15/1000
ph = 9
QFratio = c(1, 10, 2, 2)
temperature_C = 15
molType = "poly"
ph = 9




get_transition_probs = function(kcHD, kcDH, QFratio, deltaT) {
  Dfraction = QFratio[1]/(QFratio[1] + QFratio[2] + QFratio[3])
  Hfraction =(QFratio[2] + QFratio[3])/(QFratio[1] + QFratio[2] + QFratio[3])
  
  list(HD = (1 - exp(- kcHD * deltaT / pf)) * Dfraction,
       DH = (1 - exp(- kcDH * deltaT / pf)) * Hfraction)
}



get_HD_matrix = function(sequence, time_sequence, transition_probs, M = 3000, N = length(sequence)) {
  
  HDmatrix = matrix(1, N, M) #0=H; 1=D
  
  for(time in time_sequence) {
    rand_unif_matrix = matrix(runif(M * N), N, M)
    HD_zeros = (HDmatrix == 0)
    HD_ones = (HDmatrix == 1)
    HDmatrix[HD_zeros & rand_unif_matrix <= transition_probs[["HD"]]] = 1
    HDmatrix[HD_ones  & rand_unif_matrix <= transition_probs[["DH"]]] = 0
    print(time)
  }
  
  HDmatrix = t(HDmatrix)
  # consider N-term two residues and Prolines
  for (i in 1:N) {
    if (i < 3 || sequence[i] == "P") {
      #HDmatrix[,i] <- zeros(M, 1) #all must be H
      HDmatrix[, i] = rep(0, M)
    }
  }
  HDmatrix
}




get_obsDistr = function(HDmatrix, M = 3000) {
 
  Distr = rep(0, maxD + 1)
  deltaMass = rep(0, M)
  for (i in 1:M){
    deltaMass[i] = sum(HDmatrix[i, ])
    Distr[deltaMass[i] + 1] = Distr[deltaMass[i] + 1] + 1 # Distr(x) means x-1 units of mass above monoisotopic
  }
  Distr = Distr / sum(Distr) #normalization
  #do convolution with allH peaks:
  obsDistr = conv(distND, Distr)
  obsDistr = obsDistr / sum(obsDistr) # normalization
  obsDistr
}


do_simulation = function(sequence, charge, pf = 1, time_p = 15/1000, 
                         ph = 9, QFratio = c(1, 10, 2, 2), temperature_C = 15) {
  
  sequence = strsplit(sequence, "")[[1]]
  pf = rep(pf, length(sequence))
  ###data prep:
  pepinfoData = pepinfo(sequence)
  peptideMass = pepinfoData[1][[1]]
  distND = pepinfoData[2][[1]]
  maxND = pepinfoData[3][[1]]
  maxD = pepinfoData[4][[1]]
  kcHD = fbmme_hd(sequence, ph, temperature_C, 'poly', 0) # call fbmme_hd.R
  kcDH = fbmme_dh(sequence, ph, temperature_C, 'poly') # call fbmme_dh.R
  
  ###do simulation:
  M = 3000 #number of simulating molecules
  kmax = max(max(kcDH), max(kcHD))
  deltaT = 0.1/kmax   #step size of simulation time
  time_sequence = seq(0, time_p, deltaT)
  
  transition_probs = get_transition_probs(kcHD, kcDH, QFratio, deltaT)
  HDmatrix = get_HD_matrix(sequence, time_sequence, transition_probs, M)
  obsDistr = get_obsDistr(HDmatrix, M)
  
  # get MS observable peaks:
  obsPeaks = matrix(0, maxD + maxND + 1, 2)
  DM = 1.00628 #delta mass between a deuteron(2.014102 u) and a proton(1.007825 u)
  obsPeaks[1, 1] = peptideMass / charge + 1.007276 # m/z of mono; 1.007276 is the mass of proton
  obsPeaks[1, 2] = obsDistr[1]

  for (i in 2:(maxD + maxND + 1)){
    obsPeaks[i, 1] = obsPeaks[i - 1, 1] + DM / charge
    obsPeaks[i, 2] = obsDistr[i]
  }

  data.frame(
    mz = obsPeaks[, 1],
    intensity = obsPeaks[, 2]
  )
}


do_simulation(sequence, charge, pf = 1, time_p = 15/1000, 
                         ph = 9, QFratio = c(1, 10, 2, 2), temperature_C = 15)





