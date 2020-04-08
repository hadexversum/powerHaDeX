source("fbmme_hd.R")
source("pepinfo.R")


library(phonTools)
library(signal)
library(BRAIN)

inputSeq = "KITEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELAKDPRIAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTRITK"

Charge = 1
inputPF = 1
START = 1
END = NULL
Time_p = 15/1000
pH_p = 9
QFratio = c(1, 10, 2, 2)
Temp = 15
molType = "poly"
pH = 9
TempC = 15

Factors_values
get_constants(pKc_Asp, pKc_Glu, pKc_His, Ka, Kb, pH)

do_simulation = function(inputSeq, Charge, inputPF = 1, START = 1,
                         END = NULL,
                         Time_p = 15/1000, pH_p = 9, QFratio = c(1, 10, 2, 2), Temp = 15) {

  Seq = strsplit(inputSeq, "")[[1]]
  pf = inputPF

  if(is.null(END)) {
    END = length(Seq)
  }else{
    if(length(pf) == 1L) {
      pf = rep(pf, END)
    }
  }
  ###data prep:

  pepinfoData = pepinfo(Seq[START:END])
  peptideMass = pepinfoData[1]
  distND = pepinfoData[2]
  maxND = pepinfoData[3]
  maxD = pepinfoData[4]

  kcHD = fbmme_hd(inputSeq = Seq, pH_p, Temp, 'poly', 0) # call fbmme_hd.m
  kcDH = fbmme_dh(Seq, pH_p, Temp, 'poly') # call fbmme_dh.m
  Dfraction = QFratio[1]/(QFratio[1] + QFratio[2] + QFratio[3])
  Hfraction =(QFratio[2] + QFratio[3])/(QFratio[1] + QFratio[2] + QFratio[3])

  ###do simulation:
  M = 3000 #number of simulating molecules
  N = END - START + 1
  HDmatrix = matrix(1, M, N) #0=H; 1=D
  kmax = max(max(kcDH[START:END]), max(kcHD[START:END]))
  deltaT = 0.1/kmax   #step size of simulation time
  Time_seq = seq(0, Time_p, deltaT)

  for(time in Time_seq){
    ran1 = matrix(runif(M*N), M, N)
    for(i in 1:M){
      for(j in START:END){
        if(HDmatrix[i, j - START + 1] == 0){
          if(ran1[i, j - START + 1] <= (1 - exp(-kcHD[j]*deltaT/pf[j])) * Dfraction){
            HDmatrix[i, j - START + 1] = 1
          }
        }else{
          if(HDmatrix[i, j - START + 1] == 1){
            if(ran1[i, j - START + 1] <= (1 - exp(-kcDH[j]*deltaT/pf[j])) * Hfraction){
              HDmatrix[i,j-START+1] = 0
            }
          }
        }
      }
    }
    print(time)
  }

  # consider N-term two residues and Prolines
  for (i in 1:N){
    if (i < 3 || Seq[i + START - 1] == 'P'){
      #HDmatrix[,i] <- zeros(M, 1) #all must be H
      HDmatrix[,i] = rep(0, M)
    }
  }

  Distr = rep(0, maxD[[1]] + 1)
  deltaMass = rep(0, M)

  for (i in 1:M){
    deltaMass[i] <- sum(HDmatrix[i,])
    Distr[deltaMass[i] + 1] <- Distr[deltaMass[i] + 1] + 1 # Distr(x) means x-1 units of mass above monoisotopic
  }

  Distr <- Distr/sum(Distr) #normalization

  #do convolution with allH peaks:
  obsDistr <- conv(distND[[1]], Distr)
  obsDistr <- obsDistr/sum(obsDistr) # normalization

  # get MS observable peaks:
  obsPeaks <- matrix(0, maxD[[1]] + maxND[[1]] + 1, 2)
  DM <- 1.00628 #delta mass between a deuteron(2.014102 u) and a proton(1.007825 u)
  obsPeaks[1,1] <- peptideMass[[1]]/Charge + 1.007276 # m/z of mono; 1.007276 is the mass of proton
  obsPeaks[1,2] <- obsDistr[1]

  for (i in 2:(maxD[[1]] + maxND[[1]] + 1)){
    obsPeaks[i,1] <- obsPeaks[i-1, 1] + DM/Charge
    obsPeaks[i,2] <- obsDistr[i]
  }

  data.frame(
    mz = obsPeaks[, 1],
    intensity = obsPeaks[, 2]
  )
}
