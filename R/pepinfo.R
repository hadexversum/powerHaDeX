library(phonTools)
library(signal)

pepinfo <- function(subSeq){

  AAshort <- c('A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y')

  AAcarbonNum <- c(3,  3,  4,  5,  9,  2,  6,  6,  6,  6,  5,  4,  12, 5,  5,  6,  3,  4,  3,  5,  11, 9)
  AAnitrogenNum <- c(1,  1,  1,  1,  1,  1,  3,  1,  2,  1,  1,  2,   3, 1,  2,  4,  1,  1,  1,  1,   2, 1)
  AAoxygenNum <- c(1, 1, 3,  3,  1,  1,  1,  1,  1,  1,  1,  2,   3, 1,  2,  1,  2,  2,  1,  1,   1, 2)
  AAsulferNum <- c(0, 1, 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,   0, 0,  0,  0,  0,  0,  0,  0,   0, 0)

  AAmonoMass <- c(71.037110,
              103.009190,
              115.026940,
              129.042590,
              147.068410,
              57.021460,
              137.058910,
              113.084060,
              128.094960,
              113.084060,
              131.040490,
              114.042930,
              255.158290,
              97.052760,
              128.058580,
              156.101110,
              87.032030,
              101.047680,
              168.964200,
              99.068410,
              186.079310,
              163.063330)

  peptideMass <- 0
  C <- 0
  N <- 0
  O <- 0
  S <- 0

  for (i in 1:length(subSeq)){
    index <- which(subSeq[i] == AAshort)
    peptideMass <- peptideMass + AAmonoMass[index]
    C <- C + AAcarbonNum[index]
    N <- N + AAnitrogenNum[index]
    O <- O + AAoxygenNum[index]
    S <- S + AAsulferNum[index]
  }

  peptideMass <- peptideMass + (1.007825 * 2 + 15.994915) #peptide's mass is the sum of the residue masses plus the mass of water.

  ## calculate maxND:
  obsCThreshold <- 1e-3 #set threshold

  pC13 <- 0.0111 #natural richness of C13
  distC <- dbinom(seq(0, C), C, pC13)

  pN15 <- 0.00364 #natural richness of N15
  distN <- dbinom(seq(0, N), N, pN15)

  pO18 <- 0.00205 #natural richness of O18
  dist <- dbinom(seq(0,O), O, pO18)
  distO <- zeros(1, 2*O+1)

  for (i in 1:(O+1)){
    distO[i*2-1] <- dist[i]
  }

  pS34 <- 0.04293  #natural richness of S34

  if (S > 0) {
    dist <- dbinom(seq(0, S), S, pS34)
    distS <- zeros(1, 2*S+1)
    for (i in 1:(S+1)){
      distS[i*2-1] <- dist[i]
    }
  } else {
    distS <- 1
  }

  finalDist <- conv(distS, conv(distO, conv(distC, distN))) # jakiego to jest wymiaru?? wektor?

  maxND <- length(finalDist) - 1  # tego nie jestem pewna??
  for (m in 3:(maxND + 1)){
    if (finalDist[m] < obsCThreshold & finalDist[m-1] < obsCThreshold & finalDist[m-2] >= obsCThreshold){
      maxND <- m-3
      break
    }
  }

  distND <- finalDist[1:(maxND+1)]
  # calculate maxD:
  maxD <- length(subSeq) - 2 # exclude N-terminal two residues
  for (m in 3:length(subSeq)){
    if (subSeq[m] == 'P'){  #exclude Proline
    maxD <- maxD - 1
    }
  }

  return(list(peptideMass, distND, maxND, maxD))
}
