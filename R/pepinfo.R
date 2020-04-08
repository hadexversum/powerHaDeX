#sequence = strsplit("KITEGKLVIWINGDKGYNGLAEVGKK", "")[[1]]
# Internal data, left here for convenience for now:
AAshort = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "O",
            "P", "Q", "R", "S", "T", "U", "V", "W", "Y")
AAcarbonNum = c(3, 3, 4, 5, 9, 2, 6, 6, 6, 6, 5, 4, 12, 5, 5, 6, 3, 4, 3, 5, 11, 9)
AAnitrogenNum = c(1, 1, 1, 1, 1, 1, 3, 1, 2, 1, 1, 2, 3, 1, 2, 4, 1, 1, 1, 1, 2, 1)
AAoxygenNum = c(1, 1, 3, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 1, 2, 2, 1, 1, 1, 2)
AAsulferNum = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
AAmonoMass = c(71.03711, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146,
               137.05891, 113.08406, 128.09496, 113.08406, 131.04049, 114.04293,
               255.15829, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
               168.96420, 99.06841, 186.07931, 163.06333)
names(AAcarbonNum) = AAshort
names(AAnitrogenNum) = AAshort
names(AAoxygenNum) = AAshort
names(AAmonoMass) = AAshort

pepinfo = function(sequence) {
    h2_o_mass = 1.007825 * 2 + 15.994915
    obsCThreshold = 1e-3 #set threshold
    pC13 = 0.0111 #natural richness of C13
    pN15 = 0.00364 #natural richness of N15
    pO18 = 0.00205 #natural richness of O18
    pS34 = 0.04293  #natural richness of S34

    peptide_mass = sum(AAmonoMass[sequence]) + h2_o_mass
    n_carbon = sum(AAcarbonNum[sequence])
    n_nitrogen = sum(AAnitrogenNum[sequence])
    n_oxygen = sum(AAoxygenNum[sequence])
    n_sulfer = sum(AAsulferNum[sequence])

    ## calculate maxND:
    distC = dbinom(seq(0, n_carbon), n_carbon, pC13)
    distN = dbinom(seq(0, n_nitrogen), n_nitrogen, pN15)
    dist = dbinom(seq(0, n_oxygen), n_oxygen, pO18)
    distO = zeros(1, 2*n_oxygen + 1)

    for (i in 1:(n_oxygen + 1)) {
        distO[i*2 - 1] = dist[i]
    }
    if (!is.na(n_sulfer)) {
        dist = dbinom(seq(0, n_sulfer), n_sulfer, pS34)
        distS = zeros(1, 2 * sulfer+1)
        for (i in 1:(n_sulfer + 1)){
            distS[i*2 - 1] = dist[i]
        }
    } else {
        distS = 1
    }

    finalDist = sort(conv(distS, conv(distO, conv(distC, distN))),
                     decreasing = TRUE) # jakiego to jest wymiaru?? wektor?
    maxND = length(finalDist) - 1  # tego nie jestem pewna??
    for (m in 3:(maxND + 1)){
        if (finalDist[m] < obsCThreshold & finalDist[m-1] < obsCThreshold & finalDist[m - 2] >= obsCThreshold){
            maxND = m - 3
            break
        }
    }

    distND = finalDist[1:(maxND+1)]
    # calculate maxD:
    maxD = length(sequence) - 2 # exclude N-terminal two residues
    for (m in 3:length(sequence)) {
        if (sequence[m] == 'P') {  #exclude Proline
            maxD = maxD - 1
        }
    }

    return(list(peptide_mass, distND, maxND, maxD))
}

# pepinfo = function(sequence) {
#     BRAIN::useBRAIN(BRAIN::getAtomsFromSeq(sequence))
# }
