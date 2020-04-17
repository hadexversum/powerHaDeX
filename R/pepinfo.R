#' @importFrom signal conv
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
    distO = rep(0, 2*n_oxygen + 1)
    distO[seq(1, 2*n_oxygen + 1, 2)] = dist[1:(n_oxygen + 1)]

    if (!is.na(n_sulfer)) {
        dist = dbinom(seq(0, n_sulfer), n_sulfer, pS34)
        distS = rep(0, 2 * n_sulfer+1)
        distS[seq(1, 2*n_sulfer + 1, 2)] = dist[1:(n_sulfer + 1)]
    } else {
        distS = 1
    }

    finalDist = sort(conv(distS, conv(distO, conv(distC, distN))),
                     decreasing = TRUE) # jakiego to jest wymiaru?? wektor?
    maxND = length(finalDist) - 1  # tego nie jestem pewna??

    for (m in 3:(maxND + 1)){
        if (finalDist[m] < obsCThreshold & finalDist[m - 1] < obsCThreshold & finalDist[m - 2] >= obsCThreshold){
            maxND = m - 3
            break
        }
    }

    distND = finalDist[1:(maxND + 1)]
    # calculate maxD:
    maxD = length(sequence) - 2 # exclude N-terminal two residues
    maxD = maxD - sum(sequence[3:length(sequence)] == 'P')

    return(list(peptide_mass, distND, maxND, maxD))
}
