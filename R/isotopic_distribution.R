#' Approximates isotopic distribution via sampling
#' @param sequence amino acid sequence of a peptide
#' @param min_probability minimum isotopic probability that will be considered
#' @importFrom signal conv
#' @return list
#' @keywords internal
#' @export
get_approx_isotopic_distribution = function(sequence, min_probability = 1e-3) {
    h2_o_mass = 1.007825 * 2 + 15.994915
    pC13 = 0.0111
    pN15 = 0.00364
    pO18 = 0.00205
    pS34 = 0.04293

    peptide_mass = sum(AAmonoMass[sequence]) + h2_o_mass
    n_carbon = sum(AAcarbonNum[sequence])
    n_nitrogen = sum(AAnitrogenNum[sequence])
    n_oxygen = sum(AAoxygenNum[sequence])
    n_sulfer = sum(AAsulferNum[sequence])

    distC = dbinom(seq(0, n_carbon), n_carbon, pC13)
    distN = dbinom(seq(0, n_nitrogen), n_nitrogen, pN15)
    dist = dbinom(seq(0, n_oxygen), n_oxygen, pO18)
    distO = rep(0, 2*n_oxygen + 1)
    distO[seq(1, 2*n_oxygen + 1, 2)] = dist[1:(n_oxygen + 1)]

    if (!is.na(n_sulfer)) {
        dist = dbinom(seq(0, n_sulfer), n_sulfer, pS34)
        distS = rep(0, 2 * n_sulfer + 1)
        distS[seq(1, 2*n_sulfer + 1, 2)] = dist[1:(n_sulfer + 1)]
    } else {
        distS = 1
    }

    finalDist = sort(signal::conv(distS, signal::conv(distO, signal::conv(distC, distN))),
                     decreasing = TRUE)

    # maxND = length(finalDist) - 1
    #
    # for (m in 3:(maxND + 1)) {
    #     if (finalDist[m] < min_probability & finalDist[m - 1] < min_probability & finalDist[m - 2] >= min_probability) {
    #         maxND = m - 3
    #         break
    #     }
    # }
    # distND = finalDist[1:(maxND + 1)]

    maxND = sum(finalDist >= min_probability)
    distND = finalDist[1:maxND]

    maxD = length(sequence)- sum(sequence[3:length(sequence)] == 'P')

    return(list(mass = peptide_mass,
                isotopic_distribution = distND,
                max_ND = maxND,
                n_exchangeable = maxD))
}
