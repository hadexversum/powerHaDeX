

#' Peptide mass
#' @description Calculate mass of undeuterated peptide
#' @param sequence character vector of amino acid sequence of a peptide
#' @details Calculates peptide mass as a sum of amino acids' from \code{sequence}
#'  masses and H2O mass (1.007825 * 2 + 15.994915 = 18.01056).
#' @keywords internal
#' @export

calculate_peptide_mass <- function(sequence) {
    h2_o_mass = 18.01056
    sum(AAmonoMass[sequence]) + h2_o_mass
}


#' Approximates isotopic distribution
#' @description Internal function used in the simulation of theoretical spectra.
#' It calculates the isotopic distribution of an undeuterated peptide
#' that is required to get an empirical distribution.
#' @param sequence character vector of amino acid sequence of a peptide
#' @param min_probability minimum isotopic probability that will be considered
#' @details Additional file \code{sysdata.RDA} contains the maximal possible occurrence
#' of the isotopes C13, N15, O18, S34 (carbon, nitrogen, oxygen, and sulfur, respectively)
#' in the respective amino acids, and their masses. Based on that, the maximal possible
#' number of molecules of the isotopes in the sequence is calculated. Peptide mass is the
#' sum of the masses of amino acids and H2O mass - as it includes the N terminal group
#' (H) and C terminal group (OH).
#'
#' Next, the distributions of mentioned isotopes are calculated under the assumption
#' that the occurrence of ith considered isotope has a binomial distribution B(n_i, p_i)
#'  with parameters n_i (maximal possible occurrence in the sequence) and p_i
#'  (natural richness - possibility of occurrence in the universe).  For the oxygen
#'  molecules, we have to take into account that oxygen occurs in a diatomic molecule.
#'  Calculation of the sulfur distribution takes into account its rare occurrence.
#'
#'  The final isotopic distribution is computed as a convolution of obtained distributions
#'  with probabilities greater than \code{min_probability}. It is a vector of probabilities
#'  of possible monoisotopic masses. The number of exchangeable amides is computed as the
#'  length of the sequence, reduced by the number of prolines located on the third of
#'  further position.
#' @importFrom signal conv
#' @return list of elements: the mass of the peptide (\code{peptide_mass}),
#' final distribution (\code{isotopic_distribution}) of the isotopes,
#' number of significant probabilities minus one (\code{max_ND}) and
#' number of exchangeable amino acids (\code{n_exchangeable}).
#' @keywords internal
#' @export
get_approx_isotopic_distribution = function(sequence, min_probability = 1e-3) {
    pC13 = 0.0111
    pN15 = 0.00364
    pO18 = 0.00205
    pS34 = 0.04293

    peptide_mass = calculate_peptide_mass(sequence)
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

    maxND = sum(finalDist >= min_probability)
    distND = finalDist[1:maxND]

    maxD = length(sequence)- sum(sequence[3:length(sequence)] == 'P')

    return(list(mass = peptide_mass,
                isotopic_distribution = distND,
                max_ND = maxND - 1,
                n_exchangeable = maxD))
}
