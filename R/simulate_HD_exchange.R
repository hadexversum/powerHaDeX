#' Get probability of an exchange (HD and DH)
#'
#' @description Calculate probabilities of exchanges that are required to
#' simulate the exchange process.
#'
#' @param HD_rate rate of hydrogen-deuterium exchange calculated via
#' \code{\link[powerHaDeX]{get_exchange_rates}}
#' @param DH_rate rate of deuterium-hydrogen exchange (back-exchange) calculated
#' via \code{\link[powerHaDeX]{get_exchange_rates}}
#' @param time_step size of a single time step of a simulation
#' @inheritParams simulate_theoretical_spectra
#'
#' @details he process is defined as a series of steps from the time sequence,
#' and each step depends on the state in the previous one.  Therefore, the
#' probabilities of changing the state are conditional probabilities -
#' probabilities of particular state in k+1^th step given particular state in
#' k^th step. For i^th amino acid probabilities are calculated as follows
#'
#' \deqn{ P_i( H -> D ) = 1 - exp((-kcHD_i * time_step)/protection_factor)}
#' \deqn{ P_i( D -> H ) = 1 - exp((-kcDH_i * time_step)/protection_factor)}
#' \deqn{ P_i( H -> H ) = 1 - P_i( H -> D )}
#' \deqn{ P_i( D -> D ) = 1 - P_i( D -> H )}
#' where the last two equations describe the probabilities of staying in the
#' same state.
#'
#' @return a list of four vectors: vector \code{HD} for probabilities
#' \eqn{P_i( H -> D )}, vector \code{DH} for probabilities \eqn{P_i( D -> H )},
#' vector \code{HH} for probabilities \eqn{P_i( H -> H )} and vector \code{DD}
#' for probabilities \eqn{P_i( D -> D)}.
#'
#' @keywords internal
#'
get_exchange_probabilities <- function(HD_rate,
                                       DH_rate,
                                       time_step,
                                       protection_factor) {

    if(any(protection_factor == 0)) stop("Protection factor can not be 0.")

    transition_probs <- list(HD = (1 - exp(-HD_rate * time_step / protection_factor)),
                             DH = (1 - exp(-DH_rate * time_step / protection_factor)))
    transition_probs[["HH"]] <- 1 - transition_probs[["HD"]]
    transition_probs[["DD"]] <- 1 - transition_probs[["DH"]]

    transition_probs

}


#' Get a matrix of simulated exchanged hydrogens for each experiment time point
#'
#' @importFrom Rcpp evalCpp
#'
#' @description Calculate matrices of simulated exchange required for obtaining
#' empirical distribution.
#'
#' @param sequence amino acid sequence of a peptide as a character vector
#' @param transition_probs list of probabilities of exchange returned by the
#' \code{\link[powerHaDeX]{get_exchange_probabilities}} function
#' @param experiment_times numeric vector of times at which exchange will happen
#' @param times_to_record numeric vector of times for which deuteration level
#' measurement should be made
#' @param n_molecules number of peptide molecules
#'
#' @return  Matrices are stored in a list of matrices (\code{HD_matrices}) -
#' each matrix for the respective time point of the measurement \code{times}.
#'
#' @details  At each time point in the time sequence:
#'
#' change \code{H} to \code{D} with probability \eqn{P(H -> D)} in each entry of
#' the matrix from the previous iteration,
#'
#' change \code{D} to \code{H} with probability \eqn{P(D -> H)} in each entry of
#' the matrix from the previous iteration.
#'
#' @keywords internal
#'
#' @useDynLib powerHaDeX, .registration = TRUE
#'

get_HD_matrices <- function(sequence, transition_probs, experiment_times,
                            times_to_record, n_molecules = 100) {

    peptide_length <- length(sequence)

    time_intervals <- cut(experiment_times, c(0, times_to_record), right = TRUE,
                          include.lowest = TRUE)
    separated_times <- split(experiment_times, time_intervals)

    hd_matrices <- vector("list", length(times_to_record))
    HDmatrix <- matrix(0L, n_molecules, peptide_length) # 0 denotes hydrogen, 1 denotes deuterium

    for (i in 1:length(hd_matrices)) {

        if(length(separated_times[[i]]) == 0) {

            hd_matrices[[i]] <- HDmatrix

        }else {

            HDmatrix <- get_deuteration_single_timepoint(HDmatrix,
                                                         separated_times[[i]],
                                                         transition_probs[[1]],
                                                         transition_probs[[2]])
            HDmatrix[, unique(c(1, 2, which(sequence == "P")))] = 0
            hd_matrices[[i]] <- HDmatrix

        }
    }

    hd_matrices

}

#' Get a matrix of simulated exchanged hydrogens for each experiment time point
#' using markov chains
#'
#' @description Calculate matrices of simulated exchange required for obtaining
#' empirical distribution.
#'
#' @param steps_between_time_points A vector containing sum of steps between
#' times at which deuteration levels are measured.
#' @inheritParams get_HD_matrices
#'
#' @details The improvement is based on the observation that the considered
#' process is a Markov chain with transition probabilities \eqn{P(H -> D)} and
#' \eqn{P(D -> H)}, and states \code{H} and \code{D}.
#'
#' @return  Matrices are stored in a list of matrices (\code{HD_matrices}) -
#' each matrix for the respective time point of the measurement \code{times}.
#'
#' Using the distributions of process in given times of mesurements, states
#' \code{H} or \code{D} are sampled for \code{m} peptide molecules
#' (\code{n_molecules)} for each of \eqn{i = 1,..., n} amino acids and stored
#' in a \code{m x n} dimensional matrix for each of the time points of the
#' measurement given by \code{times}.
#'
#' @keywords internal
#'
get_HD_matrices_using_markov <- function(sequence, transition_probs,
                                         steps_between_time_points,
                                         n_molecules = 100) {

    peptide_length <- length(sequence)
    initial_state <- c(1, 0)

    hd_matrices <- lapply(1:length(steps_between_time_points), function(i) {

        steps <- sum(steps_between_time_points[1:i])

        HDmatrix <- sapply(1:peptide_length, function(amino_acid) {

            transition_matrix <- matrix(c(transition_probs[["HH"]][amino_acid],
                                          transition_probs[["HD"]][amino_acid],
                                          transition_probs[["DH"]][amino_acid],
                                          transition_probs[["DD"]][amino_acid]),
                                        nrow = 2,
                                        byrow = TRUE)
            probabilities <- as.vector(initial_state%*%(transition_matrix%^%steps))
            sample(c(0, 1), n_molecules, replace = TRUE, prob = probabilities)

        })

        HDmatrix[, unique(c(1, 2, which(sequence == "P")))] <- 0
        HDmatrix

    })

    hd_matrices
}


#' Get observed distribution of ions
#'
#' @importFrom signal conv
#'
#' @description Calculate isotopic probabilities (intensity).
#'
#' @param HDmatrix simulated matrix after hydrogen-deuterium-exchange
#' @param isotopic_distribution vector of isotopic probabilities of a peptide
#' @inheritParams get_iso_probs_deut
#'
#' @details The exchangeable-hydrogen distribution describing the increase of
#' the mass is obtained from the exchange matrix from
#' \code{\link[powerHaDeX]{get_HD_matrices}} or
#' \code{\link[powerHaDeX]{get_HD_matrices_using_markov}} and the number of
#' exchangeable hydrogens \code{n_exchangeable}. First, the numbers of hydrogens
#' exchanged in each molecule are calculated as sums of rows of the exchange
#' matrix. Next, a  vector of the counts is built and stored in a vector of
#' length \code{n_exchangeable} plus one (for the lack of exchange). To obtain
#' fractions counts are averaged.
#'
#' The isotopic probabilities for the deuterated peptide are computed as the
#' convolution of obtained distribution and the isotopic distribution for the
#' undeuterated peptide (\code{isotopic_distribution}) as it is a sum of those
#' variables (Claesen and Burzykowski 2017, Deconvolution-Based Approach).
#' Namely
#' \deqn{M_delta = M_{mol} - M_{mon}}
#' where \code{M_mol} is the random variable describing molecular mass,
#' \code{M_mon} is the random variable describing monoisotopic mass and
#' \code{M_delta} is the random variable describing the increase in mass.
#' @return a vector of observed isotopic distribution (\code{observed_dist}) and
#' the observed peaks for mass spectrum (observed isotopic probabilities).
#'
#' @keywords internal
#'
get_observed_iso_dist <- function(HDmatrix, isotopic_distribution, maxD) {

    Distr <- rep(0, maxD + 1)
    deltaMass <- rowSums(HDmatrix)
    Distr[sort(unique(deltaMass)) + 1] <- table(deltaMass)
    Distr <- Distr / sum(Distr)
    obsDistr <- signal::conv(isotopic_distribution, Distr)
    obsDistr <- obsDistr / sum(obsDistr)

    obsDistr

}


#' Calculate isotopic probabilities (intensity) and mass-to-charge ratio (m/z).
#'
#' @description Compute the isotopic probabilities for the deuterated peptide
#' as a convolution of the isotopic distribution for the undeuterated peptide
#' and the observed isotopic distribution after hydrogen-deuterium exchange
#' computed by \code{get_observed_iso_dist}.
#'
#' @param HD_matrices list. Simulated matrices for every time point after
#' hydrogen-deuterium-exchange. Calculated via
#' \code{\link[powerHaDeX]{get_HD_matrices}} or
#' \code{\link[powerHaDeX]{get_HD_matrices_using_markov}}
#' @param maxD length of the sequence - amount of prolines
#' @param maxND length of the isotopic distribution - 1
#' @param peptide_mass mass of the peptide + mass of H2O
#' @param isotopic_probs the isotopic distribution for the undeuterated peptide.
#' @inheritParams simulate_theoretical_spectra
#'
#' @details  The m/z values for the deuterated peptide are calculated using the
#' \code{peptide_mass}, \code{charge} and constants - deuteron mass (1.00628)
#' and proton mass (1.007276). Starting from the m/z value for the monoisotopic
#' peak, the difference between the mass of deuteron and proton divided by the
#' charge of the peptide ion is added.
#'
#' @return The output is a data frame with the variables:
#'
#' \code{Exposure} (time point of measurement consistent with given HD matrix),
#'
#' \code{Mz} - m/z values,
#'
#' \code{Intensity} - isotopic probabilities,
#'
#' \code{PH} - pH.
#'
#' @keywords internal
#'

get_iso_probs_deut <- function(HD_matrices, maxD, maxND, isotopic_probs,
                               peptide_mass, times, charge, pH) {

    isotope_dists <- lapply(1:length(times), function(ith_time) {

        do.call(rbind, lapply(1:length(charge), function(ith_charge) {

            observed_dist <- get_observed_iso_dist(HD_matrices[[ith_time]],
                                                   isotopic_probs, maxD)
            observed_peaks <- matrix(0, maxD + maxND + 1, 2)
            DM <- 1.00628
            observed_peaks[1, 1] <- peptide_mass / charge[ith_charge] + 1.007276
            observed_peaks[1, 2] <- observed_dist[1]

            for (i in 2:(maxD + maxND + 1)) {

                observed_peaks[i, 1] <- observed_peaks[i - 1, 1] + DM / charge[ith_charge]
                observed_peaks[i, 2] <- observed_dist[i]

            }

            data.frame(Exposure = times[ith_time],
                       Mz = observed_peaks[, 1],
                       Intensity = observed_peaks[, 2],
                       PH = pH,
                       Charge = charge[ith_charge]
            )
        }))
    })

    isotope_dists <- data.table::rbindlist(isotope_dists)

    isotope_dists

}

