#' Get probability of an exchange (H->D and D->H)
#' @param HD_rate rate of hydrogen-deuterium exchange
#' @param DH_rate rate of deuterium-hydrogen exchange (back-exchange)
#' @param time_step time step of a simulation
#' @param protection_factor protection factor
#' @return list
#' @keywords internal
#' @export
get_exchange_probabilities = function(HD_rate, DH_rate, time_step, protection_factor) {
    list(HD = (1 - exp(-HD_rate * time_step / protection_factor)),
         DH = (1 - exp(-DH_rate * time_step / protection_factor)))
}


#' Time times at which deuteration levels will be recorded
#'
#' For a sequence of exchange times (calculated based on exchange rates) and
#' times at which measurements should be taken,
#' get moments at which deuteration levels will be saved.
#'
#' @param exchange_times numeric vector of times at which exchange happens
#' @param experiment_times numeric vector of times at which deuteration levels
#' are be measured
#' @return numeric vector
#' @keywords internal
#' @export
get_recording_times = function(exchange_times, experiment_times) {
    unique(exchange_times[findInterval(experiment_times, exchange_times)])
}


#' Get a matrix of simulated exchanged hydrogens for each experiment time point
#' @param sequence amino acid sequence of a peptide as a character vector
#' @param transition_probs list of probabilities of exchange returned by the get_exchange_probabilities function
#' @param experiment_times numeric vector of times at which exchange will happen
#' @param times_to_record numeric vector of times for which deuteration level measurement should be made
#' @param n_molecules number of peptide molecules
#' @return list of matrices
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @useDynLib powerHDX, .registration = TRUE
#' @export
get_HD_matrices = function(sequence, transition_probs, experiment_times,
                           times_to_record, n_molecules = 100) {
    peptide_length = length(sequence)

    time_intervals <- cut(experiment_times, c(0, times_to_record), right = TRUE, include.lowest = TRUE)
    separated_times <- split(experiment_times, time_intervals)

    hd_matrices = vector("list", length(times_to_record))
    HDmatrix = matrix(0L, n_molecules, peptide_length) # 0 denotes hydrogen, 1 denotes deuterium

    for (i in 1:length(hd_matrices)) {
        if(length(separated_times[[i]]) == 0) {
            hd_matrices[[i]] = HDmatrix
        }else {
            HDmatrix = get_deuteration_single_timepoint(HDmatrix, separated_times[[i]],
                                                        transition_probs[[1]], transition_probs[[2]])
            HDmatrix[, unique(c(1, 2, which(sequence == "P")))] = 0
            hd_matrices[[i]] = HDmatrix
        }
    }
    hd_matrices
}


#' Get observed isotopic distribition after hydrogen-deuterium exchange
#' @param HDmatrix simulated matrix after hydrogen-deuterium-exchange
#' @param isotopic_distribution vector of isotopic probabilities of a peptide
#' @param maxD length of the isotopic distribution - 1
#' @importFrom signal conv
#' @keywords internal
#' @export
get_observed_iso_dist = function(HDmatrix, isotopic_distribution, maxD) {
    Distr = rep(0, maxD + 1)
    deltaMass = rowSums(HDmatrix)
    Distr[sort(unique(deltaMass)) + 1] = table(deltaMass)
    Distr = Distr / sum(Distr)
    obsDistr = signal::conv(isotopic_distribution, Distr)
    obsDistr = obsDistr / sum(obsDistr)
    obsDistr
}



