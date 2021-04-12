#' Simulate theoretical spectra of a deuterated peptide over time
#' @param sequence amino acid sequence of a peptide as a single string
#' @param charge vector of charges of the peptide ion. If NULL, oe value is sampled
#' from vector 2:6. Default NULL.
#' @param protection_factor protection factor. If a single number of provided,
#' same protection factor will be assumed for each amide. Default value: 1
#' (indicates that the exchange rate is equal to the intristic exchange rate)
#' @param times a vector of times at which deuteration levels will be measured (seconds)
#' @param pH pH of the reaction. Default to 7.5.
#' @param temperature temperature of the reaction (Celsius)
#' @param n_molecules number of peptide molecules. Default to 100.
#' @param time_step_const constant that will be dived by a maximum exchange rate
#' to obtain time step. Hydrogen-deuteration exchange may occur at each of these steps.
#' @inheritParams get_exchange_rates
#' @param min_probability smallest isotopic probability to consider
#' @param use_markov logical. If TRUE algorithm basing on Markov chain will be used.
#' If FALSE simulation provided by Zhong-Yuan Kan will be executed. Default to TRUE,
#' as it fastens the calculation
#' @return data.table
#' @export
simulate_theoretical_spectra = function(sequence, charge = NULL, protection_factor = 1,
                                        times = c(60, 600), pH = 7.5,
                                        temperature = 15, n_molecules = 100,
                                        time_step_const = 1, if_corr = 0,
                                        min_probability = 1e-4,
                                        use_markov = TRUE) {

    sequence = strsplit(sequence, "")[[1]]
    if (length(protection_factor) == 1L) {
        protection_factor = rep(protection_factor, length(sequence))
    }
    if (is.null(charge)) {
        charge = sample(2:6, 1)
    }

    peptide_iso_dist = get_approx_isotopic_distribution(sequence, min_probability)
    peptide_mass = peptide_iso_dist[[1]]
    isotopic_probs = peptide_iso_dist[[2]]
    maxND = peptide_iso_dist[[3]]
    maxD = peptide_iso_dist[[4]]
    kcHD = get_exchange_rates(sequence, "HD", pH, temperature, 'poly', if_corr)
    kcDH = get_exchange_rates(sequence, "DH", pH, temperature, 'poly')

    kmax = max(max(kcDH), max(kcHD))
    deltaT = time_step_const / kmax
    steps_between_time_points = ceiling(times/deltaT)

    if (floor(max(times)/deltaT) == 0) {
        print("There is no deuteration before given time point.")
        isotope_dists = data.frame()
    } else {
        tryCatch({
            transition_probs = get_exchange_probabilities(kcHD, kcDH, deltaT, protection_factor)

            if(use_markov) {
                HD_matrices = get_HD_matrices_using_markov(sequence, transition_probs,
                                                           steps_between_time_points,
                                                           n_molecules)
            }else {
                time_sequence = seq(0, max(times), deltaT)
                HD_matrices = get_HD_matrices(sequence, transition_probs,
                                              time_sequence, times,
                                              n_molecules)
            }

            isotope_dists = get_iso_probs_deut(HD_matrices, maxD, maxND,
                                               isotopic_probs, peptide_mass,
                                               times, charge, pH)
        })
    }
    isotope_dists[["Sequence"]] = paste0(sequence, collapse = "")
    if (length(unique(protection_factor)) == 1) {
        isotope_dists[["PF"]] = protection_factor[1]
    } else {
        isotope_dists[["PF"]] = paste(protection_factor,
                                      sep = ",", collapse = ",")
    }
    isotope_dists = isotope_dists[isotope_dists[["Intensity"]] > min_probability, ]
    data.table::as.data.table(isotope_dists)
}
