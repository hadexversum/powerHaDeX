#' Simulate theoretical spectra of a deuterated peptide over time
#' @param sequence amino acid sequence of a petpide as a single string
#' @param charge charge of the peptide ion.
#' @param protection_factor protection factor. If a single number of provided,
#' same protection factor will be assumed for each amide.
#' @param times times at which deuteration levels will be measured (seconds)
#' @param pH pH of the reaction
#' @param temperature temperature of the reaction (Celsius)
#' @param n_molecules number of peptide molecules
#' @param time_step_const constant that will be dived by a maximum exchange rate
#' to obtain time step. Hydrogen-deuteration exchange may occur at each of these steps.
#' @param if_corr correction factor for pH - 0 or 1
#' @param min_probability smallest isotopic probability to consider
#' @return data.frame
#' @export
simulate_theoretical_spectra = function(sequence, charge, protection_factor = 1,
                                        times = c(60, 600), pH = 7.5,
                                        temperature = 15, n_molecules = 100,
                                        time_step_const = 1, if_corr = 0,
                                        min_probability = 1e-4) {

    sequence = strsplit(sequence, "")[[1]]
    if (length(protection_factor) == 1L) {
        protection_factor = rep(protection_factor, length(sequence))
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
    time_sequence = seq(0, max(times), deltaT)

    if (time_sequence == 0 && length(time_sequence) == 1) {
        print("There is no deuteration before given time point.")
        isotope_dists = data.frame()
    } else {
        tryCatch({
            times_to_record = get_recording_times(time_sequence, times)
            times_to_record = setdiff(times_to_record, 0)

            transition_probs = get_exchange_probabilities(kcHD, kcDH, deltaT, protection_factor)
            HD_matrices = get_HD_matrices(sequence, transition_probs,
                                          time_sequence, times_to_record,
                                          n_molecules)

            isotope_dists = lapply(1:length(times_to_record), function(ith_time) {
                observed_dist = get_observed_iso_dist(HD_matrices[[ith_time]], isotopic_probs, maxD)
                observed_peaks = matrix(0, maxD + maxND + 1, 2)
                DM = 1.00628
                observed_peaks[1, 1] = peptide_mass / charge + 1.007276
                observed_peaks[1, 2] = observed_dist[1]

                for (i in 2:(maxD + maxND + 1)) {
                    observed_peaks[i, 1] = observed_peaks[i - 1, 1] + DM / charge
                    observed_peaks[i, 2] = observed_dist[i]
                }
                data.frame(
                    time = times[ith_time],
                    mz = observed_peaks[, 1],
                    intensity = observed_peaks[, 2],
                    pH = pH
                )
            })
            isotope_dists = do.call("rbind", isotope_dists)
        })

    }

    isotope_dists = rbind(data.frame(time = 0,
                                     mz = peptide_mass / charge + 1.007276,
                                     intensity = isotopic_probs,
                                     pH = pH),
                          isotope_dists)
    isotope_dists$sequence = paste0(sequence, collapse = "")
    isotope_dists$protection_factor = protection_factor[1] # TODO: what if it's a vector?
    isotope_dists$charge = charge
    isotope_dists = isotope_dists[isotope_dists$intensity > min_probability, ]
    isotope_dists
}
