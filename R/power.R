get_noisy_deuteration_curves = function(theoretical_spectrum, n_runs, n_replicates,
                                        mass_deviations, intensity_deviations = NULL,
                                        deuteration_deviations = NULL, 
                                        min_probability = 1e-4, relative = TRUE) {
    theoretical_spectrum = theoretical_spectrum[theoretical_spectrum$intensity > min_probability, ]
    perturbed_spectra = perturb_spectra(theoretical_spectrum, mass_deviations,
                                        intensity_deviations, n_replicates, n_runs)
    deuteration_curves = get_deuteration_curves(perturbed_spectra)
    deuteration_curves = perturb_curves(deuteration_curves, deuteration_deviations)
    deuteration_curves = get_relative_mass(deuteration_curves, relative)
    deuteration_curves
}

get_relative_mass = function(deuterations_list, relative) {
    if (relative) {
        deuterations_list = lapply(deuterations_list, function(curves) {
            by_timepoint = split(curves, curves[, c("Rep", "Sequence", "Charge", "State")])
            relative_masses = lapply(by_timepoint, function(single_timepoint) {
                single_timepoint$Mass = single_timepoint$Mass - single_timepoint$Mass[single_timepoint$Exposure == 0]
                single_timepoint
            })
            do.call("rbind", c(relative_masses, make.row.names = FALSE))
        })
    }
    deuterations_list
}


#' Calculate power for a given set of tests and protection factors
#' @param sequence amino acid sequence of a peptide (as a single string).
#' @param protection_factors vector of protection factors that will be used in simulations.
#' @param tests a list of tests that will be evaluated.
#' @param time_points time points for simulation.
#' @param n_runs number of replicate runs of the mass spectrometer that will be simulated.
#' @param n_replicates number of replicates that will be used for power calculation.
#' @param mass_deviations a scalar or a vector of magnitudes of error in m/z dimension in ppm.
#' @param intensity_deviations standard deviations for isotopic probabilities, which
#' are used in place of intensities in simulations. If `NULL`, only error on the m/z dimension
#' will be used. If scalar, the same distribution will be used for each time point.
#' If a vector of parameters is provided, the random error will have a different variance
#' at each time point.
#' @param deuteration_deviations standard deviation for normal distribution that will be
#' used to model errors in deuteration curves. If `NULL`, only error on the m/z dimension
#' will be used. If scalar, the same distribution will be used for each time point.
#' If a vector of parameters is provided, the random error will have a different variance
#' at each time point.
#' @param ... additional parameters passed to `simulate_theoretical_spectra` function.
#' @return list
#' @export
calculate_hdx_power = function(sequence, protection_factors, tests,
                               time_points, n_runs, n_replicates,
                               mass_deviations, intensity_deviations = NULL,
                               deuteration_deviations = NULL, ...) {
    theoretical_spectra = lapply(protection_factors, function(PF) {
        simulate_theoretical_spectra(
            sequence = sequence, protection_factor = PF, times = time_points, ...
        )
    })
    theoretical_spectra = do.call("rbind", c(theoretical_spectra,
                                             make.row.names = FALSE))
    perturbed_spectra = perturb_spectra(theoretical_spectra, mass_deviations,
                                        intensity_deviations, n_replicates, n_runs)
    deuteration_curves = get_deuteration_curves(perturbed_spectra)
    deuteration_curves = perturb_curves(deuteration_curves, deuteration_deviations)
    deuteration_curves
}


perturb_curves = function(deuteration_curves, deuteration_deviations) {
    if (is.null(deuteration_deviations)) {
        return(deuteration_curves)
    }
    lapply(deuteration_curves, function(curve) {
        curve$Mass = curve$Mass + rnorm(nrow(curve), sd = deuteration_deviations)
        curve
    })
}


perturb_spectra = function(spectrum, mass_deviations, intensity_deviations,
                           n_replicates, n_runs) {
    undeuterated_mass = unique(spectrum$charge[spectrum$time == 0] * (spectrum$mz[spectrum$time == 0] - 1.007276))[1]
    spectrum = split(spectrum, list(spectrum$protection_factor, spectrum$time))
    replicates = lapply(1:n_replicates, function(ith_replicate) {
        replicate = lapply(1:n_runs, function(ith_run) {
            all_times = lapply(1:length(spectrum), function(ith_time) { # Remove error in the first measurement?
                perturb_single_timepoint(spectrum[[ith_time]],
                                         mass_deviations,
                                         intensity_deviations,
                                         undeuterated_mass) # + charge?
            })
            perturbed = do.call("rbind", c(spectrum[1], all_times,
                                           make.row.names = FALSE))
            perturbed$Rep = ith_run
            perturbed
        })
        do.call("rbind", c(replicate, make.row.names = FALSE))
    })
    replicates
}


perturb_single_timepoint = function(spectrum, mass_deviation,
                                    intensity_deviation, undeuterated_mass) {
    spectrum$mz = spectrum$mz + rnorm(nrow(spectrum),
                                      sd = (mass_deviation * undeuterated_mass) / 1e6)
    if (!is.null(intensity_deviation)) {
        spectrum$intensity = spectrum$intensity + rnorm(nrow(spectrum),
                                                        sd = intensity_deviation)
    }
    spectrum
}


get_deuteration_curves = function(replicated_experiments) {
    lapply(replicated_experiments, function(ith_replicate) {
        ith_replicate$mz = ith_replicate$charge * (ith_replicate$mz - 1.007276)
        columns_names = colnames(ith_replicate)[!colnames(ith_replicate) %in% c("mz", "intensity")]
        curves = lapply(split(ith_replicate,
                              ith_replicate[, columns_names]),
                        function(x) {
                            if (nrow(x) > 0) {
                                deuteration = weighted.mean(x$mz, x$intensity)
                                result = cbind(unique(x[, columns_names]), Mass = deuteration)
                                colnames(result) = c("Exposure", "pH", "Sequence",
                                                     "State", "Charge", "Rep", "Mass")
                                result[, c("Rep", "Sequence", "Charge", "State",
                                           "Exposure", "Mass")]
                            } else {
                                data.frame()
                            }
                        })
        do.call("rbind", c(curves, make.row.names = FALSE))
    })
}
