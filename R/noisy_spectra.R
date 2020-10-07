#' Create a list of lists of noisy deuteration curves from theoretical spectra
#'
#' @param theoretical_spectra data.table of spectra created by the
#' `simulate_theoretical_spectra` function
#' @param compare_pairs if FALSE, all groups (defined by the protection factor)
#' will be considered jointly. If TRUE (default), each protection factor will be
#' considered together with the protection factor given by the `reference` parameter
#' @param reference protection factor that will be used for comparison to other
#' protection factors in
#' @param n_runs number of technical replicates to create
#' @param n_replicates number of replicates of a spectrum for power calculation
#' @param mass_deviations mass deviation in parts per million. Either a single number
#' (then the error at each time point will be the same) or a vector of the same length
#' as number of unique time points in the experiment. The error will be sampled from
#' normal distribution with standard deviation equal to `mass_deviations`*`undeuterated_mass`/1e6
#' @param intensity_deviations optional, standard deviations of random noise that
#' will be added to intensities. Either a single number (then the error at each
#' time point will be the same) or a vector of the same length as number of unique
#' time points in the experiment. The error will be sampled from normal distribution
#' with these standard deviations.
#' @param per_run_deviations optional, standard deviations of random noise that
#' will be added to deuteration curves. Either a single number (then the error at each
#' time point will be the same) or a vector of the same length as number of unique
#' time points in the experiment. The error will be sampled from normal distribution
#' with these standard deviations.
#' @param relative logical, if TRUE (default), each deuteration curve will start at 0
#'
#' @import data.table
#' @return list of lists of data.tables
#' @export
#'
get_noisy_deuteration_curves = function(theoretical_spectra,
                                        compare_pairs = TRUE, reference = NA,
                                        n_runs = 4, n_replicates = 100,
                                        mass_deviations = 50,
                                        intensity_deviations = NULL,
                                        per_run_deviations = NULL,
                                        relative = TRUE) {
    undeuterated_mass = get_undeuterated_mass(theoretical_spectra)
    spectra = get_spectra_list(theoretical_spectra, compare_pairs, reference)
    spectra = add_noise_to_spectra(spectra, n_runs, n_replicates, undeuterated_mass,
                                   mass_deviations, intensity_deviations)
    curves = get_deuteration_curves_from_spectra(spectra)
    curves = add_noise_to_curves(curves, per_run_deviations, relative)
    curves = fix_columns_names_types(curves)
    curves
}
# Create a list of spectra ----
#' Get a list of spectra
#' @inheritParams get_noisy_deuteration_curves
#' @return list of data.tables
#' @keywords internal
get_spectra_list = function(theoretical_spectra, compare_pairs, reference) {
    if (compare_pairs & (is.na(reference) |
                         !(as.integer(reference) %in% as.integer(unique(theoretical_spectra$PF))))) {
        stop("With pair comparisons, reference protection factor must be provided")
    }

    if (compare_pairs) {
        spectra_list = get_paired_spectra(theoretical_spectra, reference)
    } else {
        spectra_list = list(theoretical_spectra)
    }
    unname(spectra_list)
}
#' Get data.tables of spectra for pairs of protection factors
#' @inheritParams get_noisy_deuteration_curves
#' @return list of data.tables
#' @keywords internal
get_paired_spectra = function(theoretical_spectra, reference) {
    PF = NULL
    reference_spectrum = theoretical_spectra[abs(PF - reference) < 1e-9, ]
    theoretical_spectra = split(theoretical_spectra[abs(PF - reference) > 1e-9, ],
                                theoretical_spectra[abs(PF - reference) > 1e-9, ]$PF)
    theoretical_spectra = lapply(theoretical_spectra,
                                 function(spectrum) rbind(spectrum, reference_spectrum))
    unname(theoretical_spectra)
}
# Add noise to spectra in mass and intensity domains ----
#' Create a list of spectra with replicates from a a list of spectra
#' @param spectra list created by the `get_spectra_list` function
#' @param n_runs number of technical replicates in the experiment
#' @return list of data.tables
#' @keywords internal
make_experimental_design = function(spectra, n_runs) {
    lapply(spectra, function(spectrum) {
        data.table::rbindlist(lapply(1:n_runs, function(i) {
            spectrum$Rep = as.character(i)
            spectrum
        }))
    })
}
#' Creates spectra with technical replicates and noise
#' @param spectra list of spectra created by the `get_spectra_list` function
#' @param n_runs number of technical replicates to create
#' @param n_replicates number of replicates of a spectrum for power calculation
#' @param undeuterated_mass undeuterated mass of a peptide
#' @param mass_deviations mass deviation in parts per million. Either a single number
#' (then the error at each time point will be the same) or a vector of the same length
#' as number of unique time points in the experiment. The error will be sampled from
#' normal distribution with standard deviation equal to `mass_deviations`*`undeuterated_mass`/1e6
#' @param intensity_deviations optional, standard deviations of random noise that
#' will be added to intensities. Either a single number (then the error at each
#' time point will be the same) or a vector of the same length as number of unique
#' time points in the experiment. The error will be sampled from normal distribution
#' with these standard deviations.
#' @return list of lists of data.tables
#' @keywords internal
add_noise_to_spectra = function(spectra, n_runs, n_replicates, undeuterated_mass,
                                mass_deviations, intensity_deviations) {
    spectra = make_experimental_design(spectra, n_runs)
    make_noisy_spectra(spectra, n_replicates, undeuterated_mass,
                       mass_deviations, intensity_deviations)
}
#' Creates spectra with technical replicates and noise
#' @inheritParams add_noise_to_spectra
#' @return list of lists of data.tables
#' @keywords internal
make_noisy_spectra = function(spectra, n_replicates, undeuterated_mass,
                              mass_deviations, intensity_deviations) {
    lapply(spectra, function(spectrum) {
        lapply(1:n_replicates, function(ith_replicate) {
            add_noise_to_one_spectrum(spectrum, undeuterated_mass,
                                      mass_deviations, intensity_deviations)
        })
    })
}
#' Add noise to a single spectrum
#' @inheritParams add_noise_to_spectra
#' @return data.table
#' @keywords internal
add_noise_to_one_spectrum = function(spectrum, undeuterated_mass,
                                     mass_deviations, intensity_deviations) {
    if (length(mass_deviations) == 1) {
        mass_deviations = rep(mass_deviations, data.table::uniqueN(spectrum$Exposure))
    }
    sds = mass_deviations * undeuterated_mass / 1e6
    names(sds) = as.character(unique(spectrum$Exposure))

    spectrum = spectrum[, add_noise_to_one_timepoint(.SD, sds), by = "Exposure",
                        .SDcols = colnames(spectrum)]

    if (!is.null(intensity_deviations)) {
        if (length(intensity_deviations) == 1) {
            intensity_deviations = rep(intensity_deviations,
                                       data.table::uniqueN(spectrum$Exposure))
        }
        names(intensity_deviations) = names(sds)
        spectrum = spectrum[, add_noise_to_intensities(.SD, sds), by = "Exposure"]
    }
    spectrum
}
# Calculate deuteration curves ----
#' Create a list of deuteration curves from a list of spectra
#' @param spectra list created by `get_spectra_list` and `add_noise_to_spectra`
#' functions.
#' @return list of lists of data.tables
#' @keywords internal
#' @export
#'
get_deuteration_curves_from_spectra = function(spectra) {
    lapply(spectra, function(spectrum) {
        lapply(spectrum, function(replicate) {
            get_deuteration_curve_single_spectrum(replicate)
        })
    })
}
#' Get a deuteration curve from a single spectrum
#' @param spectrum data.table with a spectrum
#' @return data.table
#' @keywords internal
get_deuteration_curve_single_spectrum = function(spectrum) {
    Mz = Intensity = Charge = NULL

    grouping_columns = setdiff(colnames(spectrum), c("Mz", "Intensity", "Charge"))
    spectrum = spectrum[, list(Mass = weighted.mean(Charge * (Mz - 1.007276),
                                                    Intensity)),
                        by = grouping_columns]
    spectrum
}
# Noise for deuteration curves ----
#' Adds noise to deuteration curves
#' @param curves list of lists of deuteration curves
#' @param per_run_deviations vector of standard deviations of per-replicate error
#' @param relative if TRUE (default), each curve will start at 0
#' @return list of lists of data.tables
#' @keywords internal
add_noise_to_curves = function(curves, per_run_deviations, relative) {
    lapply(curves, function(curve) {
        lapply(curve, function(replicate_curve) {
            add_noise_to_single_curve(replicate_curve, per_run_deviations,
                                      relative)
        })
    })
}
#' Adds noise to a single deuteration curve
#' @param replicate_curve `data.table` with technical replicates of a deuteration curve
#' @param per_run_deviations vector of standard deviations of random error
#' @param relative if TRUE (default), each curve will start at 0
#' @return `data.table`
#' @keywords internal
add_noise_to_single_curve = function(replicate_curve, per_run_deviations, relative) {
    Exposure = Mass = NULL
    if (!is.null(per_run_deviations)) {
        if (length(per_run_deviations) == 1) {
            per_run_deviations = rep(per_run_deviations, length(unique(replicate_curve$Rep)))
        }
        names(per_run_deviations) = as.character(unique(replicate_curve$Rep))

        replicate_curve = replicate_curve[, add_noise(.SD, per_run_deviations),
                                          by = "Rep", .SDcols = colnames(replicate_curve)]
        replicate_curve = replicate_curve[, -1, with = FALSE]
    }
    if (relative) {
        grouping_columns = setdiff(colnames(replicate_curve),
                                   c("Mass", "Charge", "Exposure"))
        replicate_curve = replicate_curve[, list(Exposure, Mass = get_relative_mass(Mass, Exposure)),
                                          by = grouping_columns]
    }
    replicate_curve
}
# Smaller utility functions ----
#' Get mass of an undeuterated peptide
#' @param theoretical_spectra `data.table` returned by the `simulate_theoretical_spectra` function.
#' @return `data.table`
#' @keywords internal
get_undeuterated_mass = function(theoretical_spectra) {
    unique(theoretical_spectra$Charge * (theoretical_spectra$Mz - 1.007276))[1]
}
#' Adds noise to mass domain
#' @param curve `data.table`
#' @param standard_deviations named vector of standard deviations
#' @return `data.table`
#' @keywords internal
add_noise = function(curve, standard_deviations) {
    sd = standard_deviations[as.character(unique(curve$Rep))]
    curve$Mass = ifelse(curve$Mass == 0, 0,
                        curve$Mass + rnorm(nrow(curve), 0, sd = sd))
    curve
}
#' Adds noise to a single time point in mass domain
#' @param spectrum data.table with a single time point
#' @param standard_deviations vector of standard deviations for the random noise
#' @return data.table
#' @keywords internal
add_noise_to_one_timepoint = function(spectrum, standard_deviations) {
    sd = standard_deviations[as.character(unique(spectrum$Exposure))]
    spectrum$Mz = spectrum$Mz + rnorm(nrow(spectrum), 0, sd)
    spectrum[, !(colnames(spectrum) == "Exposure"), with = FALSE]
}
#' Adds noise to a single time point in intensity domain
#' @inheritParams add_noise_to_one_timepoint
#' @return data.table
#' @keywords internal
add_noise_to_intensities = function(spectrum, standard_deviations) {
    sd = standard_deviations[as.character(unique(spectrum$Exposure))]
    spectrum$Intensity = spectrum$Intensity + rnorm(nrow(spectrum), 0, sd)
    spectrum[, !(colnames(spectrum) == "Exposure"), with = FALSE]
}
#' Mass minus mass at zero
#' @param mass vector of masses
#' @param time vector of time points
#' @keywords internal
get_relative_mass = function(mass, time) {
    mass - mass[time == 0]
}
#' Standardize column names and types
#' @param curves list of lists of data.tables
fix_columns_names_types = function(curves) {
    Sequence = Rep = PF = Exposure = Mass = NULL
    lapply(curves, function(curve) {
        lapply(curve, function(replicate_curve) {
            replicate_curve[, list(Sequence = as.factor(Sequence),
                                   Rep = as.factor(as.character(Rep)),
                                   State = as.factor(as.character(PF)),
                                   Exposure, Mass)]
        })
    })
}
