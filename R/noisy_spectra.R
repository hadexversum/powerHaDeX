#' Replicated deuterium uptake curves
#'
#' @description This function creates a list of lists of noisy deuteration curves
#' based on theoretical spectra in order to imitate the data from the HDX experiments.
#'
#' @param theoretical_spectra a data table or a list of data tables of theoretical
#' spectra created by the function \code{\link[powerHaDeX]{simulate_theoretical_spectra}}.
#' @param compare_pairs if FALSE, all groups (defined by the protection factor)
#' will be considered jointly. If TRUE (default), each protection factor will be
#' considered together with the protection factor given by the `reference` parameter
#' @param reference protection factor that will be used for comparison to other
#' protection factors in. The function accepts either \code{NA} (for comparing all
#' protection factors), a number (for comparing with reference value of protection factor)
#' or "all" (for pairwise comparisons of all the possible combinations). Default \code{NA}.
#' @param n_replicates number of technical replicates to create
#' @param n_experiments number of replicates of an experiment for power calculation
#' @param mass_deviations mass deviation in parts per million. Either a single number
#' (then the error at each time point will be the same) or a vector of the same length
#' as number of unique time points in the experiment. The error will be sampled from
#' normal distribution with standard deviation equal to
#' \deqn{mass_deviations * undeuterated_mass/1e6}
#' Default to 50.
#' @param intensity_deviations optional, standard deviations of random noise that
#' will be added to intensities. Either a single number (then the error at each
#' time point will be the same) or a vector of the same length as number of unique
#' time points in the experiment. The error will be sampled from normal distribution
#' with these standard deviations.Default \code{NULL}.
#' @param per_run_deviations optional, standard deviations of random noise that
#' will be added to deuteration curves. Either a single number (then the error at each
#' time point will be the same) or a vector of the same length as number of unique
#' time points in the experiment. The error will be sampled from normal distribution
#' with these standard deviations. Default \code{NULL}.
#' @param relative logical, if TRUE (default), each deuteration curve will start at 0
#' (relative mass will be returned). Default \code{TRUE}.
#'
#' @return a list (for paired states when \code{compare_pairs} is \code{TRUE}) of lists
#' (repetitions of experiment for power calculations) of data tables of the variables:
#'
#' - \code{Sequence} - provided amino acid sequence
#'
#' - \code{Rep} - technical replication
#'
#' - \code{State} - provided protection factor (the theoretical - in practice unknown -
#' state of the protein)
#'
#' - \code{Exposure} - exposure time
#'
#' - \code{Mass} - mass or deuterium uptake when \code{relative} is \code{TRUE}.
#'
#' - \code{Charge} - charge
#'
#' - \code{Experimental_state} - the biological state (from the viewpoint of the
#' experimenter) provided in the case when \code{compare_pairs} is \code{TRUE}.
#'
#' @examples
#' theo_spectra_pf_100 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 100,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#'
#' theo_spectra_pf_200 <- simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                                                     charge = c(3, 5),
#'                                                     protection_factor = 200,
#'                                                     times = c(0.167, 5),
#'                                                     pH = 7.5,
#'                                                     temperature = 15,
#'                                                     n_molecules = 500,
#'                                                     time_step_const = 1,
#'                                                     use_markov = TRUE)
#'
#' theo_spectra_two_states <- rbind(theo_spectra_pf_100, theo_spectra_pf_200)
#'
#' deuteration_curves_paired_states <- get_noisy_deuteration_curves(theo_spectra_two_states,
#'                                                                  n_replicates = 4,
#'                                                                  n_experiments = 2,
#'                                                                  compare_pairs = TRUE,
#'                                                                  reference = "all")
#'
#' @export
#'
get_noisy_deuteration_curves <- function(theoretical_spectra,
                                        compare_pairs = TRUE,
                                        reference = NA,
                                        n_replicates = 4,
                                        n_experiments = 100,
                                        mass_deviations = 50,
                                        intensity_deviations = NULL,
                                        per_run_deviations = NULL,
                                        relative = TRUE) {

    undeuterated_mass <- get_undeuterated_mass(theoretical_spectra)
    spectra <- get_spectra_list(theoretical_spectra, compare_pairs, reference)
    spectra <- add_noise_to_spectra(spectra, n_replicates, n_experiments, undeuterated_mass,
                                   mass_deviations, intensity_deviations)
    curves <- get_deuteration_curves_from_spectra(spectra)
    curves <- add_noise_to_curves(curves, per_run_deviations, relative)
    curves <- fix_columns_names_types(curves)

    curves

}


#' Get a list of spectra
#'
#' @description Create a list of data tables of spectra for all states jointly
#'  or paired states.
#'
#' @inheritParams get_noisy_deuteration_curves
#'
#' @details If the parameter \code{compare_pairs} is \code{FALSE} then all the
#' provided protection factors will be considered jointly. If \code{compare_pairs} is
#' \code{TRUE}, then the parameter \code{reference} is necessary (a single number or
#' \code{"all"}). Then the data is split via the supplementary function
#' \code{\link[powerHaDeX]{get_paired_spectra}} into data tables of spectra with
#' paired biological states (the reference protection factor and the protection
#' factor of interest if provided, or all the possible pairs if \code{reference}
#' equals "all").
#'
#' @return list of data.tables containing spectra - for paired states or all states.
#'
#' @export
get_spectra_list <- function(theoretical_spectra,
                            compare_pairs = FALSE,
                            reference = NA) {

    if (compare_pairs & (is.na(reference))) stop("With pair comparisons, reference protection factor must be provided")

    if (compare_pairs) {

        spectra_list <- get_paired_spectra(theoretical_spectra, reference)
    } else {

        spectra_list <- list(theoretical_spectra)
    }

    unname(spectra_list)
}
#' Get list of data.tables of spectra for pairs of protection factors
#'
#' @description A supplementary function for \code{\link[powerHaDeX]{get_spectra_list}}.
#'
#' @inheritParams get_noisy_deuteration_curves
#'
#' @details For more details see \code{\link[powerHaDeX]{get_spectra_list}}.
#'
#' @return list of data.tables
#'
#' @keywords internal
get_paired_spectra <- function(theoretical_spectra,
                              reference = NA) {

    if(reference == "all") {

        states <- unique(theoretical_spectra[["PF"]])
        pairs <- rbind(t(combn(states, 2)), cbind(states, states))
        split_spectra <- split(theoretical_spectra, theoretical_spectra[["PF"]])

        lapply(1:nrow(pairs), function(i) {

            spectrum <- split_spectra[[as.character(pairs[i, 1])]]
            reference_spectrum <- split_spectra[[as.character(pairs[i, 2])]]

            spectrum[["Experimental_state"]] <- "A"
            reference_spectrum[["Experimental_state"]] <- "B"
            paired_spectra <- rbind(spectrum, reference_spectrum)
        })

    }else {

        if (!(as.integer(reference) %in% as.integer(unique(theoretical_spectra[["PF"]])))) {
            stop("Reference protection factor does not fit the data.")
        }

        PF <- NULL
        reference_spectrum <- theoretical_spectra[abs(PF - reference) < 1e-9, ]
        theoretical_spectra <- split(theoretical_spectra[abs(PF - reference) > 1e-9, ],
                                    theoretical_spectra[abs(PF - reference) > 1e-9, "PF"])
        theoretical_spectra <- lapply(theoretical_spectra,
                                     function(spectrum) {
                                         spectrum[["Experimental_state"]] <- "A"
                                         reference_spectrum[["Experimental_state"]] <- "B"
                                         rbind(spectrum, reference_spectrum)
                                     })

        unname(theoretical_spectra)
    }
}

# Add noise to spectra in mass and intensity domains ----

#' Create a list of spectra with replicates from a list of spectra
#'
#' @description This function is used to prepare technical replicates from the
#' provided spectra. It is a supplementary function used in
#' \code{\link[powerHaDeX]{add_noise_to_spectra}}.
#'
#' @param spectra list created by the \code{\link[powerHaDeX]{get_spectra_list}}
#'  function
#' @inheritParams get_noisy_deuteration_curves
#'
#' @return list of data.tables
#'
#' @keywords internal
make_experimental_design <- function(spectra,
                                    n_replicates = 4) {

    lapply(spectra, function(spectrum) {

        data.table::rbindlist(lapply(1:n_replicates, function(i) {

            spectrum[["Rep"]] <- as.character(i)
            spectrum

        }))
    })
}
#' Creates spectra with technical replicates and noise
#'
#' @inheritParams get_noisy_deuteration_curves
#' @inheritParams make_experimental_design
#' @param undeuterated_mass the value of mass of an undeuterated peptide computed by
#' \code{\link[powerHaDeX]{get_undeuterated_mass}}.
#'
#' @return list of lists of data.tables
#'
#' @details This function uses \code{\link[powerHaDeX]{make_experimental_design}} and
#' \code{\link[powerHaDeX]{make_noisy_spectra}}.
#'
#' @keywords internal
#'
#' @export
#'
add_noise_to_spectra <- function(spectra,
                                n_replicates = 4,
                                n_experiments = 100,
                                undeuterated_mass,
                                mass_deviations = 50,
                                intensity_deviations = NULL) {

    spectra <- make_experimental_design(spectra, n_replicates)
    make_noisy_spectra(spectra, n_experiments, undeuterated_mass,
                       mass_deviations, intensity_deviations)

}


#' Creates spectra with technical replicates and noise
#'
#' @param spectra a list of spectra obtained via \code{\link[powerHaDeX]{make_experimental_design}}.
#' @inheritParams add_noise_to_spectra
#'
#' @return list of lists of data.tables
#'
#' @details This function uses  \code{\link[powerHaDeX]{add_noise_to_one_spectrum}}.
#'
#' @keywords internal
#'
make_noisy_spectra <- function(spectra,
                              n_experiments = 100,
                              undeuterated_mass,
                              mass_deviations = 50,
                              intensity_deviations = NULL) {

    lapply(spectra, function(spectrum) {

        lapply(1:n_experiments, function(ith_replicate) {

            add_noise_to_one_spectrum(spectrum, undeuterated_mass,
                                      mass_deviations, intensity_deviations)

        })
    })
}


#' Add noise to a single spectrum
#'
#' @description This function adds noise to intensities and/or masses in a spectrum.
#'
#' @param spectrum a single spectrum from the list obtained using
#'  \code{\link[powerHaDeX]{make_experimental_design}}.
#' @inheritParams add_noise_to_spectra
#'
#' @details The deviations are calculated as follows
#' \deqn{sd = mass_deviations * undeuterated_mass/10^6}.
#' To add noise this function uses \code{\link[powerHaDeX]{add_noise_to_one_timepoint}}
#' and \code{\link[powerHaDeX]{add_noise_to_intensities}}.
#'
#' @return data.table containing a single noisy spectrum
#'
#' @keywords internal
#'
add_noise_to_one_spectrum <- function(spectrum,
                                     undeuterated_mass,
                                     mass_deviations = 50,
                                     intensity_deviations = NULL) {

    if (length(mass_deviations) == 1) {
        mass_deviations <- rep(mass_deviations, data.table::uniqueN(spectrum[["Exposure"]]))
    }

    sds <- mass_deviations * undeuterated_mass / 1e6
    names(sds) <- as.character(unique(spectrum[["Exposure"]]))

    spectrum <- spectrum[, add_noise_to_one_timepoint(.SD, sds), by <- "Exposure",
                        .SDcols = colnames(spectrum)]

    if (!is.null(intensity_deviations)) {

        if (length(intensity_deviations) == 1) {

            intensity_deviations <- rep(intensity_deviations,
                                       data.table::uniqueN(spectrum[["Exposure"]]))
        }

        names(intensity_deviations) <- names(sds)
        spectrum <- spectrum[, add_noise_to_intensities(.SD, sds), by = "Exposure"]
    }

    spectrum

}
# Calculate deuteration curves ----

#' Calculate deuteration curves
#'
#' @description Create a list of deuteration curves from a list of spectra
#'
#' @param spectra list created by \code{\link[powerHaDeX]{get_spectra_list}}
#'  and \code{\link[powerHaDeX]{add_noise_to_spectra}} functions.
#'
#' @return list of lists of data.tables
#'
#' @details From each spectrum provided in the list \code{spectra} the deuterium
#' uptake will be calculated. This function uses
#' \code{\link[powerHaDeX]{get_deuteration_curve_single_spectrum}}.
#'
#' @keywords internal
#'
#' @export
#'
get_deuteration_curves_from_spectra <- function(spectra) {

    lapply(spectra, function(spectrum) {

        lapply(spectrum, function(replicate) {

            get_deuteration_curve_single_spectrum(replicate)

        })
    })
}
#' Get a deuteration curve from a single spectrum
#'
#' @description This funtion calculates deuterium uptake from spectra.
#'
#' @param spectrum data.table with a spectrum being an element of list created
#' by \code{\link[powerHaDeX]{get_spectra_list}} and
#' \code{\link[powerHaDeX]{add_noise_to_spectra}} functions.
#'
#' @return data.table containing one deuterium uptake curve.
#'
#' @details The centroid mass value from spectrum is calculated as a weighted mean
#' from peaks based on the formulas
#' \deqn{M = charge * (Mz - p_m)}
#' and
#' \deqn{m = 1/N \sum_{k = 1}^{N} I_k M_k.}
#'
#' @keywords internal
#'
#' @export
#'
get_deuteration_curve_single_spectrum <- function(spectrum) {

    Mz <- Intensity <- Charge <- NULL

    grouping_columns <- setdiff(colnames(spectrum), c("Mz", "Intensity"))
    spectrum <- spectrum[, list(Mass = weighted.mean(Charge * (Mz - 1.007276),
                                                    Intensity)),
                        by = grouping_columns]

    spectrum

}
# Noise for deuteration curves ----

#' Adds noise to deuteration curves
#'
#' @description This function makes noisy deuteration curves for power calculation purposes.
#'
#' @param curves list of lists of deuteration curves returned by
#' \code{\link[powerHaDeX]{get_deuteration_curves_from_spectra}}.
#' @inheritParams get_noisy_deuteration_curves
#'
#' @return list of lists of data.tables
#'
#' @details This function uses \code{\link[powerHaDeX]{add_noise_to_single_curve}}
#'
#' @keywords internal
#'
add_noise_to_curves <- function(curves,
                               per_run_deviations = NULL,
                               relative = TRUE) {

    lapply(curves, function(curve) {

        lapply(curve, function(replicate_curve) {

            add_noise_to_single_curve(replicate_curve, per_run_deviations, relative)

        })
    })
}
#' Adds noise to a single deuteration curve
#'
#' @param replicate_curve `data.table` with technical replicates of a deuteration curve
#' (one element of a list obtained by \code{\link[powerHaDeX]{get_deuteration_curves_from_spectra}})
#' @inheritParams get_noisy_deuteration_curves
#'
#' @return `data.table` containing noisy replicated deuteration curves.
#'
#' @details This function uses \code{\link[powerHaDeX]{add_noise}} and
#' \code{\link[powerHaDeX]{get_relative_mass}}.
#'
#' @keywords internal
#'
#' @export
#'
add_noise_to_single_curve <- function(replicate_curve,
                                     per_run_deviations = NULL,
                                     relative = TRUE) {
    Exposure <- Mass <- Charge <- NULL

    if (!is.null(per_run_deviations)) {

        if (length(per_run_deviations) == 1) {

            per_run_deviations <- rep(per_run_deviations, length(unique(replicate_curve[["Rep"]])))
        }
        names(per_run_deviations) <- as.character(unique(replicate_curve[["Rep"]]))

        replicate_curve <- replicate_curve[, add_noise(.SD, per_run_deviations),
                                          by = "Rep", .SDcols = colnames(replicate_curve)]
        replicate_curve <- replicate_curve[, -1, with = FALSE]
    }

    if (relative) {

        grouping_columns <- setdiff(colnames(replicate_curve),
                                   c("Mass", "Charge", "Exposure"))
        replicate_curve <- replicate_curve[, list(Exposure, Charge, Mass = get_relative_mass(Mass, Exposure)),
                                          by = grouping_columns]
    }

    replicate_curve

}

# Smaller utility functions ----

#' Mass of an undeuterated peptide
#'
#' @description This function gets mass of an undeuterated peptide based on its spectrum.
#'
#' @inheritParams get_noisy_deuteration_curves
#'
#' @details For the calculations the formula below is used
#' \deqn{undeuterated_mass = charge * Mz - p_m}
#' where \code{Mz} is mass-to-charge ratio for the peaks from the provided theoretical
#' spectrum and p_m is the mass of proton equal to 1.007276.
#'
#' @return `data.table` of calculated mass value for the first peak (the smallest one)
#' as it is usually the peak corresponding to the monoisotopic mass.
#'
#' @keywords internal
#'
#' @export
#'
get_undeuterated_mass <- function(theoretical_spectra) {

    unique(theoretical_spectra[["Charge"]] * (theoretical_spectra[["Mz"]] - 1.007276))[1]

}

#' Adds noise to mass domain
#'
#' @inheritParams add_noise_to_single_curve
#' @param standard_deviations named vector of standard deviations
#'
#' @return `data.table` containing a single noisy curve (for one replication).
#'
#' @details  The noise is sampled from normal distribution with mean \code{0}
#' and standard deviation equal to \code{per_run_deviations} and added to the
#' \code{Mass} values unless they are not zeroes (there is no noise at the time \code{0}).
#'
#' @keywords internal
#'
add_noise <- function(replicate_curve, standard_deviations) {

    sd <- standard_deviations[as.character(unique(replicate_curve[["Rep"]]))]
    replicate_curve[["Mass"]] <- ifelse(replicate_curve[["Mass"]] == 0, 0,
                                       replicate_curve[["Mass"]] + rnorm(nrow(replicate_curve), 0, sd = sd))
    replicate_curve
}

#' Adds noise to a single time point in mass domain
#'
#' @inheritParams add_noise_to_one_spectrum
#' @param standard_deviations vector of standard deviations for the random noise
#'
#' @return data.table
#'
#' @details the noise is sampled from a normal distribution with mean \code{0}
#' and standard deviation equal to \code{standard_deviations} and added to \code{Mz}
#' values for the time points (based on supplied parameters for time points).
#'
#' @keywords internal
#'
add_noise_to_one_timepoint <- function(spectrum, standard_deviations) {

    sd <- standard_deviations[as.character(unique(spectrum[["Exposure"]]))]
    spectrum[["Mz"]] <- spectrum[["Mz"]] + rnorm(nrow(spectrum), 0, sd)
    spectrum[, !(colnames(spectrum) == "Exposure"), with = FALSE]

}

#' Adds noise to a single time point in intensity domain
#'
#' @inheritParams add_noise_to_one_timepoint
#'
#' @return data.table
#'
#' @details if the \code{intensity_deviations} were provided, then noise is sampled
#' from a normal distribution with mean \code{0} and standard deviation equal to
#' those deviations and added to \code{Intensity}.
#'
#' @keywords internal
#'
add_noise_to_intensities <- function(spectrum, standard_deviations) {

    spectrum[["Intensity"]] <- spectrum[["Intensity"]] + rnorm(nrow(spectrum), 0, standard_deviations)
    spectrum[, !(colnames(spectrum) == "Exposure"), with = FALSE]

}
#' Get relative mass
#'
#' @description Mass minus mass at zero is calculated as a form of deuterium uptake.
#'
#' @param mass vector of masses
#' @param time vector of time points
#'
#' @keywords internal
get_relative_mass <- function(mass, time) {

    mass - mass[time == 0]

}

#' Standardize column names and types
#'
#' @param curves list of lists of data.tables
#'
fix_columns_names_types <- function(curves) {

    Sequence <- Rep <- PF <- Exposure <- Mass <- Charge <- Experimental_state <-  NULL

    lapply(curves, function(curve) {

        lapply(curve, function(replicate_curve) {

            replicate_curve[, list(Sequence = as.factor(Sequence),
                                   Rep = as.factor(as.character(Rep)),
                                   State = as.factor(as.character(PF)),
                                   Exposure, Mass, Charge, Experimental_state)]

        })
    })
}
