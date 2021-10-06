#' Simulate theoretical spectra of a deuterated peptide over time
#'
#' @param sequence amino acid sequence of a peptide as a single string
#' @param charge vector of charges of the peptide ion. If NULL, one value is sampled
#' from vector 2:6. Default NULL.
#' @param protection_factor protection factor. If a single number of provided,
#' same protection factor will be assumed for each amide. Default value: 1
#' (indicates that the exchange rate is equal to the intristic exchange rate)
#' @param times a vector of times at which deuteration levels will be measured (seconds)
#' @param pH pH of the reaction. Default to 7.5.
#' @param temperature temperature of the reaction (Celsius)
#' @param n_molecules number of peptide molecules. Default to 100.
#' @param time_step_const time step constant. Default value: $1$. Value that
#' indicates the length of the time step of the simulation. The bigger the time
#' step, the fewer time points are simulated (the fewer iterations in case of
#' Zhong-Yuan Kan's approach).
#' @param min_probability smallest isotopic probability to consider
#' @param use_markov logical. If TRUE algorithm basing on Markov chain will be used.
#' If FALSE simulation provided by Zhong-Yuan Kan will be executed. Default to TRUE,
#' as it fastens the calculation
#' @inheritParams get_exchange_rates
#'
#' @details To the results calculated by \code{\link[powerHaDeX]{get_iso_probs_deut}}
#' is added a minimal exchange control - for time point \code{0} (directly after
#' adding a buffer). The m/z values are obtained as a ratio of the
#' \code{peptide_mass} magnified by proton mass and the peptide charge.
#' The distribution of undeuterated peptide is the
#' intensities vector.
#'
#' @return a data table of variables:
#'
#' - \code{Exposure} - time point of a measurement,
#'
#'  - \code{Mz} - mass-to-charge ratio,
#'
#'  - \code{Intensity} - isotopic probabilities larger than \code{min_probability}
#'  (the smaller ones are zeroes)
#'
#'  and the variables provided by user
#'
#'  - \code{Sequence},
#'
#'  - \code{PF},
#'
#'  - \code{Charge},
#'
#'  - \code{PH}.
#'
#' @seealso The algorithm that is used to simulate theoretical spectra is based on
#' Zhong-Yuan Kan's implementation in Matlab. The original version of codes is
#' located in the repository \url{https://github.com/kanzy/HX-MS-Simulations}
#' (as at 29.06.2020). In the \code{powerHaDeX} package can be found the Kan's
#' algorithm re-implemented in R (using Rcpp) and the accelerated implementation
#' (that uses Markov chains' properties). Moreover, the package \code{powerHaDeX}
#' allows the user to simulate spectra for more than one exposure time for both
#' (Rcpp and Markov) approaches.
#'
#' @examples
#' simulate_theoretical_spectra(sequence = "LVRKDLQN",
#'                              charge = c(3, 5),
#'                              protection_factor = 100,
#'                              times = c(0.167, 5),
#'                              pH = 7.5,
#'                              temperature = 15,
#'                              n_molecules = 500,
#'                              time_step_const = 1,
#'                              use_markov = TRUE)
#'
#' @export
simulate_theoretical_spectra <- function(sequence, charge = NULL, protection_factor = 1,
                                         times = c(60, 600), pH = 7.5,
                                         temperature = 15, n_molecules = 100,
                                         time_step_const = 1, if_corr = FALSE,
                                         min_probability = 1e-4,
                                         use_markov = TRUE) {

    sequence <- strsplit(sequence, "")[[1]]

    if (length(protection_factor) == 1L) {

        protection_factor <- rep(protection_factor, length(sequence))
    }
    if (is.null(charge)) {

        charge <- sample(2:6, 1)
    }

    peptide_iso_dist <- get_approx_isotopic_distribution(sequence, min_probability)
    peptide_mass <- peptide_iso_dist[[1]]
    isotopic_probs <- peptide_iso_dist[[2]]
    maxND <- peptide_iso_dist[[3]]
    maxD <- peptide_iso_dist[[4]]
    kcHD <- get_exchange_rates(sequence, "HD", pH, temperature, 'poly', if_corr)
    kcDH <- get_exchange_rates(sequence, "DH", pH, temperature, 'poly')

    kmax <- max(max(kcDH), max(kcHD))
    deltaT <- time_step_const / kmax
    steps_between_time_points <- ceiling(times/deltaT)

    if (floor(max(times)/deltaT) == 0) {

        print("There is no deuteration before given time point. The measurement at the control time (conventionally 0) is returned.")
        isotope_dists <- data.table::data.table()
    } else {
        tryCatch({
            transition_probs <- get_exchange_probabilities(kcHD, kcDH, deltaT, protection_factor)

            if(use_markov) {
                HD_matrices <- get_HD_matrices_using_markov(sequence, transition_probs,
                                                            steps_between_time_points,
                                                            n_molecules)
            }else {
                time_sequence <- seq(0, max(times), deltaT)
                HD_matrices <- get_HD_matrices(sequence, transition_probs,
                                               time_sequence, times,
                                               n_molecules)
            }

            isotope_dists <- get_iso_probs_deut(HD_matrices, maxD, maxND,
                                                isotopic_probs, peptide_mass,
                                                times, charge, pH)
        })
    }

    isotope_dists <- rbind(merge(data.frame(Exposure = 0,
                                            Intensity = isotopic_probs,
                                            PH = pH),
                                 data.frame(Mz = peptide_mass / charge + 1.007276,
                                            PH = pH,
                                            Exposure = 0,
                                            Charge = charge)),
                           isotope_dists)

    isotope_dists[["Sequence"]] <- paste0(sequence, collapse = "")

    if (length(unique(protection_factor)) == 1) {
        isotope_dists[["PF"]] <- protection_factor[1]
    } else {
        isotope_dists[["PF"]] <- paste(protection_factor,
                                       sep = ",", collapse = ",")
    }

    isotope_dists <- isotope_dists[isotope_dists[["Intensity"]] > min_probability, ]

    data.table::as.data.table(isotope_dists)
}
