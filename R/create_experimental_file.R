#' Create experimental file
#'
#' @description This function generates replications of mass spectra that are
#' consistent with common experimental data files
#'
#' @param peptides a data frame of sequences (\code{sequence}), \code{Protein},
#' and \code{Start}, \code{End} and parameters except \code{times}
#' that can be used for simulating mass spectra.
#' See \code{\link[powerHaDeX]{simulate_theoretical_spectra}} for more details
#' about the additional parameters.
#' @inheritParams simulate_theoretical_spectra
#' @inheritParams get_noisy_deuteration_curves
#' @param file_type ...
#'
#' @examples
#' \dontrun{
#' peptides <- data.frame(sequence = c("FPTTKTY", "LVRKDLQN"),
#'                        protection_factor = c(10, 100))
#' create_experimental_file(peptides, charge = 1:3)
#'
#' }
#'
#' @export
#'

create_experimental_file <- function(peptides,
                                     times = c(0.167, 1, 5, 25, 1440),
                                     charge,
                                     n_replicates = 3,
                                     mass_deviations = 50,
                                     intensity_deviations = NULL,
                                     file_type = "DynamX") {

    Modification <- MaxUptake <- Fragment <- Sequence <- RT <- NULL

    if (is.null(peptides[["sequence"]])) stop("You must provide at least one sequence to simulate spectra.
                                              See simulate_theoretical_spectra for more details.")
    all_params <- prepare_input_peptides(peptides)

    noisy_spectra_data <- rbindlist(lapply(1:nrow(all_params), function(ith_row) {

        print(paste("Simulation row:", ith_row))
        spectrum <- tryCatch(
            simulate_theoretical_spectra(sequence = all_params[ith_row,
                                                               "sequence"],
                                         charge = charge,
                                         protection_factor = all_params[ith_row,
                                                                        "protection_factor"],
                                         times = times,
                                         pH = all_params[ith_row,
                                                         "pH"],
                                         temperature = all_params[ith_row,
                                                                  "temperature"],
                                         n_molecules = all_params[ith_row,
                                                                  "n_molecules"],
                                         time_step_const = all_params[ith_row,
                                                                      "time_step_const"],
                                         use_markov = all_params[ith_row,
                                                                 "use_markov"]),
            error = function(e) {
                print(e)
                data.frame()})

        spectra_by_charge <- split(spectrum, f = spectrum[["Charge"]])

        rbindlist(lapply(spectra_by_charge, function(one_charge_spectrum) {

            undeuterated_mass <- get_undeuterated_mass(one_charge_spectrum)
            spectra <- get_spectra_list(one_charge_spectrum,
                                        compare_pairs = FALSE,
                                        reference = NA)
            noisy_spectra <- add_noise_to_spectra(spectra,
                                                  n_replicates = n_replicates,
                                                  n_experiments = 1,
                                                  undeuterated_mass,
                                                  mass_deviations = mass_deviations,
                                                  intensity_deviations = intensity_deviations)
            noisy_spectra[[1]][[1]][["MHP"]] <- calculate_peptide_mass(strsplit(all_params[ith_row, "sequence"], "")[[1]])
            noisy_spectra_center <- get_deuteration_curve_single_spectrum(noisy_spectra[[1]][[1]])
            noisy_spectra_center[["Experimental_state"]] <- NULL
            noisy_spectra_center[["Protein"]] <- all_params[ith_row, "Protein"]
            noisy_spectra_center[["Start"]] <- all_params[ith_row, "Start"]
            noisy_spectra_center[["End"]] <- all_params[ith_row, "End"]
            noisy_spectra_center[["PH"]] <- NULL
            colnames(noisy_spectra_center) <- c("Exposure", "z", "Sequence",
                                                "State", "File", "MHP", "Center",
                                                "Protein", "Start", "End")
            noisy_spectra_center[["State"]] <- paste0("PF_",
                                                      noisy_spectra_center[["State"]])
            noisy_spectra_center

        }))
    }))

    noisy_spectra_data[, Modification := NA]
    noisy_spectra_data[, Fragment := NA]
    noisy_spectra_data[, MaxUptake := nchar(Sequence) - 2 - lengths(regmatches(Sequence,
                                                                               gregexpr("P", Sequence)))]
    noisy_spectra_data[, RT := 0]
    DynamX_order <- c("Protein", "Start", "End", "Sequence", "Modification",
                      "Fragment", "MaxUptake", "MHP", "State", "Exposure",
                      "File", "z", "RT", "Center")
    setcolorder(noisy_spectra_data, DynamX_order)

    data.table(noisy_spectra_data)

}

#' Prepare input for \code{\link[powerHaDeX]{create_experimental_file}}
#'
#' @description Supplementary function providing appropriate input.
#' @param peptides a data frame of parameters for which
#' \code{\link[powerHaDeX]{simulate_theoretical_spectra}} will be executed.
#'
#' @return a data frame being a proper input for
#' \code{\link[powerHaDeX]{create_experimental_file}}.
#'

prepare_input_peptides <- function(peptides) {

    peptides <- add_column(peptides, "protection_factor", 1)
    peptides <- add_column(peptides, "pH",  7.5)
    peptides <- add_column(peptides, "temperature", 15)
    peptides <- add_column(peptides, "n_molecules", 100)
    peptides <- add_column(peptides, "time_step_const",  1)
    peptides <- add_column(peptides, "if_corr", 0)
    peptides <- add_column(peptides, "min_probability", 1e-4)
    peptides <- add_column(peptides, "use_markov", TRUE)
    peptides <- add_column(peptides, "Start", "")
    peptides <- add_column(peptides, "End", "")
    peptides <- add_column(peptides, "Protein", "Protein")

    peptides

}




#' Complete data frame with columns
#'
#' @description This function adds column if does not exist and fill it with
#' provided value.
#'
#' @param data a data frame of interest.
#' @param col_name a character. Name of column that should be created if it does
#' not exist
#' @param value optional. A value to fill with.
#'
#' @examples
#' \dontrun{
#' add_column(mtcars, "mpg")
#' add_column(mtcars, c("column1", "column2"))
#' }
#'

add_column <- function(data, col_name, value = NULL) {

    missing_column <- col_name[!col_name%in%names(data)]

    if(length(missing_column)!= 0){

        data[[missing_column]] <- value
    }

    data

}


