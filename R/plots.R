
#' Draw mass spectrum
#' @param spectrum data table. Result of \code{\link[powerHDX]{simulate_theoretical_spectra}}
#' @param time_point a value of exposure for which plot will be provided. If \code{NULL}
#' exposure time will be taken from \code{spectrum} but only when there will be one value.
#' Default to \code{NULL}.
#' @param ... additional arguments passing to the \code{\link[ggplot2]{theme}}.
#' @return ggplot object
#' @details This function draws mass spectrum from data obtained via
#' \code{\link[powerHDX]{simulate_theoretical_spectra}}.
#' @examples
#' \dontrun{
#' spectrum <-   simulate_theoretical_spectra(sequence = "AAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
#'               charge = 3,
#'               protection_factor = 100,
#'               times = 60,
#'               pH = 6,
#'               temperature = 15,
#'               n_molecules = 500,
#'               time_step_const = 1,
#'               use_markov = TRUE)
#' #one exposure in data so time_point can be NULL
#' plot_single_spectrum(spectrum)
#'
#' spectrum <-   simulate_theoretical_spectra(sequence = "AAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
#'               charge = 3,
#'               protection_factor = 100,
#'               times = c(60, 600),
#'               pH = 6,
#'               temperature = 15,
#'               n_molecules = 500,
#'               time_step_const = 1,
#'               use_markov = TRUE)
#' #in case of two time points one must be chosen
#' plot_single_spectrum(spectrum, 60)
#' }
#'
#' @export


plot_single_spectrum <- function(spectrum, time_point = NULL, ...) {
    Mz = Intensity = NULL

    times_from_data <- setdiff(unique(spectrum[["Exposure"]]), 0)

    if(is.null(time_point)) {
        if (length(times_from_data) > 1) stop("One time point must be chosen.")
        time_point <- times_from_data
    }
    match.arg(time_point, times_from_data)
    spectrum_time_point <- spectrum[spectrum[["Exposure"]] == time_point, ]

    p <- ggplot(spectrum_time_point, aes(x = Mz, xend = Mz, yend = Intensity, y = 0)) +
        geom_segment() +
        ylab("Intensity") +
        xlab("m/z") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5), ...) +
        ggtitle(paste("Sequence:", unique(spectrum_time_point[["Sequence"]]), ", Exposure:" , time_point, "sec"))

    return(p)
}

