
#' Draw mass spectra
#' @description Graphical visualization of mass spectra obtained using the function
#' \code{\link[powerHaDeX]{simulate_theoretical_spectra}}.
#' @param spectra data table. Result of \code{\link[powerHaDeX]{simulate_theoretical_spectra}}.
#' @param time_points vector of values of exposure times to be displayed on the plot.
#' Default \code{unique(spectra[["Exposure"]])}.
#' @param charges vector of charges to be displayed on the plot. Default \code{unique(spectra[["Charge"]])}.
#' @param control_time logical. Indicates whether the spectrum at the control time
#' (conventionally equal to 0) should be drawn.
#' @param ... additional arguments passing to the \code{\link[ggplot2]{theme}}.
#' @return ggplot object
#' @details This function draws mass spectra from data obtained via
#' \code{\link[powerHaDeX]{simulate_theoretical_spectra}}.
#' @export


plot_spectra <- function(spectra,
                         time_points = unique(spectra[["Exposure"]]),
                         charges = unique(spectra[["Charge"]]),
                         control_time = FALSE,
                         ...) {
    Mz = Intensity = NULL

    if(control_time) {
        spectra <- spectra[spectra[["Exposure"]] %in% time_points &
                               spectra[["Charge"]] %in% charges]
    }else{
        spectra <- spectra[spectra[["Exposure"]] %in% time_points &
                               spectra[["Charge"]] %in% charges &
                               spectra[["Exposure"]] != 0]
    }
    spectra[["Exposure"]] <- paste0("Exposure: ", spectra[["Exposure"]])
    spectra[["Charge"]] <- paste0("Charge: ", spectra[["Charge"]])

    p <- ggplot(spectra, aes(x = Mz, xend = Mz, yend = Intensity, y = 0)) +
        geom_segment() +
        ylab("Intensity") +
        xlab("m/z") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5), ...) +
        ggtitle(paste("Sequence:", unique(spectra[["Sequence"]]))) +
        facet_wrap(Charge~ Exposure, scales = "free")

    return(p)
}



