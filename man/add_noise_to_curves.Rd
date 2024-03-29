% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noisy_spectra.R
\name{add_noise_to_curves}
\alias{add_noise_to_curves}
\title{Adds noise to deuteration curves}
\usage{
add_noise_to_curves(curves, per_run_deviations = NULL, relative = TRUE)
}
\arguments{
\item{curves}{list of lists of deuteration curves returned by
\code{\link[powerHaDeX]{get_deuteration_curves_from_spectra}}.}

\item{per_run_deviations}{optional, standard deviations of random noise that
will be added to deuteration curves. Either a single number (then the error
at each time point will be the same) or a vector of the same length as number
of unique time points in the experiment. The error will be sampled from
normal distribution with these standard deviations. Default \code{NULL}.}

\item{relative}{logical, if TRUE (default), each deuteration curve will start
at 0 (relative mass will be returned). Default \code{TRUE}.}
}
\value{
list of lists of data.tables
}
\description{
This function makes noisy deuteration curves for power
calculation purposes.
}
\details{
This function uses
\code{\link[powerHaDeX]{add_noise_to_single_curve}}
}
\keyword{internal}
