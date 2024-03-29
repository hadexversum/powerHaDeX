% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noisy_spectra.R
\name{add_noise_to_single_curve}
\alias{add_noise_to_single_curve}
\title{Adds noise to a single deuteration curve}
\usage{
add_noise_to_single_curve(
  replicate_curve,
  per_run_deviations = NULL,
  relative = TRUE
)
}
\arguments{
\item{replicate_curve}{`data.table` with technical replicates of a
deuteration curve (one element of a list obtained by
\code{\link[powerHaDeX]{get_deuteration_curves_from_spectra}})}

\item{per_run_deviations}{optional, standard deviations of random noise that
will be added to deuteration curves. Either a single number (then the error
at each time point will be the same) or a vector of the same length as number
of unique time points in the experiment. The error will be sampled from
normal distribution with these standard deviations. Default \code{NULL}.}

\item{relative}{logical, if TRUE (default), each deuteration curve will start
at 0 (relative mass will be returned). Default \code{TRUE}.}
}
\value{
`data.table` containing noisy replicated deuteration curves.
}
\description{
Adds noise to a single deuteration curve
}
\details{
This function uses \code{\link[powerHaDeX]{add_noise}} and
\code{\link[powerHaDeX]{get_relative_mass}}.
}
\keyword{internal}
