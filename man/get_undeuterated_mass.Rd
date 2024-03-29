% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noisy_spectra.R
\name{get_undeuterated_mass}
\alias{get_undeuterated_mass}
\title{Mass of an undeuterated peptide}
\usage{
get_undeuterated_mass(theoretical_spectra)
}
\arguments{
\item{theoretical_spectra}{a data table or a list of data tables of
theoretical spectra created by the function
\code{\link[powerHaDeX]{simulate_theoretical_spectra}}.}
}
\value{
`data.table` of calculated mass value for the first peak (the
smallest one) as it is usually the peak corresponding to the monoisotopic
mass.
}
\description{
This function gets mass of an undeuterated peptide based on its
spectrum.
}
\details{
For the calculations the formula below is used
\deqn{undeuterated_mass = charge * Mz - p_m}
where \code{Mz} is mass-to-charge ratio for the peaks from the provided
theoretical spectrum and p_m is the mass of proton equal to 1.007276.
}
\keyword{internal}
