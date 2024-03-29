% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exchange_rates.R
\name{get_pkc}
\alias{get_pkc}
\title{Calculating pKc values}
\usage{
get_pkc(temp_kelvin, gas_constant, exchange = "HD")
}
\arguments{
\item{temp_kelvin}{temperature reaction in Kelvins}

\item{gas_constant}{gas constant (1/(dT * R) = 1.9858775)}

\item{exchange}{type of exchange - "HD" for hydrogen to deuterium,
"DH" for deuterium to hydrogen (back-exchange). Default "HD".}
}
\value{
The function returns a list of \code{asp}, \code{glu} and \code{his}
(\code{pKc} values corresponding to amino acids).
}
\description{
calculates supplementary constants for aspartic acid (Asp),
glutamic acid (Glu) and histidine (His). Values for mentioned amino acids are
pH and temperature dependent, in contrary to the rest amino acids with fixed
values.
}
\details{
Depending on provided \code{exchange} direction tabular values of
exponents E_{const} are assigned. For \code{Asp}, \code{Glu} and \code{His}
the \code{pKc} constants are calculated based on the energies of activation
for given amino acid and the chosen \code{exchange} direction.
}
\keyword{internal}
