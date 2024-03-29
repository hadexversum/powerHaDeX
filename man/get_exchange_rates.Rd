% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exchange_rates.R
\name{get_exchange_rates}
\alias{get_exchange_rates}
\title{Hydrogen-deuterium or back-exchange exchange rates}
\usage{
get_exchange_rates(
  sequence,
  exchange = "HD",
  pH = 9,
  temperature = 15,
  mol_type = "poly",
  if_corr = FALSE
)
}
\arguments{
\item{sequence}{peptide amino acid sequence as a character vector of amino
acids}

\item{exchange}{type of exchange - "HD" for hydrogen to deuterium,
"DH" for deuterium to hydrogen (back-exchange). Default "HD".}

\item{pH}{reaction pH}

\item{temperature}{temperature of the reaction (Celsius). Default to 15.}

\item{mol_type}{character, "poly" or "oligo". If "oligo" then the calculated
constants are multiplied by c according to the considered condition (base,
water or acid).}

\item{if_corr}{logical. PH correction indicator. Default value \code{FALSE}.
The value of pH is equal to pD. If there is correction, the pD = pH + 0.4.
(Conelly et al 1993)}
}
\value{
The obtained exchange rates are stored in vector \code{kcHD} or
\code{kcDH} according to the exchange direction. They are used to calculate
the exchange probabilities thus both \code{kcHD} and \code{kcDH} are
necessary as we take the possibility of back-exchange into account.
}
\description{
Calculate exchange rates that are required to obtain exchange
probabilities.
}
\details{
The correction of \code{pH} is taken into account for calculation
of \code{pD}:
\deqn{pD = pH + 0.4 * if_corr}
Next, the provided temperature is converted into K and the internal functions
\code{\link[powerHaDeX]{get_F_const}},
\code{\link[powerHaDeX]{get_poly_const}} and
\code{\link[powerHaDeX]{get_pkc}} are evaluated.

Using the obtained matrix of constants and provided \code{sequence}
\code{F_a} and \code{F_b} are calculated for each amino acid in the sequence,
concerning the previous and next amino acid. For the amino acids in the
middle of the sequence, the following formula is used:

\deqn{F_x = 10^{ previous_x + current_x}}

where \code{x} is either \code{a} or \code{b}, and \code{previous_x} is the
acid/base factor for a previous amino acid in the sequence, and
\code{current_x} for the amino acid it is calculated for. If the amino acid
is next to the C- or N-term, the term-effect is taken into account.

Finally, the exchange rate \code{k_c} for the amino acid is the sum of
catalysis constants for acid, base and water (Conelly et al, 1993). Namely:

\deqn{k_c = k_{acid} + k_{base} + k_{water}}
where
\deqn{k_{acid} = F_a * K_a + D * F_ta,}
\deqn{k_{base} = F_b* K_b + OD * F_tb,}
\deqn{k_{water} = F_b * K_w * F_tw}

where \code{D} and \code{OD} indicates deuterium  and deuterium oxide
concentration, \code{F_a} and \code{F_b} are values calculated specifically
for given amino acid, as described before, \code{K_a} and \code{K_b} are
values computed by \code{\link[powerHaDeX]{get_poly_const}} function, based
on the mole type, \code{F_ta}, \code{F_tb} and \code{F_tw} are values
computed by \code{\link[powerHaDeX]{get_F_const}} function.
}
\keyword{internal}
