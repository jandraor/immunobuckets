% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DGP.R
\name{simulate_infections}
\alias{simulate_infections}
\title{Simulate infections}
\usage{
simulate_infections(
  total_buckets,
  n_e,
  avg_lambda_e,
  lambda,
  start_index,
  init_weight,
  final_weight,
  serostatus,
  is_vaccinated,
  rho_v,
  rho,
  vac_buckets,
  follow_up_times,
  enrolment_year,
  weights_fu,
  enrolment_weight,
  is_vac_prob,
  switch_rho,
  vac_year
)
}
\arguments{
\item{total_buckets}{A numeric value indicating the number of serotypes.}

\item{n_e}{A}

\item{avg_lambda_e}{A number corresponding to the average force of infection
that the individual experience since becoming susceptible to the time of
enrolment.}

\item{lambda}{A numeric vector for the yearly FOI.}

\item{start_index}{A}

\item{init_weight}{A}

\item{final_weight}{A}

\item{serostatus}{A}

\item{is_vaccinated}{A}

\item{rho_v}{A}

\item{rho}{A}

\item{vac_buckets}{A numeric value indicating the number of buckets that
a vaccine fills.}

\item{follow_up_times}{A}

\item{enrolment_year}{A}

\item{weights_fu}{A}

\item{enrolment_weight}{A numeric value between 0 and 1 representing the
time elapsed from the date of the first blood draw to the year's end
(31st December).}

\item{is_vac_prob}{A}

\item{switch_rho}{A numeric value,
which can be either 1 (allows re-susceptibility) or 0 (no re-susceptibility).
Employing the buckets metaphor, this parameter determines whether buckets
have a hole of waning immunity.}

\item{vac_year}{A number.}
}
\value{
A vector of zeros and ones, where one represents infection.
}
\description{
Simulate infections
}
\examples{
simulate_infections(total_buckets    = 4,
                    n_e              = 2,
                    avg_lambda_e     = 0.1,
                    lambda           = rep(0.1, 4),
                    start_index      = 1,
                    init_weight      = 0,
                    final_weight     = 1,
                    serostatus       = 1,
                    is_vaccinated    = 1,
                    rho_v            = 0.1,
                    rho              = 0.5,
                    vac_buckets      = 1,
                    follow_up_times  = c(3, 4),
                    enrolment_year   = 2,
                    weights_fu       = c(1, 1),
                    enrolment_weight = 0,
                    is_vac_prob      = 1,
                    switch_rho       = 1,
                    vac_year         = 15)
}
