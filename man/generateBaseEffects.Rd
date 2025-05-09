% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_funcs.R
\name{generateBaseEffects}
\alias{generateBaseEffects}
\title{Generate Base Fixed Effects and Covariates for Simulation}
\usage{
generateBaseEffects(N, d, T, R, distribution = "gaussian")
}
\arguments{
\item{N}{Integer. Number of units in the panel.}

\item{d}{Integer. Number of time-invariant covariates.}

\item{T}{Integer. Number of time periods.}

\item{R}{Integer. Number of treated cohorts.}

\item{distribution}{Character. Distribution to generate covariates.
Defaults to \code{"gaussian"}. If set to \code{"uniform"}, covariates are drawn uniformly
from \eqn{[-\sqrt{3}, \sqrt{3}]}.}
}
\value{
A list containing:
\item{cohort_fe}{A matrix of cohort fixed effects (dummy variables).}
\item{time_fe}{A matrix of time fixed effects (dummy variables for periods 2 to T).}
\item{X_long}{A long-format matrix of covariates, repeated for each time period.}
\item{assignments}{A vector of counts indicating how many units fall into
the never-treated group and each of the R treated cohorts.}
\item{cohort_inds}{A list where each element contains the row indices in the
long-format matrices corresponding to the units in a specific treated cohort.}
}
\description{
Creates cohort fixed effects, time fixed effects, and a long-format covariate matrix
for simulating panel data. Covariates are drawn based on the specified distribution.
}
\keyword{internal}
