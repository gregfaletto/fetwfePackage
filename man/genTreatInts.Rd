% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/did_funcs.R
\name{genTreatInts}
\alias{genTreatInts}
\title{Generate Treatment-Covariate Interactions}
\usage{
genTreatInts(
  treat_mat_long,
  X_long,
  n_treats,
  cohort_fe,
  N,
  T,
  R,
  d,
  N_UNTREATED
)
}
\arguments{
\item{treat_mat_long}{A matrix of treatment-period dummy variables (N*T rows,
\code{n_treats} columns).}

\item{X_long}{A matrix of time-invariant covariates (N*T rows, \code{d} columns).
Each unit's \code{d} covariate values are repeated \code{T} times.}

\item{n_treats}{Integer; total number of unique treatment-period dummies
(columns in \code{treat_mat_long}).}

\item{cohort_fe}{A matrix of cohort dummy variables (N*T rows, R columns).
Used to identify units within each cohort for centering covariates.}

\item{N}{Integer; total number of unique units.}

\item{T}{Integer; total number of time periods.}

\item{R}{Integer; total number of treated cohorts.}

\item{d}{Integer; total number of covariates.}

\item{N_UNTREATED}{Integer; number of never-treated units.}
}
\value{
A matrix of treatment-covariate interactions (N\emph{T rows, d}\code{n_treats}
columns). Returns an empty matrix if \code{d=0}.
}
\description{
Creates interaction terms between treatment-period dummies (\code{treat_mat_long})
and time-invariant covariates (\code{X_long}). Covariates are centered with respect
to their cohort means before interaction.
}
\details{
If \code{d} (number of covariates) is 0, returns an empty matrix.
Otherwise, covariates \code{X_long} are first centered. For each cohort (including
the never-treated group, identified by \code{rowSums(cohort_fe) == 0}), the mean
of each covariate \emph{within that cohort} is subtracted from the covariate values
of units in that cohort.
Then, each column of \code{treat_mat_long} (a specific treatment-period dummy) is
interacted with each column of the centered covariates.
The resulting matrix \code{X_long_treat} has \code{d * n_treats} columns. Columns are
ordered such that the first \code{d} columns are interactions of \code{treat_mat_long[,1]}
with the centered covariates, the next \code{d} with \code{treat_mat_long[,2]}, and so on.
}
\keyword{internal}
