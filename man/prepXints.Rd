% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/did_funcs.R
\name{prepXints}
\alias{prepXints}
\title{Prepare Design Matrix and Initial Data for FETWFE}
\usage{
prepXints(data, time_var, unit_var, treatment, covs, response, verbose = FALSE)
}
\arguments{
\item{data}{A \code{data.frame} object representing the panel data. Each row
should be an observation of a unit at a specific time.}

\item{time_var}{A character string specifying the name of the column in \code{data}
that contains the time period variable (e.g., year). Expected to be integer.}

\item{unit_var}{A character string specifying the name of the column in \code{data}
that contains the unit identifier. Expected to be character.}

\item{treatment}{A character string specifying the name of the column in \code{data}
that contains the binary treatment indicator (0 for untreated, 1 for treated).
Treatment is assumed to be an absorbing state.}

\item{covs}{A character vector specifying the names of the columns in \code{data}
to be used as covariates. These are treated as time-invariant (values from
the first period are used).}

\item{response}{A character string specifying the name of the column in \code{data}
that contains the response variable.}

\item{verbose}{Logical. If \code{TRUE}, messages about the data processing steps
will be printed. Default is \code{FALSE}.}
}
\value{
A list containing:
\item{X_ints}{The fully constructed design matrix with all fixed effects,
covariates, treatment dummies, and their interactions.}
\item{y}{The centered response vector.}
\item{N}{The final number of unique units after processing.}
\item{T}{The number of unique time periods.}
\item{d}{The final number of covariates after processing.}
\item{p}{The total number of columns in \code{X_ints} (total parameters).}
\item{in_sample_counts}{An integer vector named with cohort identifiers
(including "Never_treated"), indicating the number of units in each cohort
within the provided \code{data}.}
\item{num_treats}{The total number of unique treatment effect parameters
(e.g., \eqn{\tau_{rt}} from the paper).}
\item{first_inds}{A numeric vector indicating the starting column index within
the block of treatment effect parameters in \code{X_ints} (and subsequently in
\code{beta_hat}) for each respective cohort's first treatment period.}
}
\description{
This function serves as a pre-processing step for \code{fetwfe}. It takes the raw
panel data and transforms it into a design matrix (\code{X_ints}) and a response
vector (\code{y}) suitable for the core estimation logic. It identifies cohorts,
processes covariates (making them time-invariant and handling missing values),
and generates dummy variables for cohorts, time periods, and treatments.
}
\details{
The function performs several key steps:
\enumerate{
\item Calls \code{idCohorts} to identify treatment cohorts, unique units, and time
periods. Units treated in the very first time period are removed.
\item Calls \code{processCovs} to handle covariates. Time-varying covariates are
replaced by their value in the first (pre-treatment) period. Covariates
with missing values in the first period or covariates that are constant
across all units are removed.
\item Calls \code{addDummies} to create dummy variables for cohorts, time periods
(excluding the first), and treatment indicators (for each cohort and
each post-treatment time period). The response variable is also centered.
\item Calls \code{genXintsData} to construct the final design matrix \code{X_ints} by
combining the dummy variables, covariates, and their interactions.
The interactions include cohort-covariate, time-covariate, and
treatment-covariate interactions.
}
Input \code{covs} are expected to be numeric or integer after factor processing.
The function also calculates the number of units (\code{N}), time periods (\code{T}),
covariates (\code{d}), the total number of parameters in the full design matrix (\code{p}),
in-sample cohort counts, the number of unique treatment terms (\code{num_treats}),
and indices of the first treatment effect for each cohort (\code{first_inds}).
}
\keyword{internal}
