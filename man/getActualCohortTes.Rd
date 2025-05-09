% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_funcs.R
\name{getActualCohortTes}
\alias{getActualCohortTes}
\title{Calculate True Cohort Average Treatment Effects from Coefficients}
\usage{
getActualCohortTes(R, first_inds, treat_inds, coefs, num_treats)
}
\arguments{
\item{R}{Integer. Number of treated cohorts.}

\item{first_inds}{Integer vector. \code{first_inds[r]} is the index (within the
block of treatment effect coefficients) of the first treatment effect for cohort \code{r}.}

\item{treat_inds}{Integer vector. Indices in the full \code{coefs} vector that
correspond to the block of all treatment effect coefficients.}

\item{coefs}{Numeric vector. The full true coefficient vector \eqn{\beta}.}

\item{num_treats}{Integer. Total number of treatment effect parameters.}
}
\value{
A numeric vector of length R, where the r-th element is the true
average treatment effect for cohort r.
}
\description{
Given a full coefficient vector and information about treatment effect indices,
this function calculates the true average treatment effect for each cohort by
averaging the relevant treatment effect coefficients.
}
\keyword{internal}
