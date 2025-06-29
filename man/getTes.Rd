% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_funcs.R
\name{getTes}
\alias{getTes}
\title{Compute True Treatment Effects}
\usage{
getTes(coefs_obj)
}
\arguments{
\item{coefs_obj}{An object of class \code{"FETWFE_coefs"} containing the coefficient vector
and simulation parameters.}
}
\value{
A named list with two elements:
\describe{
\item{att_true}{A numeric value representing the overall average treatment effect on the
treated. It is computed as the (equal-weighted) mean of the cohort-specific treatment
effects.}
\item{actual_cohort_tes}{A numeric vector containing the true cohort-specific treatment
effects, calculated by averaging the coefficients corresponding to the treatment dummies
for each cohort.}
}
}
\description{
This function extracts the true treatment effects from a full coefficient vector
as generated by \code{genCoefs()}. It calculates the overall average treatment effect on the
treated (ATT) as the equal-weighted average of the cohort-specific treatment effects, and also
returns the individual treatment effects for each treated cohort.
}
\details{
The function internally uses auxiliary routines \code{getNumTreats()}, \code{getP()},
\code{getFirstInds()}, \code{getTreatInds()}, and \code{getActualCohortTes()} to determine the
correct indices of treatment effect coefficients in \code{beta}. The overall treatment effect
is computed as the simple average of these cohort-specific effects.
}
\examples{
\dontrun{
# Generate coefficients
coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)

# Compute the true treatment effects:
te_results <- getTes(coefs)

# Overall average treatment effect on the treated:
print(te_results$att_true)

# Cohort-specific treatment effects:
print(te_results$actual_cohort_tes)
}

}
