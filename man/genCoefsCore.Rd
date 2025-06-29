% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_funcs.R
\name{genCoefsCore}
\alias{genCoefsCore}
\title{Generate Coefficient Vector for Data Generation}
\usage{
genCoefsCore(R, T, d, density, eff_size, seed = NULL)
}
\arguments{
\item{R}{Integer. The number of treated cohorts (treatment is assumed to start in periods 2 to
\code{R + 1}).}

\item{T}{Integer. The total number of time periods.}

\item{d}{Integer. The number of time-invariant covariates. If \code{d > 0}, additional terms
corresponding to covariate main effects and interactions are included in \code{beta}.}

\item{density}{Numeric in (0,1). The probability that any given entry in the initial sparse
coefficient vector \code{theta} is nonzero.}

\item{eff_size}{Numeric. The magnitude used to scale nonzero entries in \code{theta}. Each
nonzero entry is set to \code{eff_size} or \code{-eff_size} (with a 60 percent chance for a
positive value).}

\item{seed}{(Optional) Integer. Seed for reproducibility.}
}
\value{
A list with two elements:
\describe{
\item{\code{beta}}{A numeric vector representing the full coefficient vector after the inverse
fusion transform.}
\item{theta}{A numeric vector representing the coefficient vector in the transformed feature
space. \code{theta} is a sparse vector, which aligns with an assumption that deviations from the
restrictions encoded in the FETWFE model are sparse. \code{beta} is derived from
\code{theta}.}
}
}
\description{
This function generates a coefficient vector \code{beta} along with a sparse auxiliary vector
\code{theta} for simulation studies of the fused extended two-way fixed effects estimator. The
returned \code{beta} is formatted to align with the design matrix created by
\code{genRandomData()}, and is a valid input for the \code{beta} argument of that function. The
vector \code{theta} is sparse, with nonzero entries occurring with probability \code{density} and
scaled by \code{eff_size}. See the simulation studies section of Faletto (2025) for details.
}
\details{
The length of \code{beta} is given by
\deqn{p = R + (T - 1) + d + dR + d(T - 1) + \mathit{num\_treats} + (\mathit{num\_treats} \times d)}{p = R + (T - 1) + d + dR + d(T - 1) + num_treats + (num_treats * d)},
where the number of treatment parameters is defined as
\deqn{\mathit{num\_treats} = T \times R - \frac{R(R+1)}{2}}{num_treats = T * R - R(R+1)/2}.

The function operates in two steps:
\enumerate{
\item It first creates a sparse vector \code{theta} of length \eqn{p}, with nonzero entries
occurring
with probability \code{density}. Nonzero entries are set to \code{eff_size} or \code{-eff_size}
(with a 60\\% chance for a positive value).
\item The full coefficient vector \code{beta} is then computed by applying an inverse fusion
transform to \code{theta} using internal routines (e.g.,
\code{genBackwardsInvFusionTransformMat()} and \code{genInvTwoWayFusionTransformMat()}).
}
}
\examples{
\dontrun{
  # Set parameters for the coefficient generation
  R <- 3         # Number of treated cohorts
  T <- 6         # Total number of time periods
  d <- 2         # Number of covariates
  density <- 0.1 # Probability that an entry in the initial vector is nonzero
  eff_size <- 1.5  # Scaling factor for nonzero coefficients
  seed <- 789    # Seed for reproducibility

  # Generate coefficients using genCoefsCore()
  coefs_core <- genCoefsCore(R = R, T = T, d = d, density = density,
  eff_size = eff_size, seed = seed)
  beta <- coefs_core$beta
  theta <- coefs_core$theta

  # For diagnostic purposes, compute the expected length of beta.
  # The length p is defined internally as:
  #   p = R + (T - 1) + d + d*R + d*(T - 1) + num_treats + num_treats*d,
  # where num_treats = T * R - (R*(R+1))/2.
  num_treats <- T * R - (R * (R + 1)) / 2
  p_expected <- R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d

  cat("Length of beta:", length(beta), "\nExpected length:", p_expected, "\n")
}

}
\references{
Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions.
\emph{arXiv preprint arXiv:2312.05985}.
\url{https://arxiv.org/abs/2312.05985}.
}
