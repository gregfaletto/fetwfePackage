% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetwfe.R
\name{etwfe}
\alias{etwfe}
\title{Extended two-way fixed effects}
\usage{
etwfe(
  pdata,
  time_var,
  unit_var,
  treatment,
  response,
  covs = c(),
  indep_counts = NA,
  sig_eps_sq = NA,
  sig_eps_c_sq = NA,
  verbose = FALSE,
  alpha = 0.05,
  add_ridge = FALSE
)
}
\arguments{
\item{pdata}{Dataframe; the panel data set. Each row should represent an
observation of a unit at a time. Should contain columns as described below.}

\item{time_var}{Character; the name of a single column containing a variable
for the time period. This column is expected to contain integer values (for
example, years). Recommended encodings for dates include format YYYY, YYYYMM,
or YYYYMMDD, whichever is appropriate for your data.}

\item{unit_var}{Character; the name of a single column containing a variable
for each unit. This column is expected to contain character values (i.e. the
"name" of each unit).}

\item{treatment}{Character; the name of a single column containing a variable
for the treatment dummy indicator. This column is expected to contain integer
values, and in particular, should equal 0 if the unit was untreated at that
time and 1 otherwise. Treatment should be an absorbing state; that is, if
unit \code{i} is treated at time \code{t}, then it must also be treated at all times
\code{t} + 1, ..., \code{T}. Any units treated in the first time period will be removed
automatically. Please make sure yourself that at least some units remain
untreated at the final time period ("never-treated units").}

\item{response}{Character; the name of a single column containing the
response for each unit at each time. The response must be an integer or
numeric value.}

\item{covs}{(Optional.) Character; a vector containing the names of the
columns for covariates. All of these columns are expected to contain integer,
numeric, or factor values, and any categorical values will be automatically
encoded as binary indicators. If no covariates are provided, the treatment
effect estimation will proceed, but it will only be valid under unconditional
versions of the parallel trends and no anticipation assumptions. Default is c().}

\item{indep_counts}{(Optional.) Integer; a vector. If you have a sufficiently
large number of units, you can optionally randomly split your data set in
half (with \code{N} units in each data set). The data for half of the units should
go in the \code{pdata} argument provided above. For the other \code{N} units, simply
provide the counts for how many units appear in the untreated cohort plus
each of the other \code{R} cohorts in this argument \code{indep_counts}. The benefit
of doing this is that the standard error for the average treatment effect
will be (asymptotically) exact instead of conservative. The length of
\code{indep_counts} must equal 1 plus the number of treated cohorts in \code{pdata}.
All entries of \code{indep_counts} must be strictly positive (if you are concerned
that this might not work out, maybe your data set is on the small side and
it's best to just leave your full data set in \code{pdata}). The sum of all the
counts in \code{indep_counts} must match the total number of units in \code{pdata}.
Default is NA (in which case conservative standard errors will be calculated
if \code{q < 1}.)}

\item{sig_eps_sq}{(Optional.) Numeric; the variance of the row-level IID
noise assumed to apply to each observation. See Section 2 of Faletto (2025)
for details. It is best to provide this variance if it is known (for example,
if you are using simulated data). If this variance is unknown, this argument
can be omitted, and the variance will be estimated using the estimator from
Pesaran (2015, Section 26.5.1) with ridge regression. Default is NA.}

\item{sig_eps_c_sq}{(Optional.) Numeric; the variance of the unit-level IID
noise (random effects) assumed to apply to each observation. See Section 2 of
Faletto (2025) for details. It is best to provide this variance if it is
known (for example, if you are using simulated data). If this variance is
unknown, this argument can be omitted, and the variance will be estimated
using the estimator from Pesaran (2015, Section 26.5.1) with ridge
regression. Default is NA.}

\item{verbose}{Logical; if TRUE, more details on the progress of the function will
be printed as the function executes. Default is FALSE.}

\item{alpha}{Numeric; function will calculate (1 - \code{alpha}) confidence intervals
for the cohort average treatment effects that will be returned in \code{catt_df}.}

\item{add_ridge}{(Optional.) Logical; if TRUE, adds a small amount of ridge
regularization to the (untransformed) coefficients to stabilize estimation.
Default is FALSE.}
}
\value{
A named list with the following elements: \item{att_hat}{The
estimated overall average treatment effect for a randomly selected treated
unit.} \item{att_se}{A standard error for the ATT. If the Gram matrix is not
invertible, this will be NA.} \item{catt_hats}{A named vector containing the
estimated average treatment effects for each cohort.} \item{catt_ses}{A named
vector containing the (asymptotically exact) standard errors for
the estimated average treatment effects within each cohort.}
\item{cohort_probs}{A vector of the estimated probabilities of being in each
cohort conditional on being treated, which was used in calculating \code{att_hat}.
If \code{indep_counts} was provided, \code{cohort_probs} was calculated from that;
otherwise, it was calculated from the counts of units in each treated
cohort in \code{pdata}.} \item{catt_df}{A dataframe displaying the cohort names,
average treatment effects, standard errors, and \code{1 - alpha} confidence
interval bounds.} \item{beta_hat}{The full vector of estimated coefficients.}
\item{treat_inds}{The indices of \code{beta_hat} corresponding to
the treatment effects for each cohort at each time.}
\item{treat_int_inds}{The indices of \code{beta_hat} corresponding to the
interactions between the treatment effects for each cohort at each time and
the covariates.} \item{sig_eps_sq}{Either the provided \code{sig_eps_sq} or
the estimated one, if a value wasn't provided.} \item{sig_eps_c_sq}{Either
the provided \code{sig_eps_c_sq} or the estimated one, if a value wasn't
provided.} \item{X_ints}{The design matrix created containing all
interactions, time and cohort dummies, etc.} \item{y}{The vector of
responses, containing \code{nrow(X_ints)} entries.} \item{X_final}{The design
matrix after applying the change in coordinates to fit the model and also
multiplying on the left by the square root inverse of the estimated
covariance matrix for each unit.} \item{y_final}{The final response after
multiplying on the left by the square root inverse of the estimated
covariance matrix for each unit.} \item{N}{The final number of units that
were in the  data set used for estimation (after any units may have been
removed because they were treated in the first time period).} \item{T}{The
number of time periods in the final data set.} \item{R}{The final number of
treated cohorts that appear in the final data set.} \item{d}{The final number
of covariates that appear in the final data set (after any covariates may
have been removed because they contained missing values or all contained the
same value for every unit).} \item{p}{The final number of columns in the full
set of covariates used to estimate the model.}
}
\description{
Implementation of extended two-way fixed effects.
Estimates overall ATT as well as CATT (cohort average treatment effects on
the treated units).
}
\references{
Wooldridge, J. M. (2021). Two-way fixed effects, the two-way mundlak
regression, and difference-in-differences estimators.
\emph{Available at SSRN 3906345}.
\doi{10.2139/ssrn.3906345}.
}
\author{
Gregory Faletto
}
