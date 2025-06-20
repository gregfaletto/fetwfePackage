% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betwfe_core.R
\name{betwfe}
\alias{betwfe}
\title{Bridge-penalized extended two-way fixed effects}
\usage{
betwfe(
  pdata,
  time_var,
  unit_var,
  treatment,
  response,
  covs = c(),
  indep_counts = NA,
  sig_eps_sq = NA,
  sig_eps_c_sq = NA,
  lambda.max = NA,
  lambda.min = NA,
  nlambda = 100,
  q = 0.5,
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

\item{lambda.max}{(Optional.) Numeric. A penalty parameter \code{lambda} will be
selected over a grid search by BIC in order to select a single model. The
largest \code{lambda} in the grid will be \code{lambda.max}. If no \code{lambda.max} is
provided, one will be selected automatically. When \code{q <= 1}, the model
will be sparse, and ideally all of the following are true at once: the
smallest model (the one corresponding to \code{lambda.max}) selects close to 0
features, the largest model (the one corresponding to \code{lambda.min}) selects
close to \code{p} features, \code{nlambda} is large enough so that models are
considered at every feasible model size, and \code{nlambda} is small enough so
that the computation doesn't become infeasible. You may
want to manually tweak \code{lambda.max}, \code{lambda.min}, and \code{nlambda} to try
to achieve these goals, particularly if the selected model size is very
close to the model corresponding to \code{lambda.max} or \code{lambda.min}, which could
indicate that the range of \code{lambda} values was too narrow or coarse. You can
use the function outputs \code{lambda.max_model_size}, \code{lambda.min_model_size}, and
\code{lambda_star_model_size} to try to assess this. Default is NA.}

\item{lambda.min}{(Optional.) Numeric. The smallest \code{lambda} penalty
parameter that will be considered. See the description of \code{lambda.max} for
details. Default is NA.}

\item{nlambda}{(Optional.) Integer. The total number of \code{lambda} penalty
parameters that will be considered. See the description of \code{lambda.max} for
details. Default is 100.}

\item{q}{(Optional.) Numeric; determines what \code{L_q} penalty is used for the
regularization. \code{q} = 1 is the lasso, and for 0 < \code{q} < 1, it is
possible to get standard errors and confidence intervals. \code{q} = 2 is ridge
regression. See Faletto (2025) for details. Default is 0.5.}

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
unit.} \item{att_se}{If \code{q < 1}, a standard error for the ATT. If
\code{indep_counts} was provided, this standard error is asymptotically exact; if
not, it is asymptotically conservative. If \code{q >= 1}, this will be NA.}
\item{catt_hats}{A named vector containing the estimated average treatment
effects for each cohort.} \item{catt_ses}{If \code{q < 1}, a named vector
containing the (asymptotically exact, non-conservative) standard errors for
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
provided.} \item{lambda.max}{Either the provided \code{lambda.max} or the one
that was used, if a value wasn't provided. (This is returned to help with
getting a reasonable range of \code{lambda} values for grid search.)}
\item{lambda.max_model_size}{The size of the selected model corresponding
\code{lambda.max} (for \code{q <= 1}, this will be the smallest model size). As
mentioned above, for \code{q <= 1} ideally this value is close to 0.}
\item{lambda.min}{Either the provided \code{lambda.min} or the one
that was used, if a value wasn't provided.} \item{lambda.min_model_size}{The
size of the selected model corresponding to \code{lambda.min} (for \code{q <= 1}, this
will be the largest model size). As mentioned above, for \code{q <= 1} ideally
this value is close to \code{p}.}\item{lambda_star}{The value of \code{lambda} chosen
by BIC. If this value is close to \code{lambda.min} or \code{lambda.max}, that could
suggest that the range of \code{lambda} values should be expanded.}
\item{lambda_star_model_size}{The size of the model that was selected. If
this value is close to \code{lambda.max_model_size} or \code{lambda.min_model_size},
That could suggest that the range of \code{lambda} values should be expanded.}
\item{X_ints}{The design matrix created containing all
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
Implementation of extended two-way fixed effects with a bridge
penalty. Estimates overall ATT as well as CATT (cohort average treatment
effects on the treated units).
}
\examples{
set.seed(23451)

library(bacondecomp)

data(divorce)

# sig_eps_sq and sig_eps_c_sq, calculated in a separate run of `fetwfe(),
# are provided to speed up the computation of the example
res <- betwfe(
    pdata = divorce[divorce$sex == 2, ],
    time_var = "year",
    unit_var = "st",
    treatment = "changed",
    covs = c("murderrate", "lnpersinc", "afdcrolls"),
    response = "suiciderate_elast_jag",
    sig_eps_sq = 0.1025361,
    sig_eps_c_sq = 4.227651e-35,
    verbose = TRUE)

# Average treatment effect on the treated units (in percentage point
# units)
100 * res$att_hat

# Conservative 95\% confidence interval for ATT (in percentage point units)

low_att <- 100 * (res$att_hat - qnorm(1 - 0.05 / 2) * res$att_se)
high_att <- 100 * (res$att_hat + qnorm(1 - 0.05 / 2) * res$att_se)

c(low_att, high_att)

# Cohort average treatment effects and confidence intervals (in percentage
# point units)

catt_df_pct <- res$catt_df
catt_df_pct[["Estimated TE"]] <- 100 * catt_df_pct[["Estimated TE"]]
catt_df_pct[["SE"]] <- 100 * catt_df_pct[["SE"]]
catt_df_pct[["ConfIntLow"]] <- 100 * catt_df_pct[["ConfIntLow"]]
catt_df_pct[["ConfIntHigh"]] <- 100 * catt_df_pct[["ConfIntHigh"]]

catt_df_pct
}
\references{
Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions.
\emph{arXiv preprint arXiv:2312.05985}.
\url{https://arxiv.org/abs/2312.05985}.
Pesaran, M. H. . Time Series and Panel Data Econometrics. Number 9780198759980 in OUP
Catalogue. Oxford University Press, 2015. URL
\url{https://ideas.repec.org/b/oxp/obooks/9780198759980.html}.
}
\author{
Gregory Faletto
}
