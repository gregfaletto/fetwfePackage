#' @import glmnet

#' @title Fused extended two-way fixed effects
#'
#' @description Implementation of fused extended two-way fixed effects.
#' Estimates overall ATT as well as CATT (cohort average treatment effects on
#' the treated units). 
#'
#' @param pdata Dataframe; the panel data set. Each row should represent an
#' observation of a unit at a time. Should contain columns as described below.
#' @param time_var Character; the name of a single column containing a variable
#' for the time period. This column is expected to contain integer values (for
#' example, years). Recommended encodings for dates include format YYYY, YYYYMM,
#' or YYYYMMDD, whichever is appropriate for your data.
#' @param unit_var Character; the name of a single column containing a variable
#' for each unit. This column is expected to contain character values (i.e. the
#' "name" of each unit).
#' @param treatment Character; the name of a single column containing a variable
#' for the treatment dummy indicator. This column is expected to contain integer
#' values, and in particular, should equal 0 if the unit was untreated at that
#' time and 1 otherwise. Treatment should be an absorbing state; that is, if
#' unit `i` is treated at time `t`, then it must also be treated at all times
#' `t` + 1, ..., `T`. Any units treated in the first time period will be removed
#' automatically. Please make sure yourself that at least some units remain
#' untreated at the final time period ("never-treated units").
#' @param response Character; the name of a single column containing the
#' response for each unit at each time. The response must be an integer or
#' numeric value.
#' @param covs (Optional.) Character; a vector containing the names of the 
#' columns for covariates. All of these columns are expected to contain integer
#' or numeric values (so if you use categorical values, encode them using e.g.
#' binary indicators before passing the data to this function). If no covariates
#' are provided, the treatment effect estimation will proceed, but it will only
#' be valid under unconditional assumptions of parallel trends and no 
#' anticipation. Default is c().
#' @param indep_counts (Optional.) Integer; a vector. If you have a sufficiently
#' large number of units, you can optionally randomly split your data set in
#' half (with `N` units in each data set). The data for half of the units should
#' go in the `pdata` argument provided above. For the other `N` units, simply
#' provide the counts for how many units appear in the untreated cohort plus
#' each of the other `R` cohorts in this argument `indep_counts`. The benefit
#' of doing this is that the standard error for the average treatment effect
#' will be (asymptotically) exact instead of conservative. The length of
#' `indep_counts` must equal 1 plus the number of treated cohorts in `pdata`.
#' All entries of `indep_counts` must be strictly positive (if you are concerned
#' that this might not work out, maybe your data set is on the small side and
#' it's best to just leave your full data set in `pdata`). The sum of all the
#' counts in `indep_counts` must match the total number of units in `pdata`.
#' Default is NA (in which case conservative standard errors will be calculated
#' if `q < 1`.)
#' @param sig_eps_sq (Optional.) Numeric; the variance of the row-level IID
#' noise assumed to apply to each observation. See Section 2 of Faletto (2024)
#' for details. It is best to provide this variance if it is known (for example,
#' if you are using simulated data). If this variance is unknown, this argument
#' can be omitted, and the variance will be estimated using the estimator from
#' Pesaran (2015, Section 26.5.1) with ridge regression. Default is NA.
#' @param sig_eps_c_sq (Optional.) Numeric; the variance of the unit-level IID
#' noise (random effects) assumed to apply to each observation. See Section 2 of
#' Faletto (2024) for details. It is best to provide this variance if it is
#' known (for example, if you are using simulated data). If this variance is
#' unknown, this argument can be omitted, and the variance will be estimated
#' using the estimator from Pesaran (2015, Section 26.5.1) with ridge
#' regression. Default is NA.
#' @param lambda.max (Optional.) Numeric. A penalty parameter `lambda` will be
#' selected over a grid search by BIC in order to select a single model. The 
#' largest `lambda` in the grid will be `lambda.max`. If no `lambda.max` is
#' provided, one will be selected automatically. For `lambda <= 1`, the model
#' will be sparse, and ideally all of the following are true at once: the
#' smallest model (the one corresponding to `lambda.max`) selects close to 0
#' features, the largest model (the one corresponding to `lambda.min`) selects
#' close to `p` features, `nlambda` is large enough so that models are
#' considered at every feasible model size, and `nlambda` is small enough so
#' that the computation doesn't become infeasible. You may
#' want to manually tweak `lambda.max`, `lambda.min`, and `nlambda` to try
#' to achieve these goals, particularly if the selected model size is very
#' close to the model corresponding to `lambda.max` or `lambda.min`, which could
#' indicate that the range of `lambda` values was too narrow. You can use the
#' function outputs `lambda.max_model_size`, `lambda.min_model_size`, and
#' `lambda_star_model_size` to try to assess this. Default is NA.
#' @param lambda.min (Optional.) Numeric. The smallest `lambda` penalty
#' parameter that will be considered. See the description of `lambda.max` for
#' details. Default is NA.
#' @param nlambda (Optional.) Integer. The total number of `lambda` penalty
#' parameters that will be considered. See the description of `lambda.max` for
#' details. Default is 100.
#' @param q (Optional.) Numeric; determines what `L_q` penalty is used for the
#' fusion regularization. `q` = 1 is the lasso, and for 0 < `q` < 1, it is 
#' possible to get standard errors and confidence intervals. `q` = 2 is ridge
#' regression. See Faletto (2024) for details. Default is 0.5.
#' @param verbose Logical; if TRUE, more details on the progress of the function will
#' be printed as the function executes. Default is FALSE.
#' @param alpha Numeric; function will calculate (1 - `alpha`) confidence intervals
#' for the cohort average treatment effects that will be returned in `catt_df`.
#' @param add_ridge (Optional.) Logical; if TRUE, adds a small amount of ridge
#' regularization to the (untransformed) coefficients to stabilize estimation.
#' Default is FALSE.
#' @return A named list with the following elements: \item{att_hat}{The
#' estimated overall average treatment effect for a randomly selected treated
#' unit.} \item{att_se}{If `q < 1`, a standard error for the ATT. If
#' `indep_counts` was provided, this standard error is asymptotically exact; if
#' not, it is asymptotically conservative. If `q >= 1`, this will be NA.}
#' \item{catt_hats}{A named vector containing the estimated average treatment
#' effects for each cohort.} \item{catt_ses}{If `q < 1`, a named vector
#' containing the (asymptotically exact, non-conservative) standard errors for
#' the estimated average treatment effects within each cohort.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each
#' cohort conditional on being treated, which was used in calculating `att_hat`.
#' If `indep_counts` was provided, `cohort_probs` was calculated from that;
#' otherwise, it was calculated from the counts of units in each treated
#' cohort in `pdata`.} \item{catt_df}{A dataframe displaying the cohort names,
#' average treatment effects, standard errors, and `1 - alpha` confidence
#' interval bounds.} \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to
#' the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the
#' interactions between the treatment effects for each cohort at each time and
#' the covariates.} \item{sig_eps_sq}{Either the provided `sig_eps_sq` or
#' the estimated one, if a value wasn't provided.} \item{sig_eps_c_sq}{Either
#' the provided `sig_eps_c_sq` or the estimated one, if a value wasn't
#' provided.} \item{lambda.max}{Either the provided `lambda.max` or the one
#' that was used, if a value wasn't provided. (This is returned to help with
#' getting a reasonable range of `lambda` values for grid search.)}
#' \item{lambda.max_model_size}{The size of the selected model corresponding
#' `lambda.max` (for `q <= 1`, this will be the smallest model size). As
#' mentioned above, for `q <= 1` ideally this value is close to 0.}
#' \item{lambda.min}{Either the provided `lambda.min` or the one
#' that was used, if a value wasn't provided.} \item{lambda.min_model_size}{The
#' size of the selected model corresponding `lambda.min` (for `q <= 1`, this
#' will be the largest model size). As mentioned above, for `q <= 1` ideally
#' this value is close to `p`.}\item{lambda_star}{The value of `lambda` chosen
#' by BIC. If this value is close to `lambda.min` or `lambda.max`, that could
#' suggest that the range of `lambda` values should be expanded.}
#' \item{lambda_star_model_size}{The size of the model that was selected. If
#' this value is close to `lambda.max_model_size` or `lambda.min_model_size`,
#' That could suggest that the range of `lambda` values should be expanded.}
#' \item{X_ints}{The design matrix created containing all
#' interactions, time and cohort dummies, etc.} \item{y}{The vector of 
#' responses, containing `nrow(X_ints)` entries.} \item{X_final}{The design
#' matrix after applying the change in coordinates to fit the model and also
#' multiplying on the left by the square root inverse of the estimated
#' covariance matrix for each unit.} \item{y_final}{The final response after
#' multiplying on the left by the square root inverse of the estimated
#' covariance matrix for each unit.} \item{N}{The final number of units that
#' were in the  data set used for estimation (after any units may have been
#' removed because they were treated in the first time period).} \item{T}{The 
#' number of time periods in the final data set.} \item{R}{The final number of
#' treated cohorts that appear in the final data set.} \item{d}{The final number
#' of covariates that appear in the final data set (after any covariates may
#' have been removed because they contained missing values or all contained the
#' same value for every unit).} \item{p}{The final number of columns in the full
#' set of covariates used to estimate the model.}
#' @author Gregory Faletto
#' @references
#' Faletto, G (2024). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#' Pesaran, M. H. . Time Series and Panel Data Econometrics. Number 9780198759980 in OUP
#' Catalogue. Oxford University Press, 2015. URL
#' \url{https://ideas.repec.org/b/oxp/obooks/9780198759980.html}.
#' @examples
#' set.seed(23451)
#' 
#' library(bacondecomp)
#' 
#' data(divorce)
#' 
#' # sig_eps_sq and sig_eps_c_sq, calculated in a separate run of `fetwfe(),
#' # are provided to speed up the computation of the example
#' res <- fetwfe(
#'     pdata = divorce[divorce$sex == 2, ],
#'     time_var = "year",
#'     unit_var = "st",
#'     treatment = "changed",
#'     covs = c("murderrate", "lnpersinc", "afdcrolls"),
#'     response = "suiciderate_elast_jag",
#'     sig_eps_sq = 0.1025361,
#'     sig_eps_c_sq = 4.227651e-35,
#'     verbose = TRUE)
#'
#' # Average treatment effect on the treated units (in percentage point
#' # units)
#' 100 * res$att_hat
#' 
#' # Conservative 95% confidence interval for ATT (in percentage point units)
#'
#' low_att <- 100 * (res$att_hat - qnorm(1 - 0.05 / 2) * res$att_se)
#' high_att <- 100 * (res$att_hat + qnorm(1 - 0.05 / 2) * res$att_se)
#'
#' c(low_att, high_att)
#'
#' # Cohort average treatment effects and confidence intervals (in percentage
#' # point units)
#'
#' catt_df_pct <- res$catt_df
#' catt_df_pct[["Estimated TE"]] <- 100 * catt_df_pct[["Estimated TE"]]
#' catt_df_pct[["SE"]] <- 100 * catt_df_pct[["SE"]]
#' catt_df_pct[["ConfIntLow"]] <- 100 * catt_df_pct[["ConfIntLow"]]
#' catt_df_pct[["ConfIntHigh"]] <- 100 * catt_df_pct[["ConfIntHigh"]]
#'
#' catt_df_pct
#' @export
fetwfe <- function(
    pdata,
    time_var,
    unit_var,
    treatment,
    response,
    covs=c(),
    indep_counts=NA,
    sig_eps_sq=NA,
    sig_eps_c_sq=NA,
    lambda.max=NA,
    lambda.min=NA,
    nlambda=100,
    q=0.5,
    verbose=FALSE,
    alpha=0.05,
    add_ridge=FALSE
    ){

    # Check inputs
    indep_count_data_available <- checkFetwfeInputs(
        pdata=pdata,
        time_var=time_var,
        unit_var=unit_var,
        treatment=treatment,
        response=response,
        covs=covs,
        indep_counts=indep_counts,
        sig_eps_sq=sig_eps_sq,
        sig_eps_c_sq=sig_eps_c_sq,
        lambda.max=lambda.max,
        lambda.min=lambda.min,
        nlambda=nlambda,
        q=q,
        verbose=verbose,
        alpha=alpha,
        add_ridge=add_ridge
    )

    pdata <- pdata[, c(response, time_var, unit_var, treatment, covs)]

    res <- prepXints(
        data=pdata,
        time_var=time_var,
        unit_var=unit_var,
        treatment=treatment,
        covs=covs,
        response=response,
        verbose=verbose
        )

    X_ints <- res$X_ints
    y <- res$y
    N <- res$N
    T <- res$T
    d <- res$d
    p <- res$p
    in_sample_counts <- res$in_sample_counts
    num_treats <- res$num_treats
    first_inds <- res$first_inds

    rm(res)

    R <- length(in_sample_counts) - 1
    stopifnot(R >= 1)
    stopifnot(R <= T - 1)
    if(R < 2){
        stop("Only one treated cohort detected in data. Currently fetwfe only supports data sets with at least two treated cohorts.")
    }

    stopifnot(N >= R + 1)

    stopifnot(sum(in_sample_counts) == N)
    stopifnot(all(in_sample_counts >= 0))
    stopifnot(is.integer(in_sample_counts))

    if(in_sample_counts[1] == 0){
        stop("No never-treated units detected in data to fit model; estimating treatment effects is not possible")
    }
    if(length(names(in_sample_counts)) != length(in_sample_counts)){
        stop("in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)")
    }

    if(length(names(in_sample_counts)) != length(unique(names(in_sample_counts)))){
        stop("in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)")
    }

    if(indep_count_data_available){
        if(sum(indep_counts) != N){
            stop("Number of units in independent cohort count data does not equal number of units in data to be used to fit model.")
        }
        if(length(indep_counts) != length(in_sample_counts)){
            stop("Number of counts in independent counts does not match number of cohorts in data to be used to fit model.")
        }
    }

    res <- fetwfe_core(
        X_ints=X_ints,
        y=y,
        in_sample_counts=in_sample_counts,
        N=N,
        T=T,
        d=d,
        p=p,
        num_treats=num_treats,
        first_inds=first_inds,
        indep_counts=indep_counts,
        sig_eps_sq=sig_eps_sq,
        sig_eps_c_sq=sig_eps_c_sq,
        lambda.max=lambda.max,
        lambda.min=lambda.min,
        nlambda=nlambda,
        q=q,
        verbose=verbose,
        alpha=alpha,
        add_ridge=add_ridge
        )

    if(indep_count_data_available){
        if(q < 1){
            stopifnot(!is.na(res$indep_att_hat))
            stopifnot(!is.na(res$indep_att_se))
        }
        stopifnot(all(!is.na(res$indep_cohort_probs)))
        att_hat <- res$indep_att_hat
        att_se <- res$indep_att_se
        cohort_probs <- res$indep_cohort_probs
    } else{
        att_hat <- res$in_sample_att_hat
        att_se <- res$in_sample_att_se
        cohort_probs <- res$cohort_probs
    }

    return(list(
        att_hat=att_hat,
        att_se=att_se,
        catt_hats=res$catt_hats,
        catt_ses=res$catt_ses,
        cohort_probs=cohort_probs,
        catt_df=res$catt_df,
        beta_hat=res$beta_hat,
        treat_inds=res$treat_inds,
        treat_int_inds=res$treat_int_inds,
        sig_eps_sq=res$sig_eps_sq,
        sig_eps_c_sq=res$sig_eps_c_sq,
        lambda.max = res$lambda.max,
        lambda.max_model_size = res$lambda.max_model_size,
        lambda.min = res$lambda.min,
        lambda.min_model_size = res$lambda.min_model_size,
        lambda_star = res$lambda_star,
        lambda_star_model_size = res$lambda_star_model_size,
        X_ints=res$X_ints,
        y=res$y,
        X_final=res$X_final,
        y_final=res$y_final,
        N=res$N,
        T=res$T,
        R=res$R,
        d=res$d,
        p=res$p
        )
    )

}


#' Generate Random Panel Data for FETWFE Simulations
#'
#' @description
#' Generates a random panel dataset for simulation studies of the fused extended two-way fixed
#' effects (FETWFE) estimator by taking an object of class  \code{"FETWFE_coefs"} (produced by
#' \code{genCoefs()}) and using it to simulate data. The function creates a balanced panel
#' with \eqn{N} units over \eqn{T} time periods, assigns treatment status across \eqn{R}
#' treated cohorts (with equal marginal probabilities for treatment and non-treatment), and
#' constructs a design matrix along with the corresponding outcome. The covariates are
#' generated according to the specified \code{distribution}: by default, covariates are drawn
#' from a normal distribution; if \code{distribution = "uniform"}, they are drawn uniformly
#' from \eqn{[-\sqrt{3}, \sqrt{3}]}. When \eqn{d = 0} (i.e. no covariates), no
#' covariate-related columns or interactions are generated.
#'
#' @param coefs_obj An object of class \code{"FETWFE_coefs"} containing the coefficient vector
#' and simulation parameters.
#' @param N Integer. Number of units in the panel.
#' @param sig_eps_sq Numeric. Variance of the idiosyncratic (observation-level) noise.
#' @param sig_eps_c_sq Numeric. Variance of the unit-level random effects.
#' @param seed (Optional) Integer. Seed for reproducibility.
#' @param distribution Character. Distribution to generate covariates.
#'   Defaults to \code{"gaussian"}. If set to \code{"uniform"}, covariates are drawn uniformly
#'   from \eqn{[-\sqrt{3}, \sqrt{3}]}.
#'
#' @return An object of class \code{"FETWFE_simulated"}, which is a list containing:
#' \describe{
#'   \item{pdf}{A dataframe containing generated data that can be passed to \code{fetwfe()}.}
#'   \item{X}{The design matrix. When \code{gen_ints = TRUE}, \eqn{X} has \eqn{p} columns with
#'     interactions; when \code{gen_ints = FALSE}, \eqn{X} has no interactions.}
#'   \item{y}{A numeric vector of length \eqn{N \times T} containing the generated responses.}
#'   \item{covs}{A character vector containing the names of the generated features (if \eqn{d > 0}),
#'          or simply an empty vector (if \eqn{d = 0})}
#'   \item{coefs}{The coefficient vector \eqn{\beta} used for data generation.}
#'   \item{first_inds}{A vector of indices indicating the first treatment effect for each treated cohort.}
#'   \item{N_UNTREATED}{The number of never-treated units.}
#'   \item{assignments}{A vector of counts (of length \eqn{R+1}) indicating how many units fall into
#'         the never-treated group and each of the \eqn{R} treated cohorts.}
#'   \item{indep_counts}{Independent cohort assignments (for auxiliary purposes).}
#'   \item{p}{The number of columns in the design matrix \eqn{X}.}
#'   \item{N}{Number of units.}
#'   \item{T}{Number of time periods.}
#'   \item{R}{Number of treated cohorts.}
#'   \item{d}{Number of covariates.}
#'   \item{sig_eps_sq}{The idiosyncratic noise variance.}
#'   \item{sig_eps_c_sq}{The unit-level noise variance.}
#' }
#'
#' @details
#' This function extracts simulation parameters from the \code{FETWFE_coefs} object and passes them,
#' along with additional simulation parameters, to the internal function \code{simulateDataCore()}.
#' It validates that all necessary components are returned and assigns the S3 class \code{"FETWFE_simulated"}
#' to the output.
#'
#' The argument \code{distribution} controls the generation of covariates. For
#' \code{"gaussian"}, covariates are drawn from \code{rnorm}; for \code{"uniform"},
#' they are drawn from \code{runif} on the interval \eqn{[-\sqrt{3}, \sqrt{3}]}.
#'
#' When \eqn{d = 0} (i.e. no covariates), the function omits any covariate-related columns
#' and their interactions.
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)
#' }
#'
#' @export
simulateData <- function(coefs_obj, N, sig_eps_sq, sig_eps_c_sq, seed = NULL,
                         distribution = "gaussian") {
  if (!inherits(coefs_obj, "FETWFE_coefs")) {
    stop("coefs_obj must be an object of class 'FETWFE_coefs'")
  }
  if (!is.numeric(N) || length(N) != 1 || N <= 0) {
    stop("N must be a positive numeric value")
  }
  if (!is.numeric(sig_eps_sq) || length(sig_eps_sq) != 1 || sig_eps_sq <= 0) {
    stop("sig_eps_sq must be a positive numeric value")
  }
  if (!is.numeric(sig_eps_c_sq) || length(sig_eps_c_sq) != 1 || sig_eps_c_sq <= 0) {
    stop("sig_eps_c_sq must be a positive numeric value")
  }
  
  # Extract parameters from the coefs object
  R <- coefs_obj$R
  T <- coefs_obj$T
  d <- coefs_obj$d
  beta <- coefs_obj$beta
  
  sim_data <- simulateDataCore(N=N, T=T, R=R, d=d, sig_eps_sq=sig_eps_sq,
    sig_eps_c_sq=sig_eps_c_sq, beta=beta, seed = seed, gen_ints = FALSE,
    distribution = distribution)
  
  required_fields <- c("pdata", "X", "y", "covs", "time_var", "unit_var",
                       "treatment", "response", "coefs", "first_inds", "N_UNTREATED",
                       "assignments", "indep_counts", "p", "N", "T", "R", "d",
                       "sig_eps_sq", "sig_eps_c_sq")
  missing_fields <- setdiff(required_fields, names(sim_data))
  if (length(missing_fields) > 0) {
    stop(paste("simulateDataCore did not return expected components:",
               paste(missing_fields, collapse = ", ")))
  }
  
  obj <- sim_data
  class(obj) <- "FETWFE_simulated"
  return(obj)
}

#' Generate Coefficient Vector for Data Generation
#'
#' This function generates a coefficient vector \code{beta} for simulation studies of the fused
#' extended two-way fixed effects estimator. Itreturns an S3 object of class
#' \code{"FETWFE_coefs"} containing \code{beta} along with simulation parameters \code{R},
#' \code{T}, and \code{d}.
#'
#' @param R Integer. The number of treated cohorts (treatment is assumed to start in periods 2 to \code{R + 1}).
#' @param T Integer. The total number of time periods.
#' @param d Integer. The number of time-invariant covariates. If \code{d > 0}, additional terms corresponding
#'   to covariate main effects and interactions are included in \code{beta}.
#' @param density Numeric in (0,1). The probability that any given entry in the initial sparse coefficient
#'   vector \code{theta} is nonzero.
#' @param eff_size Numeric. The magnitude used to scale nonzero entries in \code{theta}. Each nonzero entry is
#'   set to \code{eff_size} or \code{-eff_size} (with a 60 percent chance for a positive value).
#'
#' @return An object of class \code{"FETWFE_coefs"}, which is a list containing:
#' \describe{
#'   \item{beta}{A numeric vector representing the full coefficient vector after the inverse fusion transform.}
#'   \item{R}{The provided number of treated cohorts.}
#'   \item{T}{The provided number of time periods.}
#'   \item{d}{The provided number of covariates.}
#' }
#'
#' @details
#' The length of \code{beta} is given by
#' \deqn{p = R + (T - 1) + d + dR + d(T - 1) + \mathit{num\_treats} + (\mathit{num\_treats} \times d)}{p = R + (T - 1) + d + dR + d(T - 1) + num_treats + (num_treats * d)},
#' where the number of treatment parameters is defined as
#' \deqn{\mathit{num\_treats} = T \times R - \frac{R(R+1)}{2}}{num_treats = T * R - R(R+1)/2}.
#'
#' The function operates in two steps:
#' \enumerate{
#'   \item It first creates a sparse vector \code{theta} of length \eqn{p}, with nonzero entries occurring
#'   with probability \code{density}. Nonzero entries are set to \code{eff_size} or \code{-eff_size} (with a 60\%
#'   chance for a positive value).
#'   \item The full coefficient vector \code{beta} is then computed by applying an inverse fusion transform to \code{theta}
#'   using internal routines (e.g., \code{genBackwardsInvFusionTransformMat()} and \code{genInvTwoWayFusionTransformMat()}).
#' }
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)
#' }
#'
#' @export
genCoefs <- function(R, T, d, density, eff_size){

    # Check that T is a numeric scalar and at least 3.
  if (!is.numeric(T) || length(T) != 1 || T < 3) {
    stop("T must be a numeric value greater than or equal to 3")
  }
  
  # Check that R is a numeric scalar and at least 2.
  if (!is.numeric(R) || length(R) != 1 || R < 2) {
    stop("R must be a numeric value greater than or equal to 2 (currently there is only support for data sets with staggered adoptions, so at least two treated cohorts)")
  }
  
  # Check that R does not exceed T - 1.
  if (R > T - 1) {
    stop("R must be less than or equal to T - 1")
  }
  
  # Check that d is a numeric scalar and is non-negative.
  if (!is.numeric(d) || length(d) != 1 || d < 0) {
    stop("d must be a non-negative numeric value")
  }
  
  # Check that density is a numeric scalar strictly between 0 and 1.
  if (!is.numeric(density) || length(density) != 1 || density <= 0 || density >= 1) {
    stop("density must be numeric and strictly between 0 and 1")
  }
  
  # Check that eff_size is numeric.
  if (!is.numeric(eff_size) || length(eff_size) != 1) {
    stop("eff_size must be a numeric value")
  }

    stopifnot(R >= 2)
    stopifnot(T >= 3)
    stopifnot(R <= T - 1)
  
  core_obj <- genCoefsCore(R, T, d, density, eff_size)
  if (is.null(core_obj$beta)) {
    stop("Internal error: genCoefsCore() did not return expected components.")
  }
  
  # Create an S3 object of class "FETWFE_coefs"
  obj <- list(beta = core_obj$beta, R = R, T = T, d = d)
  class(obj) <- "FETWFE_coefs"
  return(obj)
}

#' Run FETWFE on Simulated Data
#'
#' @description
#' This function runs the fused extended two-way fixed effects estimator (\code{fetwfe()}) on simulated data.
#' It accepts an object of class \code{"FETWFE_simulated"} (produced by \code{simulateData()}) and unpacks the
#' necessary components to pass to \code{fetwfe()}.
#'
#' @param simulated_obj An object of class \code{"FETWFE_simulated"} containing the simulated panel data and design matrix.
#' @param lambda.max (Optional.) Numeric. A penalty parameter `lambda` will be
#' selected over a grid search by BIC in order to select a single model. The 
#' largest `lambda` in the grid will be `lambda.max`. If no `lambda.max` is
#' provided, one will be selected automatically. For `lambda <= 1`, the model
#' will be sparse, and ideally all of the following are true at once: the
#' smallest model (the one corresponding to `lambda.max`) selects close to 0
#' features, the largest model (the one corresponding to `lambda.min`) selects
#' close to `p` features, `nlambda` is large enough so that models are
#' considered at every feasible model size, and `nlambda` is small enough so
#' that the computation doesn't become infeasible. You may
#' want to manually tweak `lambda.max`, `lambda.min`, and `nlambda` to try
#' to achieve these goals, particularly if the selected model size is very
#' close to the model corresponding to `lambda.max` or `lambda.min`, which could
#' indicate that the range of `lambda` values was too narrow. You can use the
#' function outputs `lambda.max_model_size`, `lambda.min_model_size`, and
#' `lambda_star_model_size` to try to assess this. Default is NA.
#' @param lambda.min (Optional.) Numeric. The smallest `lambda` penalty
#' parameter that will be considered. See the description of `lambda.max` for
#' details. Default is NA.
#' @param nlambda (Optional.) Integer. The total number of `lambda` penalty
#' parameters that will be considered. See the description of `lambda.max` for
#' details. Default is 100.
#' @param q (Optional.) Numeric; determines what `L_q` penalty is used for the
#' fusion regularization. `q` = 1 is the lasso, and for 0 < `q` < 1, it is 
#' possible to get standard errors and confidence intervals. `q` = 2 is ridge
#' regression. See Faletto (2024) for details. Default is 0.5.
#' @param verbose Logical; if TRUE, more details on the progress of the function will
#' be printed as the function executes. Default is FALSE.
#' @param alpha Numeric; function will calculate (1 - `alpha`) confidence intervals
#' for the cohort average treatment effects that will be returned in `catt_df`.
#' @param add_ridge (Optional.) Logical; if TRUE, adds a small amount of ridge
#' regularization to the (untransformed) coefficients to stabilize estimation.
#' Default is FALSE.
#' @return A named list with the following elements: \item{att_hat}{The
#' estimated overall average treatment effect for a randomly selected treated
#' unit.} \item{att_se}{If `q < 1`, a standard error for the ATT. If
#' `indep_counts` was provided, this standard error is asymptotically exact; if
#' not, it is asymptotically conservative. If `q >= 1`, this will be NA.}
#' \item{catt_hats}{A named vector containing the estimated average treatment
#' effects for each cohort.} \item{catt_ses}{If `q < 1`, a named vector
#' containing the (asymptotically exact, non-conservative) standard errors for
#' the estimated average treatment effects within each cohort.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each
#' cohort conditional on being treated, which was used in calculating `att_hat`.
#' If `indep_counts` was provided, `cohort_probs` was calculated from that;
#' otherwise, it was calculated from the counts of units in each treated
#' cohort in `pdata`.} \item{catt_df}{A dataframe displaying the cohort names,
#' average treatment effects, standard errors, and `1 - alpha` confidence
#' interval bounds.} \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to
#' the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the
#' interactions between the treatment effects for each cohort at each time and
#' the covariates.} \item{sig_eps_sq}{Either the provided `sig_eps_sq` or
#' the estimated one, if a value wasn't provided.} \item{sig_eps_c_sq}{Either
#' the provided `sig_eps_c_sq` or the estimated one, if a value wasn't
#' provided.} \item{lambda.max}{Either the provided `lambda.max` or the one
#' that was used, if a value wasn't provided. (This is returned to help with
#' getting a reasonable range of `lambda` values for grid search.)}
#' \item{lambda.max_model_size}{The size of the selected model corresponding
#' `lambda.max` (for `q <= 1`, this will be the smallest model size). As
#' mentioned above, for `q <= 1` ideally this value is close to 0.}
#' \item{lambda.min}{Either the provided `lambda.min` or the one
#' that was used, if a value wasn't provided.} \item{lambda.min_model_size}{The
#' size of the selected model corresponding `lambda.min` (for `q <= 1`, this
#' will be the largest model size). As mentioned above, for `q <= 1` ideally
#' this value is close to `p`.}\item{lambda_star}{The value of `lambda` chosen
#' by BIC. If this value is close to `lambda.min` or `lambda.max`, that could
#' suggest that the range of `lambda` values should be expanded.}
#' \item{lambda_star_model_size}{The size of the model that was selected. If
#' this value is close to `lambda.max_model_size` or `lambda.min_model_size`,
#' That could suggest that the range of `lambda` values should be expanded.}
#' \item{X_ints}{The design matrix created containing all
#' interactions, time and cohort dummies, etc.} \item{y}{The vector of 
#' responses, containing `nrow(X_ints)` entries.} \item{X_final}{The design
#' matrix after applying the change in coordinates to fit the model and also
#' multiplying on the left by the square root inverse of the estimated
#' covariance matrix for each unit.} \item{y_final}{The final response after
#' multiplying on the left by the square root inverse of the estimated
#' covariance matrix for each unit.} \item{N}{The final number of units that
#' were in the  data set used for estimation (after any units may have been
#' removed because they were treated in the first time period).} \item{T}{The 
#' number of time periods in the final data set.} \item{R}{The final number of
#' treated cohorts that appear in the final data set.} \item{d}{The final number
#' of covariates that appear in the final data set (after any covariates may
#' have been removed because they contained missing values or all contained the
#' same value for every unit).} \item{p}{The final number of columns in the full
#' set of covariates used to estimate the model.}
#'
#' @details
#' The function extracts the panel data and simulation parameters from the \code{FETWFE_simulated} object
#' and passes them to \code{fetwfe()}. The estimation output is then returned without further modification.
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)
#'
#'   result <- fetwfeWithSimulatedData(sim_data, verbose = TRUE)
#' }
#'
#' @export
fetwfeWithSimulatedData <- function(simulated_obj, lambda.max = NA, lambda.min = NA, nlambda = 100,
                                    q = 0.5, verbose = FALSE, alpha = 0.05, add_ridge = FALSE) {
  if (!inherits(simulated_obj, "FETWFE_simulated")) {
    stop("simulated_obj must be an object of class 'FETWFE_simulated'")
  }
  
  pdata <- simulated_obj$pdata
  time_var <- simulated_obj$time_var
  unit_var <- simulated_obj$unit_var
  treatment <- simulated_obj$treatment
  response <- simulated_obj$response
  covs <- simulated_obj$covs
  sig_eps_sq <- simulated_obj$sig_eps_sq
  sig_eps_c_sq <- simulated_obj$sig_eps_c_sq
  indep_counts <- simulated_obj$indep_counts
  
  res <- fetwfe(
    pdata = pdata,
    time_var = time_var,
    unit_var = unit_var,
    treatment = treatment,
    response = response,
    covs = covs,
    indep_counts = indep_counts,
    sig_eps_sq = sig_eps_sq,
    sig_eps_c_sq = sig_eps_c_sq,
    lambda.max = lambda.max,
    lambda.min = lambda.min,
    nlambda = nlambda,
    q = q,
    verbose = verbose,
    alpha = alpha,
    add_ridge = add_ridge
  )
  
  return(res)
}


#' Compute True Treatment Effects
#'
#' @description 
#' This function extracts the true treatment effects from a full coefficient vector
#' obtained by fitting a fused extended two-way fixed effects model. It is designed
#' specifically for coefficients formatted by \code{genCoefs()} and/or coefficients
#' estimted on data generated by \code{genRandomData()}. It calculates the 
#' overall average treatment effect on the treated (ATT) as the equal-weighted average 
#' of the cohort-specific treatment effects, and also returns the individual treatment 
#' effects for each treated cohort.
#'
#' @param beta Numeric vector. The full vector of estimated coefficients from the model.
#'   Its length must equal the total number of parameters specified by \code{getP(R, T, d, num_treats)}.
#' @param R Integer. The number of treated cohorts in the data. Must be at least 2.
#' @param T Integer. The number of time periods in the panel data.
#' @param d Integer. The number of covariates used in the model. Set \code{d = 0} if no covariates are used.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{att_true}{A numeric value representing the overall average treatment effect on the treated.
#'         It is computed as the (equal-weighted) mean of the cohort-specific treatment effects.}
#'   \item{actual_cohort_tes}{A numeric vector containing the true cohort-specific treatment effects,
#'         calculated by averaging the coefficients corresponding to the treatment dummies for each cohort.}
#' }
#'
#' @details 
#' The function internally uses auxiliary routines \code{getNumTreats()}, \code{getP()},
#' \code{getFirstInds()}, \code{getTreatInds()}, and \code{getActualCohortTes()} to determine the 
#' correct indices of treatment effect coefficients in \code{beta}. The overall treatment effect 
#' is computed as the simple average of these cohort-specific effects.
#'
#' @examples 
#' \dontrun{
#' # Suppose you have a fitted coefficient vector from the model:
#' beta <- your_model$beta
#'
#' # Define the number of treated cohorts, time periods, and covariates:
#' R <- 3   # e.g., three cohorts
#' T <- 5   # e.g., five time periods
#' d <- 2   # e.g., two covariates
#'
#' # Compute the true treatment effects:
#' te_results <- getTes(beta, R, T, d)
#'
#' # Overall average treatment effect on the treated:
#' print(te_results$att_true)
#'
#' # Cohort-specific treatment effects:
#' print(te_results$actual_cohort_tes)
#' }
#'
#' @export
getTes <- function(beta, R, T, d){

    num_treats <- getNumTreats(R=R, T=T)

    p <- getP(R=R, T=T, d=d, num_treats=num_treats)

    stopifnot(length(beta) == p)

    first_inds <- getFirstInds(R=R, T=T)
    treat_inds <- getTreatInds(R=R, T=T, d=d, num_treats=num_treats)

    actual_cohort_tes <- getActualCohortTes(
        R=R,
        first_inds=first_inds,
        treat_inds=treat_inds,
        coefs=beta,
        num_treats=num_treats
        )

    att_true <- as.numeric(mean(actual_cohort_tes))

    return(list(att_true=att_true, actual_cohort_tes=actual_cohort_tes))
}

