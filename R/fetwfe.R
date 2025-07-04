#' @import glmnet
#' @importFrom stats rbinom rmultinom rnorm runif sd lm coef pnorm

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
#' columns for covariates. All of these columns are expected to contain integer,
#' numeric, or factor values, and any categorical values will be automatically
#' encoded as binary indicators. If no covariates are provided, the treatment
#' effect estimation will proceed, but it will only be valid under unconditional
#' versions of the parallel trends and no anticipation assumptions. Default is c().
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
#' noise assumed to apply to each observation. See Section 2 of Faletto (2025)
#' for details. It is best to provide this variance if it is known (for example,
#' if you are using simulated data). If this variance is unknown, this argument
#' can be omitted, and the variance will be estimated using the estimator from
#' Pesaran (2015, Section 26.5.1) with ridge regression. Default is NA.
#' @param sig_eps_c_sq (Optional.) Numeric; the variance of the unit-level IID
#' noise (random effects) assumed to apply to each observation. See Section 2 of
#' Faletto (2025) for details. It is best to provide this variance if it is
#' known (for example, if you are using simulated data). If this variance is
#' unknown, this argument can be omitted, and the variance will be estimated
#' using the estimator from Pesaran (2015, Section 26.5.1) with ridge
#' regression. Default is NA.
#' @param lambda.max (Optional.) Numeric. A penalty parameter `lambda` will be
#' selected over a grid search by BIC in order to select a single model. The
#' largest `lambda` in the grid will be `lambda.max`. If no `lambda.max` is
#' provided, one will be selected automatically. When `q <= 1`, the model
#' will be sparse, and ideally all of the following are true at once: the
#' smallest model (the one corresponding to `lambda.max`) selects close to 0
#' features, the largest model (the one corresponding to `lambda.min`) selects
#' close to `p` features, `nlambda` is large enough so that models are
#' considered at every feasible model size, and `nlambda` is small enough so
#' that the computation doesn't become infeasible. You may
#' want to manually tweak `lambda.max`, `lambda.min`, and `nlambda` to try
#' to achieve these goals, particularly if the selected model size is very
#' close to the model corresponding to `lambda.max` or `lambda.min`, which could
#' indicate that the range of `lambda` values was too narrow or coarse. You can
#' use the function outputs `lambda.max_model_size`, `lambda.min_model_size`, and
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
#' regression. See Faletto (2025) for details. Default is 0.5.
#' @param verbose Logical; if TRUE, more details on the progress of the function will
#' be printed as the function executes. Default is FALSE.
#' @param alpha Numeric; function will calculate (1 - `alpha`) confidence intervals
#' for the cohort average treatment effects that will be returned in `catt_df`.
#' @param add_ridge (Optional.) Logical; if TRUE, adds a small amount of ridge
#' regularization to the (untransformed) coefficients to stabilize estimation.
#' Default is FALSE.
#' @return An object of class \code{fetwfe} containing the following elements:
#' \item{att_hat}{The estimated overall average treatment effect for a randomly selected treated unit.}
#' \item{att_se}{If `q < 1`, a standard error for the ATT. If `indep_counts` was provided, this standard error is asymptotically exact; if not, it is asymptotically conservative. If `q >= 1`, this will be NA.}
#' \item{catt_hats}{A named vector containing the estimated average treatment effects for each cohort.}
#' \item{catt_ses}{If `q < 1`, a named vector containing the (asymptotically exact, non-conservative) standard errors for the estimated average treatment effects within each cohort.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each cohort conditional on being treated, which was used in calculating `att_hat`. If `indep_counts` was provided, `cohort_probs` was calculated from that; otherwise, it was calculated from the counts of units in each treated cohort in `pdata`.}
#' \item{catt_df}{A dataframe displaying the cohort names, average treatment effects, standard errors, and `1 - alpha` confidence interval bounds.}
#' \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the interactions between the treatment effects for each cohort at each time and the covariates.}
#' \item{sig_eps_sq}{Either the provided `sig_eps_sq` or the estimated one, if a value wasn't provided.}
#' \item{sig_eps_c_sq}{Either the provided `sig_eps_c_sq` or the estimated one, if a value wasn't provided.}
#' \item{lambda.max}{Either the provided `lambda.max` or the one that was used, if a value wasn't provided. (This is returned to help with getting a reasonable range of `lambda` values for grid search.)}
#' \item{lambda.max_model_size}{The size of the selected model corresponding to `lambda.max` (for `q <= 1`, this will be the smallest model size). As mentioned above, for `q <= 1` ideally this value is close to 0.}
#' \item{lambda.min}{Either the provided `lambda.min` or the one that was used, if a value wasn't provided.}
#' \item{lambda.min_model_size}{The size of the selected model corresponding to `lambda.min` (for `q <= 1`, this will be the largest model size). As mentioned above, for `q <= 1` ideally this value is close to `p`.}
#' \item{lambda_star}{The value of `lambda` chosen by BIC. If this value is close to `lambda.min` or `lambda.max`, that could suggest that the range of `lambda` values should be expanded.}
#' \item{lambda_star_model_size}{The size of the model that was selected. If this value is close to `lambda.max_model_size` or `lambda.min_model_size`, that could suggest that the range of `lambda` values should be expanded.}
#' \item{N}{The final number of units that were in the data set used for estimation (after any units may have been removed because they were treated in the first time period).}
#' \item{T}{The number of time periods in the final data set.}
#' \item{R}{The final number of treated cohorts that appear in the final data set.}
#' \item{d}{The final number of covariates that appear in the final data set (after any covariates may have been removed because they contained missing values or all contained the same value for every unit).}
#' \item{p}{The final number of columns in the full set of covariates used to estimate the model.}
#' \item{alpha}{The alpha level used for confidence intervals.}
#' \item{internal}{A list containing internal outputs that are typically not needed for interpretation:
#'   \describe{
#'     \item{X_ints}{The design matrix created containing all interactions, time and cohort dummies, etc.}
#'     \item{y}{The vector of responses, containing `nrow(X_ints)` entries.}
#'     \item{X_final}{The design matrix after applying the change in coordinates to fit the model and also multiplying on the left by the square root inverse of the estimated covariance matrix for each unit.}
#'     \item{y_final}{The final response after multiplying on the left by the square root inverse of the estimated covariance matrix for each unit.}
#'     \item{calc_ses}{Logical indicating whether standard errors were calculated.}
#'   }
#' }
#'
#' The object has methods for \code{print()}, \code{summary()}, and \code{coef()}. By default, \code{print()} and \code{summary()} only show the essential outputs. To see internal details, use \code{print(x, show_internal = TRUE)} or \code{summary(x, show_internal = TRUE)}. The \code{coef()} method returns the vector of estimated coefficients (\code{beta_hat}).
#'
#' @author Gregory Faletto
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
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
#' # Print results with internal details
#' print(res, max_cohorts = Inf)
#'
#' @export
fetwfe <- function(
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
) {
	# Check inputs
	ret <- checkFetwfeInputs(
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

	pdata <- ret$pdata
	indep_count_data_available = ret$indep_count_data_available

	rm(ret)

	res1 <- prep_for_etwfe_core(
		pdata = pdata,
		response = response,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		covs = covs,
		verbose = verbose,
		indep_count_data_available = indep_count_data_available,
		indep_counts = indep_counts
	)

	pdata <- res1$pdata
	covs <- res1$covs
	X_ints <- res1$X_ints
	y <- res1$y
	N <- res1$N
	T <- res1$T
	d <- res1$d
	p <- res1$p
	in_sample_counts <- res1$in_sample_counts
	num_treats <- res1$num_treats
	first_inds <- res1$first_inds
	R <- res1$R

	rm(res1)

	res <- fetwfe_core(
		X_ints = X_ints,
		y = y,
		in_sample_counts = in_sample_counts,
		N = N,
		T = T,
		d = d,
		p = p,
		num_treats = num_treats,
		first_inds = first_inds,
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

	if (indep_count_data_available) {
		stopifnot(!is.na(res$indep_att_hat))

		if ((q < 1) & res$calc_ses) {
			stopifnot(!is.na(res$indep_att_se))
		}

		stopifnot(all(!is.na(res$indep_cohort_probs)))
		att_hat <- res$indep_att_hat
		att_se <- res$indep_att_se
		cohort_probs <- res$indep_cohort_probs
	} else {
		stopifnot(!is.na(res$in_sample_att_hat))

		if ((q < 1) & res$calc_ses) {
			stopifnot(!is.na(res$in_sample_att_se))
		}

		att_hat <- res$in_sample_att_hat
		att_se <- res$in_sample_att_se
		cohort_probs <- res$cohort_probs
	}

	# Create the main output list with essential results
	out <- list(
		att_hat = att_hat,
		att_se = att_se,
		catt_hats = res$catt_hats,
		catt_ses = res$catt_ses,
		cohort_probs = cohort_probs,
		catt_df = res$catt_df,
		beta_hat = res$beta_hat,
		treat_inds = res$treat_inds,
		treat_int_inds = res$treat_int_inds,
		sig_eps_sq = res$sig_eps_sq,
		sig_eps_c_sq = res$sig_eps_c_sq,
		lambda.max = res$lambda.max,
		lambda.max_model_size = res$lambda.max_model_size,
		lambda.min = res$lambda.min,
		lambda.min_model_size = res$lambda.min_model_size,
		lambda_star = res$lambda_star,
		lambda_star_model_size = res$lambda_star_model_size,
		N = res$N,
		T = res$T,
		R = res$R,
		d = res$d,
		p = res$p,
		alpha = alpha
	)

	# Add internal outputs in a separate list
	out$internal <- list(
		X_ints = res$X_ints,
		y = res$y,
		X_final = res$X_final,
		y_final = res$y_final,
		calc_ses = res$calc_ses
	)

	# Add the fetwfe class
	class(out) <- "fetwfe"

	return(out)
}


#' Run FETWFE on Simulated Data
#'
#' @description
#' This function runs the fused extended two-way fixed effects estimator (\code{fetwfe()}) on
#' simulated data. It is simply a wrapper for \code{fetwfe()}: it accepts an object of class
#' \code{"FETWFE_simulated"} (produced by \code{simulateData()}) and unpacks the necessary
#' components to pass to \code{fetwfe()}. So the outputs match \code{fetwfe()}, and the needed inputs
#' match their counterparts in \code{fetwfe()}.
#'
#' @param simulated_obj An object of class \code{"FETWFE_simulated"} containing the simulated panel
#' data and design matrix.
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
#' regression. See Faletto (2025) for details. Default is 0.5.
#' @param verbose Logical; if TRUE, more details on the progress of the function will
#' be printed as the function executes. Default is FALSE.
#' @param alpha Numeric; function will calculate (1 - `alpha`) confidence intervals
#' for the cohort average treatment effects that will be returned in `catt_df`.
#' @param add_ridge (Optional.) Logical; if TRUE, adds a small amount of ridge
#' regularization to the (untransformed) coefficients to stabilize estimation.
#' Default is FALSE.
#' @return An object of class \code{fetwfe} containing the following elements:
#' \item{att_hat}{The estimated overall average treatment effect for a randomly selected treated unit.}
#' \item{att_se}{If `q < 1`, a standard error for the ATT. If `indep_counts` was provided, this standard error is asymptotically exact; if not, it is asymptotically conservative. If `q >= 1`, this will be NA.}
#' \item{catt_hats}{A named vector containing the estimated average treatment effects for each cohort.}
#' \item{catt_ses}{If `q < 1`, a named vector containing the (asymptotically exact, non-conservative) standard errors for the estimated average treatment effects within each cohort.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each cohort conditional on being treated, which was used in calculating `att_hat`. If `indep_counts` was provided, `cohort_probs` was calculated from that; otherwise, it was calculated from the counts of units in each treated cohort in `pdata`.}
#' \item{catt_df}{A dataframe displaying the cohort names, average treatment effects, standard errors, and `1 - alpha` confidence interval bounds.}
#' \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the interactions between the treatment effects for each cohort at each time and the covariates.}
#' \item{sig_eps_sq}{Either the provided `sig_eps_sq` or the estimated one, if a value wasn't provided.}
#' \item{sig_eps_c_sq}{Either the provided `sig_eps_c_sq` or the estimated one, if a value wasn't provided.}
#' \item{lambda.max}{Either the provided `lambda.max` or the one that was used, if a value wasn't provided. (This is returned to help with getting a reasonable range of `lambda` values for grid search.)}
#' \item{lambda.max_model_size}{The size of the selected model corresponding to `lambda.max` (for `q <= 1`, this will be the smallest model size). As mentioned above, for `q <= 1` ideally this value is close to 0.}
#' \item{lambda.min}{Either the provided `lambda.min` or the one that was used, if a value wasn't provided.}
#' \item{lambda.min_model_size}{The size of the selected model corresponding to `lambda.min` (for `q <= 1`, this will be the largest model size). As mentioned above, for `q <= 1` ideally this value is close to `p`.}
#' \item{lambda_star}{The value of `lambda` chosen by BIC. If this value is close to `lambda.min` or `lambda.max`, that could suggest that the range of `lambda` values should be expanded.}
#' \item{lambda_star_model_size}{The size of the model that was selected. If this value is close to `lambda.max_model_size` or `lambda.min_model_size`, that could suggest that the range of `lambda` values should be expanded.}
#' \item{N}{The final number of units that were in the data set used for estimation (after any units may have been removed because they were treated in the first time period).}
#' \item{T}{The number of time periods in the final data set.}
#' \item{R}{The final number of treated cohorts that appear in the final data set.}
#' \item{d}{The final number of covariates that appear in the final data set (after any covariates may have been removed because they contained missing values or all contained the same value for every unit).}
#' \item{p}{The final number of columns in the full set of covariates used to estimate the model.}
#' \item{alpha}{The alpha level used for confidence intervals.}
#' \item{internal}{A list containing internal outputs that are typically not needed for interpretation:
#'   \describe{
#'     \item{X_ints}{The design matrix created containing all interactions, time and cohort dummies, etc.}
#'     \item{y}{The vector of responses, containing `nrow(X_ints)` entries.}
#'     \item{X_final}{The design matrix after applying the change in coordinates to fit the model and also multiplying on the left by the square root inverse of the estimated covariance matrix for each unit.}
#'     \item{y_final}{The final response after multiplying on the left by the square root inverse of the estimated covariance matrix for each unit.}
#'     \item{calc_ses}{Logical indicating whether standard errors were calculated.}
#'   }
#' }
#'
#' The object has methods for \code{print()}, \code{summary()}, and \code{coef()}. By default, \code{print()} and \code{summary()} only show the essential outputs. To see internal details, use \code{print(x, show_internal = TRUE)} or \code{summary(x, show_internal = TRUE)}. The \code{coef()} method returns the vector of estimated coefficients (\code{beta_hat}).
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
#'
#'   result <- fetwfeWithSimulatedData(sim_data)
#' }
#'
#' @export
fetwfeWithSimulatedData <- function(
	simulated_obj,
	lambda.max = NA,
	lambda.min = NA,
	nlambda = 100,
	q = 0.5,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE
) {
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


#' @title Extended two-way fixed effects
#'
#' @description Implementation of extended two-way fixed effects.
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
#' columns for covariates. All of these columns are expected to contain integer,
#' numeric, or factor values, and any categorical values will be automatically
#' encoded as binary indicators. If no covariates are provided, the treatment
#' effect estimation will proceed, but it will only be valid under unconditional
#' versions of the parallel trends and no anticipation assumptions. Default is c().
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
#' noise assumed to apply to each observation. See Section 2 of Faletto (2025)
#' for details. It is best to provide this variance if it is known (for example,
#' if you are using simulated data). If this variance is unknown, this argument
#' can be omitted, and the variance will be estimated using the estimator from
#' Pesaran (2015, Section 26.5.1) with ridge regression. Default is NA.
#' @param sig_eps_c_sq (Optional.) Numeric; the variance of the unit-level IID
#' noise (random effects) assumed to apply to each observation. See Section 2 of
#' Faletto (2025) for details. It is best to provide this variance if it is
#' known (for example, if you are using simulated data). If this variance is
#' unknown, this argument can be omitted, and the variance will be estimated
#' using the estimator from Pesaran (2015, Section 26.5.1) with ridge
#' regression. Default is NA.
#' @param verbose Logical; if TRUE, more details on the progress of the function will
#' be printed as the function executes. Default is FALSE.
#' @param alpha Numeric; function will calculate (1 - `alpha`) confidence intervals
#' for the cohort average treatment effects that will be returned in `catt_df`.
#' @param add_ridge (Optional.) Logical; if TRUE, adds a small amount of ridge
#' regularization to the (untransformed) coefficients to stabilize estimation.
#' Default is FALSE.
#' @return A named list with the following elements: \item{att_hat}{The
#' estimated overall average treatment effect for a randomly selected treated
#' unit.} \item{att_se}{A standard error for the ATT. If the Gram matrix is not
#' invertible, this will be NA.} \item{catt_hats}{A named vector containing the
#' estimated average treatment effects for each cohort.} \item{catt_ses}{A named
#' vector containing the (asymptotically exact) standard errors for
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
#' provided.} \item{X_ints}{The design matrix created containing all
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
#' Wooldridge, J. M. (2021). Two-way fixed effects, the two-way mundlak
#' regression, and difference-in-differences estimators.
#' \emph{Available at SSRN 3906345}.
#' \doi{10.2139/ssrn.3906345}.
#' @export
etwfe <- function(
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
) {
	# Check inputs
	ret <- checkEtwfeInputs(
		pdata = pdata,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		response = response,
		covs = covs,
		indep_counts = indep_counts,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge
	)

	pdata <- ret$pdata
	indep_count_data_available = ret$indep_count_data_available

	rm(ret)

	res1 <- prep_for_etwfe_core(
		pdata = pdata,
		response = response,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		covs = covs,
		verbose = verbose,
		indep_count_data_available = indep_count_data_available,
		indep_counts = indep_counts
	)

	pdata <- res1$pdata
	covs <- res1$covs
	X_ints <- res1$X_ints
	y <- res1$y
	N <- res1$N
	T <- res1$T
	d <- res1$d
	p <- res1$p
	in_sample_counts <- res1$in_sample_counts
	num_treats <- res1$num_treats
	first_inds <- res1$first_inds
	R <- res1$R

	rm(res1)

	warning_flag <- FALSE

	for (r in 1:(R + 1)) {
		if (in_sample_counts[r] < d + 1) {
			if (add_ridge) {
				warning_flag <- TRUE
			} else {
				stop(
					"At least one cohort contains fewer than d + 1 units. The design matrix is rank-deficient. Calculating standard errors will not be possible, and estimating treatment effects is only possible using add_ridge = TRUE."
				)
			}
		}
	}

	if (warning_flag) {
		warning(
			"At least one cohort contains fewer than d + 1 units. The design matrix is rank-deficient. Calculating standard errors will not be possible, and estimating treatment effects is only possible using add_ridge = TRUE."
		)
	}

	res <- etwfe_core(
		X_ints = X_ints,
		y = y,
		in_sample_counts = in_sample_counts,
		N = N,
		T = T,
		d = d,
		p = p,
		num_treats = num_treats,
		first_inds = first_inds,
		indep_counts = indep_counts,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge
	)

	if (indep_count_data_available) {
		stopifnot(!is.na(res$indep_att_hat))
		stopifnot(!is.na(res$indep_att_se))
		stopifnot(all(!is.na(res$indep_cohort_probs)))
		att_hat <- res$indep_att_hat
		att_se <- res$indep_att_se
		cohort_probs <- res$indep_cohort_probs
	} else {
		att_hat <- res$in_sample_att_hat
		att_se <- res$in_sample_att_se
		cohort_probs <- res$cohort_probs
	}

	# Create the main output list with essential results
	out <- list(
		att_hat = att_hat,
		att_se = att_se,
		catt_hats = res$catt_hats,
		catt_ses = res$catt_ses,
		cohort_probs = cohort_probs,
		catt_df = res$catt_df,
		beta_hat = res$beta_hat,
		treat_inds = res$treat_inds,
		treat_int_inds = res$treat_int_inds,
		sig_eps_sq = res$sig_eps_sq,
		sig_eps_c_sq = res$sig_eps_c_sq,
		X_ints = res$X_ints,
		y = res$y,
		X_final = res$X_final,
		y_final = res$y_final,
		N = res$N,
		T = res$T,
		R = res$R,
		d = res$d,
		p = res$p,
		alpha = alpha,
		calc_ses = res$calc_ses
	)

	# Add the etwfe class
	class(out) <- "etwfe"

	return(out)
}
