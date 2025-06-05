#' @title Bridge-penalized extended two-way fixed effects
#'
#' @description Implementation of extended two-way fixed effects with a bridge
#' penalty. Estimates overall ATT as well as CATT (cohort average treatment
#' effects on the treated units).
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
#' regularization. `q` = 1 is the lasso, and for 0 < `q` < 1, it is
#' possible to get standard errors and confidence intervals. `q` = 2 is ridge
#' regression. See Faletto (2025) for details. Default is 0.5.
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
#' size of the selected model corresponding to `lambda.min` (for `q <= 1`, this
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
#' res <- betwfe(
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
betwfe <- function(
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

	res <- betwfe_core(
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
	return(list(
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
		X_ints = res$X_ints,
		y = res$y,
		X_final = res$X_final,
		y_final = res$y_final,
		N = res$N,
		T = res$T,
		R = res$R,
		d = res$d,
		p = res$p,
		calc_ses = res$calc_ses
	))
}

#' Run BETWFE on Simulated Data
#'
#' @description
#' This function runs the bridge-penalized extended two-way fixed effects estimator (\code{betwfe()}) on
#' simulated data. It is simply a wrapper for \code{betwfe()}: it accepts an object of class
#' \code{"FETWFE_simulated"} (produced by \code{simulateData()}) and unpacks the necessary
#' components to pass to \code{betwfe()}. So the outputs match \code{betwfe()}, and the needed inputs
#' match their counterparts in \code{betwfe()}.
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
#' size of the selected model corresponding to `lambda.min` (for `q <= 1`, this
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
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
#'
#'   result <- betwfeWithSimulatedData(sim_data)
#' }
#'
#' @export
betwfeWithSimulatedData <- function(
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

	res <- betwfe(
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


#' Core Estimation Logic for Bridge-Penalized Extended Two-Way Fixed Effects
#'
#' @description
#' This function implements the core estimation steps of the BETWFE methodology.
#' It takes a pre-processed design matrix and response, handles variance
#' components, performs bridge regression, selects the optimal penalty via BIC,
#' and calculates treatment effects and their standard errors.
#'
#' @param X_ints The design matrix with all fixed effects, covariates, treatment
#'   dummies, and their interactions, as produced by `prepXints`.
#' @param y The centered response vector, as produced by `prepXints`.
#' @param in_sample_counts An integer vector named with cohort identifiers
#'   (including "Never_treated"), indicating the number of units in each cohort
#'   within the data used for estimation.
#' @param N The number of unique units.
#' @param T The number of unique time periods.
#' @param d The number of covariates.
#' @param p The total number of columns in `X_ints` (total parameters).
#' @param num_treats The total number of unique treatment effect parameters.
#' @param first_inds A numeric vector indicating the starting column index for
#'   each cohort's first treatment effect within the treatment effect block.
#' @param indep_counts (Optional) An integer vector of counts for how many units
#'   appear in the untreated cohort plus each of the other `R` cohorts, derived
#'   from an independent dataset. Used for asymptotically exact standard errors for
#'   the ATT. Default is `NA`.
#' @param sig_eps_sq (Optional) Numeric; the known variance of the observation-level
#'   IID noise. If `NA`, it will be estimated. Default is `NA`.
#' @param sig_eps_c_sq (Optional) Numeric; the known variance of the unit-level IID
#'   noise (random effects). If `NA`, it will be estimated. Default is `NA`.
#' @param lambda.max (Optional) Numeric; the maximum `lambda` penalty parameter for
#'   the bridge regression grid search. If `NA`, `grpreg` selects it. Default is `NA`.
#' @param lambda.min (Optional) Numeric; the minimum `lambda` penalty parameter.
#'   If `NA`, `grpreg` selects it. Default is `NA`.
#' @param nlambda (Optional) Integer; the number of `lambda` values in the grid.
#'   Default is 100.
#' @param q (Optional) Numeric; the power of the Lq penalty for fusion regularization
#'   (0 < q <= 2). `q=0.5` is default, `q=1` is lasso, `q=2` is ridge.
#'   Default is 0.5.
#' @param verbose Logical; if `TRUE`, prints progress messages. Default is `FALSE`.
#' @param alpha Numeric; significance level for confidence intervals (e.g., 0.05 for
#'   95% CIs). Default is 0.05.
#' @param add_ridge (Optional) Logical; if `TRUE`, adds a small L2 penalty to
#'   the coefficients to stabilize estimation. Default is `FALSE`.
#'
#' @details
#' The function executes the following main steps:
#' \enumerate{
#'   \item **Input Checks:** Validates the provided parameters.
#'   \item **Variance Component Handling:**
#'     \itemize{
#'       \item If `sig_eps_sq` or `sig_eps_c_sq` are `NA`, `estOmegaSqrtInv` is
#'         called to estimate them from the data using a fixed-effects ridge
#'         regression.
#'       \item Constructs the covariance matrix `Omega` and its inverse square
#'         root `Omega_sqrt_inv`.
#'       \item Pre-multiplies `y` and `X_mod` by `sqrt(sig_eps_sq) * Omega_sqrt_inv`
#'         (via Kronecker product) to obtain `y_final` and `X_final`, effectively
#'         performing a GLS transformation.
#'     }
#'   \item **Optional Ridge Penalty:** If `add_ridge` is `TRUE`, `X_final_scaled`
#'     (scaled version of `X_final`) and `y_final` are augmented to add an L2
#'     penalty.
#'   \item **Cohort Probabilities:** Calculates cohort membership probabilities
#'     conditional on being treated, using `in_sample_counts` and `indep_counts`
#'     if available.
#'   \item **Bridge Regression:** Fits a bridge regression model using
#'     `grpreg::gBridge` on `X_final_scaled` and `y_final` with the specified `q`
#'     and lambda sequence.
#'   \item **Coefficient Selection (BIC):** Calls `getBetaBIC` to select the
#'     optimal `lambda` using BIC and retrieve the corresponding estimated
#'     coefficients.
#'   \item **Handle Zero-Feature Case:** If BIC selects a model with zero features,
#'     treatment effects are set to zero.
#'   \item **Treatment Effect Calculation:**
#'     \itemize{
#'       \item Extracts cohort-specific average treatment effects (CATTs) from
#'         `beta_hat`.
#'       \item Calls `getCohortATTsFinal` to calculate CATT point estimates,
#'         standard errors (if `q < 1`), and confidence intervals. This involves
#'         computing the Gram matrix and related quantities.
#'     }
#'   \item **Overall ATT Calculation:** Calls `getTeResultsOLS` to calculate the
#'     overall average treatment effect on the treated (ATT) and its standard
#'     error, using both in-sample probabilities and independent probabilities
#'     if `indep_counts` were provided.
#' }
#' The standard errors for CATTs are asymptotically exact. For ATT, if
#' `indep_counts` are provided, the SE is asymptotically exact; otherwise, it's
#' asymptotically conservative (if `q < 1`).
#'
#' @return A list containing detailed estimation results:
#'   \item{in_sample_att_hat}{Estimated overall ATT using in-sample cohort probabilities.}
#'   \item{in_sample_att_se}{Standard error for `in_sample_att_hat`.}
#'   \item{in_sample_att_se_no_prob}{SE for `in_sample_att_hat` ignoring variability from estimating cohort probabilities.}
#'   \item{indep_att_hat}{Estimated overall ATT using `indep_counts` cohort probabilities (NA if `indep_counts` not provided).}
#'   \item{indep_att_se}{Standard error for `indep_att_hat` (NA if not applicable).}
#'   \item{catt_hats}{A named vector of estimated CATTs for each cohort.}
#'   \item{catt_ses}{A named vector of SEs for `catt_hats` (NA if `q >= 1`).}
#'   \item{catt_df}{A data.frame summarizing CATTs, SEs, and confidence intervals.}
#'   \item{beta_hat}{The vector of estimated coefficients.}
#'   \item{treat_inds}{Indices in `beta_hat` corresponding to base treatment effects.}
#'   \item{treat_int_inds}{Indices in `beta_hat` corresponding to treatment-covariate interactions.}
#'   \item{cohort_probs}{Estimated cohort probabilities conditional on being treated, from `in_sample_counts`.}
#'   \item{indep_cohort_probs}{Estimated cohort probabilities from `indep_counts` (NA if not provided).}
#'   \item{sig_eps_sq}{The (possibly estimated) variance of observation-level noise.}
#'   \item{sig_eps_c_sq}{The (possibly estimated) variance of unit-level random effects.}
#'   \item{lambda.max}{The maximum lambda value used in `grpreg`.}
#'   \item{lambda.max_model_size}{Model size for `lambda.max`.}
#'   \item{lambda.min}{The minimum lambda value used in `grpreg`.}
#'   \item{lambda.min_model_size}{Model size for `lambda.min`.}
#'   \item{lambda_star}{The lambda value selected by BIC.}
#'   \item{lambda_star_model_size}{Model size for `lambda_star`.}
#'   \item{X_ints}{The original input design matrix from `prepXints`.}
#'   \item{y}{The original input centered response vector from `prepXints`.}
#'   \item{X_final}{The design matrix after GLS weighting.}
#'   \item{y_final}{The response vector after GLS weighting.}
#'   \item{N, T, R, d, p}{Dimensions used in estimation.}
#' @keywords internal
#' @noRd
betwfe_core <- function(
	X_ints,
	y,
	in_sample_counts,
	N,
	T,
	d,
	p,
	num_treats,
	first_inds,
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
	ret <- check_etwfe_core_inputs(
		in_sample_counts = in_sample_counts,
		N = N,
		T = T,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		indep_counts = indep_counts,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge
	)

	R <- ret$R
	c_names <- ret$c_names
	indep_count_data_available <- ret$indep_count_data_available

	rm(ret)

	if (any(!is.na(lambda.max))) {
		stopifnot(is.numeric(lambda.max) | is.integer(lambda.max))
		stopifnot(length(lambda.max) == 1)
		stopifnot(lambda.max >= 0)
	}

	if (any(!is.na(lambda.min))) {
		stopifnot(is.numeric(lambda.min) | is.integer(lambda.min))
		stopifnot(length(lambda.min) == 1)
		stopifnot(lambda.min >= 0)
		if (any(!is.na(lambda.max))) {
			stopifnot(lambda.max >= lambda.min)
		}
	}

	stopifnot(is.numeric(q) | is.integer(q))
	stopifnot(length(q) == 1)
	stopifnot(q > 0)
	stopifnot(q <= 2)

	res <- prep_for_etwfe_regresion(
		verbose = verbose,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		y = y,
		X_ints = X_ints,
		X_mod = X_ints, # Don't transform matrix
		N = N,
		T = T,
		R = R,
		d = d,
		p = p,
		num_treats = num_treats,
		add_ridge = add_ridge,
		first_inds = first_inds,
		in_sample_counts = in_sample_counts,
		indep_count_data_available = indep_count_data_available,
		indep_counts = indep_counts
	)

	X_final_scaled <- res$X_final_scaled
	y_final <- res$y_final
	scale_center <- res$scale_center
	scale_scale <- res$scale_scale
	cohort_probs <- res$cohort_probs
	cohort_probs_overall <- res$cohort_probs_overall
	indep_cohort_probs <- res$indep_cohort_probs
	indep_cohort_probs_overall <- res$indep_cohort_probs_overall
	X_final <- res$X_final
	lambda_ridge <- res$lambda_ridge
	sig_eps_sq <- res$sig_eps_sq
	sig_eps_c_sq <- res$sig_eps_c_sq

	rm(res)

	#
	#
	# Step 4: estimate bridge regression and extract fitted coefficients
	#
	#

	# Estimate bridge regression
	if (verbose) {
		message("Estimating bridge regression...")
		t0 <- Sys.time()
	}

	if (!is.na(lambda.max) & !is.na(lambda.min)) {
		fit <- grpreg::gBridge(
			X = X_final_scaled,
			y = y_final,
			gamma = q,
			lambda.max = lambda.max,
			lambda.min = lambda.min,
			nlambda = nlambda
		)
	} else if (!is.na(lambda.max)) {
		fit <- grpreg::gBridge(
			X = X_final_scaled,
			y = y_final,
			gamma = q,
			lambda.max = lambda.max,
			nlambda = nlambda
		)
	} else if (!is.na(lambda.min)) {
		fit <- grpreg::gBridge(
			X = X_final_scaled,
			y = y_final,
			gamma = q,
			lambda.min = lambda.min,
			nlambda = nlambda
		)
	} else {
		fit <- grpreg::gBridge(
			X = X_final_scaled,
			y = y_final,
			gamma = q,
			nlambda = nlambda
		)
	}

	if (verbose) {
		message("Done! Time for estimation:")
		message(Sys.time() - t0)
	}

	# For diagnostics later, store largest and smallest lambda, as well as
	# corresponding smallest and largest model sizes, to return.
	lambda.max <- max(fit$lambda)
	lambda.max_model_size <- sum(fit$beta[, ncol(fit$beta)] != 0)

	lambda.min <- min(fit$lambda)
	lambda.min_model_size <- sum(fit$beta[, 1] != 0)

	# Select a single set of fitted coefficients by using BIC to choose among
	# the penalties that were fitted
	res <- getBetaBIC(
		fit = fit,
		N = N,
		T = T,
		p = p,
		X_mod = X_ints, # Pass untransformed matrix
		y = y,
		scale_center = scale_center,
		scale_scale = scale_scale
	)

	beta_hat <- res$theta_hat # This includes intercept
	lambda_star_ind <- res$lambda_star_ind
	lambda_star_model_size <- res$lambda_star_model_size

	lambda_star <- fit$lambda[lambda_star_ind]

	# c_names <- names(in_sample_counts)[2:(R + 1)] # Moved definition up
	stopifnot(length(c_names) == R)

	# Indices corresponding to base treatment effects
	treat_inds <- getTreatInds(R = R, T = T, d = d, num_treats = num_treats)

	if (d > 0) {
		stopifnot(max(treat_inds) + 1 <= p)
		stopifnot(
			max(treat_inds) == R + T - 1 + d + R * d + (T - 1) * d + num_treats
		)

		treat_int_inds <- (max(treat_inds) + 1):p

		stopifnot(length(treat_int_inds) == num_treats * d)
	} else {
		stopifnot(max(treat_inds) <= p)
		stopifnot(max(treat_inds) == R + T - 1 + num_treats)

		treat_int_inds <- c()
	}

	# Handle edge case where no features are selected (model_size includes intercept)
	if (lambda_star_model_size <= 1 && all(beta_hat[2:(p + 1)] == 0)) {
		# only intercept might be non-zero
		if (verbose) {
			message(
				"No features selected (or only intercept); all treatment effects estimated to be 0."
			)
		}

		if (q < 1) {
			ret_se <- 0
		} else {
			ret_se <- NA
		}

		catt_df_to_ret <- data.frame(
			Cohort = c_names,
			`Estimated TE` = rep(0, R),
			SE = rep(ret_se, R),
			ConfIntLow = rep(ret_se, R),
			ConfIntHigh = rep(ret_se, R),
			check.names = FALSE
		)

		return(list(
			in_sample_att_hat = 0,
			in_sample_att_se = ret_se,
			in_sample_att_se_no_prob = ret_se,
			indep_att_hat = 0,
			indep_att_se = ret_se,
			catt_hats = setNames(rep(0, R), c_names),
			catt_ses = setNames(rep(ret_se, R), c_names),
			catt_df = catt_df_to_ret,
			beta_hat = beta_hat[2:(p + 1)], # Slopes are all zero
			treat_inds = treat_inds,
			treat_int_inds = treat_int_inds,
			cohort_probs = cohort_probs,
			indep_cohort_probs = indep_cohort_probs,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
			lambda.max = lambda.max,
			lambda.max_model_size = lambda.max_model_size,
			lambda.min = lambda.min,
			lambda.min_model_size = lambda.min_model_size,
			lambda_star = lambda_star,
			lambda_star_model_size = lambda_star_model_size,
			X_ints = X_ints,
			y = y,
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T,
			R = R,
			d = d,
			p = p,
			calc_ses = q < 1
		))
	}

	# estimated coefficients
	beta_hat_slopes <- beta_hat[2:(p + 1)]

	# Indices of selected features
	sel_feat_inds <- which(beta_hat_slopes != 0)

	sel_treat_inds <- sel_feat_inds[sel_feat_inds %in% treat_inds]

	stopifnot(length(sel_treat_inds) == length(unique(sel_treat_inds)))
	stopifnot(length(sel_treat_inds) <= length(sel_feat_inds))
	stopifnot(length(sel_treat_inds) <= length(treat_inds))
	stopifnot(is.integer(sel_treat_inds) | is.numeric(sel_treat_inds))

	# 1. Get beta_hat_slopes[treat_inds] -> these are the treatment coefficients
	# 2. Find which of these are non-zero: `which(beta_hat_slopes[treat_inds] != 0)` -> these are
	# indices *within* the treat_inds block. Let's call this `sel_treat_inds_relative_to_block`.
	# `sel_treat_inds` itself contains the global indices of selected treatment features.
	# So `beta_hat_slopes[sel_treat_inds]` are the non-zero treatment coefs.

	beta_hat_treat_block = beta_hat_slopes[treat_inds]
	sel_treat_inds_shifted <- which(beta_hat_treat_block != 0) # these are 1 to num_treats

	stopifnot(all(sel_treat_inds_shifted >= 1))
	stopifnot(all(sel_treat_inds_shifted <= num_treats))

	# Handle edge case where no treatment features selected
	if (length(sel_treat_inds_shifted) == 0) {
		if (verbose) {
			message(
				"No treatment features selected; all treatment effects estimated to be 0."
			)
		}

		if (q < 1) {
			ret_se <- 0
		} else {
			ret_se <- NA
		}

		catt_df_to_ret <- data.frame(
			Cohort = c_names,
			`Estimated TE` = rep(0, R),
			SE = rep(ret_se, R),
			ConfIntLow = rep(ret_se, R),
			ConfIntHigh = rep(ret_se, R),
			check.names = FALSE
		)

		if (add_ridge) {
			lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)
			beta_hat_slopes <- beta_hat_slopes * (1 + lambda_ridge)
		}

		return(list(
			in_sample_att_hat = 0,
			in_sample_att_se = ret_se,
			in_sample_att_se_no_prob = ret_se,
			indep_att_hat = 0,
			indep_att_se = ret_se,
			catt_hats = setNames(rep(0, R), c_names),
			catt_ses = setNames(rep(ret_se, R), c_names),
			catt_df = catt_df_to_ret,
			beta_hat = beta_hat_slopes, # Untransformed slopes
			treat_inds = treat_inds,
			treat_int_inds = treat_int_inds,
			cohort_probs = cohort_probs,
			indep_cohort_probs = indep_cohort_probs,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
			lambda.max = lambda.max,
			lambda.max_model_size = lambda.max_model_size,
			lambda.min = lambda.min,
			lambda.min_model_size = lambda.min_model_size,
			lambda_star = lambda_star,
			lambda_star_model_size = lambda_star_model_size,
			X_ints = X_ints,
			y = y,
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T,
			R = R,
			d = d,
			p = p,
			calc_ses = q < 1
		))
	}

	beta_hat <- beta_hat_slopes

	# If using ridge regularization, multiply the "naive" estimated coefficients
	# by 1 + lambda_ridge, similar to suggestion in original elastic net paper.
	if (add_ridge) {
		lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)

		beta_hat <- beta_hat * (1 + lambda_ridge)
	}

	# Get actual estimated treatment effects
	tes <- beta_hat[treat_inds]

	stopifnot(length(tes) <= num_treats)

	stopifnot(all(beta_hat_slopes[treat_inds][sel_treat_inds_shifted] != 0))
	stopifnot(all(
		beta_hat_slopes[treat_inds][setdiff(
			1:num_treats,
			sel_treat_inds_shifted
		)] ==
			0
	))

	stopifnot(length(first_inds) == R)
	stopifnot(max(first_inds) <= num_treats)

	stopifnot(length(sel_feat_inds) > 0) # sel_feat_inds are indices in beta_hat_slopes
	stopifnot(length(sel_treat_inds_shifted) > 0) # sel_treat_inds_shifted are indices within the treat_inds block of beta_hat_slopes

	#
	#
	# Step 6: calculate cohort-specific treatment effects and standard
	# errors
	#
	#

	res <- getCohortATTsFinal(
		X_final = X_final, # This is X_mod * GLS_transform_matrix
		sel_feat_inds = sel_feat_inds, # Indices of non-zero elements in beta_hat_slopes
		treat_inds = treat_inds, # Global indices for treatment effects
		num_treats = num_treats,
		first_inds = first_inds,
		sel_treat_inds_shifted = sel_treat_inds_shifted, # Indices (1 to num_treats) of non-zero treat. coefs.
		c_names = c_names,
		tes = tes, # Treatment effect estimates (beta_hat[treat_inds])
		sig_eps_sq = sig_eps_sq,
		R = R,
		N = N,
		T = T,
		fused = FALSE,
		calc_ses = q < 1,
		p = p, # Total number of original parameters (columns in X_ints)
		alpha = alpha
	)

	cohort_te_df <- res$cohort_te_df
	cohort_tes <- res$cohort_tes
	cohort_te_ses <- res$cohort_te_ses
	psi_mat <- res$psi_mat
	gram_inv <- res$gram_inv
	calc_ses <- res$calc_ses

	rm(res)

	#
	#
	# Step 7: calculate overall average treatment effect on treated units
	#
	#

	# Get overal estimated ATT!
	in_sample_te_results <- getTeResultsOLS(
		sig_eps_sq = sig_eps_sq,
		N = N,
		T = T,
		R = R,
		num_treats = num_treats,
		cohort_tes = cohort_tes, # CATTs (point estimates)
		cohort_probs = cohort_probs, # In-sample pi_r | treated
		psi_mat = psi_mat,
		gram_inv = gram_inv,
		tes = tes[sel_treat_inds_shifted], # Untransformed treatment effect estimates beta_hat[treat_inds]
		cohort_probs_overall = cohort_probs_overall, # In-sample pi_r (unconditional on treated)
		first_inds = first_inds,
		calc_ses = calc_ses,
		indep_probs = FALSE
	)

	in_sample_att_hat <- in_sample_te_results$att_hat
	in_sample_att_se <- in_sample_te_results$att_te_se
	in_sample_att_se_no_prob <- in_sample_te_results$att_te_se_no_prob

	if ((q < 1) & calc_ses) {
		stopifnot(!is.na(in_sample_att_se))
	}

	if (indep_count_data_available) {
		stopifnot(nrow(psi_mat) == length(tes[sel_treat_inds_shifted]))

		indep_te_results <- getTeResultsOLS(
			sig_eps_sq = sig_eps_sq,
			N = N,
			T = T,
			R = R,
			num_treats = num_treats,
			cohort_tes = cohort_tes,
			cohort_probs = indep_cohort_probs, # indep pi_r | treated
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			tes = tes[sel_treat_inds_shifted],
			cohort_probs_overall = indep_cohort_probs_overall, # indep pi_r (unconditional)
			first_inds = first_inds,
			calc_ses = calc_ses,
			indep_probs = TRUE
		)

		indep_att_hat <- indep_te_results$att_hat
		indep_att_se <- indep_te_results$att_te_se
	} else {
		indep_att_hat <- NA
		indep_att_se <- NA
	}

	return(list(
		in_sample_att_hat = in_sample_att_hat,
		in_sample_att_se = in_sample_att_se,
		in_sample_att_se_no_prob = in_sample_att_se_no_prob,
		indep_att_hat = indep_att_hat,
		indep_att_se = indep_att_se,
		catt_hats = cohort_tes, # Already named if applicable from getCohortATTsFinal
		catt_ses = cohort_te_ses, # Already named if applicable
		catt_df = cohort_te_df,
		beta_hat = beta_hat, # Untransformed slopes
		treat_inds = treat_inds,
		treat_int_inds = treat_int_inds,
		cohort_probs = cohort_probs,
		indep_cohort_probs = indep_cohort_probs,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		lambda.max = lambda.max,
		lambda.max_model_size = lambda.max_model_size,
		lambda.min = lambda.min,
		lambda.min_model_size = lambda.min_model_size,
		lambda_star = lambda_star,
		lambda_star_model_size = lambda_star_model_size,
		X_ints = X_ints,
		y = y,
		X_final = X_final,
		y_final = y_final,
		N = N,
		T = T,
		R = R,
		d = d,
		p = p,
		calc_ses = calc_ses
	))
}
