#' @title Two-way fixed effects with covariates and separate treatment effects
#' for each cohort
#'
#' @description **WARNING: This function should NOT be used for estimation. It
#' is a biased estimator of treatment effects.** Implementation of two-way fixed
#' effects with covariates and separate treatment effects for each cohort.
#' Estimates overall ATT as well as CATT (cohort average treatment effects on
#' the treated units). It is implemented only for the sake of the simulation
#' studies in Faletto (2025). This estimator is only unbiased under the
#' assumptions that treatment effects are homogeneous across covariates and are
#' identical within cohorts across all times since treatment.
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
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#' @export
twfeCovs <- function(
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
					"At least one cohort contains fewer than d + 1 units. The design matrix may be rank-deficient. Calculating standard errors may not be possible, and estimating treatment effects may only be possible using add_ridge = TRUE."
				)
			}
		}
	}

	if (warning_flag) {
		warning(
			"At least one cohort contains fewer than d + 1 units. The design matrix may be rank-deficient. Calculating standard errors may not be possible, and estimating treatment effects may only be possible using add_ridge = TRUE."
		)
	}

	res <- twfeCovs_core(
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


#' Run twfeCovs on Simulated Data
#'
#' @description
#' This function runs the bridge-penalized extended two-way fixed effects estimator (\code{twfeCovs()}) on
#' simulated data. It is simply a wrapper for \code{twfeCovs()}: it accepts an object of class
#' \code{"FETWFE_simulated"} (produced by \code{simulateData()}) and unpacks the necessary
#' components to pass to \code{twfeCovs()}. So the outputs match \code{twfeCovs()}, and the needed inputs
#' match their counterparts in \code{twfeCovs()}.
#'
#' @param simulated_obj An object of class \code{"FETWFE_simulated"} containing the simulated panel
#' data and design matrix.
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
#' provided.}
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
#'   result <- twfeCovsWithSimulatedData(sim_data)
#' }
#'
#' @export
twfeCovsWithSimulatedData <- function(
	simulated_obj,
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

	res <- twfeCovs(
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

	return(res)
}


#' Core Estimation Logic for twfeCovs
#'
#' @description
#' This function implements the core estimation steps of the twfeCovs methodology.
#' It takes a pre-processed design matrix and response, handles variance components, performs
#' ordinary least squares regression, and calculates treatment effects and their standard errors.
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
#' @param verbose Logical; if `TRUE`, prints progress messages. Default is `FALSE`.
#' @param alpha Numeric; significance level for confidence intervals (e.g., 0.05 for
#'   95% CIs). Default is 0.05.
#' @param add_ridge (Optional.) Logical; if TRUE, adds a small amount of ridge
#'   regularization to the (untransformed) coefficients to stabilize estimation.
#'   Default is FALSE.
#'
#' @details
#' The function executes the following main steps:
#' \enumerate{
#'   \item **Input Checks:** Validates the provided parameters.
#'   \item **Coordinate Transformation:** Calls `transformXintImproved` to transform
#'     `X_ints` into `X_mod`. This transformation allows a standard bridge
#'     regression penalty on `X_mod` to achieve the desired fusion penalties
#'     on the original coefficients.
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
#'     penalty on the *original* (untransformed) coefficient scale. This involves
#'     using `genFullInvFusionTransformMat` to get the inverse of the overall
#'     fusion transformation matrix.
#'   \item **Cohort Probabilities:** Calculates cohort membership probabilities
#'     conditional on being treated, using `in_sample_counts` and `indep_counts`
#'     if available.
#'   \item **Bridge Regression:** Fits a bridge regression model using
#'     `grpreg::gBridge` on `X_final_scaled` and `y_final` with the specified `q`
#'     and lambda sequence.
#'   \item **Coefficient Selection (BIC):** Calls `getBetaBIC` to select the
#'     optimal `lambda` using BIC and retrieve the corresponding estimated
#'     coefficients (`theta_hat` in the transformed space).
#'   \item **Handle Zero-Feature Case:** If BIC selects a model with zero features,
#'     treatment effects are set to zero.
#'   \item **Coefficient Untransformation:** Calls `untransformCoefImproved` to
#'     transform `theta_hat` back to the original coefficient space, yielding
#'     `beta_hat`. If `add_ridge` was true, `beta_hat` is scaled.
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
#'   \item{theta_hat}{The vector of estimated coefficients in the *transformed* (fused) space, including the intercept as the first element.}
#'   \item{beta_hat}{The vector of estimated coefficients in the *original* space (after untransforming `theta_hat`, excluding intercept).}
#'   \item{treat_inds}{Indices in `beta_hat` corresponding to base treatment effects.}
#'   \item{treat_int_inds}{Indices in `beta_hat` corresponding to treatment-covariate interactions.}
#'   \item{cohort_probs}{Estimated cohort probabilities conditional on being treated, from `in_sample_counts`.}
#'   \item{indep_cohort_probs}{Estimated cohort probabilities from `indep_counts` (NA if not provided).}
#'   \item{sig_eps_sq}{The (possibly estimated) variance of observation-level noise.}
#'   \item{sig_eps_c_sq}{The (possibly estimated) variance of unit-level random effects.}
#'   \item{X_ints}{The original input design matrix from `prepXints`.}
#'   \item{y}{The original input centered response vector from `prepXints`.}
#'   \item{X_final}{The design matrix after fusion transformation and GLS weighting.}
#'   \item{y_final}{The response vector after GLS weighting.}
#'   \item{N, T, R, d, p}{Dimensions used in estimation.}
#' @keywords internal
#' @noRd
twfeCovs_core <- function(
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

	stopifnot(length(c_names) == R)

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
		indep_counts = indep_counts,
		is_fetwfe = FALSE,
		is_twfe_covs = TRUE
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
	# Step 4: estimate OLS regression and extract fitted coefficients
	#
	#

	p_short <- R + T - 1 + d + R
	treat_inds_short <- (R + T - 1 + d + 1):p_short
	first_inds <- 1:R

	df <- data.frame(y = y_final, X_final_scaled)

	stopifnot(all(!is.na(df)))
	stopifnot("y" %in% colnames(df))
	stopifnot(ncol(df) == p_short + 1)

	t0 <- Sys.time()

	# Response already centered; no intercept needed
	fit <- lm(y ~ . + 0, df)

	stopifnot(length(coef(fit)) == p_short)
	stopifnot(length(scale_scale) == p_short)

	beta_hat_slopes <- coef(fit) / scale_scale

	if (add_ridge) {
		lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)
		beta_hat_slopes <- beta_hat_slopes * (1 + lambda_ridge)
	}

	stopifnot(length(beta_hat_slopes) == p_short)
	stopifnot(all(!is.na(beta_hat_slopes)))

	stopifnot(max(treat_inds_short) == p_short)

	treat_int_inds <- c()

	stopifnot(length(treat_inds_short) == R)

	# Get actual estimated treatment effects (in original, untransformed space)
	tes <- beta_hat_slopes[treat_inds_short]

	stopifnot(all(!is.na(tes)))

	stopifnot(length(tes) == R)

	stopifnot(length(first_inds) == R)
	stopifnot(max(first_inds) <= R)

	#
	#
	# Step 6: calculate cohort-specific treatment effects and standard
	# errors
	#
	#

	res <- getCohortATTsFinalOLS(
		X_final = X_final, # This is X_mod * GLS_transform_matrix
		treat_inds = treat_inds_short, # Global indices for treatment effects
		num_treats = R,
		first_inds = first_inds,
		c_names = c_names,
		tes = tes, # Treatment effect estimates (beta_hat_slopes[treat_inds])
		sig_eps_sq = sig_eps_sq,
		R = R,
		N = N,
		T = T,
		p = p_short, # Total number of original parameters (columns in X_ints)
		alpha = alpha
	)

	cohort_te_df <- res$cohort_te_df
	cohort_tes <- res$cohort_tes
	cohort_te_ses <- res$cohort_te_ses
	psi_mat <- res$psi_mat
	gram_inv <- res$gram_inv
	calc_ses <- res$calc_ses

	rm(res)

	stopifnot(nrow(psi_mat) == R)
	stopifnot(ncol(psi_mat) == R)

	#
	#
	# Step 7: calculate overall average treatment effect on treated units
	#
	#

	# Get overal estimated ATT!
	stopifnot(length(tes) == R)
	stopifnot(nrow(psi_mat) == length(tes))

	in_sample_te_results <- getTeResultsOLS(
		sig_eps_sq = sig_eps_sq,
		N = N,
		T = T,
		R = R,
		num_treats = R,
		cohort_tes = cohort_tes, # CATTs (point estimates)
		cohort_probs = cohort_probs, # In-sample pi_r | treated
		psi_mat = psi_mat,
		gram_inv = gram_inv,
		tes = tes, # Untransformed treatment effect estimates beta_hat[treat_inds]
		cohort_probs_overall = cohort_probs_overall, # In-sample pi_r (unconditional on treated)
		first_inds = first_inds,
		calc_ses = calc_ses,
		indep_probs = FALSE
	)

	in_sample_att_hat <- in_sample_te_results$att_hat
	in_sample_att_se <- in_sample_te_results$att_te_se
	in_sample_att_se_no_prob <- in_sample_te_results$att_te_se_no_prob

	if (indep_count_data_available) {
		indep_te_results <- getTeResultsOLS(
			sig_eps_sq = sig_eps_sq,
			N = N,
			T = T,
			R = R,
			num_treats = R,
			cohort_tes = cohort_tes,
			cohort_probs = indep_cohort_probs, # indep pi_r | treated
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			tes = tes,
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
		beta_hat = beta_hat_slopes, # Untransformed slopes
		treat_inds = treat_inds_short,
		treat_int_inds = treat_int_inds,
		cohort_probs = cohort_probs,
		indep_cohort_probs = indep_cohort_probs,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		X_ints = X_ints,
		y = y,
		X_final = X_final,
		y_final = y_final,
		N = N,
		T = T,
		R = R,
		d = d,
		p = p_short,
		calc_ses = calc_ses
	))
}
