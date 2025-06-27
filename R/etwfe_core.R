#' Run ETWFE on Simulated Data
#'
#' @description
#' This function runs the extended two-way fixed effects estimator (\code{etwfe()}) on
#' simulated data. It is simply a wrapper for \code{etwfe()}: it accepts an object of class
#' \code{"FETWFE_simulated"} (produced by \code{simulateData()}) and unpacks the necessary
#' components to pass to \code{etwfe()}. So the outputs match \code{etwfe()}, and the needed inputs
#' match their counterparts in \code{etwfe()}.
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
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
#'
#'   result <- etwfeWithSimulatedData(sim_data)
#' }
#'
#' @export
etwfeWithSimulatedData <- function(
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

	res <- etwfe(
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


# checkEtwfeInputs
#' @title Check Inputs for the main `fetwfe` function
#' @description Validates the inputs provided to the main `fetwfe` function,
#'   ensuring they meet type, dimension, and content requirements. Stops
#'   execution with an error message if any check fails.
#' @param pdata Dataframe; the panel data set.
#' @param time_var Character; name of the time variable column.
#' @param unit_var Character; name of the unit variable column.
#' @param treatment Character; name of the treatment indicator column.
#' @param response Character; name of the response variable column.
#' @param covs Character vector; names of covariate columns. Default `c()`.
#' @param indep_counts Integer vector or NA; counts for independent cohort data.
#'   Default `NA`.
#' @param sig_eps_sq Numeric or NA; variance of idiosyncratic noise. Default `NA`.
#' @param sig_eps_c_sq Numeric or NA; variance of unit-level random effects.
#'   Default `NA`.
#' @param lambda.max Numeric or NA; maximum lambda for `gBridge`. Default `NA`.
#' @param lambda.min Numeric or NA; minimum lambda for `gBridge`. Default `NA`.
#' @param nlambda Integer; number of lambdas for `gBridge`. Default `100`.
#' @param q Numeric; Lq penalty exponent for `gBridge`. Default `0.5`.
#' @param verbose Logical; if TRUE, print progress. Default `FALSE`.
#' @param alpha Numeric; significance level for confidence intervals. Default `0.05`.
#' @param add_ridge Logical; if TRUE, add small ridge penalty. Default `FALSE`.
#' @return Logical `indep_count_data_available`, which is `TRUE` if valid
#'   `indep_counts` were provided, `FALSE` otherwise.
#' @details This function performs a series of `stopifnot` checks on each
#'   parameter. For example:
#'   - `pdata` must be a dataframe with at least 4 rows.
#'   - `time_var`, `unit_var`, `treatment`, `response` must be single characters,
#'     present in `pdata`, and the corresponding columns must have the correct
#'     type (e.g., integer for time, character for unit, 0/1 integer for treatment).
#'   - `covs` if provided, must be characters, present in `pdata`, and columns
#'     must be numeric, integer, or factor.
#'   - `indep_counts` if provided, must be positive integers.
#'   - `sig_eps_sq`, `sig_eps_c_sq` if provided, must be non-negative numerics.
#'   - `lambda.max`, `lambda.min` if provided, must be valid numerics
#'     (`lambda.max > lambda.min >= 0`).
#'   - `q` must be in `(0, 2]`.
#'   - `alpha` must be in `(0, 1)`.
#'   Issues a warning if `alpha > 0.5`.
#' @keywords internal
#' @noRd
checkEtwfeInputs <- function(
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
	stopifnot(is.data.frame(pdata))
	# Check if pdata is a tibble; if so, convert to a dataframe
	if ("tbl_df" %in% class(pdata)) {
		pdata <- as.data.frame(pdata)
	}
	stopifnot(nrow(pdata) >= 4) # bare minimum, 2 units at 2 times

	stopifnot(is.character(time_var))
	stopifnot(length(time_var) == 1)
	stopifnot(time_var %in% colnames(pdata))
	stopifnot(is.integer(pdata[[time_var]]))

	stopifnot(is.character(unit_var))
	stopifnot(length(unit_var) == 1)
	stopifnot(unit_var %in% colnames(pdata))
	stopifnot(is.character(pdata[[unit_var]]))

	stopifnot(is.character(treatment))
	stopifnot(length(treatment) == 1)
	stopifnot(treatment %in% colnames(pdata))
	stopifnot(is.integer(pdata[[treatment]]))
	stopifnot(all(pdata[, treatment] %in% c(0, 1)))

	if (length(covs) > 0) {
		stopifnot(is.character(covs))
		stopifnot(all(covs %in% colnames(pdata)))
		for (cov in covs) {
			stopifnot(
				is.numeric(pdata[[cov]]) |
					is.integer(pdata[[cov]]) |
					is.factor(pdata[[cov]])
			)
		}
	}

	stopifnot(is.character(response))
	stopifnot(length(response) == 1)
	stopifnot(response %in% colnames(pdata))
	stopifnot(is.numeric(pdata[[response]]) | is.integer(pdata[[response]]))

	indep_count_data_available <- FALSE
	if (any(!is.na(indep_counts))) {
		stopifnot(is.integer(indep_counts))
		if (any(indep_counts <= 0)) {
			stop(
				"At least one cohort in the independent count data has 0 members"
			)
		}
		indep_count_data_available <- TRUE
	}

	if (any(!is.na(sig_eps_sq))) {
		stopifnot(is.numeric(sig_eps_sq) | is.integer(sig_eps_sq))
		stopifnot(length(sig_eps_sq) == 1)
		stopifnot(sig_eps_sq >= 0)
	}

	if (any(!is.na(sig_eps_c_sq))) {
		stopifnot(is.numeric(sig_eps_c_sq) | is.integer(sig_eps_c_sq))
		stopifnot(length(sig_eps_c_sq) == 1)
		stopifnot(sig_eps_c_sq >= 0)
	}

	stopifnot(is.logical(verbose))
	stopifnot(length(verbose) == 1)

	stopifnot(is.numeric(alpha))
	stopifnot(length(alpha) == 1)
	stopifnot(alpha > 0)
	stopifnot(alpha < 1)
	if (alpha > 0.5) {
		warning(
			"Provided alpha > 0.5; are you sure you didn't mean to enter a smaller alpha? The confidence level will be 1 - alpha."
		)
	}

	stopifnot(is.logical(add_ridge))
	stopifnot(length(add_ridge) == 1)

	return(list(
		pdata = pdata,
		indep_count_data_available = indep_count_data_available
	))
}

#' Core Estimation Logic for Extended Two-Way Fixed Effects
#'
#' @description
#' This function implements the core estimation steps of the ETWFE methodology.
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
#'   the untransformed coefficients to stabilize estimation. Default is `FALSE`.
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
#'   \item{lambda.max}{The maximum lambda value used in `grpreg`.}
#'   \item{lambda.max_model_size}{Model size for `lambda.max`.}
#'   \item{lambda.min}{The minimum lambda value used in `grpreg`.}
#'   \item{lambda.min_model_size}{Model size for `lambda.min`.}
#'   \item{lambda_star}{The lambda value selected by BIC.}
#'   \item{lambda_star_model_size}{Model size for `lambda_star`.}
#'   \item{X_ints}{The original input design matrix from `prepXints`.}
#'   \item{y}{The original input centered response vector from `prepXints`.}
#'   \item{X_final}{The design matrix after fusion transformation and GLS weighting.}
#'   \item{y_final}{The response vector after GLS weighting.}
#'   \item{N, T, R, d, p}{Dimensions used in estimation.}
#' @keywords internal
#' @noRd
etwfe_core <- function(
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
		is_fetwfe = FALSE
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

	df <- data.frame(y = y_final, X_final_scaled)

	stopifnot(all(!is.na(df)))
	stopifnot("y" %in% colnames(df))

	t0 <- Sys.time()

	# Response already centered; no intercept needed
	fit <- lm(y ~ . + 0, df)

	beta_hat_slopes <- coef(fit) / scale_scale

	stopifnot(length(beta_hat_slopes) == p)
	stopifnot(all(!is.na(beta_hat_slopes)))

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

	stopifnot(length(treat_inds) == num_treats)

	# If using ridge regularization, multiply the "naive" estimated coefficients
	# by 1 + lambda_ridge, similar to suggestion in original elastic net paper.
	if (add_ridge) {
		lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)
		beta_hat_slopes <- beta_hat_slopes * (1 + lambda_ridge)
		stopifnot(all(!is.na(beta_hat_slopes)))
	}

	# Get actual estimated treatment effects (in original, untransformed space)
	tes <- beta_hat_slopes[treat_inds]

	stopifnot(all(!is.na(tes)))

	stopifnot(length(tes) == num_treats)

	stopifnot(length(first_inds) == R)
	stopifnot(max(first_inds) <= num_treats)

	#
	#
	# Step 6: calculate cohort-specific treatment effects and standard
	# errors
	#
	#

	res <- getCohortATTsFinalOLS(
		X_final = X_final, # This is X_mod * GLS_transform_matrix
		treat_inds = treat_inds, # Global indices for treatment effects
		num_treats = num_treats,
		first_inds = first_inds,
		c_names = c_names,
		tes = tes, # Treatment effect estimates (beta_hat_slopes[treat_inds])
		sig_eps_sq = sig_eps_sq,
		R = R,
		N = N,
		T = T,
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

	stopifnot(nrow(psi_mat) == num_treats)
	stopifnot(ncol(psi_mat) == R)

	#
	#
	# Step 7: calculate overall average treatment effect on treated units
	#
	#

	# Get overal estimated ATT!
	stopifnot(length(tes) == num_treats)
	stopifnot(nrow(psi_mat) == length(tes))

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
			num_treats = num_treats,
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
		treat_inds = treat_inds,
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
		p = p,
		calc_ses = calc_ses
	))
}


#' Prepare Data & Design Matrix for ETWFE/TWFE Workflows
#'
#' @description
#' A helper that converts raw **panel data** into the core objects required by
#' `etwfe_core()`, `twfeCovs_core()`, and related fitting routines.  The function
#' (i) keeps only the relevant variables, (ii) one-hot-encodes or otherwise
#' expands any *factor* covariates, (iii) builds the stacked design matrix
#' `X_ints` and centred response `y` via [`prepXints()`], and (iv) performs a
#' battery of **sanity checks** on cohort counts before estimation proceeds.
#'
#' @param pdata A **data.frame** in long (unit x time) format containing at
#'   least the columns named in `response`, `time_var`, `unit_var`,
#'   `treatment`, and `covs`.
#' @param response Character; name of the response variable column.
#' @param time_var Character; name of the time variable column.
#' @param unit_var Character; name of the unit variable column.
#' @param treatment Character; name of the treatment indicator column.
#' @param covs Character vector of additional covariate column names.  Factor
#'   covariates are expanded to dummies by `processFactors()`.  May be empty.
#' @param verbose Logical; print progress/timing information?
#' @param indep_count_data_available Logical; `TRUE` when *independent* cohort
#'   counts are supplied in `indep_counts` and should be validated.
#' @param indep_counts Optional integer vector with length
#'   `1 + R` (never-treated plus `R` treated cohorts) giving cohort sizes in an
#'   independent sample.  Only used when
#'   `indep_count_data_available = TRUE`.
#'
#' @return A named **list** ready to be passed into the "_core" estimators:
#' \describe{
#'   \item{pdata}{The (possibly modified) data frame after factor processing.}
#'   \item{covs}{Updated character vector of covariate names after dummy-expansion.}
#'   \item{X_ints}{Design matrix with unit FE, time FE, covariates,
#'     treatment dummies and their interactions (dimensions \(N T \times p\)).}
#'   \item{y}{Centred response vector of length \(N T\).}
#'   \item{N, T}{Integers - number of unique units and time periods.}
#'   \item{d}{Integer - number of *raw* covariates after processing.}
#'   \item{p}{Integer - total number of columns in `X_ints`.}
#'   \item{in_sample_counts}{Named integer vector of length `1 + R` with cohort
#'     sizes in the estimation sample (first entry = never-treated).}
#'   \item{num_treats}{Integer - total number of base treatment-effect
#'     parameters in `X_ints`.}
#'   \item{first_inds}{Integer vector (length `R`) giving the first column
#'     index of each cohort's treatment-effect block inside `X_ints`.}
#'   \item{R}{Integer - number of treated cohorts detected.}
#' }
#'
#' @details
#' The routine executes the following steps:
#' \enumerate{
#'   \item **Column subset** - drops all variables not explicitly required.
#'   \item **Factor processing** - calls `processFactors()` so that every
#'     categorical covariate becomes a set of \(0/1\) dummies.
#'   \item **Design-matrix construction** - hands the cleaned data to
#'     `prepXints()`, which adds unit and time dummies, interactions, etc.
#'   \item **Cohort diagnostics** - verifies that
#'     \eqn{R \ge 2}, at least one never-treated unit exists, cohort names are
#'     unique, and (if provided) `indep_counts` are dimensionally consistent.
#' }
#'
#' Any violation of the checks triggers an informative `stop()` message so that
#' higher-level functions fail fast.
#'
#' @keywords internal
#' @noRd
prep_for_etwfe_core <- function(
	pdata,
	response,
	time_var,
	unit_var,
	treatment,
	covs,
	verbose,
	indep_count_data_available,
	indep_counts
) {
	# Subset pdata to include only the key columns
	pdata <- pdata[, c(response, time_var, unit_var, treatment, covs)]

	# Process any factor covariates:
	if (length(covs) > 0) {
		pf_res <- processFactors(pdata, covs)
		pdata <- pf_res$pdata
		covs <- pf_res$covs
	}

	# Proceed to generate design matrix and other objects
	res <- prepXints(
		data = pdata,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		covs = covs,
		response = response,
		verbose = verbose
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
	if (R < 2) {
		stop(
			"Only one treated cohort detected in data. Currently fetwfe and etwfe only support data sets with at least two treated cohorts."
		)
	}
	stopifnot(N >= R + 1)
	stopifnot(sum(in_sample_counts) == N)
	stopifnot(all(in_sample_counts >= 0))
	stopifnot(is.integer(in_sample_counts))
	if (in_sample_counts[1] == 0) {
		stop(
			"No never-treated units detected in data to fit model; estimating treatment effects is not possible"
		)
	}
	if (length(names(in_sample_counts)) != length(in_sample_counts)) {
		stop(
			"in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)"
		)
	}
	if (
		length(names(in_sample_counts)) !=
			length(unique(names(in_sample_counts)))
	) {
		stop(
			"in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)"
		)
	}
	if (indep_count_data_available) {
		if (sum(indep_counts) != N) {
			stop(
				"Number of units in independent cohort count data does not equal number of units in data to be used to fit model."
			)
		}
		if (length(indep_counts) != length(in_sample_counts)) {
			stop(
				"Number of counts in independent counts does not match number of cohorts in data to be used to fit model."
			)
		}
	}

	return(list(
		pdata = pdata,
		covs = covs,
		X_ints = X_ints,
		y = y,
		N = N,
		T = T,
		d = d,
		p = p,
		in_sample_counts = in_sample_counts,
		num_treats = num_treats,
		first_inds = first_inds,
		R = R
	))
}

# processFactors
#' @title Process Factor Covariates into Dummy Variables
#' @description Identifies factor columns within the specified covariates and
#'   converts them into dummy variables. Other covariate types are left unchanged.
#' @param pdata Dataframe; the panel data set. This dataframe will be modified
#'   in place by removing original factor columns and adding new dummy columns.
#' @param covs Character vector; names of covariate columns in `pdata` to process.
#' @return A list containing:
#'   \item{pdata}{The modified dataframe with factor covariates replaced by
#'     their dummy variable representations.}
#'   \item{covs}{The updated character vector of covariate names, where original
#'     factor names are replaced by the names of the newly created dummy variables.}
#' @details For each covariate name in `covs`:
#'   - If the corresponding column in `pdata` is a factor:
#'     - `stats::model.matrix(~ pdata[[v]] - 1)` is used to create dummy variables
#'       (without an intercept column).
#'     - If the factor has more than one level, the first dummy column created
#'       (corresponding to the first level) is dropped to serve as the baseline/reference category.
#'     - New dummy column names are generated by pasting the original factor name
#'       with the level name (e.g., "factorVar_levelName").
#'     - The original factor column is removed from `pdata`, and the new dummy
#'       columns are `cbind`ed to `pdata`.
#'     - The original factor name in the `covs` list is replaced by the new dummy names.
#'   - If the column is not a factor, it (and its name in `covs`) is unchanged.
#' @examples
#'   # df <- data.frame(y = 1:4,
#'   #                  group = factor(c("A", "B", "A", "C")),
#'   #                  x1 = 11:14)
#'   # cov_names <- c("group", "x1")
#'   # result <- processFactors(df, cov_names)
#'   # print(result$pdata)
#'   # print(result$covs)
#'   # Expected output for pdata would have 'group_B', 'group_C' (or similar)
#'   # and 'x1'. 'group' column would be removed.
#'   # 'covs' would be c("group_B", "group_C", "x1").
#' @keywords internal
#' @noRd
processFactors <- function(pdata, covs) {
	new_covs <- c()
	# Loop over each variable in covs
	for (v in covs) {
		if (is.factor(pdata[[v]])) {
			# Create dummy variables from the factor.
			# The model.matrix() call produces an intercept and dummies; we drop the intercept
			dummies <- stats::model.matrix(~ pdata[[v]] - 1)
			# If there is more than one level, drop the first column to use it as baseline.
			if (ncol(dummies) > 1) {
				dummies <- dummies[, -1, drop = FALSE]
			}
			# Rename the dummy columns: for example, if v = "group", new names will be "group_level2", etc.
			dummy_names <- paste(v, colnames(dummies), sep = "_")
			colnames(dummies) <- dummy_names
			# Remove the original factor column from pdata
			pdata[[v]] <- NULL
			# Bind the new dummy columns
			pdata <- cbind(pdata, dummies)
			# Record the new dummy variable names
			new_covs <- c(new_covs, dummy_names)
		} else {
			# Leave non-factor columns unchanged.
			new_covs <- c(new_covs, v)
		}
	}
	return(list(pdata = pdata, covs = new_covs))
}

#' @title Prepare Design Matrix and Initial Data for FETWFE
#'
#' @description
#' This function serves as a pre-processing step for `fetwfe`. It takes the raw
#' panel data and transforms it into a design matrix (`X_ints`) and a response
#' vector (`y`) suitable for the core estimation logic. It identifies cohorts,
#' processes covariates (making them time-invariant and handling missing values),
#' and generates dummy variables for cohorts, time periods, and treatments.
#'
#' @param data A `data.frame` object representing the panel data. Each row
#'   should be an observation of a unit at a specific time.
#' @param time_var A character string specifying the name of the column in `data`
#'   that contains the time period variable (e.g., year). Expected to be integer.
#' @param unit_var A character string specifying the name of the column in `data`
#'   that contains the unit identifier. Expected to be character.
#' @param treatment A character string specifying the name of the column in `data`
#'   that contains the binary treatment indicator (0 for untreated, 1 for treated).
#'   Treatment is assumed to be an absorbing state.
#' @param covs A character vector specifying the names of the columns in `data`
#'   to be used as covariates. These are treated as time-invariant (values from
#'   the first period are used).
#' @param response A character string specifying the name of the column in `data`
#'   that contains the response variable.
#' @param verbose Logical. If `TRUE`, messages about the data processing steps
#'   will be printed. Default is `FALSE`.
#'
#' @details
#' The function performs several key steps:
#' \enumerate{
#'   \item Calls `idCohorts` to identify treatment cohorts, unique units, and time
#'     periods. Units treated in the very first time period are removed.
#'   \item Calls `processCovs` to handle covariates. Time-varying covariates are
#'     replaced by their value in the first (pre-treatment) period. Covariates
#'     with missing values in the first period or covariates that are constant
#'     across all units are removed.
#'   \item Calls `addDummies` to create dummy variables for cohorts, time periods
#'     (excluding the first), and treatment indicators (for each cohort and
#'     each post-treatment time period). The response variable is also centered.
#'   \item Calls `genXintsData` to construct the final design matrix `X_ints` by
#'     combining the dummy variables, covariates, and their interactions.
#'     The interactions include cohort-covariate, time-covariate, and
#'     treatment-covariate interactions.
#' }
#' Input `covs` are expected to be numeric or integer after factor processing.
#' The function also calculates the number of units (`N`), time periods (`T`),
#' covariates (`d`), the total number of parameters in the full design matrix (`p`),
#' in-sample cohort counts, the number of unique treatment terms (`num_treats`),
#' and indices of the first treatment effect for each cohort (`first_inds`).
#'
#' @return A list containing:
#'   \item{X_ints}{The fully constructed design matrix with all fixed effects,
#'     covariates, treatment dummies, and their interactions.}
#'   \item{y}{The centered response vector.}
#'   \item{N}{The final number of unique units after processing.}
#'   \item{T}{The number of unique time periods.}
#'   \item{d}{The final number of covariates after processing.}
#'   \item{p}{The total number of columns in `X_ints` (total parameters).}
#'   \item{in_sample_counts}{An integer vector named with cohort identifiers
#'     (including "Never_treated"), indicating the number of units in each cohort
#'     within the provided `data`.}
#'   \item{num_treats}{The total number of unique treatment effect parameters
#'     (e.g., \eqn{\tau_{rt}} from the paper).}
#'   \item{first_inds}{A numeric vector indicating the starting column index within
#'     the block of treatment effect parameters in `X_ints` (and subsequently in
#'     `beta_hat`) for each respective cohort's first treatment period.}
#' @keywords internal
#' @noRd
prepXints <- function(
	data,
	time_var,
	unit_var,
	treatment,
	covs,
	response,
	verbose = FALSE
) {
	# Check inputs
	stopifnot(is.data.frame(data))
	stopifnot(nrow(data) >= 4) # bare minimum, 2 units at 2 times

	stopifnot(is.character(time_var))
	stopifnot(length(time_var) == 1)
	stopifnot(time_var %in% colnames(data))
	stopifnot(is.integer(data[[time_var]]))

	stopifnot(is.character(unit_var))
	stopifnot(length(unit_var) == 1)
	stopifnot(unit_var %in% colnames(data))
	stopifnot(is.character(data[[unit_var]]))

	stopifnot(is.character(treatment))
	stopifnot(length(treatment) == 1)
	stopifnot(treatment %in% colnames(data))
	stopifnot(is.integer(data[[treatment]]))
	stopifnot(all(data[, treatment] %in% c(0, 1)))

	if (length(covs) > 0) {
		stopifnot(is.character(covs))
		stopifnot(all(covs %in% colnames(data)))
		for (cov in covs) {
			stopifnot(is.numeric(data[[cov]]) | is.integer(data[[cov]]))
		}
	}

	stopifnot(is.character(response))
	stopifnot(length(response) == 1)
	stopifnot(response %in% colnames(data))
	stopifnot(is.numeric(data[[response]]) | is.integer(data[[response]]))

	stopifnot(is.logical(verbose))
	stopifnot(length(verbose) == 1)

	# Identify cohorts, units, and times; remove units that were treated in
	# first period; remove treatment variable after identifying cohorts
	ret <- idCohorts(
		df = data,
		time_var = time_var,
		unit_var = unit_var,
		treat_var = treatment,
		covs = covs
	)

	data <- ret$df
	cohorts <- ret$cohorts
	units <- ret$units
	times <- ret$times

	rm(ret)

	N <- length(units)
	T <- length(times)
	R <- length(cohorts)

	stopifnot(R >= 1)
	stopifnot(T >= 2)
	stopifnot(R <= T - 1)

	N_treated <- length(unlist(cohorts))

	stopifnot(N_treated <= N)

	N_UNTREATED <- N - N_treated

	# Get in-sample counts of units in each cohort
	in_sample_counts <- rep(as.numeric(NA), R + 1)

	for (r in 1:R) {
		in_sample_counts[r + 1] <- length(cohorts[[r]])
	}

	in_sample_counts[1] <- N - N_treated

	stopifnot(all(in_sample_counts == round(in_sample_counts)))
	in_sample_counts <- as.integer(in_sample_counts)
	stopifnot(all(!is.na(in_sample_counts)))
	stopifnot(sum(in_sample_counts) == N)

	if (in_sample_counts[1] == 0) {
		warning(
			"No untreated units found. It will not be possible to estimate treatment effects on this data set."
		)
	}

	names(in_sample_counts) <- c("Never_treated", names(cohorts))

	# Replace time-varying covariates with value in first (pre-treated) period;
	# remove covariates with any missing values in this period
	ret <- processCovs(
		df = data,
		units = units,
		unit_var = unit_var,
		times = times,
		time_var = time_var,
		covs = covs,
		resp_var = response,
		T = T,
		verbose = verbose
	)

	data <- ret$df
	covs <- ret$covs

	rm(ret)

	d <- length(covs)

	stopifnot(all(!is.na(data)))
	stopifnot(nrow(data) == N * T)

	# Add cohort dummies, time dummies, and treatment variables. (Names of
	# cohorts must be same as times of treatment.) Also, center response.
	ret <- addDummies(
		df = data,
		cohorts = cohorts,
		times = times,
		N = N,
		T = T,
		unit_var = unit_var,
		time_var = time_var,
		resp_var = response,
		n_cohorts = R
	)

	y <- ret$y # response
	cohort_treat_names <- ret$cohort_treat_names # List of names of treatment
	# dummies for each cohort
	time_var_names <- ret$time_var_names # Names of time dummies
	cohort_vars <- ret$cohort_vars # Names of cohort dummies
	first_inds <- ret$first_inds # Among blocks of treatment variables, indices
	# of first treatment effects corresponding to each cohort
	cohort_var_mat <- ret$cohort_var_mat # Cohort dummies; R columns total
	time_var_mat <- ret$time_var_mat # Time dummies; T - 1 columns total
	treat_var_mat <- ret$treat_var_mat # Treatment dummies

	num_treats <- ncol(treat_var_mat)

	stopifnot(num_treats == length(unlist(cohort_treat_names)))
	stopifnot(all(covs %in% colnames(data)))

	stopifnot(is.numeric(ncol(treat_var_mat)) | is.integer(ncol(treat_var_mat)))
	stopifnot(ncol(treat_var_mat) >= 1)

	covariate_mat <- as.matrix(data[, covs, drop = FALSE]) # Covariates

	rm(ret)
	rm(data)

	#### Create matrix with all interactions
	p <- getP(R = R, T = T, d = d, num_treats = num_treats)

	stopifnot(ncol(cohort_var_mat) == R)
	stopifnot(is.matrix(covariate_mat))

	X_ints <- genXintsData(
		cohort_fe = cohort_var_mat,
		time_fe = time_var_mat,
		X_long = covariate_mat,
		treat_mat_long = treat_var_mat,
		N = N,
		R = R,
		T = T,
		d = d,
		N_UNTREATED = N_UNTREATED,
		p = p
	)

	return(list(
		X_ints = X_ints,
		y = y,
		N = N,
		T = T,
		d = d,
		p = p,
		in_sample_counts = in_sample_counts,
		num_treats = num_treats,
		first_inds = first_inds
	))
}

###### Replace time-varying variables with pre-treatment values for all times;

#' Process Covariates for FETWFE
#'
#' @description
#' This function processes covariates for use in the FETWFE model. It ensures
#' covariates are time-invariant by using their first-period values. It also
#' handles missing covariate values and removes covariates that are constant
#' across all units.
#'
#' @param df A `data.frame` containing the panel data. It should have been
#'   partially processed by `idCohorts` (e.g., first-period treated units removed).
#' @param units A character vector of unique unit identifiers to process.
#' @param unit_var Character string; name of the unit identifier column in `df`.
#' @param times A numeric or integer vector of unique, sorted time periods.
#' @param time_var Character string; name of the time variable column in `df`.
#' @param covs A character vector of covariate names to process.
#' @param resp_var Character string; name of the response variable column in `df`.
#' @param T Integer; the total number of unique time periods.
#' @param verbose Logical. If `TRUE`, prints messages about covariate processing.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item If no covariates are provided (`length(covs) == 0`), it returns the
#'     dataframe ordered by unit and time, containing only response, time, and
#'     unit variables.
#'   \item For each unit, it checks for missing covariate values in the first
#'     time period (`times[1]`). Covariates with any `NA` in the first period
#'     for any unit are removed from the `covs` list, and a warning is issued.
#'   \item Covariates that have the same value across all unit-time observations
#'     are removed, and a warning is issued.
#'   \item The dataframe `df` is subsetted to keep only the response, time, unit,
#'     and valid covariate columns.
#'   \item For any remaining time-varying covariates, their values for all time
#'     periods of a unit are replaced by that unit's value from the first time
#'     period (`times[1]`).
#'   \item Finally, `df` is sorted by unit and then by time.
#' }
#'
#' @return A list containing:
#'   \item{df}{The processed `data.frame` with covariates made time-invariant
#'     and problematic ones removed. Rows are sorted by unit and time.}
#'   \item{covs}{The updated character vector of valid covariate names remaining
#'     after processing.}
#' @keywords internal
#' @noRd
processCovs <- function(
	df,
	units,
	unit_var,
	times,
	time_var,
	covs,
	resp_var,
	T,
	verbose = FALSE
) {
	# Always check that every unit has exactly T observations.
	for (s in units) {
		df_s <- df[df[, unit_var] == s, ]
		if (nrow(df_s) != T) {
			stop(paste("Unit", s, "does not have exactly", T, "observations."))
		}
	}

	# If no covariates are provided, simply return the ordered data frame.
	if (length(covs) == 0) {
		if (verbose) {
			message("No covariates provided; skipping covariate processing.")
		}
		# # TODO: change to consider
		# # Ensure df contains only necessary columns if covs is empty.
		# # This selection should ideally happen once, consistently.
		# # prepXints selects columns from the original pdata:
		# # pdata <- pdata[, c(response, time_var, unit_var, treatment, covs)]
		# # Then idCohorts gets this. Then processCovs.
		# # If covs is empty, df will contain response, time_var, unit_var.
		# df_ordered <- df[order(df[, unit_var], df[, time_var], decreasing=FALSE), ]
		# return(list(df = df_ordered, covs = covs))

		df <- df[, c(resp_var, time_var, unit_var)]
		df <- df[order(df[, unit_var], df[, time_var], decreasing = FALSE), ]
		return(list(df = df, covs = covs))
	}
	d_orig <- length(covs) # Original number of covariates passed

	# For each unit, ensure that the first period has non-missing covariate values.
	# Remove any covariates with missing values in the first period for *any* unit.
	covs_to_keep <- character()
	first_time_val <- times[1]

	for (cov_name in covs) {
		is_valid_cov <- TRUE
		for (s in units) {
			val_first_period <- df[
				(df[, unit_var] == s) & (df[, time_var] == first_time_val),
				cov_name
			]
			if (length(val_first_period) != 1 || is.na(val_first_period)) {
				if (verbose)
					message(
						"Covariate '",
						cov_name,
						"' has NA or missing first period value for unit '",
						s,
						"'. Removing covariate."
					)
				is_valid_cov <- FALSE
				break
			}
		}
		if (is_valid_cov) {
			covs_to_keep <- c(covs_to_keep, cov_name)
		}
	}

	if (length(covs_to_keep) < d_orig && length(covs_to_keep) > 0) {
		removed_covs_na <- setdiff(covs, covs_to_keep)
		warning(paste(
			length(removed_covs_na),
			"covariate(s) were removed because they contained missing values in the first time period for at least one unit: ",
			paste(removed_covs_na, collapse = ", ")
		))
	} else if (length(covs_to_keep) == 0 && d_orig > 0) {
		warning(
			"All covariates were removed due to missing values in the first period for at least one unit."
		)
	}
	covs <- covs_to_keep
	d <- length(covs)

	# Remove covariates that are constant across units (after making them time-invariant based on first period).
	if (d > 0) {
		covs_to_remove_const <- character()
		# First, make all covariates time-invariant based on first period value
		df_temp_const_check <- df # Operate on a temporary copy
		for (s in units) {
			first_period_rows_s_idx <- which(
				(df_temp_const_check[, unit_var] == s) &
					(df_temp_const_check[, time_var] == first_time_val)
			)
			stopifnot(length(first_period_rows_s_idx) == 1) # Should be one row for first period
			cov_values_s_first_period <- df_temp_const_check[
				first_period_rows_s_idx,
				covs,
				drop = FALSE
			]

			for (t_idx in seq_along(times)) {
				current_rows_s_t_idx <- which(
					(df_temp_const_check[, unit_var] == s) &
						(df_temp_const_check[, time_var] == times[t_idx])
				)
				stopifnot(length(current_rows_s_t_idx) == 1)
				df_temp_const_check[
					current_rows_s_t_idx,
					covs
				] <- cov_values_s_first_period
			}
		}

		for (cov_name in covs) {
			# Check uniqueness on the first-period values for each unit
			first_period_vals_for_cov <- sapply(units, function(u) {
				df_temp_const_check[
					(df_temp_const_check[, unit_var] == u) &
						(df_temp_const_check[, time_var] == first_time_val),
					cov_name
				]
			})
			if (length(unique(first_period_vals_for_cov)) == 1) {
				if (verbose) {
					message(
						"Removing covariate because all units have the same first-period value: ",
						cov_name
					)
				}
				covs_to_remove_const <- c(covs_to_remove_const, cov_name)
			}
		}
		covs <- covs[!(covs %in% covs_to_remove_const)]

		if (length(covs_to_remove_const) > 0 && length(covs) == 0) {
			warning(
				"All remaining covariates were removed because they were constant across units (based on first-period values). Continuing with no covariates."
			)
		} else if (length(covs_to_remove_const) > 0) {
			warning(paste(
				length(covs_to_remove_const),
				"covariate(s) were removed because all units had the same first-period value: ",
				paste(covs_to_remove_const, collapse = ", ")
			))
		}
		d <- length(covs)
	}

	# Keep only the needed columns.
	df <- df[, c(resp_var, time_var, unit_var, covs), drop = FALSE] # Ensure covs can be empty

	# For any time-varying covariates, replace all values with the first-period value.
	# This was partially done for the constant check, ensure it's finalized.
	if (d > 0) {
		# Only if there are covariates left
		if (verbose) {
			message(
				"Finalizing: Replacing time-varying covariate values with first-period values..."
			)
		}
		for (s in units) {
			first_period_rows_s_idx <- which(
				(df[, unit_var] == s) & (df[, time_var] == first_time_val)
			)
			# This assumes first_time_val is the relevant one (times[1])
			covs_s_first_period_vals <- df[
				first_period_rows_s_idx,
				covs,
				drop = FALSE
			]

			for (t_val in times) {
				ind_s_t <- (df[, unit_var] == s) & (df[, time_var] == t_val)
				stopifnot(sum(ind_s_t) == 1)
				df[ind_s_t, covs] <- covs_s_first_period_vals
			}
		}
	}

	# Sort rows: first T rows should be the observations for the first unit, and
	# so on
	df <- df[order(df[, unit_var], df[, time_var], decreasing = FALSE), ]

	return(list(df = df, covs = covs))
}

######################### Add cohort dummies and treatment variables

#-------------------------------------------------------------------------------
# Dummy Variable and Interaction Matrix Generation
#-------------------------------------------------------------------------------

#' Generate Treatment Variables for Real Data
#'
#' @description
#' Creates a matrix of dummy variables for a specific cohort. This includes one
#' column indicating membership in the cohort and additional columns for each
#' time period during which the cohort is treated.
#'
#' @param cohort_name Character string; the name for the cohort dummy variable column.
#' @param c_t_names Character vector; names for the treatment period dummy
#'   variable columns (e.g., "c1_t2", "c1_t3").
#' @param N Integer; total number of unique units in the panel.
#' @param T Integer; total number of time periods in the panel.
#' @param n_treated_times Integer; the number of time periods for which this
#'   specific cohort is treated.
#' @param unit_vars Factor or character vector from the main data, indicating the
#'   unit for each observation (length N*T).
#' @param time_vars Numeric or integer vector from the main data, indicating the
#'   time period for each observation (length N*T).
#' @param cohort Character vector; identifiers of units belonging to this specific
#'   cohort.
#' @param treated_times Numeric or integer vector; the specific time periods during
#'   which this cohort is treated.
#'
#' @return A matrix with `N*T` rows. The first column is the cohort dummy named
#'   `cohort_name`. Subsequent `n_treated_times` columns are treatment-period
#'   dummies named according to `c_t_names`.
#' @keywords internal
#' @noRd
genTreatVarsRealData <- function(
	cohort_name,
	c_t_names,
	N,
	T,
	n_treated_times,
	unit_vars,
	time_vars,
	cohort,
	treated_times
) {
	# Create matrix of variables for cohort indicator and cohort/treatment
	# variables to append to df
	treat_vars <- matrix(0L, nrow = N * T, ncol = 1 + n_treated_times) # Initialize with 0L

	colnames(treat_vars) <- c(cohort_name, c_t_names)

	# Rows corresponding to units in the current cohort
	c_i_inds_units <- unit_vars %in% cohort # Logical vector of length N*T

	# Add cohort dummy: 1 for all observations of units in this cohort
	treat_vars[c_i_inds_units, 1] <- 1L

	# Add treatment period dummies
	if (n_treated_times > 0) {
		# build a logical index of length N*T for all cohort-time pairs at once
		treat_combo <- interaction(unit_vars, time_vars, drop = TRUE)
		for (j in seq_len(n_treated_times)) {
			treat_vars[
				(unit_vars %in% cohort) & (time_vars == treated_times[j]),
				j + 1L
			] <- 1L
		}
	}

	return(treat_vars)
}

# Names of cohorts must be the same as time of treatment

#' Add Dummy Variables to Dataframe
#'
#' @description
#' Augments the processed dataframe with cohort dummies, time dummies (for all
#' but the first period), and treatment-period dummies for each cohort.
#' It also centers the response variable.
#'
#' @param df A `data.frame` that has been processed by `idCohorts` and `processCovs`.
#' @param cohorts A list where names are treatment adoption times and values are
#'   character vectors of unit IDs in that cohort.
#' @param times A numeric or integer vector of unique, sorted time periods.
#' @param N Integer; total number of unique units.
#' @param T Integer; total number of time periods.
#' @param unit_var Character string; name of the unit identifier column in `df`.
#' @param time_var Character string; name of the time variable column in `df`.
#' @param resp_var Character string; name of the response variable column in `df`.
#' @param n_cohorts Integer; the number of treated cohorts (R from the paper).
#'
#' @details
#' For each cohort:
#' \itemize{
#'   \item Identifies the time periods during which it's treated.
#'   \item Calls `genTreatVarsRealData` to create cohort and treatment-period dummies.
#'   \item These dummies are added as columns to a copy of `df` (though the modified
#'     `df` itself is not returned, its column structure informs outputs).
#'   \item Matrices `cohort_var_mat` (for cohort indicators) and `treat_var_mat`
#'     (for treatment-period indicators) are constructed.
#' }
#' Time dummies are created for `times[2]` through `times[T]`.
#' The response variable `df[[resp_var]]` is centered by subtracting its mean.
#'
#' @return A list containing:
#'   \item{time_var_mat}{Matrix of time dummies (`N*T` rows, `T-1` columns).}
#'   \item{cohort_var_mat}{Matrix of cohort dummies (`N*T` rows, `n_cohorts` columns).}
#'   \item{treat_var_mat}{Matrix of treatment-period dummies (`N*T` rows, `num_treats` columns).}
#'   \item{y}{The centered response vector.}
#'   \item{cohort_treat_names}{A list (names are cohort adoption times) where each
#'     element is a character vector of names for the treatment-period dummy
#'     variables for that cohort (e.g., "c1_t2", "c1_t3").}
#'   \item{time_var_names}{Character vector of names for the time dummy variables
#'     (e.g., "t_2", "t_3").}
#'   \item{cohort_vars}{Character vector of names for the cohort dummy variables
#'     (e.g., "c_1", "c_2").}
#'   \item{first_inds}{A numeric vector. For each cohort, this gives the starting
#'     column index within `treat_var_mat` that corresponds to its first
#'     treatment-period dummy.}
#' @keywords internal
#' @noRd
addDummies <- function(
	df,
	cohorts,
	times,
	N,
	T,
	unit_var,
	time_var,
	resp_var,
	n_cohorts
) {
	# Total number of treated times for all cohorts (later, this will be the
	# total number of treatment effects to estimate)
	num_treats <- 0

	# A list of names of indicator variables for cohorts
	cohort_vars <- character()

	# Matrix of cohort variable indicators
	cohort_var_mat <- matrix(as.integer(NA), nrow = N * T, ncol = n_cohorts)
	# It's better to initialize with 0L if these are counts/indicators
	# cohort_var_mat <- matrix(0L, nrow = N*T, ncol = n_cohorts)
	colnames(cohort_var_mat) <- rep(NA_character_, n_cohorts)

	# Matrix of treatment indicators (just creating one column now to
	# initalize; will delete first column later)
	treat_var_mat <- matrix(as.integer(NA), nrow = N * T, ncol = 1) # Placeholder
	# treat_var_mat <- matrix(0L, nrow = N*T, ncol = 0) # Better start

	# A list of names of variables for cohort/time treatments
	cohort_treat_names <- list()

	# Indices of first treatment time for each cohort within the treat_var_mat block
	current_first_ind_val <- 1
	first_inds <- integer(n_cohorts) # pre-allocate

	# Ensure cohorts list is not empty if n_cohorts > 0
	if (n_cohorts > 0 && length(cohorts) == 0) {
		stop("n_cohorts > 0 but the cohorts list is empty in addDummies.")
	}
	if (n_cohorts == 0 && length(cohorts) > 0) {
		stop("n_cohorts == 0 but the cohorts list is not empty in addDummies.")
	}

	for (i in 1:n_cohorts) {
		# Time of first treatment for this cohort
		y1_treat_i <- as.integer(names(cohorts)[i])
		stopifnot(y1_treat_i <= max(times))
		# Cohort must start treatment strictly after the first period
		stopifnot(y1_treat_i > times[1])

		treated_times_i <- times[times >= y1_treat_i] # Actual time values

		# How many treated times are there?
		n_treated_times_for_cohort_i <- length(treated_times_i)
		stopifnot(n_treated_times_for_cohort_i <= T)

		if (n_treated_times_for_cohort_i == 0 && length(cohorts[[i]]) > 0) {
			# This means a cohort is defined but has no post-treatment periods within 'times'
			# This shouldn't happen if cohorts are defined by times >= times[2]
			# and times contains at least two periods.
			warning(paste(
				"Cohort",
				names(cohorts)[i],
				"has no treated time periods within the observed times."
			))
		}

		first_inds[i] <- current_first_ind_val
		num_treats <- num_treats + n_treated_times_for_cohort_i # This is total across all cohorts so far
		current_first_ind_val <- current_first_ind_val +
			n_treated_times_for_cohort_i

		# Cohort/treatment time variable names
		c_i_names <- paste0("c", i, "_t", treated_times_i) # More descriptive names
		stopifnot(length(c_i_names) == n_treated_times_for_cohort_i)

		cohort_treat_names[[as.character(y1_treat_i)]] <- c_i_names # Name list element by cohort start time

		current_cohort_var_name <- paste0("cohort_", y1_treat_i) # Name cohort by its start time
		cohort_vars <- c(cohort_vars, current_cohort_var_name)
		stopifnot(length(cohort_vars) == i)

		# Units belonging to the current cohort
		current_cohort_units <- cohorts[[i]]
		if (
			length(current_cohort_units) == 0 &&
				n_treated_times_for_cohort_i > 0
		) {
			stop("A defined cohort from names(cohorts) has no units.")
			# This shouldn't happen if idCohorts filters empty cohorts, unless n_cohorts was passed independently.
			# warning(paste("Cohort starting at time", y1_treat_i, "has no units."))
			# If a cohort has no units, its dummies will be all zero.
		}

		treat_vars_i <- genTreatVarsRealData(
			cohort_name = current_cohort_var_name,
			c_t_names = c_i_names,
			N = N,
			T = T,
			n_treated_times = n_treated_times_for_cohort_i,
			unit_vars = df[, unit_var],
			time_vars = df[, time_var],
			cohort = current_cohort_units,
			treated_times = treated_times_i
		)

		# First column is cohort dummy and its name is current_cohort_var_name.
		# Remaining columns are treatment dummies with names c_i_names.

		stopifnot(is.na(cohort_var_mat[, i])) # If initialized with NA
		# if (any(cohort_var_mat[,i] != 0L)) stop("cohort_var_mat not 0L before assignment") # If initialized with 0L
		stopifnot(ncol(treat_vars_i) == n_treated_times_for_cohort_i + 1)

		cohort_var_mat[, i] <- treat_vars_i[, 1]
		colnames(cohort_var_mat)[i] <- current_cohort_var_name

		if (n_treated_times_for_cohort_i > 0) {
			current_treat_dummies <- treat_vars_i[,
				2:(n_treated_times_for_cohort_i + 1),
				drop = FALSE
			]
			if (ncol(treat_var_mat) == 1 && all(is.na(treat_var_mat[, 1]))) {
				# First time, replace placeholder
				treat_var_mat <- current_treat_dummies
			} else {
				treat_var_mat <- cbind(treat_var_mat, current_treat_dummies)
			}
		}
	}

	# After loop, num_treats holds the grand total number of treatment-period dummies
	# current_first_ind_val is now num_treats + 1
	stopifnot(all(first_inds %in% 1:num_treats) || n_cohorts == 0) # first_inds can be empty if n_cohorts=0
	stopifnot(length(first_inds) == n_cohorts)

	stopifnot(length(cohort_treat_names) == n_cohorts)

	if (n_cohorts == 0 || num_treats == 0) {
		# If no cohorts or no treatment effects
		treat_var_mat <- matrix(0L, nrow = N * T, ncol = 0)
	} else {
		stopifnot(
			is.numeric(ncol(treat_var_mat)) | is.integer(ncol(treat_var_mat))
		)
		stopifnot(ncol(treat_var_mat) == num_treats)
	}

	# Add time dummies for all but first time
	time_var_mat <- matrix(0L, nrow = N * T, ncol = T - 1)
	if (T > 1) {
		# Only if there's more than one time period
		time_var_names <- paste0("t_", times[2:T])
		colnames(time_var_mat) <- time_var_names
		for (t_idx in 2:T) {
			# Iterate from the second time period
			actual_time_value <- times[t_idx]
			c_t_inds <- df[, time_var] %in% actual_time_value # Use actual time value
			stopifnot(length(c_t_inds) == N * T)
			stopifnot(sum(c_t_inds) == N) # Each unit observed at this time

			# Add time dummy; column is t_idx - 1 because times[1] is baseline
			time_var_mat[c_t_inds, t_idx - 1] <- 1L
		}
		stopifnot(all(colSums(time_var_mat) == N))
		stopifnot(length(time_var_names) == T - 1)
		stopifnot(ncol(time_var_mat) == T - 1)
	} else {
		# T=1 case
		time_var_names <- character(0)
		# time_var_mat remains 0-col matrix
	}

	stopifnot(
		length(unlist(cohort_treat_names)) == num_treats || num_treats == 0
	)
	stopifnot(ncol(cohort_var_mat) == n_cohorts)

	stopifnot(ncol(treat_var_mat) == num_treats)
	# stopifnot(ncol(treat_var_mat) >= 1) # Can be 0 if num_treats is 0

	stopifnot(length(cohort_vars) == n_cohorts)

	# Center response
	y <- df[, resp_var] - mean(df[, resp_var])

	return(list(
		time_var_mat = time_var_mat,
		cohort_var_mat = cohort_var_mat,
		treat_var_mat = treat_var_mat,
		y = y,
		cohort_treat_names = cohort_treat_names,
		time_var_names = time_var_names,
		cohort_vars = cohort_vars,
		first_inds = first_inds
	))
}

#' Generate Fixed Effect Interactions with Covariates
#'
#' @description
#' Creates interaction terms between covariates (`X_long`) and cohort fixed
#' effects (`cohort_fe`) and time fixed effects (`time_fe`).
#'
#' @param X_long A matrix of time-invariant covariates, with `N*T` rows and `d`
#'   columns. Each unit's `d` covariate values are repeated `T` times.
#' @param cohort_fe A matrix of cohort dummy variables (N*T rows, R columns).
#' @param time_fe A matrix of time dummy variables (N*T rows, T-1 columns).
#' @param N Integer; total number of unique units.
#' @param T Integer; total number of time periods.
#' @param R Integer; total number of treated cohorts.
#' @param d Integer; total number of covariates.
#'
#' @details
#' If `d` (number of covariates) is 0, returns empty matrices.
#' Otherwise:
#' - `X_long_cohort`: Interaction of each covariate with each cohort dummy.
#'   The resulting matrix has `R*d` columns. Columns are ordered such that
#'   the first `d` columns are interactions of `X_long` with `cohort_fe[,1]`,
#'   the next `d` columns with `cohort_fe[,2]`, and so on.
#' - `X_long_time`: Interaction of each covariate with each time dummy.
#'   The resulting matrix has `(T-1)*d` columns, ordered similarly.
#'
#' @return A list containing:
#'   \item{X_long_cohort}{Matrix of cohort-covariate interactions (N*T rows, R*d columns).}
#'   \item{X_long_time}{Matrix of time-covariate interactions (N*T rows, (T-1)*d columns).}
#' @keywords internal
#' @noRd
generateFEInts <- function(X_long, cohort_fe, time_fe, N, T, R, d) {
	# If no covariates are present, return empty matrices.
	if (d == 0) {
		return(list(
			X_long_cohort = matrix(nrow = N * T, ncol = 0),
			X_long_time = matrix(nrow = N * T, ncol = 0)
		))
	}

	# Interact with cohort effects
	X_long_cohort <- matrix(as.numeric(NA), nrow = N * T, ncol = R * d)

	stopifnot(ncol(cohort_fe) == R)
	stopifnot(nrow(cohort_fe) == N * T)
	stopifnot(ncol(time_fe) == (T - 1) || T == 1) # time_fe can be 0-col if T=1
	stopifnot(ncol(X_long) == d)
	stopifnot(is.matrix(X_long))
	stopifnot(!is.data.frame(X_long))

	for (r_idx in 1:R) {
		# iterate R times for R cohorts
		# Notice that these are arranged one cohort at a time, interacted with
		# all covariates
		first_col_r <- (r_idx - 1) * d + 1
		last_col_r <- r_idx * d

		stopifnot(last_col_r - first_col_r + 1 == d)
		stopifnot(all(is.na(X_long_cohort[, first_col_r:last_col_r])))

		# Element-wise multiplication of cohort_fe[,r_idx] (a vector) with each column of X_long
		# This creates d columns for the r_idx-th cohort
		interaction_block_r <- cohort_fe[, r_idx] * X_long
		stopifnot(ncol(interaction_block_r) == d)

		X_long_cohort[, first_col_r:last_col_r] <- interaction_block_r
	}

	stopifnot(all(!is.na(X_long_cohort)))
	stopifnot(nrow(X_long_cohort) == N * T)
	stopifnot(ncol(X_long_cohort) == R * d)

	# Interact with time effects
	if (T > 1) {
		# Only if there are time fixed effects (T-1 > 0)
		X_long_time <- matrix(as.numeric(NA), nrow = N * T, ncol = (T - 1) * d)
		for (t_idx in 1:(T - 1)) {
			# iterate T-1 times for T-1 time dummies
			# Notice that these are arranged one time at a time, interacted with all
			# covariates
			first_col_t <- (t_idx - 1) * d + 1
			last_col_t <- t_idx * d

			stopifnot(last_col_t - first_col_t + 1 == d)
			stopifnot(all(is.na(X_long_time[, first_col_t:last_col_t])))

			interaction_block_t <- time_fe[, t_idx] * X_long
			stopifnot(ncol(interaction_block_t) == d)

			X_long_time[, first_col_t:last_col_t] <- interaction_block_t
		}
		stopifnot(all(!is.na(X_long_time)))
		stopifnot(nrow(X_long_time) == N * T)
		stopifnot(ncol(X_long_time) == (T - 1) * d)
	} else {
		# T=1 case, no time fixed effects
		X_long_time <- matrix(nrow = N * T, ncol = 0)
	}

	return(list(X_long_cohort = X_long_cohort, X_long_time = X_long_time))
}

#' Generate Treatment-Covariate Interactions
#'
#' @description
#' Creates interaction terms between treatment-period dummies (`treat_mat_long`)
#' and time-invariant covariates (`X_long`). Covariates are centered with respect
#' to their cohort means before interaction.
#'
#' @param treat_mat_long A matrix of treatment-period dummy variables (N*T rows,
#'   `n_treats` columns).
#' @param X_long A matrix of time-invariant covariates (N*T rows, `d` columns).
#'   Each unit's `d` covariate values are repeated `T` times.
#' @param n_treats Integer; total number of unique treatment-period dummies
#'   (columns in `treat_mat_long`).
#' @param cohort_fe A matrix of cohort dummy variables (N*T rows, R columns).
#'   Used to identify units within each cohort for centering covariates.
#' @param N Integer; total number of unique units.
#' @param T Integer; total number of time periods.
#' @param R Integer; total number of treated cohorts.
#' @param d Integer; total number of covariates.
#' @param N_UNTREATED Integer; number of never-treated units.
#'
#' @details
#' If `d` (number of covariates) is 0, returns an empty matrix.
#' Otherwise, covariates `X_long` are first centered. For each cohort (including
#' the never-treated group, identified by `rowSums(cohort_fe) == 0`), the mean
#' of each covariate *within that cohort* is subtracted from the covariate values
#' of units in that cohort.
#' Then, each column of `treat_mat_long` (a specific treatment-period dummy) is
#' interacted with each column of the centered covariates.
#' The resulting matrix `X_long_treat` has `d * n_treats` columns. Columns are
#' ordered such that the first `d` columns are interactions of `treat_mat_long[,1]`
#' with the centered covariates, the next `d` with `treat_mat_long[,2]`, and so on.
#'
#' @return A matrix of treatment-covariate interactions (N*T rows, d*`n_treats`
#'   columns). Returns an empty matrix if `d=0`.
#' @keywords internal
#' @noRd
genTreatInts <- function(
	treat_mat_long,
	X_long,
	n_treats,
	cohort_fe,
	N,
	T,
	R,
	d,
	N_UNTREATED
) {
	# If there are no covariates, return an empty matrix.
	if (d == 0) {
		return(matrix(nrow = N * T, ncol = 0))
	}

	stopifnot(ncol(cohort_fe) == R)
	stopifnot(ncol(treat_mat_long) == n_treats)

	stopifnot(is.numeric(d) || is.integer(d))
	stopifnot(d >= 1)

	stopifnot(is.numeric(n_treats) || is.integer(n_treats))
	stopifnot(n_treats >= 1)

	stopifnot(is.numeric(N) || is.integer(N))
	stopifnot(N >= 2)

	stopifnot(is.numeric(T) || is.integer(T))
	stopifnot(T >= 2)
	# Interact with covariates
	X_long_treat <- matrix(as.numeric(NA), nrow = N * T, ncol = d * n_treats)

	# First need to center covariates with respect to cohort means
	X_long_centered <- matrix(as.numeric(NA), nrow = N * T, ncol = d)

	# Never treated group
	unt_row_inds <- which(rowSums(cohort_fe) == 0)

	# Assume balanced panel
	stopifnot(length(unt_row_inds) == N_UNTREATED * T)

	stopifnot(all(is.na(X_long_centered[unt_row_inds, ])))

	X_long_centered[unt_row_inds, ] <- scale(
		X_long[unt_row_inds, , drop = FALSE],
		center = TRUE,
		scale = FALSE
	)

	for (r in 1:R) {
		# Notice that these are arranged column-wise one cohort at a time,
		# interacted with all treatment effects

		cohort_r_inds <- which(cohort_fe[, r] == 1)

		new_mat <- scale(
			X_long[cohort_r_inds, , drop = FALSE],
			center = TRUE,
			scale = FALSE
		)

		stopifnot(all(!is.na(new_mat)))
		stopifnot(all(is.na(X_long_centered[cohort_r_inds, ])))

		X_long_centered[cohort_r_inds, ] <- new_mat
	}
	# Final index should correspond to last row
	stopifnot(all(!is.na(X_long_centered)))

	for (i in 1:n_treats) {
		first_col_i <- (i - 1) * d + 1
		last_col_i <- i * d

		stopifnot(all(is.na(X_long_treat[, first_col_i:last_col_i])))

		X_long_treat[, first_col_i:last_col_i] <- treat_mat_long[, i] *
			X_long_centered
	}

	stopifnot(all(!is.na(X_long_treat)))
	stopifnot(nrow(X_long_treat) == N * T)
	stopifnot(ncol(X_long_treat) == n_treats * d)

	return(X_long_treat)
}

#' Generate Full Design Matrix with Interactions (`X_ints`)
#'
#' @description
#' Constructs the complete design matrix `X_ints` used in the ETWFE model.
#' This matrix includes cohort fixed effects, time fixed effects, covariates,
#' treatment-period dummies, and all their relevant two-way interactions
#' (cohort-covariate, time-covariate, treatment-covariate).
#'
#' @param cohort_fe Matrix of cohort dummy variables (N*T rows, R columns).
#' @param time_fe Matrix of time dummy variables (N*T rows, T-1 columns).
#' @param X_long Matrix of time-invariant covariates (N*T rows, `d` columns).
#' @param treat_mat_long Matrix of treatment-period dummy variables (N*T rows,
#'   `num_treats` columns).
#' @param N Integer; total number of unique units.
#' @param R Integer; total number of treated cohorts.
#' @param T Integer; total number of time periods.
#' @param d Integer; total number of covariates.
#' @param N_UNTREATED Integer; number of never-treated units.
#' @param p Integer; the expected total number of columns in the final `X_ints` matrix.
#'
#' @details
#' The matrix `X_ints` is constructed by column-binding the following blocks in order:
#' \enumerate{
#'   \item `cohort_fe` (`R` columns)
#'   \item `time_fe` (`T-1` columns, or 0 if `T=1`)
#'   \item `X_long` (`d` columns, or 0 if `d=0`)
#'   \item Cohort-covariate interactions from `generateFEInts` (`R*d` columns)
#'   \item Time-covariate interactions from `generateFEInts` (`(T-1)*d` columns)
#'   \item `treat_mat_long` (`num_treats` columns)
#'   \item Treatment-covariate interactions from `genTreatInts` (`d*num_treats` columns)
#' }
#' The function checks if the number of columns in the resulting matrix matches `p`.
#'
#' @return The complete design matrix `X_ints` (`N*T` rows, `p` columns).
#' @keywords internal
#' @noRd
genXintsData <- function(
	cohort_fe,
	time_fe,
	X_long,
	treat_mat_long,
	N,
	R,
	T,
	d,
	N_UNTREATED,
	p
) {
	stopifnot(is.matrix(X_long))
	# X_long may be an empty matrix if d == 0.
	X_int <- cbind(cohort_fe, time_fe, X_long)
	stopifnot(ncol(X_int) == R + T - 1 + d)

	# Generate interactions of X with time and cohort fixed effects
	# Note: columns of X_long_cohort are arranged in blocks of size d for one
	# cohort at a time (that is, the first d columns are the interactions of
	# all d features with the indicator variables for the first block, and
	# so on). Similarly, the columns of X_long_time are arranged in T - 1
	# blocks of size d.
	stopifnot(ncol(cohort_fe) == R)
	res <- generateFEInts(X_long, cohort_fe, time_fe, N, T, R, d)

	X_long_cohort <- res$X_long_cohort
	X_long_time <- res$X_long_time

	X_int <- cbind(X_int, X_long_cohort, X_long_time, treat_mat_long)

	rm(res)

	stopifnot(
		is.integer(ncol(treat_mat_long)) | is.numeric(ncol(treat_mat_long))
	)
	stopifnot(ncol(treat_mat_long) >= 1)

	stopifnot(ncol(cohort_fe) == R)

	# Generate interactions between treatment effects and X (if any)
	X_long_treat <- genTreatInts(
		treat_mat_long = treat_mat_long,
		X_long = X_long,
		n_treats = ncol(treat_mat_long),
		cohort_fe = cohort_fe,
		N = N,
		T = T,
		R = R,
		d = d,
		N_UNTREATED = N_UNTREATED
	)

	# The first R columns are the cohort fixed effects in order
	# Next T are time fixed effects in order
	# Next d are X
	# Next d*R are cohort effects interacted with X. (The first d of these
	# are first cohort effects interacted with X, and so on until the Rth
	# cohort.)
	# Next d*T are time effects interacted with X (similarly to the above,
	# the first d are first time effects interacted with X, and so on until
	# the Tth time).
	# Next num_treats = R*T - R*(R + 1)/2 columns are base treatment effects
	# (for each cohort and time)
	# Finally, the next num_treats*d are interactions of all of these
	# treatment effects over time--first d are the first column of
	# treat_mat_long interacted with all of the columns of X, and so on.
	X_int <- cbind(X_int, X_long_treat)

	if (ncol(X_int) != p) {
		stop(paste(
			"ncol(X_int) =",
			ncol(X_int),
			", p =",
			p,
			", R =",
			R,
			", T =",
			T,
			", d =",
			d
		))
	}

	stopifnot(ncol(X_int) == p)
	stopifnot(nrow(X_int) == N * T)

	return(X_int)
}
