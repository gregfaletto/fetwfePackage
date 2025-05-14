#' @import glmnet
#' @importFrom Matrix bdiag
#' @importFrom expm sqrtm
#' @importFrom stats qnorm predict coef model.matrix setNames

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

	covariate_mat <- as.matrix(data[, covs]) # Covariates

	rm(ret)
	rm(data)

	#### Create matrix with all interactions
	p <- getP(R = R, T = T, d = d, num_treats = num_treats)

	stopifnot(ncol(cohort_var_mat) == R)

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

#' Core Estimation Logic for Fused Extended Two-Way Fixed Effects
#'
#' @description
#' This function implements the core estimation steps of the FETWFE methodology.
#' It takes a pre-processed design matrix and response, applies transformations
#' for fusion penalties, handles variance components, performs bridge regression,
#' selects the optimal penalty via BIC, and calculates treatment effects and
#' their standard errors.
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
#'   \item **Overall ATT Calculation:** Calls `getTeResults2` to calculate the
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
fetwfe_core <- function(
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
	# ... (rest of the function code)
	# NOTE: The existing stopifnot calls provide good internal checks.
	# Consider if any should be user-facing `stop()` messages if triggered by bad inputs
	# that somehow bypass `checkFetwfeInputs`.

	# ... (Function body as provided)

	# Ensure all paths through the function that return have the documented list structure.
	# For instance, the early returns when `lambda_star_model_size == 0` or
	# `length(sel_treat_inds_shifted) == 0` should also return all documented fields,
	# filling with NA or appropriate defaults where necessary.
	# The current code for these early returns seems to do this.

	# Example for an early return structure (conceptual):
	# if (some_early_exit_condition) {
	#   return(list(
	#     in_sample_att_hat = 0, in_sample_att_se = ifelse(q < 1, 0, NA), ...,
	#     N = N, T = T, R = R, d = d, p = p # Ensure all documented fields are present
	#   ))
	# }
	R <- length(in_sample_counts) - 1 # This needs to be defined before the early exits
	# if c_names relies on it.
	# c_names is used in catt_df_to_ret.
	c_names <- names(in_sample_counts)[2:(R + 1)]

	# (Actual function body as provided in the original file)
	# ...
	# Check inputs
	#
	#

	stopifnot(N >= 2) # bare minimum, 2 units at 2 times

	stopifnot(T >= 2) # bare minimum, 2 units at 2 times

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

	stopifnot(sum(in_sample_counts) == N)
	stopifnot(all(in_sample_counts >= 0))
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

	# R <- length(in_sample_counts) - 1 # Moved up
	stopifnot(R >= 1)
	stopifnot(R <= T - 1)

	indep_count_data_available <- FALSE
	if (any(!is.na(indep_counts))) {
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
		if (any(indep_counts <= 0)) {
			stop(
				"At least one cohort in the independent count data has 0 members"
			)
		}
		indep_count_data_available <- TRUE
	}

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

	stopifnot(is.logical(verbose))
	stopifnot(length(verbose) == 1)

	stopifnot(is.numeric(alpha))
	stopifnot(length(alpha) == 1)
	stopifnot(alpha > 0)
	stopifnot(alpha < 1)

	stopifnot(is.logical(add_ridge))
	stopifnot(length(add_ridge) == 1)

	#
	#
	# Step 1: change coordinates of data so that regular bridge regression
	# penalty applied to transformed dataframe results in FETWFE penalty applied
	# to original data
	#
	#

	if (verbose) {
		message("Transforming matrix...")
	}

	# Transform matrix (change of coordinates so that fitting regular bridge
	# regression results in FETWFE fusion penalties)
	X_mod <- transformXintImproved(
		X_ints,
		N = N,
		T = T,
		R = R,
		d = d,
		num_treats = num_treats,
		first_inds = first_inds
	)

	#
	#
	# Step 2: get (known or estimated) covariance matrix within observations
	# (due to umit-level random effects) and pre-multiply X and y by
	# inverse square root matrix
	#
	#

	if (verbose) {
		message("Getting omega sqrt inverse estimate...")
		t0 <- Sys.time()
	}

	if (is.na(sig_eps_sq) | is.na(sig_eps_c_sq)) {
		# Get omega_sqrt_inv matrix to multiply y and X_mod by on the left
		omega_res <- estOmegaSqrtInv(
			y,
			X_ints,
			N = N,
			T = T,
			p = p
		)

		sig_eps_sq <- omega_res$sig_eps_sq
		sig_eps_c_sq <- omega_res$sig_eps_c_sq

		rm(omega_res)

		if (verbose) {
			message("Done! Time to estimate noise variances:")
			message(Sys.time() - t0)
			t0 <- Sys.time()
		}
	}

	stopifnot(!is.na(sig_eps_sq) & !is.na(sig_eps_c_sq))

	Omega <- diag(rep(sig_eps_sq, T)) + matrix(sig_eps_c_sq, T, T)

	Omega_sqrt_inv <- expm::sqrtm(solve(Omega))

	if (verbose) {
		message("Time to get sqrt inverse matrix:")
		message(Sys.time() - t0)
	}

	y_final <- kronecker(diag(N), sqrt(sig_eps_sq) * Omega_sqrt_inv) %*% y
	X_final <- kronecker(diag(N), sqrt(sig_eps_sq) * Omega_sqrt_inv) %*% X_mod

	#
	#
	# Optional: if using ridge regularization on untransformed coefficients,
	# add those rows now
	#
	#

	X_final_scaled <- my_scale(X_final)
	scale_center <- attr(X_final_scaled, "scaled:center")
	scale_scale <- attr(X_final_scaled, "scaled:scale")

	if (add_ridge) {
		# Add rows to X_final. First need to get D^{-1}:
		D_inverse <- genFullInvFusionTransformMat(
			first_inds = first_inds,
			T = T,
			R = R,
			d = d,
			num_treats = num_treats
		)

		stopifnot(ncol(D_inverse) == ncol(X_final))
		stopifnot(ncol(D_inverse) == ncol(X_final_scaled))

		# Now add rows
		lambda_ridge <- 0.00001 *
			(sig_eps_sq + sig_eps_c_sq) *
			sqrt(p / (N * T))

		X_final_scaled <- rbind(X_final_scaled, sqrt(lambda_ridge) * D_inverse)
		y_final <- c(y_final, rep(0, nrow(D_inverse)))

		stopifnot(nrow(X_final_scaled) == length(y_final))
		stopifnot(nrow(X_final_scaled) == N * T + p)
	}

	#
	#
	# Step 3: get cohort-specific sample proportions (estimated treatment
	# probabilities)
	#
	#

	cohort_probs <- in_sample_counts[2:(R + 1)] /
		sum(in_sample_counts[2:(R + 1)])

	stopifnot(all(!is.na(cohort_probs)))
	stopifnot(all(cohort_probs >= 0))
	stopifnot(all(cohort_probs <= 1))
	stopifnot(length(cohort_probs) == R)
	stopifnot(abs(sum(cohort_probs) - 1) < 10^(-6))

	cohort_probs_overall <- in_sample_counts[2:(R + 1)] / N

	stopifnot(
		abs(1 - sum(cohort_probs_overall) - in_sample_counts[1] / N) < 10^(-6)
	)

	if (indep_count_data_available) {
		indep_cohort_probs <- indep_counts[2:(R + 1)] /
			sum(indep_counts[2:(R + 1)])

		stopifnot(all(!is.na(indep_cohort_probs)))
		stopifnot(all(indep_cohort_probs >= 0))
		stopifnot(all(indep_cohort_probs <= 1))
		stopifnot(length(indep_cohort_probs) == R)
		stopifnot(abs(sum(indep_cohort_probs) - 1) < 10^(-6))

		indep_cohort_probs_overall <- indep_counts[2:(R + 1)] / N

		stopifnot(
			abs(
				1 -
					sum(
						indep_cohort_probs_overall
					) -
					indep_counts[1] / N
			) <
				10^(-6)
		)
	} else {
		indep_cohort_probs <- NA
		indep_cohort_probs_overall <- NA
	}

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
		fit,
		N = N,
		T = T,
		p = p,
		X_mod = X_mod,
		y = y,
		scale_center = scale_center,
		scale_scale = scale_scale
	)

	theta_hat <- res$theta_hat # This includes intercept
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
	if (lambda_star_model_size <= 1 && all(theta_hat[2:(p + 1)] == 0)) {
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
			theta_hat = theta_hat, # Includes intercept
			beta_hat = rep(0, p), # Slopes are all zero
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
			p = p
		))
	}

	# intercept
	# eta_hat <- theta_hat[1] # theta_hat from getBetaBIC already has intercept as first element

	# estimated coefficients (slopes in transformed space)
	theta_hat_slopes <- theta_hat[2:(p + 1)]

	# Indices of selected features in transformed feature space (among slopes)
	sel_feat_inds <- which(theta_hat_slopes != 0)

	sel_treat_inds <- sel_feat_inds[sel_feat_inds %in% treat_inds]

	stopifnot(length(sel_treat_inds) == length(unique(sel_treat_inds)))
	stopifnot(length(sel_treat_inds) <= length(sel_feat_inds))
	stopifnot(length(sel_treat_inds) <= length(treat_inds))
	stopifnot(is.integer(sel_treat_inds) | is.numeric(sel_treat_inds))

	# Shift sel_treat_inds to be relative to the start of the treat_inds block
	# This seems to be what sel_treat_inds_shifted intends
	# The current sel_treat_inds_shifted calculation appears complex and might be error-prone.
	# A simpler way:
	# 1. Get theta_hat_slopes[treat_inds] -> these are the transformed treatment coefficients
	# 2. Find which of these are non-zero: `which(theta_hat_slopes[treat_inds] != 0)` -> these are indices *within* the treat_inds block.
	# Let's call this `sel_treat_inds_relative_to_block`.
	# `sel_treat_inds` itself contains the global indices of selected treatment features.
	# So `theta_hat_slopes[sel_treat_inds]` are the non-zero transformed treatment coefs.

	theta_hat_treat_block_transformed = theta_hat_slopes[treat_inds]
	sel_treat_inds_shifted <- which(theta_hat_treat_block_transformed != 0) # these are 1 to num_treats

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

		# Need to untransform theta_hat_slopes to get beta_hat for consistency
		beta_hat_early_exit <- untransformCoefImproved(
			theta_hat_slopes, # Pass slopes only
			first_inds,
			T = T,
			R = R,
			p = p,
			d = d,
			num_treats = num_treats
		)
		if (add_ridge) {
			beta_hat_early_exit <- beta_hat_early_exit * (1 + lambda_ridge)
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
			theta_hat = theta_hat, # Full theta_hat with intercept
			beta_hat = beta_hat_early_exit, # Untransformed slopes
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
			p = p
		))
	}

	#
	#
	# Step 5: transform estimated coefficients back to original feature
	# space
	#
	#

	beta_hat <- untransformCoefImproved(
		theta_hat_slopes, # Pass slopes only
		first_inds,
		T = T,
		R = R,
		p = p,
		d = d,
		num_treats = num_treats
	)

	# If using ridge regularization, multiply the "naive" estimated coefficients
	# by 1 + lambda_ridge, similar to suggestion in original elastic net paper.
	if (add_ridge) {
		beta_hat <- beta_hat * (1 + lambda_ridge)
	}

	# Get actual estimated treatment effects (in original, untransformed space)
	tes <- beta_hat[treat_inds]

	stopifnot(length(tes) == num_treats)
	# Checks based on transformed coefficients (theta_hat_slopes)
	stopifnot(all(theta_hat_slopes[treat_inds][sel_treat_inds_shifted] != 0))
	stopifnot(all(
		theta_hat_slopes[treat_inds][setdiff(
			1:num_treats,
			sel_treat_inds_shifted
		)] ==
			0
	))

	stopifnot(length(first_inds) == R)
	stopifnot(max(first_inds) <= num_treats)

	stopifnot(length(sel_feat_inds) > 0) # sel_feat_inds are indices in theta_hat_slopes
	stopifnot(length(sel_treat_inds_shifted) > 0) # sel_treat_inds_shifted are indices within the treat_inds block of theta_hat_slopes

	#
	#
	# Step 6: calculate cohort-specific treatment effects and standard
	# errors
	#
	#

	res <- getCohortATTsFinal(
		X_final = X_final, # This is X_mod * GLS_transform_matrix
		sel_feat_inds = sel_feat_inds, # Indices of non-zero elements in theta_hat_slopes
		treat_inds = treat_inds, # Global indices for treatment effects
		num_treats = num_treats,
		first_inds = first_inds,
		sel_treat_inds_shifted = sel_treat_inds_shifted, # Indices (1 to num_treats) of non-zero transformed treat. coefs.
		c_names = c_names,
		tes = tes, # Untransformed treatment effect estimates (beta_hat[treat_inds])
		sig_eps_sq = sig_eps_sq,
		R = R,
		N = N,
		T = T,
		fused = TRUE, # This parameter might be redundant if this function is only for fused
		calc_ses = q < 1,
		p = p, # Total number of original parameters (columns in X_ints)
		alpha = alpha
	)

	cohort_te_df <- res$cohort_te_df
	cohort_tes <- res$cohort_tes
	cohort_te_ses <- res$cohort_te_ses
	psi_mat <- res$psi_mat
	gram_inv <- res$gram_inv
	d_inv_treat_sel <- res$d_inv_treat_sel
	calc_ses <- res$calc_ses

	rm(res)

	if (calc_ses) {
		stopifnot(nrow(d_inv_treat_sel) == num_treats)
		stopifnot(ncol(d_inv_treat_sel) == length(sel_treat_inds_shifted))
	}

	#
	#
	# Step 7: calculate overall average treatment effect on treated units
	#
	#

	# Get overal estimated ATT!
	# theta_hat_treat_sel needs to be the selected non-zero *transformed* treatment coefficients
	# sel_treat_inds contains global indices of selected transformed features that are treatment effects
	theta_hat_treat_sel_for_att <- theta_hat_slopes[sel_treat_inds]

	in_sample_te_results <- getTeResults2(
		sig_eps_sq = sig_eps_sq,
		N = N,
		T = T,
		R = R,
		num_treats = num_treats,
		cohort_tes = cohort_tes, # CATTs (point estimates)
		cohort_probs = cohort_probs, # In-sample pi_r | treated
		psi_mat = psi_mat,
		gram_inv = gram_inv,
		sel_treat_inds_shifted = sel_treat_inds_shifted,
		tes = tes, # Untransformed treatment effect estimates beta_hat[treat_inds]
		d_inv_treat_sel = d_inv_treat_sel,
		cohort_probs_overall = cohort_probs_overall, # In-sample pi_r (unconditional on treated)
		first_inds = first_inds,
		theta_hat_treat_sel = theta_hat_treat_sel_for_att, # Selected non-zero transformed treat coefs
		calc_ses = calc_ses,
		indep_probs = FALSE
	)

	in_sample_att_hat <- in_sample_te_results$att_hat
	in_sample_att_se <- in_sample_te_results$att_te_se
	in_sample_att_se_no_prob <- in_sample_te_results$att_te_se_no_prob

	if (indep_count_data_available) {
		indep_te_results <- getTeResults2(
			sig_eps_sq = sig_eps_sq,
			N = N,
			T = T,
			R = R,
			num_treats = num_treats,
			cohort_tes = cohort_tes,
			cohort_probs = indep_cohort_probs, # indep pi_r | treated
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			tes = tes,
			d_inv_treat_sel = d_inv_treat_sel,
			cohort_probs_overall = indep_cohort_probs_overall, # indep pi_r (unconditional)
			first_inds = first_inds,
			theta_hat_treat_sel = theta_hat_treat_sel_for_att,
			calc_ses = calc_ses,
			indep_probs = TRUE
		)
		indep_att_hat <- indep_te_results$att_hat
		indep_att_se <- indep_te_results$att_te_se
		# indep_att_se_no_prob <- indep_te_results$att_te_se_no_prob # This was commented out
	} else {
		indep_att_hat <- NA
		indep_att_se <- NA
		# indep_att_se_no_prob <- NA # Keep commented for consistency
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
		theta_hat = theta_hat, # Full theta_hat (with intercept)
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
		p = p
	))
}


#-------------------------------------------------------------------------------
# Helper Functions for Data Processing
#-------------------------------------------------------------------------------

#' Identify Cohorts and Basic Panel Structure
#'
#' @description
#' This function processes the input panel data to identify treatment cohorts,
#' unique units, and time periods. It also removes units that were treated in the
#' very first time period, as these units cannot be used for identifying
#' pre-treatment trends for themselves. The treatment variable is removed after
#' cohorts are identified.
#'
#' @param df A `data.frame` object representing the panel data. Contains columns for response, time
#' period variable (a single categorical or numeric/integer variable indicating the time period
#' for each observation), unit variable (a single categorical variable
#' indicating which unit each observation is), treatment (a binary variable;
#' 1 if unit is treated at that time and 0 if not; treatment must be an
#' absorbing state in this model); covariates (which are fixed over time)
#' @param time_var Character string; name of the time variable column.
#' @param unit_var Character string; name of the unit identifier column.
#' @param treat_var Character string; name of the treatment indicator column.
#' @param covs Character vector; names of covariate columns (used for subsetting `df`).
#'
#' @details
#' The function iterates through each unit to determine its cohort assignment based
#' on the first period it is observed as treated. It assumes treatment is an
#' absorbing state.
#' Units treated in the first time period (`times[1]`) are removed from the dataset
#' and a warning is issued.
#' Cohorts are stored in a list, where names are the treatment adoption times and
#' values are character vectors of unit identifiers belonging to that cohort.
#' The list of cohorts is ordered by adoption time, and the first time period is
#' ensured to have an empty cohort (representing never-treated or baseline).
#'
#' @return A list containing:
#'   \item{df}{The processed `data.frame` with first-period-treated units removed
#'     and the `treat_var` column dropped. Contains only `time_var`, `unit_var`,
#'     `response` (implicitly, as it's not dropped), and `covs`.}
#'   \item{cohorts}{A list where each element is a character vector of unit IDs
#'     belonging to a specific treatment cohort. The names of the list elements
#'     are the times of treatment adoption. The first cohort (earliest adoption
#'     time after the first period) is listed first.}
#'   \item{units}{A character vector of unique unit identifiers remaining after
#'     processing.}
#'   \item{times}{A numeric or integer vector of unique, sorted time periods
#'     present in the original data.}
#' @keywords internal
#' @noRd
idCohorts <- function(df, time_var, unit_var, treat_var, covs) {
	stopifnot(time_var %in% colnames(df))
	stopifnot(unit_var %in% colnames(df))
	stopifnot(treat_var %in% colnames(df))
	# stopifnot(all(covs %in% colnames(df))) # covs might be empty, or not all present if some were
	# factors

	# Form design matrix
	units <- unique(df[, unit_var])
	N <- length(units)
	times <- sort(unique(df[, time_var]))
	T <- length(times)

	# Variable to identify cohorts
	cohorts <- list()
	for (t in 1:T) {
		cohorts[[t]] <- character()
	}
	names(cohorts) <- times

	for (s in units) {
		df_s <- df[df[, unit_var] == s, ]
		# Assume this is a balanced panel
		if (nrow(df_s) != T) {
			stop(paste(
				"Panel does not appear to be balanced (unit",
				s,
				"does not have exactly T observations for T =",
				T
			))
		}

		if (any(df_s[, treat_var] == 1)) {
			# Identify first year of treatment (and cohort)
			# Ensure df_s is sorted by time before finding min
			df_s_sorted <- df_s[order(df_s[, time_var]), ]
			treat_year_s_ind_in_sorted_times <- min(which(
				df_s_sorted[, treat_var] == 1
			))
			actual_treat_time <- df_s_sorted[
				treat_year_s_ind_in_sorted_times,
				time_var
			]

			# Make sure treatment is absorbing
			# Check from the actual_treat_time onwards in the original df for unit s
			original_unit_times_from_treatment <- df[
				(df[, unit_var] == s) & (df[, time_var] >= actual_treat_time),
				treat_var
			]
			if (any(original_unit_times_from_treatment != 1)) {
				stop(paste(
					"Treatment does not appear to be an absorbing state for unit",
					s
				))
			}
			cohorts[[as.character(actual_treat_time)]] <- c(
				cohorts[[as.character(actual_treat_time)]],
				s
			)
		}
	}

	stopifnot(length(unlist(cohorts)) <= N)
	# Keep only cohorts that actually have units
	cohorts <- cohorts[lengths(cohorts) > 0]

	# Need at least one untreated period, so have to omit units that were
	# treated in the very first time period
	first_time_val_char <- as.character(times[1])
	if (first_time_val_char %in% names(cohorts)) {
		first_year_cohort_units <- cohorts[[first_time_val_char]]
		if (length(first_year_cohort_units) > 0) {
			df <- df[!(df[, unit_var] %in% first_year_cohort_units), ]
			units <- unique(df[, unit_var]) # Update units
			if (length(first_year_cohort_units) > 0) {
				# N is original N
				warning(paste(
					length(first_year_cohort_units),
					"units were removed because they were treated in the first time period:",
					paste(first_year_cohort_units, collapse = ", ")
				))
			}
		}
		cohorts[[first_time_val_char]] <- NULL # Remove this cohort entry
	}

	if (length(units) == 0) {
		stop(
			"All units were treated in the first time period or no units remain after filtering; estimating treatment effects is not possible"
		)
	}

	N <- length(units) # Update N to reflect remaining units

	# Treatment no longer needed
	df <- df[, colnames(df) != treat_var]

	# Make sure there is an empty cohort for the first time
	cohorts[[as.character(times[1])]] <- character()

	# Order cohorts in order of times
	if (length(cohorts) > 0) {
		cohorts <- cohorts[order(as.numeric(names(cohorts)))]
	}

	# If after removing first-period treated, there are no treated cohorts left:
	if (length(cohorts) == 0) {
		# This implies all treated units were treated in the first period.
		# The code in prepXints would later error if R < 2.
		# This function's main job is to return the filtered df and cohort structure.
		stop("all units appear to have been treated in the first period")
	}

	# This should have been the first cohort
	stopifnot(length(cohorts[[1]]) == 0)

	cohorts <- cohorts[-1]

	stopifnot(all(lengths(cohorts) >= 1))
	stopifnot(length(cohorts) <= T)
	stopifnot(length(unlist(cohorts)) <= N)

	return(list(df = df, cohorts = cohorts, units = units, times = times))
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
			# df_s <- df[df[, unit_var] == s, ] # Not needed if iterating by index
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
	# ... (Function body as provided)
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
	if (n_cohorts > 0 && length(cohorts) != n_cohorts) {
		# This might happen if idCohorts filters out all treated units.
		# prepXints already checks R >= 1 (where R is length(cohorts))
		# and R < 2.
		# This stop might be too strict if R is derived from length(cohorts) earlier.
		# stop(paste("n_cohorts is", n_cohorts, "but length(cohorts) is", length(cohorts)))
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
		# stopifnot(n_treated_times_for_cohort_i == max(times) - y1_treat_i + 1) # This might fail if times are not consecutive integers

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
			# A defined cohort from names(cohorts) has no units
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
		# df <- data.frame(df, treat_vars_i) # Original code modifies df, but not strictly needed if matrices are goal
	}

	# After loop, num_treats holds the grand total number of treatment-period dummies
	# current_first_ind_val is now num_treats + 1
	stopifnot(all(first_inds %in% 1:num_treats) || n_cohorts == 0) # first_inds can be empty if n_cohorts=0
	stopifnot(length(first_inds) == n_cohorts)

	stopifnot(length(cohort_treat_names) == n_cohorts)
	# names(cohort_treat_names) <- names(cohorts) # Already done inside loop if using y1_treat_i

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

#-------------------------------------------------------------------------------
# Transformation Matrices and Coefficient Transformations
#-------------------------------------------------------------------------------

#' Get Indices of First Treatment Effects for Each Cohort
#'
#' @description
#' Calculates the starting indices of treatment effect parameters for each cohort
#' within a concatenated block of all treatment effects. This is used for constructing
#' fusion penalty matrices.
#'
#' @param R Integer; the number of treated cohorts.
#' @param T Integer; the total number of time periods.
#'
#' @details
#' Assumes treatment for cohort `r` (where `r` is 1-indexed for calculation,
#' corresponding to actual adoption times `times[r+1]`) starts at `times[r+1]`
#' and continues until `times[T]`.
#' The number of treatment effects for cohort `r` is `T - (r+1) + 1 = T - r`.
#' The `num_treats` is the sum of these counts.
#' The formula for the starting index of the `r`-th cohort's treatment effects
#' (1-indexed `r` from 1 to `R`):
#' `f_inds[r] = 1 + sum_{k=1}^{r-1} (T - k) = 1 + (r-1)T - (r-1)r/2`.
#'
#' Gemini check:
#' The paper's formula `1 + (r - 1)*(2*T - r)/2` seems to be for a slightly
#' different definition of `num_treats` or indexing. This function should align
#' with how `treat_var_mat` is constructed in `addDummies` and `num_treats`
#' is calculated there.
#' Let's re-evaluate the formula based on typical DiD setups:
#' Cohort 1 (starts at time 2) has T-1 effects.
#' Cohort 2 (starts at time 3) has T-2 effects.
#' ...
#' Cohort R (starts at time R+1) has T-R effects.
#' `first_inds[1] = 1`
#' `first_inds[2] = (T-1) + 1`
#' `first_inds[3] = (T-1) + (T-2) + 1`
#' So, `first_inds[r] = 1 + sum_{j=1}^{r-1} (T-j)`.
#'
#' The current code uses: `f_inds[r] <- 1 + (r - 1)*(2*T - r)/2`.
#' Let's test this formula:
#' r=1: 1
#' r=2: 1 + 1*(2T-2)/2 = 1 + T - 1 = T
#' r=3: 1 + 2*(2T-3)/2 = 1 + 2T - 3 = 2T - 2
#' This implies the number of effects for cohort 1 is T-1 (indices 1 to T-1).
#' Number of effects for cohort 2 is (2T-2) - T = T-2.
#' This seems to match the sequence T-1, T-2, ..., T-R.
#' Total number of treats (`num_treats`) from `getNumTreats(R,T)` is `T*R - R*(R+1)/2`.
#' This is `sum_{j=1 to R} (T-j)`. This is consistent.
#'
#' @return An integer vector of length `R`, where `f_inds[r]` is the 1-based
#'   starting index of the `r`-th cohort's treatment effects in the combined block.
#' @keywords internal
#' @noRd
getFirstInds <- function(R, T) {
	# Let's identify the indices of the first treatment effects for each cohort.
	# The first one is index 1, then the second one is index (T - 1) + 1 = T,
	# then the third one is (T - 1) + (T - 2) + 1 = 2*T - 2. In general, for
	# r > 1 the rth one will occur at index
	#
	# (T - 1) + (T - 2) + ... + (T - (r - 1)) + 1
	# = 1 + (r - 1)*(T - 1 + T - r + 1)/2
	# = 1 + (r - 1)*(2*T - r)/2.
	#
	# (Looks like the formula works for r = 1 too.)

	n_treats <- getNumTreats(R = R, T = T)

	f_inds <- integer(R)
	if (R == 0) return(f_inds) # No cohorts, no first_inds

	for (r in 1:R) {
		f_inds[r] <- 1 + (r - 1) * (2 * T - r) / 2
	}
	stopifnot(all(f_inds <= n_treats))
	stopifnot(f_inds[1] == 1)
	# Last cohort has T - R treatment effects to estimate. So last first_ind
	# should be at position num_treats - (T - R) + 1 = num_treats - T + R + 1.
	stopifnot(f_inds[R] == n_treats - T + R + 1)

	# Additional checks from Gemini below

	stopifnot(all(f_inds <= n_treats) || R == 0)
	if (R > 0) stopifnot(f_inds[1] == 1)

	# Last first_ind: f_inds[R] = 1 + sum_{j=1}^{R-1} (T-j)
	# Total effects = sum_{j=1}^{R} (T-j) = n_treats
	# Last block of effects for cohort R has T-R effects.
	# So f_inds[R] should be n_treats - (T-R) + 1.
	if (R > 0) {
		stopifnot(f_inds[R] == n_treats - (T - R) + 1)
	}

	# The original formula was: `1 + (r - 1)*(2*T - r)/2`
	# Let's re-check the paper's formula logic, it implies number of effects for cohort `k` (1-indexed) is `T-k`.
	# sum_{j=1}^{r-1} (T-j) = (r-1)T - (r-1)r/2. So f_inds[r] = 1 + (r-1)T - r(r-1)/2.
	# (r-1)*(2T-r)/2 = (r-1)T - r(r-1)/2. Yes, they are the same.

	# Using the paper's original loop for safety, assuming it's correct with getNumTreats
	f_inds_paper <- integer(R)
	if (R > 0) {
		for (r_val in 1:R) {
			# r_val is 1-indexed cohort number
			f_inds_paper[r_val] <- 1 + (r_val - 1) * (2 * T - r_val) / 2 # This is not sum of T-k, but related to cumulative sum for triangular numbers.
			# This formula is sum_{k=0}^{r-2} (T-1-k) if T-1 is max effects
			# Or sum_{j=0}^{r-1-1} ( (T-1) - j )
			# Let's use the direct sum formulation which is clearer.
		}
		stopifnot(f_inds == f_inds_paper) # Verify my derivation matches paper's formula
	}

	return(f_inds)
}

#' Transform Design Matrix for Fusion Penalties (`transformXintImproved`)
#'
#' @description Applies a series of transformations to the input design matrix
#'   `X_int`. This transformation is a change of coordinates such that applying a
#'   standard bridge regression penalty to the transformed matrix results in the
#'   desired FETWFE fusion penalties on the original parameterization.
#' Transforms the design matrix `X_int` (containing original fixed effects,
#' covariates, interactions, and treatment dummies) into `X_mod`. This
#' transformation is a change of coordinates such that applying a standard
#' Lq penalty to coefficients of `X_mod` is equivalent to applying the
#' desired fusion penalties to coefficients of `X_int`.
#' This corresponds to multiplying `X_int` by `D_N^{-1}` where `theta = D_N beta`.
#' @param X_int Numeric matrix; the initial design matrix containing cohort fixed
#'   effects, time fixed effects, covariates, treatment indicators, and all
#'   their interactions.
#' @param N Integer; the total number of unique units in the panel.
#' @param T Integer; the total number of time periods in the panel.
#' @param R Integer; the total number of treated cohorts (excluding the
#'   never-treated group).
#' @param d Integer; the number of time-invariant covariates.
#' @param num_treats Integer; the total number of unique treatment effect
#'   parameters (i.e., number of cohort-time treatment indicators).
#' @param first_inds Integer vector; (Optional) A vector indicating the index of
#'   the first treatment effect for each of the `R` cohorts within the block of
#'   `num_treats` treatment effect parameters. If NA (default), it's computed
#'   internally using `getFirstInds(R, T)`.
#' @return A numeric matrix `X_mod` of the same dimensions as `X_int`,
#'   representing the transformed design matrix.
#' @details The function sequentially transforms blocks of columns in `X_int`
#'   corresponding to:
#'   \itemize{
#'     \item Cohort fixed effects
#'     \item Time fixed effects
#'     \item Covariate main effects (copied directly if `d > 0`)
#'     \item Interactions between covariates and cohort fixed effects (if `d > 0`)
#'     \item Interactions between covariates and time fixed effects (if `d > 0`)
#'     \item Base treatment effects
#'     \item Interactions between covariates and treatment effects (if `d > 0`)
#'   }
#'   Transformations for fusion penalization are applied using helper functions
#'   like `genBackwardsInvFusionTransformMat()` and `genInvTwoWayFusionTransformMat()`.
#' @keywords internal
#' @noRd
transformXintImproved <- function(
	X_int,
	N,
	T,
	R,
	d,
	num_treats,
	first_inds = NA
) {
	p <- getP(R = R, T = T, d = d, num_treats = num_treats)
	stopifnot(p == ncol(X_int))
	X_mod <- matrix(as.numeric(NA), nrow = N * T, ncol = p)
	stopifnot(nrow(X_int) == N * T)

	# Transform cohort fixed effects
	X_mod[, 1:R] <- X_int[, 1:R] %*% genBackwardsInvFusionTransformMat(R)

	# Transform time fixed effects
	X_mod[, (R + 1):(R + T - 1)] <- X_int[, (R + 1):(R + T - 1)] %*%
		genBackwardsInvFusionTransformMat(T - 1)

	# Copy X (the main covariate block; may be empty when d==0)
	if (d > 0) {
		stopifnot(all(is.na(X_mod[, (R + T - 1 + 1):(R + T - 1 + d)])))
		X_mod[, (R + T - 1 + 1):(R + T - 1 + d)] <- X_int[,
			(R + T - 1 + 1):(R + T - 1 + d)
		]

		stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d)])))
		stopifnot(all(is.na(X_mod[, (R + T - 1 + d + 1):p])))
	}

	# For cohort effects interacted with X: we have d*R columns to deal with.
	# For each individual feature, this will be handled using
	# genTransformedMatFusion.
	if (any(is.na(first_inds))) {
		first_inds <- getFirstInds(R = R, T = T)
	}

	if (d > 0) {
		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and cohort
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on R times
			feat_1 <- R + T - 1 + d + j
			feat_R <- R + T - 1 + d + (R - 1) * d + j
			feat_inds_j <- seq(feat_1, feat_R, by = d)
			stopifnot(length(feat_inds_j) == R)

			stopifnot(all(is.na(X_mod[, feat_inds_j])))

			X_mod[, feat_inds_j] <- X_int[, feat_inds_j] %*%
				genBackwardsInvFusionTransformMat(R)
		}
		stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d + R * d)])))
		stopifnot(all(is.na(X_mod[, (R + T - 1 + d + R * d + 1):p])))

		# Similar for time effects interacted with X

		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and time
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on T - 1 times
			feat_1 <- R + T - 1 + d + R * d + j
			feat_T_minus_1 <- R + T - 1 + d + R * d + (T - 2) * d + j
			feat_inds_j <- seq(feat_1, feat_T_minus_1, by = d)
			stopifnot(length(feat_inds_j) == T - 1)

			stopifnot(all(is.na(X_mod[, feat_inds_j])))
			X_mod[, feat_inds_j] <- X_int[, feat_inds_j] %*%
				genBackwardsInvFusionTransformMat(T - 1)
		}
		stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d + R * d + (T - 1) * d)])))
		stopifnot(all(is.na(X_mod[,
			(R + T - 1 + d + R * d + (T - 1) * d + 1):p
		])))
	}

	# Now base treatment effects. For each cohort, will penalize base term, then
	# fuse remaining terms toward it. Also, for each cohort, will penalize base
	# treatment effect of this cohort to base of previous cohort. New function
	# genTransformedMatTwoWayFusion does this.

	feat_inds <- (R + T - 1 + d + R * d + (T - 1) * d + 1):(R +
		T -
		1 +
		d +
		R * d +
		(T - 1) * d +
		num_treats)

	# Now ready to generate the appropriate transformed matrix
	stopifnot(all(is.na(X_mod[, feat_inds])))

	X_mod[, feat_inds] <- X_int[, feat_inds] %*%
		genInvTwoWayFusionTransformMat(num_treats, first_inds, R)

	# stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d + R*d + (T - 1)*d + num_treats)])))
	# stopifnot(all(is.na(X_mod[, (R + T - 1 + d + R*d + (T - 1)*d + num_treats + 1):p])))

	if (d > 0) {
		# Lastly, penalize interactions between each treatment effect and each feature.
		# Feature-wise, we can do this with genTransformedMatTwoWayFusion, in the same
		# way that we did for previous interactions with X.
		for (j in 1:d) {
			# Recall that we have arranged the last d*num_Feats features in X_int
			# as follows: the first d are the first column of treat_mat_long interacted
			# with all of the columns of X, and so on. So, the columns that interact
			# the jth feature with all of the treatment effects are columns j, j + 1*d,
			# j + 2*d, ..., j + (num_treats - 1)*d.
			inds_j <- seq(j, j + (num_treats - 1) * d, by = d)
			stopifnot(length(inds_j) == num_treats)
			inds_j <- inds_j + R + T - 1 + d + R * d + (T - 1) * d + num_treats

			# Now ready to generate the appropriate transformed matrix
			stopifnot(all(is.na(X_mod[, inds_j])))

			X_mod[, inds_j] <- X_int[, inds_j] %*%
				genInvTwoWayFusionTransformMat(num_treats, first_inds, R)

			stopifnot(all(!is.na(X_mod[, inds_j])))
		}
	}

	stopifnot(all(!is.na(X_mod)))
	stopifnot(ncol(X_mod) == p)
	stopifnot(nrow(X_mod) == N * T)

	return(X_mod)
}

# untransformCoefImproved
#' @title Untransform Estimated Coefficients Back to Original Scale
#' @description Reverses the transformation applied by `transformXintImproved()`
#'   to a vector of estimated coefficients. This converts coefficients estimated
#'   in the transformed (fused) space back to the original parameterization of
#'   the ETWFE model.
#' @param beta_hat_mod Numeric vector; the estimated coefficients from the
#'   penalized regression on the transformed design matrix `X_mod`. Length `p`.
#' @param T Integer; the total number of time periods in the panel.
#' @param R Integer; the total number of treated cohorts.
#' @param p Integer; the total number of parameters (columns in the original
#'   design matrix `X_int`).
#' @param d Integer; the number of time-invariant covariates.
#' @param num_treats Integer; the total number of unique treatment effect
#'   parameters.
#' @param first_inds Integer vector; (Optional) A vector indicating the index of
#'   the first treatment effect for each of the `R` cohorts. If NA (default),
#'   it's computed internally.
#' @return A numeric vector `beta_hat` of length `p`, representing the
#'   coefficients in the original, untransformed feature space.
#' @details This function mirrors `transformXintImproved()` by applying the
#'   inverse of the fusion transformations to the respective blocks of
#'   `beta_hat_mod`. It uses helper functions such as
#'   `genBackwardsInvFusionTransformMat()` and `genInvTwoWayFusionTransformMat()`
#'   (which are their own inverses in this context, acting as the forward
#'   transformation from the sparse basis to the original basis).
#' @keywords internal
#' @noRd
untransformCoefImproved <- function(
	beta_hat_mod,
	T,
	R,
	p,
	d,
	num_treats,
	first_inds = NA
) {
	stopifnot(length(beta_hat_mod) == p)
	beta_hat <- rep(as.numeric(NA), p)

	if (any(is.na(first_inds))) {
		first_inds <- getFirstInds(R = R, T = T)
	}

	# First handle R cohort fixed effects effects
	beta_hat[1:R] <- genBackwardsInvFusionTransformMat(R) %*% beta_hat_mod[1:R]

	stopifnot(all(!is.na(beta_hat[1:R])))
	stopifnot(all(is.na(beta_hat[(R + 1):p])))

	# Next, T - 1 time fixed effects
	beta_hat[(R + 1):(R + T - 1)] <- genBackwardsInvFusionTransformMat(
		T - 1
	) %*%
		beta_hat_mod[(R + 1):(R + T - 1)]

	stopifnot(all(!is.na(beta_hat[1:(R + T - 1)])))
	stopifnot(all(is.na(beta_hat[(R + T):p])))

	# Coefficients for X (if any)
	if (d > 0) {
		beta_hat[(R + T):(R + T - 1 + d)] <- beta_hat_mod[
			(R + T):(R + T - 1 + d)
		]

		stopifnot(all(!is.na(beta_hat[1:(R + T - 1 + d)])))
		stopifnot(all(is.na(beta_hat[(R + T + d):p])))

		# Next, coefficients for cohort effects interacted with X. For each individual
		# feature, this will be handled using untransformVecFusion.
		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and cohort
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on R times
			feat_1 <- R + T - 1 + d + j
			feat_R <- R + T - 1 + d + (R - 1) * d + j
			feat_inds_j <- seq(feat_1, feat_R, by = d)
			stopifnot(length(feat_inds_j) == R)

			stopifnot(all(is.na(beta_hat[feat_inds_j])))

			beta_hat[feat_inds_j] <- genBackwardsInvFusionTransformMat(R) %*%
				beta_hat_mod[feat_inds_j]
			stopifnot(all(!is.na(beta_hat[feat_inds_j])))
		}
		stopifnot(all(!is.na(beta_hat[1:(R + T - 1 + d + R * d)])))
		stopifnot(all(is.na(beta_hat[(R + T - 1 + d + R * d + 1):p])))

		# Similar for time effects interacted with X

		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and time
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on T - 1 times
			feat_1 <- R + T - 1 + d + R * d + j
			feat_T_minus_1 <- R + T - 1 + d + R * d + (T - 2) * d + j
			feat_inds_j <- seq(feat_1, feat_T_minus_1, by = d)
			stopifnot(length(feat_inds_j) == T - 1)

			stopifnot(all(is.na(beta_hat[feat_inds_j])))

			beta_hat[feat_inds_j] <- genBackwardsInvFusionTransformMat(
				T - 1
			) %*%
				beta_hat_mod[feat_inds_j]
			stopifnot(all(!is.na(beta_hat[feat_inds_j])))
		}
		stopifnot(all(
			!is.na(beta_hat[1:(R + T - 1 + d + R * d + (T - 1) * d)])
		))
		stopifnot(all(is.na(beta_hat[
			(R + T - 1 + d + R * d + (T - 1) * d + 1):p
		])))
	}

	# Now base treatment effects.

	feat_inds <- (R + T - 1 + d + R * d + (T - 1) * d + 1):(R +
		T -
		1 +
		d +
		R * d +
		(T - 1) * d +
		num_treats)

	stopifnot(all(is.na(beta_hat[feat_inds])))

	beta_hat[feat_inds] <- genInvTwoWayFusionTransformMat(
		num_treats,
		first_inds,
		R
	) %*%
		beta_hat_mod[feat_inds]

	if (d > 0) {
		stopifnot(all(
			!is.na(beta_hat[
				1:(R + T - 1 + d + R * d + (T - 1) * d + num_treats)
			])
		))
		stopifnot(all(is.na(beta_hat[
			(R + T - 1 + d + R * d + (T - 1) * d + num_treats + 1):p
		])))
		# Lastly, interactions between each treatment effect and each feature.
		# Feature-wise, we can do this with untransformTwoWayFusionCoefs, in the same
		# way that we did for previous interactions with X.
		for (j in 1:d) {
			# Recall that we have arranged the last d*num_Feats features in X_int
			# as follows: the first d are the first column of treat_mat_long interacted
			# with all of the columns of X, and so on. So, the columns that interact
			# the jth feature with all of the treatment effects are columns j, j + 1*d,
			# j + 2*d, ..., j + (num_treats - 1)*d.
			inds_j <- seq(j, j + (num_treats - 1) * d, by = d)
			stopifnot(length(inds_j) == num_treats)
			inds_j <- inds_j + R + T - 1 + d + R * d + (T - 1) * d + num_treats

			# Now ready to untransform the estimated coefficients
			stopifnot(all(is.na(beta_hat[inds_j])))

			beta_hat[inds_j] <- genInvTwoWayFusionTransformMat(
				num_treats,
				first_inds,
				R
			) %*%
				beta_hat_mod[inds_j]

			stopifnot(all(!is.na(beta_hat[inds_j])))
		}
	}

	stopifnot(all(!is.na(beta_hat)))

	return(beta_hat)
}

# sse_bridge
#' @title Calculate Sum of Squared Errors for Bridge Regression
#' @description Computes the sum of squared errors (SSE) for a given set of
#'   intercept and slope coefficients from a bridge regression model, relative
#'   to the observed responses and the (potentially transformed) design matrix.
#'   The result is then averaged by dividing by the total number of observations (N*T).
#' @param eta_hat Numeric scalar; the estimated intercept term.
#' @param beta_hat Numeric vector; the estimated slope coefficients. Its length
#'   should match the number of columns in `X_mod`.
#' @param y Numeric vector; the observed response variable, of length `N*T`.
#' @param X_mod Numeric matrix; the design matrix (possibly transformed, e.g.,
#'   for FETWFE) used in the regression. It has `N*T` rows.
#' @param N Integer; the total number of unique units.
#' @param T Integer; the total number of time periods.
#' @return A numeric scalar representing the mean squared error (MSE), i.e.,
#'   SSE / (N*T).
#' @keywords internal
#' @noRd
sse_bridge <- function(eta_hat, beta_hat, y, X_mod, N, T) {
	stopifnot(length(eta_hat) == 1)
	stopifnot(length(beta_hat) == ncol(X_mod))

	y_hat <- X_mod %*% beta_hat + eta_hat
	stopifnot(length(y_hat) == N * T)
	stopifnot(all(!is.na(y_hat)))

	ret <- sum((y - y_hat)^2) / (N * T)

	stopifnot(!is.na(ret))
	stopifnot(ret >= 0)

	return(ret)
}


# genBackwardsFusionTransformMat
#' @title Generate Backward Fusion Transformation Matrix
#' @description Creates a square transformation matrix `D` of size `n_vars` x
#'   `n_vars`. When pre-multiplied by a coefficient vector `beta`, `D %*% beta`
#'   yields a transformed vector `theta` where `theta_i = beta_i - beta_{i+1}`
#'   for `i < n_vars`, and `theta_{n_vars} = beta_{n_vars}`. This is used to
#'   penalize coefficients towards the *next* coefficient in sequence, and the
#'   last coefficient directly.
#' @param n_vars Integer; the number of variables (coefficients) in the block
#'   to be transformed. This will be the dimension of the output matrix.
#' @return A numeric matrix of dimension `n_vars` x `n_vars`.
#' @details The resulting matrix `D` has 1s on the main diagonal. For each row
#'   `i` (from 1 to `n_vars - 1`), it has a -1 at column `i+1`. All other
#'   elements are 0.
#' @examples
#'   genBackwardsFusionTransformMat(3)
#'   # Output:
#'   #      [,1] [,2] [,3]
#'   # [1,]    1   -1    0
#'   # [2,]    0    1   -1
#'   # [3,]    0    0    1
#' @keywords internal
#' @noRd
genBackwardsFusionTransformMat <- function(n_vars) {
	# Generates D matrix in relation theta = D beta, where D beta is what
	# we want to penalize (for a single set of coefficients where we want to
	# penalize last coefficient directly and penalize remaining coefficints
	# towards the next coefficient)
	D <- matrix(0, n_vars, n_vars)

	for (i in 1:n_vars) {
		for (j in 1:n_vars) {
			if (i == j) {
				D[i, j] <- 1
			}
			if (j == i + 1) {
				D[i, j] <- -1
			}
		}
	}

	return(D)
}

# genBackwardsInvFusionTransformMat
#' @title Generate Inverse of Backward Fusion Transformation Matrix
#' @description Creates the inverse of the transformation matrix `D` generated by
#'   `genBackwardsFusionTransformMat(n_vars)`. If `theta = D %*% beta`, then
#'   `beta = D_inv %*% theta`, where `D_inv` is the matrix returned by this
#'   function. This matrix effectively transforms coefficients from a space where
#'   differences are penalized back to the original coefficient space.
#' @param n_vars Integer; the number of variables (coefficients), determining
#'   the dimension of the output matrix.
#' @return A numeric matrix `D_inv` of dimension `n_vars` x `n_vars`.
#' @details The resulting matrix `D_inv` is an upper triangular matrix with all
#'   elements on and above the main diagonal equal to 1, and all elements below
#'   the main diagonal equal to 0.
#' @examples
#'   genBackwardsInvFusionTransformMat(3)
#'   # Output:
#'   #      [,1] [,2] [,3]
#'   # [1,]    1    1    1
#'   # [2,]    0    1    1
#'   # [3,]    0    0    1
#' @keywords internal
#' @noRd
genBackwardsInvFusionTransformMat <- function(n_vars) {
	# Generates inverse of D matrix in relation theta = D beta, where D beta is
	# what we want to penalize (for a single set of coefficients where we want
	# to penalize last coefficient directly and penalize remaining coefficints
	# towards the next coefficient)
	D_inv <- matrix(0, n_vars, n_vars)

	diag(D_inv) <- 1

	D_inv[upper.tri(D_inv)] <- 1

	stopifnot(nrow(D_inv) == n_vars)
	stopifnot(ncol(D_inv) == n_vars)

	return(D_inv)
}

# genInvFusionTransformMat
#' @title Generate Inverse of Forward Fusion Transformation Matrix
#' @description Creates the inverse of a "forward" fusion transformation matrix.
#'   A forward fusion matrix `D_forward` (not explicitly generated here) would
#'   transform `beta` to `theta` such that `theta_1 = beta_1` and
#'   `theta_i = beta_i - beta_{i-1}` for `i > 1`. This function returns
#'   `D_forward_inv`. If `theta = D_forward %*% beta`, then
#'   `beta = D_forward_inv %*% theta`. This is used when the first coefficient
#'   in a sequence is penalized directly, and subsequent coefficients are
#'   penalized towards the *previous* coefficient.
#' @param n_vars Integer; the number of variables (coefficients), determining
#'   the dimension of the output matrix.
#' @return A numeric matrix `D_forward_inv` of dimension `n_vars` x `n_vars`.
#' @details The resulting matrix is a lower triangular matrix with all elements
#'   on and below the main diagonal equal to 1, and all elements above the main
#'   diagonal equal to 0.
#' @examples
#'   genInvFusionTransformMat(3)
#'   # Output:
#'   #      [,1] [,2] [,3]
#'   # [1,]    1    0    0
#'   # [2,]    1    1    0
#'   # [3,]    1    1    1
#' @keywords internal
#' @noRd
genInvFusionTransformMat <- function(n_vars) {
	# Generates inverse of D matrix in relation theta = D beta, where D beta is
	# what we want to penalize (for a single set of coefficients where we want
	# to penalize first coefficient directly and penalize remaining coefficints
	# towards the previous coefficient)
	D_inv <- matrix(0, n_vars, n_vars)

	diag(D_inv) <- 1

	D_inv[lower.tri(D_inv)] <- 1

	return(D_inv)
}


# estOmegaSqrtInv
#' @title Estimate Noise Variance Components and Omega Matrix
#' @description Estimates the idiosyncratic error variance (`sig_eps_sq`) and
#'   the unit-level random effect variance (`sig_eps_c_sq`) using the method
#'   described by Pesaran (2015, Section 26.5.1) with ridge regression.
#'   This function is called when these variances are not provided by the user.
#' @param y Numeric vector; the observed response variable, length `N*T`.
#' @param X_ints Numeric matrix; the design matrix including all fixed effects,
#'   covariates, treatment dummies, and interactions. `N*T` rows.
#' @param N Integer; the total number of unique units.
#' @param T Integer; the total number of time periods.
#' @param p Integer; the number of columns in `X_ints`.
#' @return A list containing two named elements:
#'   \item{sig_eps_sq}{Estimated variance of the idiosyncratic error term.}
#'   \item{sig_eps_c_sq}{Estimated variance of the unit-level random effect.}
#' @details The function first demeans `y` and `X_ints` within each unit (fixed
#'   effects transformation). Then, it fits a ridge regression (`alpha=0` in
#'   `glmnet::cv.glmnet`) to the demeaned data. Residuals from this model are
#'   used to estimate `sigma_hat_sq` (which corresponds to `sig_eps_sq`).
#'   The estimated unit-specific intercepts (`alpha_hat`) are then used to
#'   estimate `sigma_c_sq_hat` (corresponding to `sig_eps_c_sq`).
#' @references Pesaran, M. H. (2015). Time Series and Panel Data Econometrics.
#'   Oxford University Press.
#' @keywords internal
#' @noRd
estOmegaSqrtInv <- function(y, X_ints, N, T, p) {
	if (N * (T - 1) - p <= 0) {
		stop("Not enough units available to estimate the noise variance.")
	}
	stopifnot(N > 1)

	# Estimate standard deviations
	y_fe <- y
	X_ints_fe <- X_ints
	for (i in 1:N) {
		inds_i <- ((i - 1) * T + 1):(i * T)
		stopifnot(length(inds_i) == T)
		y_fe[inds_i] <- y[inds_i] - mean(y[inds_i])
		X_ints[inds_i, ] <- X_ints[inds_i, ] - colMeans(X_ints[inds_i, ])
	}

	lin_mod_fe <- glmnet::cv.glmnet(x = X_ints_fe, y = y_fe, alpha = 0)

	# Get residuals
	y_hat_fe <- stats::predict(lin_mod_fe, s = "lambda.min", newx = X_ints_fe)

	# Get coefficients
	beta_hat_fe <- stats::coef(lin_mod_fe, s = "lambda.min")[2:(p + 1)]

	resids <- y_fe - y_hat_fe

	tau <- rep(1, T)
	M <- diag(rep(1, T)) - outer(tau, tau) / T

	sigma_hat_sq <- 0
	alpha_hat <- rep(as.numeric(NA), N)

	for (i in 1:N) {
		inds_i <- ((i - 1) * T + 1):(i * T)
		sigma_hat_sq <- sigma_hat_sq +
			as.numeric(resids[inds_i] %*% M %*% resids[inds_i])
		alpha_hat[i] <- mean(y_hat_fe[inds_i]) -
			colMeans(X_ints_fe[inds_i, ] %*% beta_hat_fe)
	}

	stopifnot(all(!is.na(alpha_hat)))

	sigma_hat_sq <- sigma_hat_sq / (N * (T - 1) - p)

	sigma_c_sq_hat <- sum((alpha_hat - mean(alpha_hat))^2) / (N - 1)

	return(list(sig_eps_sq = sigma_hat_sq, sig_eps_c_sq = sigma_c_sq_hat))
}

# getSecondVarTermDataApp
#' @title Calculate Second Variance Term for ATT Standard Error (Data Application)
#' @description Computes the second component of the variance for the Average
#'   Treatment Effect on the Treated (ATT). This component accounts for the
#'   variability due to the estimation of cohort membership probabilities. This
#'   version seems tailored for contexts where cohort probabilities are
#'   estimated from the data sample.
#' @param cohort_probs Numeric vector; estimated probabilities of belonging to
#'   each treated cohort, conditional on being treated. Length `R`.
#' @param psi_mat Numeric matrix; a matrix where each column `r` is the `psi_r`
#'   vector used in calculating the ATT for cohort `r`. Dimensions:
#'   `length(sel_treat_inds_shifted)` x `R`.
#' @param sel_treat_inds_shifted Integer vector; indices of the selected
#'   treatment effects within the `num_treats` block, shifted to start from 1.
#' @param tes Numeric vector; the estimated treatment effects for all
#'   `num_treats` possible cohort-time combinations.
#' @param d_inv_treat_sel Numeric matrix; the relevant block of the inverse
#'   two-way fusion transformation matrix corresponding to selected treatment
#'   effects. Dimensions: `num_treats` (or fewer if selection occurs) x
#'   `length(sel_treat_inds_shifted)`.
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort (P(W=r)). Length `R`.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort within the `num_treats` block.
#' @param theta_hat_treat_sel Numeric vector; estimated coefficients in the
#'   transformed (fused) space, corresponding only to the selected treatment
#'   effects.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param R Integer; total number of treated cohorts.
#' @return A numeric scalar representing the second variance component for the
#'   ATT.
#' @details This function calculates `Sigma_pi_hat`, the covariance matrix of
#'   the cohort assignment indicators, and a Jacobian matrix. These are then
#'   combined with `theta_hat_treat_sel` to compute the variance term as
#'   `T * t(theta_hat_treat_sel) %*% t(jacobian_mat) %*% Sigma_pi_hat %*% jacobian_mat %*%
#'   theta_hat_treat_sel / (N * T)`. The construction of the Jacobian involves averaging parts of
#'   `d_inv_treat_sel` corresponding to different cohorts.
#' @keywords internal
#' @noRd
getSecondVarTermDataApp <- function(
	cohort_probs,
	psi_mat,
	sel_treat_inds_shifted,
	tes,
	d_inv_treat_sel,
	cohort_probs_overall,
	first_inds,
	theta_hat_treat_sel,
	num_treats,
	N,
	T,
	R
) {
	stopifnot(ncol(d_inv_treat_sel) == length(sel_treat_inds_shifted))
	stopifnot(length(theta_hat_treat_sel) == length(sel_treat_inds_shifted))
	Sigma_pi_hat <- -outer(
		cohort_probs_overall[1:(R)],
		cohort_probs_overall[1:(R)]
	)
	diag(Sigma_pi_hat) <- cohort_probs_overall[1:(R)] *
		(1 - cohort_probs_overall[1:(R)])

	stopifnot(nrow(Sigma_pi_hat) == R)
	stopifnot(ncol(Sigma_pi_hat) == R)

	# Jacobian
	jacobian_mat <- matrix(
		as.numeric(NA),
		nrow = R,
		ncol = length(sel_treat_inds_shifted)
	)

	sel_inds <- list()

	for (r in 1:R) {
		first_ind_r <- first_inds[r]

		if (r < R) {
			last_ind_r <- first_inds[r + 1] - 1
		} else {
			last_ind_r <- num_treats
		}
		stopifnot(last_ind_r >= first_ind_r)
		sel_inds[[r]] <- first_ind_r:last_ind_r
		if (r > 1) {
			stopifnot(min(sel_inds[[r]]) > max(sel_inds[[r - 1]]))
			stopifnot(length(sel_inds[[r]]) < length(sel_inds[[r - 1]]))
		}
	}
	stopifnot(all.equal(unlist(sel_inds), 1:num_treats))

	for (r in 1:R) {
		cons_r <- (sum(cohort_probs_overall) -
			cohort_probs_overall[r]) /
			sum(cohort_probs_overall)^2

		if (length(sel_treat_inds_shifted) > 1) {
			jacobian_mat[r, ] <- cons_r *
				colMeans(d_inv_treat_sel[sel_inds[[r]], , drop = FALSE])
		} else {
			jacobian_mat[r, ] <- cons_r *
				mean(d_inv_treat_sel[sel_inds[[r]], , drop = FALSE])
		}
		for (r_double_prime in setdiff(1:R, r)) {
			cons_r_double_prime <- (sum(cohort_probs_overall) -
				cohort_probs_overall[r_double_prime]) /
				sum(cohort_probs_overall)^2

			if (length(sel_treat_inds_shifted) > 1) {
				jacobian_mat[r, ] <- jacobian_mat[r, ] -
					cons_r_double_prime *
						colMeans(d_inv_treat_sel[
							sel_inds[[r_double_prime]],
							,
							drop = FALSE
						])
			} else {
				jacobian_mat[r, ] <- jacobian_mat[r, ] -
					cons_r_double_prime *
						mean(d_inv_treat_sel[
							sel_inds[[r_double_prime]],
							,
							drop = FALSE
						])
			}
		}
	}

	stopifnot(all(!is.na(jacobian_mat)))

	att_var_2 <- T *
		as.numeric(
			t(theta_hat_treat_sel) %*%
				t(jacobian_mat) %*%
				Sigma_pi_hat %*%
				jacobian_mat %*%
				theta_hat_treat_sel
		) /
		(N * T)

	return(att_var_2)
}

# getCohortATTsFinal
#' @title Calculate Cohort-Specific ATTs and Standard Errors
#' @description Computes the Average Treatment Effect on the Treated (ATT) for
#'   each cohort, along with their standard errors and confidence intervals if
#'   requested and feasible.
#' @param X_final Numeric matrix; the final design matrix, potentially
#'   transformed by `Omega_sqrt_inv` and the fusion transformation.
#' @param sel_feat_inds Integer vector; indices of all features selected by the
#'   penalized regression in the transformed space.
#' @param treat_inds Integer vector; indices in the original (untransformed)
#'   coefficient vector that correspond to the base treatment effects.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort within the block of `num_treats` treatment effect parameters.
#' @param sel_treat_inds_shifted Integer vector; indices of selected treatment
#'   effects within the `num_treats` block (shifted to start from 1).
#' @param c_names Character vector; names of the `R` treated cohorts.
#' @param tes Numeric vector; estimated treatment effects in the original
#'   parameterization for all `num_treats` possible cohort-time combinations.
#' @param sig_eps_sq Numeric scalar; variance of the idiosyncratic error term.
#' @param R Integer; total number of treated cohorts.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param fused Logical; if `TRUE`, assumes fusion penalization was used,
#'   affecting how standard errors and related matrices are computed.
#' @param calc_ses Logical; if `TRUE`, attempts to calculate standard errors.
#'   This is typically `TRUE` if `q < 1`.
#' @param p Integer; total number of parameters in the model.
#' @param alpha Numeric scalar; significance level for confidence intervals
#'   (e.g., 0.05 for 95% CIs).
#' @return A list containing:
#'   \item{cohort_te_df}{Dataframe with cohort names, estimated ATTs, SEs, and
#'     confidence interval bounds.}
#'   \item{cohort_tes}{Named numeric vector of estimated ATTs for each cohort.}
#'   \item{cohort_te_ses}{Named numeric vector of standard errors for cohort ATTs.}
#'   \item{psi_mat}{Matrix used in SE calculation for overall ATT.}
#'   \item{gram_inv}{(Potentially NA) Inverse of the Gram matrix for selected
#'     features, used in SE calculation.}
#'   \item{d_inv_treat_sel}{(If `fused=TRUE`) Relevant block of the inverse
#'     fusion matrix for selected treatment effects.}
#'   \item{calc_ses}{Logical, indicating if SEs were actually calculated.}
#' @details The function first computes the Gram matrix inverse (`gram_inv`) if
#'   `calc_ses` is `TRUE`. Then, for each cohort `r`, it calculates the average
#'   of the relevant `tes`. If SEs are calculated, it uses `getPsiRFused` or
#'   `getPsiRUnfused` to get a `psi_r` vector, which is then used with
#'   `gram_inv` to find the standard error for that cohort's ATT.
#' @keywords internal
#' @noRd
getCohortATTsFinal <- function(
	X_final,
	sel_feat_inds,
	treat_inds,
	num_treats,
	first_inds,
	sel_treat_inds_shifted,
	c_names,
	tes,
	sig_eps_sq,
	R,
	N,
	T,
	fused,
	calc_ses,
	p,
	alpha = 0.05
) {
	stopifnot(max(sel_treat_inds_shifted) <= num_treats)
	stopifnot(min(sel_treat_inds_shifted) >= 1)
	stopifnot(length(tes) == num_treats)
	stopifnot(all(!is.na(tes)))

	stopifnot(nrow(X_final) == N * T)
	X_to_pass <- X_final

	# Start by getting Gram matrix needed for standard errors
	if (calc_ses) {
		res <- getGramInv(
			N = N,
			T = T,
			X_final = X_to_pass,
			sel_feat_inds = sel_feat_inds,
			treat_inds = treat_inds,
			num_treats = num_treats,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			calc_ses = calc_ses
		)

		gram_inv <- res$gram_inv
		calc_ses <- res$calc_ses
	} else {
		gram_inv <- NA
	}

	if (fused) {
		# Get the parts of D_inv that have to do with treatment effects
		d_inv_treat <- genInvTwoWayFusionTransformMat(num_treats, first_inds, R)
	}

	# First, each cohort
	cohort_tes <- rep(as.numeric(NA), R)
	cohort_te_ses <- rep(as.numeric(NA), R)

	psi_mat <- matrix(0, length(sel_treat_inds_shifted), R)

	d_inv_treat_sel <- matrix(
		0,
		nrow = 0,
		ncol = length(sel_treat_inds_shifted)
	)

	for (r in 1:R) {
		# Get indices corresponding to rth treatment
		first_ind_r <- first_inds[r]
		if (r < R) {
			last_ind_r <- first_inds[r + 1] - 1
		} else {
			last_ind_r <- num_treats
		}

		stopifnot(last_ind_r >= first_ind_r)
		stopifnot(all(first_ind_r:last_ind_r %in% 1:num_treats))

		cohort_tes[r] <- mean(tes[first_ind_r:last_ind_r])

		if (calc_ses) {
			# Calculate standard errors

			if (fused) {
				res_r <- getPsiRFused(
					first_ind_r,
					last_ind_r,
					sel_treat_inds_shifted,
					d_inv_treat
				)

				psi_r <- res_r$psi_r

				stopifnot(
					nrow(res_r$d_inv_treat_sel) ==
						last_ind_r -
							first_ind_r +
							1
				)
				stopifnot(
					ncol(res_r$d_inv_treat_sel) ==
						length(sel_treat_inds_shifted)
				)

				stopifnot(is.matrix(res_r$d_inv_treat_sel))

				d_inv_treat_sel <- rbind(d_inv_treat_sel, res_r$d_inv_treat_sel)

				if (nrow(d_inv_treat_sel) != last_ind_r) {
					err_mes <- paste(
						"nrow(d_inv_treat_sel) == last_ind_r is not TRUE. ",
						"nrow(d_inv_treat_sel): ",
						nrow(d_inv_treat_sel),
						". num_treats: ",
						num_treats,
						". R: ",
						R,
						". first_inds: ",
						paste(first_inds, collapse = ", "),
						". r: ",
						r,
						". first_ind_r: ",
						first_ind_r,
						". last_ind_r: ",
						last_ind_r,
						". nrow(res_r$d_inv_treat_sel):",
						nrow(res_r$d_inv_treat_sel)
					)
					stop(err_mes)
				}

				rm(res_r)
			} else {
				psi_r <- getPsiRUnfused(
					first_ind_r,
					last_ind_r,
					sel_treat_inds_shifted,
					gram_inv
				)
			}

			stopifnot(length(psi_r) == length(sel_treat_inds_shifted))

			psi_mat[, r] <- psi_r
			# Get standard errors

			cohort_te_ses[r] <- sqrt(
				sig_eps_sq *
					as.numeric(
						t(psi_r) %*%
							gram_inv %*%
							psi_r
					) /
					(N * T)
			)
		}
	}

	if (fused & calc_ses) {
		if (nrow(d_inv_treat_sel) != num_treats) {
			err_mes <- paste(
				"nrow(d_inv_treat_sel) == num_treats is not TRUE. ",
				"nrow(d_inv_treat_sel): ",
				nrow(d_inv_treat_sel),
				". num_treats: ",
				num_treats,
				". R: ",
				R,
				". first_inds: ",
				paste(first_inds, collapse = ", "),
				"."
			)
			stop(err_mes)
		}
	}

	stopifnot(length(c_names) == R)
	stopifnot(length(cohort_tes) == R)

	if (calc_ses & all(!is.na(gram_inv))) {
		stopifnot(length(cohort_te_ses) == R)

		cohort_te_df <- data.frame(
			c_names,
			cohort_tes,
			cohort_te_ses,
			cohort_tes - stats::qnorm(1 - alpha / 2) * cohort_te_ses,
			cohort_tes + stats::qnorm(1 - alpha / 2) * cohort_te_ses
		)

		names(cohort_te_ses) <- c_names
		names(cohort_tes) <- c_names
	} else {
		cohort_te_df <- data.frame(
			c_names,
			cohort_tes,
			rep(NA, R),
			rep(NA, R),
			rep(NA, R)
		)
	}

	colnames(cohort_te_df) <- c(
		"Cohort",
		"Estimated TE",
		"SE",
		"ConfIntLow",
		"ConfIntHigh"
	)

	if (fused) {
		stopifnot(is.matrix(d_inv_treat_sel))
		ret <- list(
			cohort_te_df = cohort_te_df,
			cohort_tes = cohort_tes,
			cohort_te_ses = cohort_te_ses,
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			d_inv_treat_sel = d_inv_treat_sel,
			calc_ses = calc_ses
		)
	} else {
		ret <- list(
			cohort_te_df = cohort_te_df,
			cohort_tes = cohort_tes,
			cohort_te_ses = cohort_te_ses,
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			calc_ses = calc_ses
		)
	}
	return(ret)
}

# getPsiRUnfused
#' @title Calculate Psi Vector for Cohort ATT (Unfused Case)
#' @description Computes the `psi_r` vector for a specific cohort `r` when no
#'   fusion penalization is applied to the treatment effects (or when calculating
#'   SEs as if it were an OLS on selected variables). This vector is used in
#'   standard error calculations for the cohort's Average Treatment Effect on
#'   the Treated (ATT).
#' @param first_ind_r Integer; the index of the first treatment effect for
#'   cohort `r` within the `num_treats` block of treatment effects.
#' @param last_ind_r Integer; the index of the last treatment effect for
#'   cohort `r` within the `num_treats` block.
#' @param sel_treat_inds_shifted Integer vector; indices of all selected
#'   treatment effects within the `num_treats` block, shifted to start from 1.
#' @param gram_inv Numeric matrix; the inverse of the Gram matrix for the
#'   selected treatment effect features. (Note: This parameter seems unused in the
#'   current function body provided, `psi_r` is constructed based on `sel_treat_inds_shifted` only).
#' @return A numeric vector `psi_r` of length equal to
#'   `length(sel_treat_inds_shifted)`. It contains weights (typically 1/k for
#'   k selected effects in cohort r, 0 otherwise) to average the selected
#'   treatment effect coefficients for cohort `r`.
#' @details The function identifies which of the `sel_treat_inds_shifted` fall
#'   within the range `[first_ind_r, last_ind_r]`. For these identified
#'   indices in `psi_r`, it assigns a value of 1. `psi_r` is then normalized by
#'   dividing by its sum, effectively creating an averaging vector for the
#'   selected treatment effects belonging to cohort `r`. If no treatment effects
#'   for cohort `r` were selected, `psi_r` will be a zero vector.
#' @keywords internal
#' @noRd
getPsiRUnfused <- function(
	first_ind_r,
	last_ind_r,
	sel_treat_inds_shifted,
	gram_inv
) {
	which_inds_ir <- sel_treat_inds_shifted %in% (first_ind_r:last_ind_r)

	psi_r <- rep(0, length(sel_treat_inds_shifted))

	if (sum(which_inds_ir) > 0) {
		inds_r <- which(which_inds_ir)

		stopifnot(is.integer(inds_r) | is.numeric(inds_r))
		stopifnot(identical(inds_r, as.integer(round(inds_r))))
		stopifnot(length(inds_r) >= 1)
		stopifnot(length(inds_r) == length(unique(inds_r)))
		stopifnot(length(inds_r) <= length(sel_treat_inds_shifted))
		stopifnot(all(inds_r %in% 1:length(sel_treat_inds_shifted)))

		stopifnot(max(inds_r) <= nrow(gram_inv))
		stopifnot(max(inds_r) <= ncol(gram_inv))
		stopifnot(min(inds_r) >= 0)

		psi_r[inds_r] <- 1

		stopifnot(sum(psi_r) > 0)

		psi_r <- psi_r / sum(psi_r)
	}

	return(psi_r)
}

# genInvTwoWayFusionTransformMat
#' @title Generate Inverse Two-Way Fusion Transformation Matrix
#' @description Creates a square transformation matrix `D_inv` of size `n_vars` x
#'   `n_vars` used for the two-way fusion penalization of treatment effects.
#'   This matrix transforms coefficients from a sparse "fused" basis back to
#'   the original parameterization of treatment effects. The penalization
#'   structure involves fusing effects within cohorts over time and fusing the
#'   initial effect of each cohort to the initial effect of the previous cohort.
#' @param n_vars Integer; the total number of treatment effect parameters,
#'   determining the dimension of the output matrix.
#' @param first_inds Integer vector; a vector of length `R` where `first_inds[j]`
#'   is the index (1-based) of the first treatment effect parameter corresponding
#'   to the j-th treated cohort within the block of `n_vars` parameters.
#' @param R Integer; the total number of treated cohorts.
#' @return A numeric matrix `D_inv` of dimension `n_vars` x `n_vars`.
#' @details The matrix is constructed such that:
#'   \itemize{
#'     \item It has 1s on the main diagonal.
#'     \item For each cohort `j` (from 1 to `R-1`):
#'       \itemize{
#'         \item `D[first_inds[j]:n_vars, first_inds[j]] <- 1`
#'         \item `D[first_inds[j+1]:n_vars, first_inds[j+1]] <- 1`
#'       }
#'     \item Within each cohort `j`, for effects from `first_inds[j]+1` up to
#'       the one before `first_inds[j+1]`, say effect `k`,
#'       `D[k, (first_inds[j]+1):k] <- 1`. This creates a lower triangular
#'       structure for within-cohort fusion.
#'   }
#'   This structure implies that an original coefficient is a sum of certain
#'   transformed (theta) coefficients.
#' @keywords internal
#' @noRd
genInvTwoWayFusionTransformMat <- function(n_vars, first_inds, R) {
	stopifnot(length(n_vars) == 1)
	stopifnot(length(first_inds) == R)
	D <- matrix(0, n_vars, n_vars)

	diag(D) <- 1

	if (R < 2) {
		stop(
			"Only one treated cohort detected in data. Currently fetwfe only supports data sets with at least two treated cohorts."
		)
	}

	for (j in 1:(R - 1)) {
		index_j <- first_inds[j]
		next_index <- first_inds[j + 1]

		D[index_j:n_vars, index_j] <- 1
		D[next_index:n_vars, next_index] <- 1

		if (index_j + 1 <= next_index - 1) {
			for (k in (index_j + 1):(next_index - 1)) {
				D[k, (index_j + 1):k] <- 1
			}
		}
	}

	return(D)
}

# getBetaBIC
#' @title Select Optimal Coefficients using BIC from gBridge Fit
#' @description From a `gBridge` fit object (which contains solutions for a
#'   path of lambda penalties), this function selects the optimal set of
#'   coefficients based on the Bayesian Information Criterion (BIC). It also
#'   returns the chosen lambda index and the size of the selected model.
#'   Coefficients are returned on their original scale.
#' @param fit A `gBridge` fit object, typically the output from `grpreg::gBridge()`.
#' @param N Integer; the total number of unique units.
#' @param T Integer; the total number of time periods.
#' @param p Integer; the total number of predictor variables (excluding intercept)
#'   in the model matrix `X_mod`.
#' @param X_mod Numeric matrix; the design matrix (potentially transformed for
#'   FETWFE, and **not** yet GLStransformed or scaled/centered by `my_scale`) that was used to generate `y`.
#'   It's used here to calculate SSE on the original scale of `y`.
#' @param y Numeric vector; the original response variable (before GLS transform and centering)
#'   used to fit the model. Length `N*T`.
#' @param scale_center Numeric vector; the centering values used to scale `X_mod`
#'   before fitting `gBridge`. Length `p`.
#' @param scale_scale Numeric vector; the scaling values used to scale `X_mod`
#'   before fitting `gBridge`. Length `p`.
#' @return A list containing:
#'   \item{theta_hat}{Numeric vector of length `p+1`. The selected coefficients
#'     (including intercept at `theta_hat[1]`) on their original data scale.}
#'   \item{lambda_star_ind}{Integer; the index of the lambda value in `fit$lambda`
#'     that resulted in the best BIC.}
#'   \item{lambda_star_model_size}{Integer; the number of non-zero coefficients
#'     (excluding intercept) in the selected model.}
#' @details The function iterates through each lambda in `fit$lambda`. For each:
#'   1. It extracts the intercept (`eta_s`) and slopes (`beta_s`) on the scaled data.
#'   2. It converts these coefficients back to the original data scale using
#'      `scale_center` and `scale_scale`.
#'   3. It calculates the Sum of Squared Errors (SSE) using `sse_bridge()` with
#'      the original-scale coefficients, original `y`, and `X_mod`.
#'   4. It computes the BIC value: `N*T*log(SSE/(N*T)) + s*log(N*T)`, where `s`
#'      is the number of non-zero coefficients (including intercept).
#'   The set of coefficients corresponding to the minimum BIC is chosen. If multiple
#'   lambdas yield the same minimum BIC, the one resulting in the smallest model
#'   size (fewest non-zero coefficients) is selected.
#'   The final returned `theta_hat` also has its slopes and intercept adjusted back to the original scale.
#' @keywords internal
#' @noRd
getBetaBIC <- function(fit, N, T, p, X_mod, y, scale_center, scale_scale) {
	stopifnot(length(y) == N * T)
	n_lambda <- ncol(fit$beta)
	BICs <- rep(as.numeric(NA), n_lambda)
	model_sizes <- rep(as.integer(NA), n_lambda)

	stopifnot(nrow(fit$beta) == p + 1)

	for (k in 1:n_lambda) {
		## --- extract coefficients on the scaled data -------------
		eta_s <- fit$beta[1, k] # intercept (scaled space)
		beta_s <- fit$beta[2:(p + 1), k] # slopes    (scaled space)

		## --- convert to original scale ---------------------------
		beta_hat_k <- beta_s / scale_scale
		eta_hat_k <- eta_s - sum(scale_center * beta_hat_k)

		# Residual sum of squares
		mse_hat <- sse_bridge(
			eta_hat_k,
			beta_hat_k,
			y = y,
			X_mod = X_mod,
			N = N,
			T = T
		)
		# Number of fitted coefficients
		s <- sum(fit$beta[, k] != 0)
		model_sizes[k] <- s

		stopifnot(is.na(BICs[k]))
		BICs[k] <- N * T * log(mse_hat) + s * log(N * T)
	}

	lambda_star_ind <- which(BICs == min(BICs))
	if (length(lambda_star_ind) == 1) {
		lambda_star_final_ind <- lambda_star_ind
		theta_hat <- fit$beta[, lambda_star_final_ind]
	} else {
		# Choose smallest model size among models with equal BIC
		model_sizes_star <- model_sizes[lambda_star_ind]
		min_model_size_ind <- which(model_sizes_star == min(model_sizes_star))
		lambda_star_final_ind <- lambda_star_ind[min_model_size_ind][1]
		stopifnot(length(lambda_star_final_ind) == 1)
		theta_hat <- fit$beta[, lambda_star_final_ind]
	}
	stopifnot(length(lambda_star_final_ind) == 1)
	stopifnot(length(theta_hat) == p + 1)
	stopifnot(all(!is.na(theta_hat)))

	#
	# Rescale coefficients back to the original scale.
	# The coefficient vector theta_hat is of length (p+1), with the first entry
	# as the intercept.
	# For predictors: original coefficient = beta_scaled / scale_j.
	# The intercept is adjusted as: intercept_original =
	# intercept_scaled - sum(center_j * (beta_scaled/scale_j)).
	#
	adjusted_theta_hat <- theta_hat
	if (length(scale_scale) != p) {
		stop("Length of scale_scale does not match number of predictors (p).")
	}
	for (j in 2:(p + 1)) {
		adjusted_theta_hat[j] <- theta_hat[j] / scale_scale[j - 1]
	}
	adjusted_theta_hat[1] <- theta_hat[1] -
		sum(scale_center * (theta_hat[2:(p + 1)] / scale_scale))

	return(list(
		theta_hat = adjusted_theta_hat,
		lambda_star_ind = lambda_star_final_ind,
		lambda_star_model_size = model_sizes[lambda_star_final_ind]
	))
}

# getGramInv
#' @title Compute Inverse of Gram Matrix for Selected Features
#' @description Calculates the inverse of the Gram matrix formed by the selected
#'   features from the (potentially transformed) design matrix `X_final`. This
#'   is a key component in calculating standard errors for the estimated
#'   coefficients.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param X_final Numeric matrix; the final design matrix (e.g., after GLS and
#'   fusion transformations). Dimensions: `N*T` x `p_model`, where `p_model` is
#'   the total number of columns in this matrix.
#' @param sel_feat_inds Integer vector; indices of the features selected by the
#'   penalized regression, corresponding to columns in `X_final`.
#' @param treat_inds Integer vector; original indices (before selection) of the
#'   base treatment effect parameters. Used to subset `sel_feat_inds` to get
#'   only selected *treatment* features for the final Gram matrix.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param sel_treat_inds_shifted Integer vector; indices of the selected
#'   treatment effects within the `num_treats` block, shifted to start from 1.
#'   Used to check dimensions.
#' @param calc_ses Logical; if `FALSE`, the function may return `NA` for
#'   `gram_inv`.
#' @return A list containing:
#'   \item{gram_inv}{The inverse of the Gram matrix corresponding to the
#'     *selected treatment effect features*. Returns `NA` if `calc_ses` is
#'     `FALSE` or if the Gram matrix is found to be singular.}
#'   \item{calc_ses}{Logical, potentially modified to `FALSE` if the Gram
#'     matrix is singular.}
#' @details
#'   1. Subsets `X_final` to include only columns specified by `sel_feat_inds`.
#'   2. Centers these selected columns.
#'   3. Computes the Gram matrix: `(1/(N*T)) * t(X_sel_centered) %*% X_sel_centered`.
#'   4. Checks if the minimum eigenvalue of the Gram matrix is too small (close to zero).
#'      If so, it issues a warning and sets `calc_ses` to `FALSE`, returning `NA`
#'      for `gram_inv`.
#'   5. Otherwise, it computes the inverse of the Gram matrix.
#'   6. It then subsets this inverse Gram matrix to retain only the rows/columns
#'      that correspond to the *selected treatment effects* (identified via
#'      `sel_feat_inds` and `treat_inds`).
#' @keywords internal
#' @noRd
getGramInv <- function(
	N,
	T,
	X_final,
	sel_feat_inds,
	treat_inds,
	num_treats,
	sel_treat_inds_shifted,
	calc_ses
) {
	stopifnot(nrow(X_final) == N * T)
	X_sel <- X_final[, sel_feat_inds, drop = FALSE]
	X_sel_centered <- scale(X_sel, center = TRUE, scale = FALSE)

	gram <- 1 / (N * T) * (t(X_sel_centered) %*% X_sel_centered)

	stopifnot(nrow(gram) == length(sel_feat_inds))
	stopifnot(ncol(gram) == length(sel_feat_inds))

	min_gram_eigen <- min(
		eigen(gram, symmetric = TRUE, only.values = TRUE)$values
	)

	if (min_gram_eigen < 10^(-12)) {
		warning(
			"Gram matrix corresponding to selected features is not invertible. Assumptions needed for inference are not satisfied. Standard errors will not be calculated."
		)
		return(list(gram_inv = NA, calc_ses = FALSE))
	}

	gram_inv <- solve(gram)

	# Get only the parts of gram_inv that have to do with treatment effects
	sel_treat_inds <- sel_feat_inds %in% treat_inds

	stopifnot(is.logical(sel_treat_inds))
	stopifnot(sum(sel_treat_inds) <= length(sel_feat_inds))
	stopifnot(length(sel_feat_inds) == length(sel_feat_inds))

	gram_inv <- gram_inv[sel_treat_inds, sel_treat_inds]

	stopifnot(nrow(gram_inv) <= num_treats)
	stopifnot(nrow(gram_inv) <= length(sel_feat_inds))
	stopifnot(nrow(gram_inv) == ncol(gram_inv))
	stopifnot(nrow(gram_inv) == length(sel_treat_inds_shifted))

	return(list(gram_inv = gram_inv, calc_ses = calc_ses))
}

# getPsiRFused
#' @title Calculate Psi Vector and D-inverse Block for Cohort ATT (Fused Case)
#' @description Computes the `psi_r` vector and the relevant block of the
#'   inverse fusion transformation matrix (`d_inv_treat_sel`) for a specific
#'   cohort `r`. These are used in standard error calculations for the cohort's
#'   Average Treatment Effect on the Treated (ATT) when fusion penalization
#'   has been applied.
#' @param first_ind_r Integer; the index of the first treatment effect parameter
#'   for cohort `r` within the original `num_treats` block (1-based).
#' @param last_ind_r Integer; the index of the last treatment effect parameter
#'   for cohort `r` within the original `num_treats` block (1-based).
#' @param sel_treat_inds_shifted Integer vector; indices (1-based) of the
#'   treatment effects that were selected by the model, relative to the start
#'   of the `num_treats` block. E.g., if original indices 5, 7 were selected
#'   from a block starting at index 1, this would be c(5, 7).
#' @param d_inv_treat Numeric matrix; the full inverse two-way fusion
#'   transformation matrix for all `num_treats` treatment effects. Dimensions:
#'   `num_treats` x `num_treats`.
#' @return A list containing:
#'   \item{psi_r}{Numeric vector. It's the column means of the sub-matrix of
#'     `d_inv_treat` corresponding to rows `first_ind_r:last_ind_r` and columns
#'     specified by `sel_treat_inds_shifted`. If `first_ind_r == last_ind_r`,
#'     it's just that specific row of `d_inv_treat` (subsetted by selected columns).}
#'   \item{d_inv_treat_sel}{Numeric matrix. The sub-matrix of `d_inv_treat`
#'     with rows `first_ind_r:last_ind_r` and columns corresponding to
#'     `sel_treat_inds_shifted`.}
#' @details `psi_r` effectively averages the rows of `d_inv_treat` (that correspond
#'   to cohort `r`'s treatment effects) for the columns that were actually
#'   selected by the model. `d_inv_treat_sel` is this specific block of the
#'   `d_inv_treat` matrix.
#' @keywords internal
#' @noRd
getPsiRFused <- function(
	first_ind_r,
	last_ind_r,
	sel_treat_inds_shifted,
	d_inv_treat
) {
	stopifnot(length(sel_treat_inds_shifted) >= 0)
	stopifnot(last_ind_r >= first_ind_r)
	# Get psi vector: the part of D inverse that we need to look at is the
	# block corresponding to the treatment effect estimates, which is the
	# num_treats x num_treats matrix yielded by
	# genInvTwoWayFusionTransformMat(num_treats, first_inds).

	# Correct rows of matrix

	if (last_ind_r > first_ind_r) {
		if (length(sel_treat_inds_shifted) > 1) {
			psi_r <- colMeans(d_inv_treat[
				first_ind_r:last_ind_r,
				sel_treat_inds_shifted
			])
		} else {
			psi_r <- mean(d_inv_treat[
				first_ind_r:last_ind_r,
				sel_treat_inds_shifted
			])
		}
		if (length(sel_treat_inds_shifted) == 1) {
			# Need to coerce this object to be a matrix with one column so it
			# works smoothly with rbind() later
			d_inv_treat_sel <- matrix(
				d_inv_treat[first_ind_r:last_ind_r, sel_treat_inds_shifted],
				ncol = 1
			)
		} else {
			d_inv_treat_sel <- d_inv_treat[
				first_ind_r:last_ind_r,
				sel_treat_inds_shifted
			]
		}
	} else {
		psi_r <- d_inv_treat[first_ind_r:last_ind_r, sel_treat_inds_shifted]
		# Since first_ind_r and last_ind_r are the same, need to coerce this
		# object to be a matrix with one row so that it works smoothly with
		# rbind() later
		d_inv_treat_sel <- matrix(
			d_inv_treat[first_ind_r:last_ind_r, sel_treat_inds_shifted],
			nrow = 1
		)
	}

	stopifnot(is.matrix(d_inv_treat_sel))

	return(list(psi_r = psi_r, d_inv_treat_sel = d_inv_treat_sel))
}

# getTeResults2
#' @title Calculate Overall ATT and its Standard Error (Version 2)
#' @description Computes the overall Average Treatment Effect on the Treated
#'   (ATT) by taking a weighted average of cohort-specific ATTs. It also
#'   calculates the standard error for this overall ATT, potentially considering
#'   variability from estimated cohort probabilities.
#' @param sig_eps_sq Numeric scalar; variance of the idiosyncratic error term.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param R Integer; total number of treated cohorts.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param cohort_tes Numeric vector; estimated ATTs for each of the `R` cohorts.
#' @param cohort_probs Numeric vector; weights for each cohort, typically
#'   estimated probabilities of belonging to cohort `r` conditional on being
#'   treated. Length `R`. Sums to 1.
#' @param psi_mat Numeric matrix; matrix where column `r` is `psi_r` (from
#'   `getCohortATTsFinal`). Dimensions: `length(sel_treat_inds_shifted)` x `R`.
#' @param gram_inv Numeric matrix; inverse of the Gram matrix for selected
#'   treatment effect features.
#' @param sel_treat_inds_shifted Integer vector; indices of selected treatment
#'   effects within the `num_treats` block (shifted to start from 1).
#' @param tes Numeric vector; all `num_treats` estimated treatment effects
#'   (original parameterization).
#' @param d_inv_treat_sel Numeric matrix; block of the inverse fusion matrix for
#'   selected treatment effects.
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort P(W=r). Length `R`.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort.
#' @param theta_hat_treat_sel Numeric vector; estimated coefficients in
#'   transformed (fused) space for selected treatment effects.
#' @param calc_ses Logical; if `TRUE`, calculate standard errors.
#' @param indep_probs Logical; if `TRUE`, assumes `cohort_probs` (and
#'   `cohort_probs_overall`) were estimated from an independent sample, leading
#'   to a different SE formula (sum of variances) compared to when they are
#'   estimated from the same sample (conservative SE including a covariance term).
#' @return A list containing:
#'   \item{att_hat}{Numeric scalar; the estimated overall ATT.}
#'   \item{att_te_se}{Numeric scalar; the standard error for `att_hat`. NA if
#'     `calc_ses` is `FALSE`.}
#'   \item{att_te_se_no_prob}{Numeric scalar; standard error for `att_hat`
#'     ignoring variability from estimating cohort probabilities (i.e., only
#'     `att_var_1`). NA if `calc_ses` is `FALSE`.}
#' @details The overall ATT (`att_hat`) is `cohort_tes %*% cohort_probs`.
#'   If `calc_ses` is `TRUE`:
#'   - `att_var_1` (variance from `theta_hat` estimation) is computed using
#'     `psi_att = psi_mat %*% cohort_probs` and `gram_inv`.
#'   - `att_var_2` (variance from cohort probability estimation) is computed by
#'     calling `getSecondVarTermDataApp`.
#'   - `att_te_se` is `sqrt(att_var_1 + att_var_2)` if `indep_probs` is `TRUE`,
#'     otherwise it's a conservative SE: `sqrt(att_var_1 + att_var_2 + 2*sqrt(att_var_1 * att_var_2))`.
#'   - `att_te_se_no_prob` is `sqrt(att_var_1)`.
#' @keywords internal
#' @noRd
getTeResults2 <- function(
	# model,
	sig_eps_sq,
	N,
	T,
	R,
	num_treats,
	cohort_tes,
	cohort_probs,
	psi_mat,
	gram_inv,
	sel_treat_inds_shifted,
	tes,
	d_inv_treat_sel,
	cohort_probs_overall,
	first_inds,
	theta_hat_treat_sel,
	calc_ses,
	indep_probs = FALSE
) {
	att_hat <- as.numeric(cohort_tes %*% cohort_probs)

	if (calc_ses) {
		# Get ATT standard error
		# first variance term: convergence of theta
		psi_att <- psi_mat %*% cohort_probs

		att_var_1 <- sig_eps_sq *
			as.numeric(t(psi_att) %*% gram_inv %*% psi_att) /
			(N * T)

		# Second variance term: convergence of cohort membership probabilities
		att_var_2 <- getSecondVarTermDataApp(
			cohort_probs = cohort_probs,
			psi_mat = psi_mat,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			tes = tes,
			d_inv_treat_sel = d_inv_treat_sel,
			cohort_probs_overall = cohort_probs_overall,
			first_inds = first_inds,
			theta_hat_treat_sel = theta_hat_treat_sel,
			num_treats = num_treats,
			N = N,
			T = T,
			R = R
		)

		if (indep_probs) {
			att_te_se <- sqrt(att_var_1 + att_var_2)
		} else {
			att_te_se <- sqrt(
				att_var_1 +
					att_var_2 +
					2 *
						sqrt(
							att_var_1 * att_var_2
						)
			)
		}

		att_te_se_no_prob <- sqrt(att_var_1)
	} else {
		att_te_se <- NA
		att_te_se_no_prob <- NA
	}

	return(list(
		att_hat = att_hat,
		att_te_se = att_te_se,
		att_te_se_no_prob = att_te_se_no_prob
	))
}

# getNumTreats
#' @title Calculate the Total Number of Base Treatment Effect Parameters
#' @description Computes the total number of unique treatment dummy variables
#'   (and thus base treatment effect parameters `tau_rt`) in a staggered
#'   adoption setting.
#' @param R Integer; the number of treated cohorts. Treatment is assumed to
#'   start in periods 2 to `R+1`.
#' @param T Integer; the total number of time periods.
#' @return An integer representing the total number of treatment effect
#'   parameters (`num_treats`).
#' @details The formula used is `num_treats = T * R - (R * (R + 1)) / 2`.
#'   This corresponds to summing the number of post-treatment periods for each
#'   cohort:
#'   Cohort 1 (starts period 2): T-1 effects
#'   Cohort 2 (starts period 3): T-2 effects
#'   ...
#'   Cohort R (starts period R+1): T-R effects
#'   Summing these gives `R*T - (1 + 2 + ... + R) = R*T - R(R+1)/2`.
#' @keywords internal
#' @noRd
getNumTreats <- function(R, T) {
	return(T * R - (R * (R + 1)) / 2)
}

# getTreatInds
#' @title Get Indices of Base Treatment Effect Parameters
#' @description Determines the column indices in the full design matrix `X_ints` (or
#'   corresponding coefficient vector `beta`) that correspond to the base
#'   treatment effect parameters (`tau_rt`).
#' @param R Integer; the number of treated cohorts.
#' @param T Integer; the number of time periods.
#' @param d Integer; the number of time-invariant covariates.
#' @param num_treats Integer; the total number of base treatment effect
#'   parameters, typically calculated by `getNumTreats(R, T)`.
#' @return An integer vector containing the indices for the `num_treats` base
#'   treatment effect parameters.
#' @details The full design matrix `X_ints` is structured with several blocks of
#'   variables:
#'   1. Cohort fixed effects (`R` columns)
#'   2. Time fixed effects (`T-1` columns)
#'   3. Covariate main effects (`d` columns, if `d > 0`)
#'   4. Cohort-Covariate interactions (`d*R` columns, if `d > 0`)
#'   5. Time-Covariate interactions (`d*(T-1)` columns, if `d > 0`)
#'   The base treatment effects form the next block. This function calculates the
#'   total number of columns in the preceding blocks (`base_cols`) and returns
#'   `base_cols + 1` through `base_cols + num_treats`.
#' @keywords internal
#' @noRd
getTreatInds <- function(R, T, d, num_treats) {
	base_cols <- if (d > 0) {
		R + (T - 1) + d + d * R + d * (T - 1)
	} else {
		R + (T - 1)
	}

	treat_inds <- seq(from = base_cols + 1, length.out = num_treats)

	stopifnot(length(treat_inds) == num_treats)
	if (d > 0) {
		stopifnot(
			max(treat_inds) == R + T - 1 + d + R * d + (T - 1) * d + num_treats
		)
	} else {
		stopifnot(max(treat_inds) == R + T - 1 + num_treats)
	}

	return(treat_inds)
}

# getP
#' @title Calculate Total Number of Parameters in the Full Design Matrix
#' @description Computes `p`, the total number of columns (parameters) in the
#'   full design matrix `X_ints` used by the FETWFE model, including all fixed
#'   effects, covariates, treatment dummies, and all their interactions.
#' @param R Integer; the number of treated cohorts.
#' @param T Integer; the total number of time periods.
#' @param d Integer; the number of time-invariant covariates.
#' @param num_treats Integer; the total number of base treatment effect
#'   parameters (e.g., from `getNumTreats(R,T)`).
#' @return An integer `p` representing the total number of parameters.
#' @details The total number of parameters `p` is the sum of:
#'   - Cohort fixed effects: `R`
#'   - Time fixed effects: `T-1`
#'   - Covariate main effects: `d`
#'   - Cohort-Covariate interactions: `d*R`
#'   - Time-Covariate interactions: `d*(T-1)`
#'   - Base treatment effects: `num_treats`
#'   - Treatment-Covariate interactions: `num_treats*d`
#'   The formula is `p = R + (T-1) + d + d*R + d*(T-1) + num_treats + num_treats*d`.
#'   If `d=0`, terms involving `d` become zero.
#' @keywords internal
#' @noRd
getP <- function(R, T, d, num_treats) {
	return(R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d)
}

# checkFetwfeInputs
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
checkFetwfeInputs <- function(
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

	if (any(!is.na(lambda.max))) {
		stopifnot(is.numeric(lambda.max) | is.integer(lambda.max))
		stopifnot(length(lambda.max) == 1)
		stopifnot(lambda.max > 0)
	}

	if (any(!is.na(lambda.min))) {
		stopifnot(is.numeric(lambda.min) | is.integer(lambda.min))
		stopifnot(length(lambda.min) == 1)
		stopifnot(lambda.min >= 0)
		if (any(!is.na(lambda.max))) {
			stopifnot(lambda.max > lambda.min)
		}
	}

	stopifnot(is.numeric(q) | is.integer(q))
	stopifnot(length(q) == 1)
	stopifnot(q > 0)
	stopifnot(q <= 2)

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

#' Generate the full \eqn{D^{-1}} transformation matrix.
#'
#' @param first_inds A vector of indices corresponding to the first treatment effect for each
#' treated cohort.
#' @param T Total number of time periods.
#' @param R Number of treated cohorts.
#' @param d Number of covariates (time–invariant).
#' @param num_treats Total number of base treatment effect parameters.
#'
#' @return A matrix of dimension p x p where p = R + (T - 1) + d + d*R + d*(T - 1) + num_treats +
#' d*num_treats.
#' @noRd
genFullInvFusionTransformMat <- function(first_inds, T, R, d, num_treats) {
	# Load required package for block diagonal concatenation.
	if (!requireNamespace("Matrix", quietly = TRUE)) {
		stop("The 'Matrix' package is required but not installed.")
	}

	# Block 1: Cohort fixed effects block, size R x R.
	block1 <- genBackwardsInvFusionTransformMat(R)

	# Block 2: Time fixed effects block, size (T - 1) x (T - 1).
	block2 <- genBackwardsInvFusionTransformMat(T - 1)

	# Block 3: Covariate main effects, identity of dimension d.
	block3 <- if (d > 0) diag(d) else NULL

	# Block 4: Cohort-X interactions: I_d \otimes genBackwardsInvFusionTransformMat(R)
	block4 <- if (d > 0)
		kronecker(diag(d), genBackwardsInvFusionTransformMat(R)) else NULL

	# Block 5: Time-X interactions: I_d \otimes genBackwardsInvFusionTransformMat(T - 1)
	block5 <- if (d > 0)
		kronecker(diag(d), genBackwardsInvFusionTransformMat(T - 1)) else NULL

	# Block 6: Base treatment effects: genInvTwoWayFusionTransformMat(num_treats, first_inds, R)
	block6 <- genInvTwoWayFusionTransformMat(num_treats, first_inds, R)

	# Block 7: Treatment-X interactions: I_d \otimes genInvTwoWayFusionTransformMat(num_treats, first_inds, R)
	block7 <- if (d > 0)
		kronecker(
			diag(d),
			genInvTwoWayFusionTransformMat(num_treats, first_inds, R)
		) else NULL

	# Combine blocks into a block-diagonal matrix.
	# Use Matrix::bdiag which returns a sparse matrix; convert to dense if needed.
	blocks <- list(block1, block2)
	if (!is.null(block3)) blocks <- c(blocks, list(block3))
	if (!is.null(block4)) blocks <- c(blocks, list(block4))
	if (!is.null(block5)) blocks <- c(blocks, list(block5))
	blocks <- c(blocks, list(block6))
	if (!is.null(block7)) blocks <- c(blocks, list(block7))

	full_D_inv <- as.matrix(Matrix::bdiag(blocks))

	p <- getP(R = R, T = T, d = d, num_treats = num_treats)

	stopifnot(nrow(full_D_inv) == p)
	stopifnot(ncol(full_D_inv) == p)

	return(full_D_inv)
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

# my_scale
#' @title Custom Scaling Function Handling Zero-Variance Columns
#' @description Centers and scales the columns of a numeric matrix. This function
#'   is similar to `base::scale()` but explicitly handles columns with zero
#'   variance by setting their scale factor to 1, thus avoiding division by zero
#'   and NaN results.
#' @param x Numeric matrix; the matrix whose columns are to be scaled.
#' @return A numeric matrix of the same dimensions as `x`, with columns centered
#'   and scaled. Attributes "scaled:center" and "scaled:scale" are attached,
#'   containing the means and standard deviations (or 1 for zero-variance
#'   columns) used for the transformation.
#' @details
#'   1. Computes column means (`ctr`) and column standard deviations (`sds`) of `x`.
#'   2. Identifies columns where `sds == 0`.
#'   3. For these zero-variance columns, the scaling factor in `sds2` is set to 1.
#'      For other columns, `sds2` is the same as `sds`. The centering values (`ctr2`)
#'      are the original column means.
#'   4. Columns of `x` are first centered using `ctr2` and then scaled (divided)
#'      by `sds2`.
#' @examples
#'\dontrun{
#'   mat <- matrix(c(1, 1, 1, 2, 3, 4), ncol = 2)
#'   # mat is:
#'   #      [,1] [,2]
#'   # [1,]    1    2
#'   # [2,]    1    3
#'   # [3,]    1    4
#'   # Column 1 has zero variance.
#'   scaled_mat <- my_scale(mat)
#'   print(scaled_mat)
#'   # attr(scaled_mat, "scaled:center") # Should be c(1, 3)
#'   # attr(scaled_mat, "scaled:scale")  # Should be c(1, 1)
#'   # Expected output for scaled_mat:
#'   #      [,1] [,2]
#'   # [1,]    0   -1
#'   # [2,]    0    0
#'   # [3,]    0    1
#'}
#' @keywords internal
#' @noRd
my_scale <- function(x) {
	# Compute column means and standard deviations
	ctr <- colMeans(x)
	sds <- apply(x, 2, sd)

	# Identify zero-variance columns
	zero_sd <- (sds == 0)

	# For zero-variance columns, set scale=1 to avoid dividing by 0
	ctr2 <- ctr
	sds2 <- sds
	sds2[zero_sd] <- 1

	# Center and scale
	scaled <- sweep(x, 2, ctr2, FUN = "-")
	scaled <- sweep(scaled, 2, sds2, FUN = "/")

	# Attach attributes so behavior mimics base::scale()
	attr(scaled, "scaled:center") <- ctr2
	attr(scaled, "scaled:scale") <- sds2

	return(scaled)
}
