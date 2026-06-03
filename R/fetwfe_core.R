# checkFetwfeInputs
#' @title Check Inputs for the main `fetwfe` function
#' @description Validates the inputs provided to the main `fetwfe` function
#'   (and, via the shared validator wiring in `.run_estimator_input_prep()`,
#'   `betwfe()`). Stops execution with a multi-line error message listing
#'   every malformed input on a single call. Delegates the etwfe-shared
#'   per-arg checks to `.collect_etwfe_input_violations()` and appends
#'   the fetwfe/betwfe-specific bridge-penalty (`lambda.max`,
#'   `lambda.min`, `q`) checks before stopping.
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
#' @param q Numeric; Lq penalty exponent for `gBridge`. Default `0.5`.
#' @param verbose Logical; if TRUE, print progress. Default `FALSE`.
#' @param alpha Numeric; significance level for confidence intervals. Default `0.05`.
#' @param add_ridge Logical; if TRUE, add small ridge penalty. Default `FALSE`.
#' @param lambda_selection Character scalar; either `"cv"` (10-fold
#'   cross-validation, the v1.13.0+ default) or `"bic"` (BIC over the
#'   `grpreg`-generated lambda path). Default `"cv"`.
#' @param cv_folds Integer; number of CV folds when `lambda_selection
#'   = "cv"`. Default `10L`.
#' @param cv_seed Integer or NULL; the seed to pass to `set.seed()`
#'   before calling `grpreg::cv.grpreg()` so CV fold assignment is
#'   reproducible. `NULL` defaults to `as.integer(N * T)` (clipped to
#'   `.Machine$integer.max` with a warning on overflow). Default
#'   `NULL`.
#' @return A list with two elements:
#'   - `pdata`: the (possibly tibble-coerced) panel data frame.
#'   - `indep_count_data_available`: logical; `TRUE` if a valid
#'     `indep_counts` argument was provided.
#' @details The error-message UX (header `"Invalid inputs:"` followed by
#'   one bullet per malformed arg) is shared with `checkEtwfeInputs()`
#'   so that, when the only failing arg is etwfe-shared, the four entry
#'   points (`fetwfe`, `etwfe`, `betwfe`, `twfeCovs`) produce
#'   byte-identical messages (a property the PR #103 snapshot tests
#'   assert). The fetwfe-specific bridge-penalty violations append
#'   AFTER the etwfe-shared ones, so when only `q`/`lambda.*` are
#'   malformed, `fetwfe()` and `betwfe()` agree with each other
#'   verbatim.
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
	q = 0.5,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE,
	lambda_selection = "cv",
	cv_folds = 10L,
	cv_seed = NULL
) {
	res <- .collect_etwfe_input_violations(
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
	violations <- res$violations
	pdata <- res$pdata
	indep_count_data_available <- res$indep_count_data_available

	# --- lambda.max ------------------------------------------------------
	if (any(!is.na(lambda.max))) {
		if (!(is.numeric(lambda.max) || is.integer(lambda.max))) {
			violations <- c(
				violations,
				sprintf(
					"lambda.max must be numeric or integer; got %s",
					paste(class(lambda.max), collapse = "/")
				)
			)
		} else if (length(lambda.max) != 1) {
			violations <- c(
				violations,
				sprintf(
					"lambda.max must be a single value; got length %d",
					length(lambda.max)
				)
			)
		} else if (!(lambda.max > 0)) {
			violations <- c(
				violations,
				sprintf(
					"lambda.max must be > 0; got %s",
					format(lambda.max)
				)
			)
		}
	}

	# --- lambda.min ------------------------------------------------------
	if (any(!is.na(lambda.min))) {
		if (!(is.numeric(lambda.min) || is.integer(lambda.min))) {
			violations <- c(
				violations,
				sprintf(
					"lambda.min must be numeric or integer; got %s",
					paste(class(lambda.min), collapse = "/")
				)
			)
		} else if (length(lambda.min) != 1) {
			violations <- c(
				violations,
				sprintf(
					"lambda.min must be a single value; got length %d",
					length(lambda.min)
				)
			)
		} else if (lambda.min < 0) {
			violations <- c(
				violations,
				sprintf(
					"lambda.min must be >= 0; got %s",
					format(lambda.min)
				)
			)
		} else if (
			# Both well-formed; check ordering only when lambda.max also
			# survived its own checks (i.e., is a valid positive scalar).
			any(!is.na(lambda.max)) &&
				is.numeric(lambda.max) &&
				length(lambda.max) == 1 &&
				lambda.max > 0 &&
				!(lambda.max > lambda.min)
		) {
			violations <- c(
				violations,
				sprintf(
					"lambda.max must be > lambda.min; got lambda.max = %s, lambda.min = %s",
					format(lambda.max),
					format(lambda.min)
				)
			)
		}
	}

	# --- q ---------------------------------------------------------------
	if (!(is.numeric(q) || is.integer(q)) || length(q) != 1) {
		violations <- c(
			violations,
			sprintf(
				"q must be a single numeric; got %s (length %d)",
				paste(class(q), collapse = "/"),
				length(q)
			)
		)
	} else if (!(q > 0 && q <= 2)) {
		violations <- c(
			violations,
			sprintf(
				"q must be in (0, 2]; got %s",
				format(q)
			)
		)
	}

	# --- lambda_selection ------------------------------------------------
	if (
		!is.character(lambda_selection) ||
			length(lambda_selection) != 1 ||
			!(lambda_selection %in% c("cv", "bic"))
	) {
		violations <- c(
			violations,
			sprintf(
				"lambda_selection must be a single character, either \"cv\" or \"bic\"; got %s",
				format(lambda_selection)
			)
		)
	}

	# --- cv_folds --------------------------------------------------------
	if (
		!(is.numeric(cv_folds) || is.integer(cv_folds)) ||
			length(cv_folds) != 1
	) {
		violations <- c(
			violations,
			sprintf(
				"cv_folds must be a single integer; got %s (length %d)",
				paste(class(cv_folds), collapse = "/"),
				length(cv_folds)
			)
		)
	} else if (
		!is.finite(cv_folds) || cv_folds < 2 || cv_folds != round(cv_folds)
	) {
		violations <- c(
			violations,
			sprintf(
				"cv_folds must be an integer >= 2; got %s",
				format(cv_folds)
			)
		)
	}

	# --- cv_seed ---------------------------------------------------------
	if (!is.null(cv_seed)) {
		if (
			!(is.numeric(cv_seed) || is.integer(cv_seed)) ||
				length(cv_seed) != 1
		) {
			violations <- c(
				violations,
				sprintf(
					"cv_seed must be NULL or a single integer; got %s (length %d)",
					paste(class(cv_seed), collapse = "/"),
					length(cv_seed)
				)
			)
		} else if (!is.finite(cv_seed) || cv_seed != round(cv_seed)) {
			violations <- c(
				violations,
				sprintf(
					"cv_seed must be an integer (or NULL); got %s",
					format(cv_seed)
				)
			)
		}
	}

	if (length(violations) > 0) {
		stop(.format_input_violations(violations), call. = FALSE)
	}
	list(
		pdata = pdata,
		indep_count_data_available = indep_count_data_available
	)
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
#'   appear in the untreated cohort plus each of the other `G` cohorts, derived
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
#' @param se_type Character; the standard-error type, one of `"default"`
#'   (tight Gaussian variance under (Psi-IF), Theorem (c$'$)),
#'   `"conservative"` (Cauchy-Schwarz upper bound from Theorem (c) for
#'   non-(Psi-IF) propensity estimators), or `"cluster"` (experimental
#'   unit-clustered Liang-Zeger sandwich SE on the bridge-selected
#'   support). See the exported wrapper `fetwfe()` for details. Default
#'   is `"default"`.
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
#'   \item{catt_df}{A data frame (with S3 class `c("catt_df", "data.frame")`) summarizing CATTs (`cohort`, `estimate`, `se`, `ci_low`, `ci_high`, `p_value`) and a `selected` logical flag (`TRUE` when the bridge penalty left the cohort's CATT nonzero). The `catt_df` S3 class makes `[[` / `$` / `[` access on the pre-1.11.0 Title-Case column names `stop()` with a migration message pointing to the new name.}
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
#'   \item{N, T, G, d, p}{Dimensions used in estimation.}
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
	add_ridge = FALSE,
	se_type = "default",
	lambda_selection = "cv",
	cv_folds = 10L,
	cv_seed = NULL
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)
	lambda_selection <- match.arg(lambda_selection, c("cv", "bic"))
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

	G <- ret$G
	c_names <- ret$c_names
	indep_count_data_available <- ret$indep_count_data_available

	rm(ret)

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
		X_int = X_ints,
		N = N,
		T = T,
		G = G,
		d = d,
		num_treats = num_treats,
		first_inds = first_inds
	)

	res <- prep_for_etwfe_regression(
		verbose = verbose,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		y = y,
		X_ints = X_ints,
		X_mod = X_mod,
		N = N,
		T = T,
		G = G,
		d = d,
		p = p,
		num_treats = num_treats,
		add_ridge = add_ridge,
		first_inds = first_inds,
		in_sample_counts = in_sample_counts,
		indep_count_data_available = indep_count_data_available,
		indep_counts = indep_counts,
		is_fetwfe = TRUE
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

	# Dispatch on `lambda_selection` via the shared CV/BIC helper. Both
	# paths produce the same return contract; FETWFE-specific input is
	# the transformed `X_mod` handed to `getBetaBIC()` for SSE.
	bridge_sel <- .dispatch_bridge_selection(
		lambda_selection = lambda_selection,
		X_final_scaled = X_final_scaled,
		y_final = y_final,
		q = q,
		lambda.max = lambda.max,
		lambda.min = lambda.min,
		nlambda = nlambda,
		cv_folds = cv_folds,
		cv_seed = cv_seed,
		N = N,
		T = T,
		p = p,
		X_mod_bic = X_mod,
		y_bic = y,
		scale_center = scale_center,
		scale_scale = scale_scale,
		verbose = verbose
	)
	theta_hat <- bridge_sel$theta_hat
	lambda_star_ind <- bridge_sel$lambda_star_ind
	lambda_star_model_size <- bridge_sel$lambda_star_model_size
	fit <- bridge_sel$fit
	lambda.max <- bridge_sel$lambda.max
	lambda.min <- bridge_sel$lambda.min
	lambda.max_model_size <- bridge_sel$lambda.max_model_size
	lambda.min_model_size <- bridge_sel$lambda.min_model_size
	cv_seed_used <- bridge_sel$cv_seed_used

	lambda_star <- fit$lambda[lambda_star_ind]

	# c_names <- names(in_sample_counts)[2:(G + 1)] # Moved definition up
	stopifnot(length(c_names) == G)

	# Indices corresponding to base treatment effects
	ti <- .compute_treat_inds(
		G = G,
		T = T,
		d = d,
		num_treats = num_treats,
		p = p
	)
	treat_inds <- ti$treat_inds
	treat_int_inds <- ti$treat_int_inds

	# Handle edge case where no features are selected (model_size includes intercept)
	if (lambda_star_model_size <= 1 && all(theta_hat[2:(p + 1)] == 0)) {
		# Only the intercept might be non-zero. Delegate to the shared
		# helper (`.build_selected_out_result()` in `R/core_funcs.R`) that
		# also serves the no-treatment branch below and the two BETWFE
		# early-exits. FETWFE blocks pass `theta_hat` (with intercept)
		# and set `include_theta = TRUE` to insert it between `catt_df`
		# and `beta_hat`.
		return(.build_selected_out_result(
			message_text = "No features selected (or only intercept); all treatment effects estimated to be 0.",
			verbose = verbose,
			G = G,
			c_names = c_names,
			q = q,
			beta_hat = rep(0, p), # Slopes are all zero
			treat_inds = treat_inds,
			treat_int_inds = treat_int_inds,
			cohort_probs = cohort_probs,
			indep_cohort_probs = indep_cohort_probs,
			cohort_probs_overall = cohort_probs_overall,
			indep_cohort_probs_overall = indep_cohort_probs_overall,
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
			d = d,
			p = p,
			theta_hat = theta_hat, # Includes intercept
			include_theta = TRUE,
			cv_seed_used = cv_seed_used
		))
	}

	# estimated coefficients (slopes in transformed space)
	theta_hat_slopes <- theta_hat[2:(p + 1)]

	# Indices of selected features in transformed feature space (among slopes)
	sel_feat_inds <- which(theta_hat_slopes != 0)

	sel_treat_inds <- sel_feat_inds[sel_feat_inds %in% treat_inds]

	stopifnot(length(sel_treat_inds) == length(unique(sel_treat_inds)))
	stopifnot(length(sel_treat_inds) <= length(sel_feat_inds))
	stopifnot(length(sel_treat_inds) <= length(treat_inds))
	stopifnot(is.integer(sel_treat_inds) | is.numeric(sel_treat_inds))

	# 1. Get theta_hat_slopes[treat_inds] -> these are the transformed treatment coefficients
	# 2. Find which of these are non-zero: `which(theta_hat_slopes[treat_inds] != 0)` -> these are
	# indices *within* the treat_inds block. Let's call this `sel_treat_inds_relative_to_block`.
	# `sel_treat_inds` itself contains the global indices of selected treatment features.
	# So `theta_hat_slopes[sel_treat_inds]` are the non-zero transformed treatment coefs.

	theta_hat_treat_block_transformed = theta_hat_slopes[treat_inds]
	sel_treat_inds_shifted <- which(theta_hat_treat_block_transformed != 0) # these are 1 to num_treats

	stopifnot(all(sel_treat_inds_shifted >= 1))
	stopifnot(all(sel_treat_inds_shifted <= num_treats))

	# Handle edge case where no treatment features selected
	if (length(sel_treat_inds_shifted) == 0) {
		# Untransform `theta_hat_slopes` back to original (beta) basis for
		# the helper's contract — `untransformCoefImproved()` is the
		# FETWFE-specific reparameterization (theta -> beta) and stays at
		# the caller (the helper is parameterization-agnostic).
		beta_hat_early_exit <- untransformCoefImproved(
			beta_hat_mod = theta_hat_slopes, # Pass slopes only
			first_inds = first_inds,
			T = T,
			G = G,
			p = p,
			d = d,
			num_treats = num_treats
		)
		# Apply ridge adjustment locally before early-exit return; doesn't
		# affect later code paths (they re-compute `beta_hat` separately
		# and run the non-early-exit `add_ridge` scaling at line ~1083).
		# Plan D3.
		if (add_ridge) {
			lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)
			beta_hat_early_exit <- beta_hat_early_exit * (1 + lambda_ridge)
		}

		return(.build_selected_out_result(
			message_text = "No treatment features selected; all treatment effects estimated to be 0.",
			verbose = verbose,
			G = G,
			c_names = c_names,
			q = q,
			beta_hat = beta_hat_early_exit, # Untransformed slopes
			treat_inds = treat_inds,
			treat_int_inds = treat_int_inds,
			cohort_probs = cohort_probs,
			indep_cohort_probs = indep_cohort_probs,
			cohort_probs_overall = cohort_probs_overall,
			indep_cohort_probs_overall = indep_cohort_probs_overall,
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
			d = d,
			p = p,
			theta_hat = theta_hat, # Full theta_hat with intercept
			include_theta = TRUE,
			cv_seed_used = cv_seed_used
		))
	}

	#
	#
	# Step 5: transform estimated coefficients back to original feature
	# space
	#
	#

	beta_hat <- untransformCoefImproved(
		beta_hat_mod = theta_hat_slopes, # Pass slopes only
		first_inds = first_inds,
		T = T,
		G = G,
		p = p,
		d = d,
		num_treats = num_treats
	)

	# If using ridge regularization, multiply the "naive" estimated coefficients
	# by 1 + lambda_ridge, similar to suggestion in original elastic net paper.
	if (add_ridge) {
		lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)
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

	stopifnot(length(first_inds) == G)
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
		G = G,
		N = N,
		T = T,
		fused = TRUE,
		calc_ses = q < 1,
		include_selected = TRUE,
		alpha = alpha,
		se_type = se_type,
		y_final = y_final
	)

	cohort_te_df <- res$cohort_te_df
	cohort_tes <- res$cohort_tes
	cohort_te_ses <- res$cohort_te_ses
	psi_mat <- res$psi_mat
	gram_inv <- res$gram_inv
	d_inv_treat_sel <- res$d_inv_treat_sel
	calc_ses <- res$calc_ses
	sandwich_full <- res$sandwich_full
	treat_block_mask <- res$treat_block_mask

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
		G = G,
		num_treats = num_treats,
		cohort_tes = cohort_tes, # CATTs (point estimates)
		cohort_probs = cohort_probs, # In-sample pi_g | treated
		psi_mat = psi_mat,
		gram_inv = gram_inv,
		sel_treat_inds_shifted = sel_treat_inds_shifted,
		d_inv_treat_sel = d_inv_treat_sel,
		cohort_probs_overall = cohort_probs_overall, # In-sample pi_g (unconditional on treated)
		first_inds = first_inds,
		theta_hat_treat_sel = theta_hat_treat_sel_for_att, # Selected non-zero transformed treat coefs
		calc_ses = calc_ses,
		indep_probs = FALSE,
		se_type = se_type,
		sandwich_full = sandwich_full,
		treat_block_mask = treat_block_mask
	)

	in_sample_att_hat <- in_sample_te_results$att_hat
	in_sample_att_se <- in_sample_te_results$att_te_se
	in_sample_att_se_no_prob <- in_sample_te_results$att_te_se_no_prob
	in_sample_att_var_1 <- in_sample_te_results$att_var_1
	in_sample_att_var_2 <- in_sample_te_results$att_var_2

	if ((q < 1) & calc_ses) {
		stopifnot(!is.na(in_sample_att_se))
	}

	if (indep_count_data_available) {
		indep_te_results <- getTeResults2(
			sig_eps_sq = sig_eps_sq,
			N = N,
			T = T,
			G = G,
			num_treats = num_treats,
			cohort_tes = cohort_tes,
			cohort_probs = indep_cohort_probs, # indep pi_g | treated
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			d_inv_treat_sel = d_inv_treat_sel,
			cohort_probs_overall = indep_cohort_probs_overall, # indep pi_g (unconditional)
			first_inds = first_inds,
			theta_hat_treat_sel = theta_hat_treat_sel_for_att,
			calc_ses = calc_ses,
			indep_probs = TRUE,
			se_type = se_type,
			sandwich_full = sandwich_full,
			treat_block_mask = treat_block_mask
		)
		indep_att_hat <- indep_te_results$att_hat
		indep_att_se <- indep_te_results$att_te_se
		indep_att_var_1 <- indep_te_results$att_var_1
		indep_att_var_2 <- indep_te_results$att_var_2
	} else {
		indep_att_hat <- NA
		indep_att_se <- NA
		indep_att_var_1 <- NA
		indep_att_var_2 <- NA
	}

	return(list(
		in_sample_att_hat = in_sample_att_hat,
		in_sample_att_se = in_sample_att_se,
		in_sample_att_se_no_prob = in_sample_att_se_no_prob,
		in_sample_att_var_1 = in_sample_att_var_1,
		in_sample_att_var_2 = in_sample_att_var_2,
		indep_att_hat = indep_att_hat,
		indep_att_se = indep_att_se,
		indep_att_var_1 = indep_att_var_1,
		indep_att_var_2 = indep_att_var_2,
		catt_hats = cohort_tes, # Already named if applicable from getCohortATTsFinal
		catt_ses = cohort_te_ses, # Already named if applicable
		catt_df = cohort_te_df,
		theta_hat = theta_hat, # Full theta_hat (with intercept)
		beta_hat = beta_hat, # Untransformed slopes
		treat_inds = treat_inds,
		treat_int_inds = treat_int_inds,
		cohort_probs = cohort_probs,
		indep_cohort_probs = indep_cohort_probs,
		cohort_probs_overall = cohort_probs_overall,
		indep_cohort_probs_overall = indep_cohort_probs_overall,
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
		G = G,
		d = d,
		p = p,
		calc_ses = calc_ses,
		# v1.13.0 (#164): CV-path provenance. cv_seed_used is the integer
		# seed actually fed to set.seed() before cv.grpreg (NA_integer_
		# under BIC).
		cv_seed_used = cv_seed_used
	))
}
