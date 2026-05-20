# checkEtwfeInputs
#' @title Check Inputs for the main `etwfe` / `betwfe` / `twfeCovs` functions
#' @description Validates the inputs provided to the OLS-style entry points
#'   (`etwfe`, `betwfe`, `twfeCovs`), ensuring they meet type, dimension,
#'   and content requirements. Stops execution with a multi-line error
#'   message listing every malformed input on a single call (instead of
#'   one round-trip per arg). Delegates the per-arg checks to
#'   `.collect_etwfe_input_violations()`.
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
#' @param verbose Logical; if TRUE, print progress. Default `FALSE`.
#' @param alpha Numeric; significance level for confidence intervals. Default `0.05`.
#' @param add_ridge Logical; if TRUE, add small ridge penalty. Default `FALSE`.
#' @return A list with two elements:
#'   - `pdata`: the (possibly tibble-coerced) panel data frame.
#'   - `indep_count_data_available`: logical; `TRUE` if a valid
#'     `indep_counts` argument was provided.
#' @details All per-argument checks live in
#'   `.collect_etwfe_input_violations()`; this thin wrapper renders the
#'   collected violations into a multi-line `stop()` message. The
#'   user-facing wording was upgraded from the legacy
#'   `is.character(time_var) is not TRUE` style to named, actionable
#'   lines naming the arg, the expected shape, and (when feasible) the
#'   received shape (GitHub #84). The message header is
#'   `"Invalid inputs:"` — deliberately caller-agnostic so error text
#'   stays byte-identical across the four entry points (locked by the
#'   PR #103 snapshot tests).
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
	if (length(res$violations) > 0) {
		stop(.format_input_violations(res$violations), call. = FALSE)
	}
	list(
		pdata = res$pdata,
		indep_count_data_available = res$indep_count_data_available
	)
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
#' @param verbose Logical; if `TRUE`, prints progress messages. Default is `FALSE`.
#' @param alpha Numeric; significance level for confidence intervals (e.g., 0.05 for
#'   95% CIs). Default is 0.05.
#' @param add_ridge (Optional) Logical; if `TRUE`, adds a small L2 penalty to
#'   the untransformed coefficients to stabilize estimation. Default is `FALSE`.
#'
#' @details
#' The function executes the following main steps:
#' \enumerate{
#'   \item **Input Checks:** Validates the provided parameters via
#'     `check_etwfe_core_inputs()`.
#'   \item **Design preparation + GLS weighting:** Calls
#'     `prep_for_etwfe_regression(X_mod = X_ints, is_fetwfe = FALSE)`. This
#'     estimates the variance components via REML (`estOmegaSqrtInv()`) when
#'     `sig_eps_sq` or `sig_eps_c_sq` are `NA`, then GLS-weights the design and
#'     response by `sqrt(sig_eps_sq) * Omega_sqrt_inv` via a Kronecker product.
#'     No fusion transformation is applied (`X_mod = X_ints`); the OLS path
#'     fits on the original coefficient parameterization. Also computes cohort
#'     membership probabilities from `in_sample_counts` and (if provided)
#'     `indep_counts`.
#'   \item **Optional ridge penalty:** If `add_ridge = TRUE`, `X_final_scaled`
#'     and `y_final` are augmented to add a small L2 penalty on the
#'     coefficients. The fitted coefficients are then multiplied by
#'     `(1 + lambda_ridge)` per the elastic-net adjustment.
#'   \item **OLS regression:** Fits `lm(y ~ . + 0, df)` on the GLS-weighted
#'     design. No bridge penalty, no `lambda` grid, no BIC selection step.
#'     The coefficient vector is in the original (untransformed) parameter
#'     space and is returned as `beta_hat`.
#'   \item **Cohort-specific treatment effects:** Calls
#'     `getCohortATTsFinal()` (with `fused = FALSE`, `include_selected =
#'     FALSE`) to compute per-cohort point estimates, standard errors (when
#'     the Gram matrix is invertible), and confidence intervals.
#'   \item **Overall ATT:** Calls `getTeResultsOLS()` for the overall ATT
#'     point estimate and SE, using both in-sample probabilities and
#'     independent probabilities when `indep_counts` is provided.
#' }
#' The standard errors for CATTs are asymptotically exact. For ATT, if
#' `indep_counts` are provided, the SE is asymptotically exact; otherwise, it
#' is asymptotically conservative.
#'
#' @return A list containing detailed estimation results:
#'   \item{in_sample_att_hat}{Estimated overall ATT using in-sample cohort probabilities.}
#'   \item{in_sample_att_se}{Standard error for `in_sample_att_hat`.}
#'   \item{in_sample_att_se_no_prob}{SE for `in_sample_att_hat` ignoring variability from estimating cohort probabilities.}
#'   \item{indep_att_hat}{Estimated overall ATT using `indep_counts` cohort probabilities (NA if `indep_counts` not provided).}
#'   \item{indep_att_se}{Standard error for `indep_att_hat` (NA if not applicable).}
#'   \item{catt_hats}{A named vector of estimated CATTs for each cohort.}
#'   \item{catt_ses}{A named vector of SEs for `catt_hats` (NA when the Gram matrix is not invertible).}
#'   \item{catt_df}{A data.frame summarizing CATTs, SEs, and confidence intervals.}
#'   \item{beta_hat}{The vector of estimated coefficients in the *original*
#'     space (no bridge fusion transformation; etwfe is pure OLS).}
#'   \item{treat_inds}{Indices in `beta_hat` corresponding to base treatment effects.}
#'   \item{treat_int_inds}{Indices in `beta_hat` corresponding to treatment-covariate interactions.}
#'   \item{cohort_probs}{Estimated cohort probabilities conditional on being treated, from `in_sample_counts`.}
#'   \item{indep_cohort_probs}{Estimated cohort probabilities from `indep_counts` (NA if not provided).}
#'   \item{sig_eps_sq}{The (possibly estimated) variance of observation-level noise.}
#'   \item{sig_eps_c_sq}{The (possibly estimated) variance of unit-level random effects.}
#'   \item{X_ints}{The original input design matrix from `prepXints`.}
#'   \item{y}{The original input centered response vector from `prepXints`.}
#'   \item{X_final}{The design matrix after GLS weighting (no fusion
#'     transformation for `etwfe_core`).}
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
	add_ridge = FALSE,
	se_type = "default"
) {
	se_type <- match.arg(se_type, c("default", "cluster"))
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

	res <- prep_for_etwfe_regression(
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
	# OLS regression and fitted-coefficient extraction
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
	ti <- .compute_treat_inds(
		R = R,
		T = T,
		d = d,
		num_treats = num_treats,
		p = p
	)
	treat_inds <- ti$treat_inds
	treat_int_inds <- ti$treat_int_inds

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
	# Calculate cohort-specific treatment effects and standard
	# errors
	#
	#

	res <- getCohortATTsFinal(
		X_final = X_final, # This is X_mod * GLS_transform_matrix
		sel_feat_inds = NULL, # OLS path: no penalty selection occurred
		treat_inds = treat_inds, # Global indices for treatment effects
		num_treats = num_treats,
		first_inds = first_inds,
		sel_treat_inds_shifted = seq_len(num_treats),
		c_names = c_names,
		tes = tes, # Treatment effect estimates (beta_hat_slopes[treat_inds])
		sig_eps_sq = sig_eps_sq,
		R = R,
		N = N,
		T = T,
		fused = FALSE,
		calc_ses = TRUE,
		include_selected = FALSE, # etwfe has no bridge selection
		alpha = alpha,
		se_type = se_type,
		y_final = y_final
	)

	cohort_te_df <- res$cohort_te_df
	cohort_tes <- res$cohort_tes
	cohort_te_ses <- res$cohort_te_ses
	psi_mat <- res$psi_mat
	gram_inv <- res$gram_inv
	calc_ses <- res$calc_ses
	sandwich_full <- res$sandwich_full
	treat_block_mask <- res$treat_block_mask

	rm(res)

	stopifnot(nrow(psi_mat) == num_treats)
	stopifnot(ncol(psi_mat) == R)

	#
	#
	# Calculate overall average treatment effect on treated units
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
		indep_probs = FALSE,
		se_type = se_type,
		sandwich_full = sandwich_full,
		treat_block_mask = treat_block_mask
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
			indep_probs = TRUE,
			se_type = se_type,
			sandwich_full = sandwich_full,
			treat_block_mask = treat_block_mask
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
		cohort_probs_overall = cohort_probs_overall,
		indep_cohort_probs_overall = indep_cohort_probs_overall,
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
