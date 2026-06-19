# Result assembly for the ETWFE / bridge fitting core. Split out of
# R/core_funcs.R for issue #186 (no behavior change). Assembles the final
# selected-coefficient output object returned by the core.

#' @title Build the "no features selected" early-exit return list shared by
#'   `fetwfe_core()` and `betwfe_core()`
#'
#' @description
#' Internal helper that consolidates four near-duplicate early-exit blocks
#' that previously lived inline in `R/betwfe_core.R` and `R/fetwfe_core.R`
#' (BETWFE intercept-only, BETWFE no-treatment-selected, FETWFE
#' intercept-only, FETWFE no-treatment-selected). All four blocks construct
#' the same 33-/34-field return list with all-zero CATTs / ATTs and the same
#' `data.frame` of CATT placeholders; they differ only in the `beta_hat`
#' value, whether `theta_hat` is in the return, and the verbose-message
#' text. Caller computes the final `beta_hat` (including any
#' `add_ridge` adjustment and, for FETWFE, the
#' `untransformCoefImproved()` transformation back from `theta`) and
#' passes it in; helper builds the return list with the contractually-fixed
#' field ordering. Resolves drift surface that produced 9 stale `q < 1`
#' references in v1.9.5 (#56).
#'
#' @param message_text Character scalar; emitted via `message()` when
#'   `verbose` is `TRUE`. Per-block string (see callers).
#' @param verbose Logical; whether to emit `message_text`.
#' @param G,N,T,d,p Integers; problem dimensions.
#' @param c_names Character vector of length `G`; cohort labels for the
#'   `catt_df_to_ret` rows.
#' @param q Numeric; bridge regression exponent. With `gls`, drives the
#'   SE-validity gate `ses_valid <- (q < 1) && gls` that shapes `ret_se` /
#'   `ret_var` / `calc_ses`.
#' @param gls Logical; whether the fit used GLS whitening. A valid (zero)
#'   degenerate SE needs `(q < 1) && gls` (mirrors the normal path's `calc_ses`);
#'   a `gls = FALSE` fit reports `calc_ses = FALSE` with NA SEs (#304). Defaults
#'   to `TRUE` so the BETWFE / OLS callers (which have no `gls`) keep `q < 1`.
#' @param beta_hat Numeric vector of length `p`; the final slope estimates
#'   (caller-computed). For BETWFE intercept-only this is
#'   `beta_hat[2:(p+1)]` (all zero); for BETWFE no-treatment this is
#'   `beta_hat_slopes` after caller-applied `add_ridge` adjustment; for
#'   FETWFE intercept-only this is `rep(0, p)`; for FETWFE no-treatment
#'   this is `untransformCoefImproved(theta_hat_slopes)` after
#'   caller-applied `add_ridge` adjustment.
#' @param treat_inds,treat_int_inds Integer vectors; treatment-effect column
#'   indices.
#' @param cohort_probs,indep_cohort_probs Numeric vectors; cohort-probability
#'   weights computed by `prep_for_etwfe_regression()`.
#' @param cohort_probs_overall,indep_cohort_probs_overall Numeric scalars;
#'   overall cohort-probability weights.
#' @param sig_eps_sq,sig_eps_c_sq Numeric scalars; variance components.
#' @param lambda.max,lambda.min,lambda_star Numeric scalars; bridge penalty
#'   path endpoints and the BIC-selected lambda.
#' @param lambda.max_model_size,lambda.min_model_size,lambda_star_model_size
#'   Integer scalars; model sizes at each lambda.
#' @param X_ints Numeric matrix; original-basis design.
#' @param y Numeric vector; centered response.
#' @param X_final,y_final Numeric matrix/vector; GLS-whitened design and
#'   response actually fed to the bridge regression.
#' @param theta_hat Optional numeric vector of length `p + 1` (includes
#'   intercept). For FETWFE blocks only; BETWFE blocks omit (default
#'   `NULL`).
#' @param include_theta Logical; if `TRUE`, inserts `theta_hat` into the
#'   return list between `catt_df` and `beta_hat` (FETWFE field order).
#'   Default `FALSE` (BETWFE field order: no `theta_hat`).
#'
#' @return A 33-field (BETWFE) or 34-field (FETWFE) list with the same
#'   field ordering as the pre-refactor inline blocks. The ordering is
#'   load-bearing: downstream code in `betwfe()` / `fetwfe()` and the S3
#'   class machinery in `R/*_class.R` indexes the result by name, but the
#'   snapshot tests in `tests/testthat/_snaps/print-method-snapshot.md` and
#'   the constructors in `R/class_helpers.R` are sensitive to the *names()*
#'   ordering.
#'
#' @keywords internal
#' @noRd
.build_selected_out_result <- function(
	message_text,
	verbose,
	G,
	c_names,
	q,
	gls = TRUE,
	beta_hat,
	treat_inds,
	treat_int_inds,
	cohort_probs,
	indep_cohort_probs,
	cohort_probs_overall,
	indep_cohort_probs_overall,
	sig_eps_sq,
	sig_eps_c_sq,
	lambda.max,
	lambda.max_model_size,
	lambda.min,
	lambda.min_model_size,
	lambda_star,
	lambda_star_model_size,
	X_ints,
	y,
	X_final,
	y_final,
	N,
	T,
	d,
	p,
	theta_hat = NULL,
	include_theta = FALSE,
	cv_seed_used = NA_integer_
) {
	if (verbose) {
		message(message_text)
	}

	# A valid (zero) SE for this degenerate fit needs the SAME condition the normal
	# path uses (`fetwfe_core.R`: `calc_ses = (q < 1) && gls`): bridge selection
	# (q < 1) AND GLS whitening (gls). A `gls = FALSE` fit has no oracle SE, so its
	# degenerate SE is NA -- and `calc_ses`, `ret_se`, and `ret_var` must agree, or
	# the C1 SE-consistency validator (`att_se` must be NA when `calc_ses = FALSE`)
	# trips at fit construction. (#304 reconciliation: previously keyed on `q < 1`
	# alone, so a gls = FALSE degenerate fit wrongly reported calc_ses = TRUE.)
	ses_valid <- (q < 1) && gls
	ret_se <- if (ses_valid) 0 else NA

	catt_df_to_ret <- data.frame(
		cohort = c_names,
		estimate = rep(0, G),
		se = rep(ret_se, G),
		ci_low = rep(ret_se, G),
		ci_high = rep(ret_se, G),
		p_value = rep(NA_real_, G),
		selected = rep(FALSE, G),
		stringsAsFactors = FALSE
	)
	# Attach the `catt_df` S3 class so [[ / $ / [ accessors fire the
	# helpful-error layer (R/catt_df_class.R) when users hit the old
	# Title-Case column names from versions <= 1.10.0.
	class(catt_df_to_ret) <- c("catt_df", "data.frame")

	# Build the FULL union list with every possible field, then drop
	# the fields that don't belong to this block via name-vector
	# selection (per PR #92's `.summary_estimator_output()` pattern).
	# The `keep` vector encodes the contractual field ordering:
	# `theta_hat` slots between `catt_df` and `beta_hat` for FETWFE
	# blocks (3 + 4); BETWFE blocks (1 + 2) omit it.
	# Variance-component slots (issue #141/#146): when the bridge zeros
	# out every coefficient, `att_var_1` / `att_var_2` are zero (mirroring
	# `att_se = 0` under `q < 1`) or NA (mirroring `att_se = NA` under
	# `q >= 1`).
	ret_var <- if (ses_valid) 0 else NA

	out <- list(
		in_sample_att_hat = 0,
		in_sample_att_se = ret_se,
		in_sample_att_se_no_prob = ret_se,
		in_sample_att_var_1 = ret_var,
		in_sample_att_var_2 = ret_var,
		indep_att_hat = 0,
		indep_att_se = ret_se,
		indep_att_var_1 = ret_var,
		indep_att_var_2 = ret_var,
		catt_hats = setNames(rep(0, G), c_names),
		catt_ses = setNames(rep(ret_se, G), c_names),
		catt_df = catt_df_to_ret,
		theta_hat = theta_hat,
		beta_hat = beta_hat,
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
		calc_ses = ses_valid,
		# v1.13.0 (#164): CV-path provenance. NA_integer_ when reached
		# via the BIC path or an OLS-only early-exit (etwfe / twfeCovs).
		cv_seed_used = cv_seed_used
	)

	keep <- c(
		"in_sample_att_hat",
		"in_sample_att_se",
		"in_sample_att_se_no_prob",
		"in_sample_att_var_1",
		"in_sample_att_var_2",
		"indep_att_hat",
		"indep_att_se",
		"indep_att_var_1",
		"indep_att_var_2",
		"catt_hats",
		"catt_ses",
		"catt_df",
		if (include_theta) "theta_hat",
		"beta_hat",
		"treat_inds",
		"treat_int_inds",
		"cohort_probs",
		"indep_cohort_probs",
		"cohort_probs_overall",
		"indep_cohort_probs_overall",
		"sig_eps_sq",
		"sig_eps_c_sq",
		"lambda.max",
		"lambda.max_model_size",
		"lambda.min",
		"lambda.min_model_size",
		"lambda_star",
		"lambda_star_model_size",
		"X_ints",
		"y",
		"X_final",
		"y_final",
		"N",
		"T",
		"G",
		"d",
		"p",
		"calc_ses",
		"cv_seed_used"
	)
	out[keep]
}
