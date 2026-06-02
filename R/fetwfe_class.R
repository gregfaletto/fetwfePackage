#' @title Fused Extended Two-Way Fixed Effects Output Class
#' @description S3 class for the output of \code{fetwfe()}.
#' @name fetwfe-class
NULL

#—--------------------------------------------------------------------
# coef() method (unchanged)
#—--------------------------------------------------------------------
#' @export
coef.fetwfe <- function(object, ...) {
	.check_for_coef(object)
	object$beta_hat
}

#—--------------------------------------------------------------------
# print() / summary() / print.summary() methods. Method bodies live in
# R/class_helpers.R as `.print_estimator_output()`,
# `.summary_estimator_output()`, `.print_summary_estimator_output()`;
# each per-class method here is a thin wrapper that pre-resolves the
# per-class options (header text, gating flags, path-extractor
# functions, max_cohorts / order_by) and delegates. Issue #77 step 2.
#—--------------------------------------------------------------------
#' @export
print.fetwfe <- function(
	x,
	max_cohorts = getOption("fetwfe.max_cohorts", 10),
	order_by = c("cohort", "estimate", "abs_estimate", "pvalue", "none"),
	show_internal = FALSE,
	max_event_times = getOption("fetwfe.max_event_times", 10),
	...
) {
	.print_estimator_output(
		x,
		header = "Fused Extended Two-Way Fixed Effects Results\n===========================================\n\n",
		show_att_selected = TRUE,
		show_lambda = TRUE,
		X_ints_path = function(x) x$internal$X_ints,
		y_path = function(x) x$internal$y,
		calc_ses_path = function(x) x$internal$calc_ses,
		max_cohorts = max_cohorts,
		order_by = match.arg(order_by),
		show_internal = show_internal,
		max_event_times = max_event_times,
		...
	)
}

#' @export
summary.fetwfe <- function(object, full_catt = FALSE, ...) {
	.summary_estimator_output(
		object,
		output_class = "summary.fetwfe",
		include_att_selected = TRUE,
		include_lambda = TRUE,
		full_catt = full_catt
	)
}

#' @export
print.summary.fetwfe <- function(x, ...) {
	.print_summary_estimator_output(
		x,
		header = "Summary of Fused Extended Two-Way Fixed Effects\n================================================\n\n",
		show_att_selected = TRUE,
		show_lambda = TRUE
	)
}

#-------------------------------------------------------------------------------
# Constructor validator (#85). Asserts the documented invariants for a
# `fetwfe`-classed object every time one is constructed. Called from the
# bottom of `fetwfe()` (R/fetwfe.R) before `class(out) <- "fetwfe"`.
#
# The list of expected slots is the source-of-truth for what a well-formed
# `fetwfe` object looks like; the doc-slot-parity test (#70) cross-checks
# that this list matches the rendered @return docs.
#-------------------------------------------------------------------------------

.EXPECTED_SLOTS_FETWFE <- c(
	"att_hat",
	"att_se",
	"att_p_value",
	"att_selected",
	"catt_hats",
	"catt_ses",
	"cohort_probs",
	"cohort_probs_overall",
	"catt_df",
	"beta_hat",
	"treat_inds",
	"treat_int_inds",
	"sig_eps_sq",
	"sig_eps_c_sq",
	"lambda.max",
	"lambda.max_model_size",
	"lambda.min",
	"lambda.min_model_size",
	"lambda_star",
	"lambda_star_model_size",
	"lambda_selection",
	"cv_folds",
	"cv_seed",
	"N",
	"T",
	"G",
	"R",
	"d",
	"p",
	"alpha",
	"calc_ses",
	"indep_counts_used",
	"se_type",
	"y_mean",
	"response_col_name",
	"time_var",
	"unit_var",
	"treatment",
	"covs",
	"ci_type",
	"internal"
)

.EXPECTED_INTERNAL_SLOTS_FETWFE <- c(
	"X_ints",
	"y",
	"X_final",
	"y_final",
	"theta_hat",
	"calc_ses",
	"variance_components",
	"first_year"
)

#' @title Validate a `fetwfe`-classed object's contracts
#' @description
#' Asserts the cross-slot invariants of a well-formed `fetwfe` object.
#' Stops with a structured error message on violation. Called from the
#' bottom of `fetwfe()` (R/fetwfe.R) before class assignment. Also callable
#' externally (`fetwfe:::.validate_fetwfe(x)`) for use by method-entry
#' preconditions (#86).
#' @param x An object to validate as `fetwfe`-shaped.
#' @return `invisible(x)` if all contracts hold; `stop()`s otherwise.
#' @keywords internal
#' @noRd
.validate_fetwfe <- function(x) {
	cls <- "fetwfe"
	.stop_if_missing_slots(x, .EXPECTED_SLOTS_FETWFE, cls)
	.stop_if_missing_slots(
		x$internal,
		.EXPECTED_INTERNAL_SLOTS_FETWFE,
		cls,
		where = "internal"
	)
	.check_type_sanity(x, cls, has_alpha = TRUE, has_att_selected = TRUE)
	.check_se_consistency(x, calc_ses_path = "internal$calc_ses", cls)
	.check_selection_consistency(x, cls)
	.check_p_value_na(x, cls)
	.check_catt_df_shape(x, cls)
	.check_ci_band_width(x, cls)
	.check_cohort_probs(x, cls)
	.check_lambda_monotonicity(x, cls)
	# C6 dimensions (internal-nested for fetwfe)
	.assert_contract(
		length(x$beta_hat) == x$p,
		"C6 length(beta_hat) == p",
		cls
	)
	.assert_contract(
		length(x$internal$y) == x$N * x$T,
		"C6 length(internal$y) == N * T",
		cls
	)
	.assert_contract(
		nrow(x$internal$X_ints) == x$N * x$T,
		"C6 nrow(internal$X_ints) == N * T",
		cls
	)
	invisible(x)
}
