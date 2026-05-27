#' @title Bridge-Penalized Extended Two-Way Fixed Effects Output Class
#' @description S3 class for the output of \code{betwfe()}.
#' @name betwfe-class
NULL

#—--------------------------------------------------------------------
# coef() method
#—--------------------------------------------------------------------
#' @export
coef.betwfe <- function(object, ...) {
	.check_for_coef(object)
	object$beta_hat
}

#—--------------------------------------------------------------------
# print() / summary() / print.summary() methods. Method bodies live in
# R/class_helpers.R; each per-class method here is a thin wrapper.
# See R/fetwfe_class.R for the design rationale (issue #77 step 2).
#
# `betwfe` is the hybrid: has `att_selected` and lambda.* slots like
# fetwfe (so `show_att_selected = TRUE`, `show_lambda = TRUE`), but
# `X_ints` / `y` / `calc_ses` live at the TOP level like etwfe (so the
# three `_path` extractor functions read top-level fields, NOT the
# `$internal`-nested shape that fetwfe uses). Asymmetry documented at
# .plans/feat-etwfe-betwfe-vignette-17/PLAN.md.
#—--------------------------------------------------------------------
#' @export
print.betwfe <- function(
	x,
	max_cohorts = getOption("betwfe.max_cohorts", 10),
	order_by = c("cohort", "estimate", "abs_estimate", "pvalue", "none"),
	show_internal = FALSE,
	max_event_times = getOption("betwfe.max_event_times", 10),
	...
) {
	.print_estimator_output(
		x,
		header = "Bridge-Penalized Extended Two-Way Fixed Effects Results\n=======================================================\n\n",
		show_att_selected = TRUE,
		show_lambda = TRUE,
		X_ints_path = function(x) x$X_ints,
		y_path = function(x) x$y,
		calc_ses_path = function(x) x$calc_ses,
		max_cohorts = max_cohorts,
		order_by = match.arg(order_by),
		show_internal = show_internal,
		max_event_times = max_event_times,
		...
	)
}

#' @export
summary.betwfe <- function(object, full_catt = FALSE, ...) {
	.summary_estimator_output(
		object,
		output_class = "summary.betwfe",
		include_att_selected = TRUE,
		include_lambda = TRUE,
		full_catt = full_catt
	)
}

#' @export
print.summary.betwfe <- function(x, ...) {
	.print_summary_estimator_output(
		x,
		header = "Summary of Bridge-Penalized Extended Two-Way Fixed Effects\n==========================================================\n\n",
		show_att_selected = TRUE,
		show_lambda = TRUE
	)
}

#-------------------------------------------------------------------------------
# Constructor validator (#85). See R/fetwfe_class.R for the design rationale.
# BETWFE is the hybrid: has `att_selected` and lambda.* slots like fetwfe,
# but X_ints/y/X_final/y_final/calc_ses are at TOP level like etwfe.
#-------------------------------------------------------------------------------

.EXPECTED_SLOTS_BETWFE <- c(
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
	"X_ints",
	"y",
	"X_final",
	"y_final",
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
	# Phase 8 (#33): see fetwfe_class.R for description.
	"cohort_means_external",
	"internal"
)

# Inner slots under `$internal` (#144). Parallel to
# `.EXPECTED_INTERNAL_SLOTS_FETWFE` in `R/fetwfe_class.R`. The OLS-family
# `$internal` omits `theta_hat` (bridge-selection-only, no analog here).
# The first five sub-slots (`X_ints`, `y`, `X_final`, `y_final`,
# `calc_ses`) are also duplicated at top level for backward compat per
# the pure-additive strategy adopted in PR #154; `variance_components`
# and `first_year` live only under `$internal`. Both top-level and
# `$internal` are canonical and kept in sync (#180).
.EXPECTED_INTERNAL_SLOTS_BETWFE <- c(
	"X_ints",
	"y",
	"X_final",
	"y_final",
	"calc_ses",
	"variance_components",
	"first_year"
)

#' @title Validate a `betwfe`-classed object's contracts
#' @keywords internal
#' @noRd
.validate_betwfe <- function(x) {
	cls <- "betwfe"
	.stop_if_missing_slots(x, .EXPECTED_SLOTS_BETWFE, cls)
	.stop_if_missing_slots(
		x$internal,
		.EXPECTED_INTERNAL_SLOTS_BETWFE,
		cls,
		where = "internal"
	)
	.check_type_sanity(x, cls, has_alpha = TRUE, has_att_selected = TRUE)
	.check_se_consistency(x, calc_ses_path = "calc_ses", cls)
	.check_selection_consistency(x, cls)
	.check_p_value_na(x, cls)
	.check_catt_df_shape(x, cls)
	.check_ci_band_width(x, cls)
	.check_cohort_probs(x, cls)
	.check_lambda_monotonicity(x, cls)
	.assert_contract(
		length(x$beta_hat) == x$p,
		"C6 length(beta_hat) == p",
		cls
	)
	.assert_contract(
		length(x$y) == x$N * x$T,
		"C6 length(y) == N * T",
		cls
	)
	.assert_contract(
		nrow(x$X_ints) == x$N * x$T,
		"C6 nrow(X_ints) == N * T",
		cls
	)
	.assert_contract(
		is.logical(x$calc_ses) && length(x$calc_ses) == 1L,
		"C8 calc_ses is length-1 logical",
		cls
	)
	invisible(x)
}
