#' @title Extended Two-Way Fixed Effects Output Class
#' @description S3 class for the output of \code{etwfe()}.
#' @name etwfe-class
NULL

#----------------------------------------------------------------------
# coef() method
#----------------------------------------------------------------------
#' @export
coef.etwfe <- function(object, ...) {
	.check_for_coef(object)
	object$beta_hat
}

#----------------------------------------------------------------------
# print() / summary() / print.summary() methods. Method bodies live in
# R/class_helpers.R; each per-class method here is a thin wrapper.
# See R/fetwfe_class.R for the design rationale (issue #77 step 2).
#----------------------------------------------------------------------
#' @export
print.etwfe <- function(
	x,
	max_cohorts = getOption("etwfe.max_cohorts", 10),
	order_by = c("cohort", "estimate", "abs_estimate", "pvalue", "none"),
	show_internal = FALSE,
	max_event_times = getOption("etwfe.max_event_times", 10),
	...
) {
	.print_estimator_output(
		x,
		header = "Extended Two-Way Fixed Effects Results\n=====================================\n\n",
		show_att_selected = FALSE,
		show_lambda = FALSE,
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
summary.etwfe <- function(object, full_catt = FALSE, ...) {
	.summary_estimator_output(
		object,
		output_class = "summary.etwfe",
		include_att_selected = FALSE,
		include_lambda = FALSE,
		full_catt = full_catt
	)
}

#' @export
print.summary.etwfe <- function(x, ...) {
	.print_summary_estimator_output(
		x,
		header = "Summary of Extended Two-Way Fixed Effects\n========================================\n\n",
		show_att_selected = FALSE,
		show_lambda = FALSE
	)
}

#-------------------------------------------------------------------------------
# Constructor validator (#85). See R/fetwfe_class.R for the design rationale.
# ETWFE differs from FETWFE in: no `att_selected`, no `internal` sublist
# (X_ints/y/X_final/y_final/calc_ses are top-level), no lambda.* slots
# (pure OLS -- no bridge regularization).
#-------------------------------------------------------------------------------

.EXPECTED_SLOTS_ETWFE <- c(
	"att_hat",
	"att_se",
	"att_p_value",
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
	"X_ints",
	"y",
	"X_final",
	"y_final",
	"N",
	"T",
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
	# Phase 8 (#33): see fetwfe_class.R for description.
	"cohort_means_external",
	"internal"
)

# Inner slots under `$internal` (#144). Parallel to
# `.EXPECTED_INTERNAL_SLOTS_FETWFE` in `R/fetwfe_class.R`. The OLS-family
# `$internal` omits `theta_hat` (bridge-selection-only, no analog here).
# Same values are also duplicated at top level for backward compat per
# the pure-additive strategy adopted in PR #154; the canonical path is
# `$internal` going forward.
.EXPECTED_INTERNAL_SLOTS_ETWFE <- c(
	"X_ints",
	"y",
	"X_final",
	"y_final",
	"calc_ses"
)

#' @title Validate an `etwfe`-classed object's contracts
#' @keywords internal
#' @noRd
.validate_etwfe <- function(x) {
	cls <- "etwfe"
	.stop_if_missing_slots(x, .EXPECTED_SLOTS_ETWFE, cls)
	.stop_if_missing_slots(
		x$internal,
		.EXPECTED_INTERNAL_SLOTS_ETWFE,
		cls,
		where = "internal"
	)
	.check_type_sanity(x, cls, has_alpha = TRUE, has_att_selected = FALSE)
	.check_se_consistency(x, calc_ses_path = "calc_ses", cls)
	.check_p_value_na(x, cls)
	.check_catt_df_shape(x, cls)
	.check_cohort_probs(x, cls)
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
