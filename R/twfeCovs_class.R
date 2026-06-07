#' @title TWFE-With-Covariates Output Class
#' @description S3 class for the output of \code{twfeCovs()}. Carries the same
#'   styled \code{print} / \code{summary} / \code{coef} surface as the three
#'   sibling estimators, plus \code{tidy} / \code{glance} (broom) and
#'   \code{simultaneousCIs}. \code{plot} and \code{augment} are intentionally
#'   not provided: \code{twfeCovs()} estimates one pooled effect per cohort, so
#'   it has no per-(cohort, time) / event-study basis to plot, and its
#'   coefficient vector is in a reduced basis that \code{augment()}'s
#'   fitted-value path does not match (#58). Both raise an informative error.
#' @name twfeCovs-class
NULL

#----------------------------------------------------------------------
# coef() method
#----------------------------------------------------------------------
#' @export
coef.twfeCovs <- function(object, ...) {
	.check_for_coef(object)
	object$beta_hat
}

#----------------------------------------------------------------------
# print() / summary() / print.summary() methods. Thin wrappers over the
# shared helpers in R/class_helpers.R (parallel to R/etwfe_class.R), with
# `include_event_study = FALSE`: twfeCovs has one pooled effect per cohort,
# so there is no event-study section (#58; upgrades the pre-#76 bare list
# dump to the styled sibling-estimator output).
#----------------------------------------------------------------------
#' @export
print.twfeCovs <- function(
	x,
	max_cohorts = getOption("twfeCovs.max_cohorts", 10),
	order_by = c("cohort", "estimate", "abs_estimate", "pvalue", "none"),
	show_internal = FALSE,
	max_event_times = getOption("twfeCovs.max_event_times", 10),
	...
) {
	.print_estimator_output(
		x,
		header = "TWFE (with covariates) Results\n==============================\n\n",
		show_att_selected = FALSE,
		show_lambda = FALSE,
		X_ints_path = function(x) x$X_ints,
		y_path = function(x) x$y,
		calc_ses_path = function(x) x$calc_ses,
		max_cohorts = max_cohorts,
		order_by = match.arg(order_by),
		show_internal = show_internal,
		max_event_times = max_event_times,
		include_event_study = FALSE,
		...
	)
}

#' @export
summary.twfeCovs <- function(object, full_catt = FALSE, ...) {
	.summary_estimator_output(
		object,
		output_class = "summary.twfeCovs",
		include_att_selected = FALSE,
		include_lambda = FALSE,
		full_catt = full_catt,
		include_event_study = FALSE
	)
}

#' @export
print.summary.twfeCovs <- function(x, ...) {
	.print_summary_estimator_output(
		x,
		header = "Summary of TWFE (with covariates)\n=================================\n\n",
		show_att_selected = FALSE,
		show_lambda = FALSE
	)
}

#-------------------------------------------------------------------------------
# Constructor validator (#85). See R/fetwfe_class.R for the design rationale.
# twfeCovs carries an `alpha` slot like etwfe (#204): the user's `alpha` is
# threaded into the default simultaneous catt_df band and read by the C10
# band-width validator.
#-------------------------------------------------------------------------------

.EXPECTED_SLOTS_TWFECOVS <- c(
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
	"G",
	"R",
	"d",
	"p",
	"calc_ses",
	"indep_counts_used",
	"se_type",
	"alpha",
	"y_mean",
	"response_col_name",
	"time_var",
	"unit_var",
	"treatment",
	"covs",
	"ci_type",
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
.EXPECTED_INTERNAL_SLOTS_TWFECOVS <- c(
	"X_ints",
	"y",
	"X_final",
	"y_final",
	"calc_ses",
	"variance_components",
	"first_year"
)

#' @title Validate a `twfeCovs`-shaped object's contracts
#' @keywords internal
#' @noRd
.validate_twfeCovs <- function(x) {
	cls <- "twfeCovs"
	.stop_if_missing_slots(x, .EXPECTED_SLOTS_TWFECOVS, cls)
	.stop_if_missing_slots(
		x$internal,
		.EXPECTED_INTERNAL_SLOTS_TWFECOVS,
		cls,
		where = "internal"
	)
	.check_type_sanity(x, cls, has_alpha = TRUE, has_att_selected = FALSE)
	.check_se_consistency(x, calc_ses_path = "calc_ses", cls)
	.check_p_value_na(x, cls)
	.check_catt_df_shape(x, cls)
	.check_ci_band_width(x, cls)
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
