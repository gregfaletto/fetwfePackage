#' @title TWFE-With-Covariates Output Class
#' @description S3 class for the output of \code{twfeCovs()}. Minimal
#'   surface (coef + a bare print that preserves the pre-#76 behavior
#'   of just dumping the list); a full styled \code{print} / \code{summary}
#'   like the three sibling estimators is a separate follow-up.
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
# print() method
#
# Minimum-viable per #76 Item 1: `print(unclass(x))` preserves today's
# behavior (a bare list dump) so the class addition is a pure additive
# contract change. A full styled print is queued as a separate follow-up
# (#76 explicitly framed this as "a full S3 class is a separate
# follow-up").
#----------------------------------------------------------------------
#' @export
print.twfeCovs <- function(x, ...) {
	print(unclass(x))
}

#-------------------------------------------------------------------------------
# Constructor validator (#85). See R/fetwfe_class.R for the design rationale.
# twfeCovs differs from etwfe in: no `alpha` slot (twfeCovs has no inference
# output).
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
	"R",
	"d",
	"p",
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
	.check_type_sanity(x, cls, has_alpha = FALSE, has_att_selected = FALSE)
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
