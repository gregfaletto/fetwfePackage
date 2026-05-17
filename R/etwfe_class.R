#' @title Extended Two-Way Fixed Effects Output Class
#' @description S3 class for the output of \code{etwfe()}.
#' @name etwfe-class
NULL

#----------------------------------------------------------------------
# coef() method
#----------------------------------------------------------------------
#' @export
coef.etwfe <- function(object, ...) object$beta_hat

#----------------------------------------------------------------------
# print() method for etwfe objects
#----------------------------------------------------------------------
#' @export
print.etwfe <- function(
	x,
	max_cohorts = getOption("etwfe.max_cohorts", 10),
	order_by = c("cohort", "estimate", "abs_estimate", "pvalue", "none"),
	show_internal = FALSE,
	...
) {
	order_by <- match.arg(order_by)

	cat(
		"Extended Two-Way Fixed Effects Results\n",
		"=====================================\n\n",
		sep = ""
	)

	## Overall ATT
	ci_pct <- 100 * (1 - x$alpha)
	ci_low <- x$att_hat - qnorm(1 - x$alpha / 2) * x$att_se
	ci_high <- x$att_hat + qnorm(1 - x$alpha / 2) * x$att_se
	cat(sprintf(
		"Overall Average Treatment Effect (ATT):\n  Estimate:   %.4f\n",
		x$att_hat
	))
	if (identical(x$se_type, "cluster")) {
		cat(sprintf(
			"  Std. Error (cluster-robust): %.4f\n",
			x$att_se
		))
	} else {
		cat(sprintf("  Std. Error: %.4f\n", x$att_se))
	}
	if (!is.null(x$att_p_value) && !is.na(x$att_p_value)) {
		cat(sprintf("  P-value:    %.4g\n", x$att_p_value))
	} else {
		cat("  P-value:    NA\n")
	}
	cat(sprintf(
		"  %.0f%% CI:    [%.4f, %.4f]\n\n",
		ci_pct,
		ci_low,
		ci_high
	))

	## Cohort effects
	catt_df <- .truncate_catt(x$catt_df, max_cohorts, order_by)
	cat("Cohort Average Treatment Effects (CATT):\n")
	.print_catt_tbl(catt_df)
	if (isTRUE(attr(catt_df, "truncated"))) {
		cat(sprintf(
			"  ... and %d more cohorts.\n",
			attr(catt_df, "n_discarded")
		))
	}
	cat("\n")

	## Model info
	cat("Model Details:\n")
	cat(sprintf("  Units (N)           : %d\n", x$N))
	cat(sprintf("  Time periods (T)    : %d\n", x$T))
	cat(sprintf("  Treated cohorts (R) : %d\n", x$R))
	cat(sprintf("  Covariates (d)      : %d\n", x$d))
	cat(sprintf("  Features (p)        : %d\n", x$p))

	if (show_internal) {
		cat("\nInternal Details:\n")
		cat(
			"  X dims     :",
			paste(dim(x$X_ints), collapse = " x "),
			"\n"
		)
		cat("  y length   :", length(x$y), "\n")
		cat("  SEs computed:", x$calc_ses, "\n")
	}

	invisible(x)
}

#----------------------------------------------------------------------
# summary()
#----------------------------------------------------------------------
#' @export
summary.etwfe <- function(object, full_catt = FALSE, ...) {
	list(
		att = c(
			estimate = object$att_hat,
			se = object$att_se,
			p_value = object$att_p_value
		),
		catt = if (full_catt) {
			object$catt_df
		} else {
			.truncate_catt(object$catt_df, max_cohorts = 20)
		},
		model_info = list(
			N = object$N,
			T = object$T,
			R = object$R,
			d = object$d,
			p = object$p,
			sig_eps_sq = object$sig_eps_sq,
			sig_eps_c_sq = object$sig_eps_c_sq
		),
		alpha = object$alpha,
		se_type = object$se_type
	) |>
		structure(class = "summary.etwfe")
}

#' @export
print.summary.etwfe <- function(x, ...) {
	cat(
		"Summary of Extended Two-Way Fixed Effects\n",
		"========================================\n\n",
		sep = ""
	)
	ci_pct <- 100 * (1 - x$alpha)
	ci_low <- x$att["estimate"] - qnorm(1 - x$alpha / 2) * x$att["se"]
	ci_high <- x$att["estimate"] + qnorm(1 - x$alpha / 2) * x$att["se"]
	p_val <- x$att["p_value"]
	p_str <- if (is.na(p_val)) "NA" else sprintf("%.4g", p_val)
	se_label <- if (identical(x$se_type, "cluster")) {
		"SE (cluster-robust)"
	} else {
		"SE"
	}
	cat(sprintf(
		"Overall ATT: %.4f  (%s = %.4f, p = %s, %.0f%% CI = [%.4f, %.4f])\n\n",
		x$att["estimate"],
		se_label,
		x$att["se"],
		p_str,
		ci_pct,
		ci_low,
		ci_high
	))

	cat("CATT (preview):\n")
	.print_catt_tbl(x$catt)
	if (isTRUE(attr(x$catt, "truncated"))) {
		cat(sprintf("  ... + %d more cohorts.\n", attr(x$catt, "n_discarded")))
	}
	cat("\n")

	## Model info
	cat("Model Details:\n")
	cat(sprintf("  Units (N)           : %d\n", x$model_info$N))
	cat(sprintf("  Time periods (T)    : %d\n", x$model_info$T))
	cat(sprintf("  Treated cohorts (R) : %d\n", x$model_info$R))
	cat(sprintf("  Covariates (d)      : %d\n", x$model_info$d))
	cat(sprintf("  Features (p)        : %d\n", x$model_info$p))

	invisible(x)
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
	"covs"
)

#' @title Validate an `etwfe`-classed object's contracts
#' @keywords internal
#' @noRd
.validate_etwfe <- function(x) {
	cls <- "etwfe"
	.stop_if_missing_slots(x, .EXPECTED_SLOTS_ETWFE, cls)
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
