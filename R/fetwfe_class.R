#' @title Fused Extended Two-Way Fixed Effects Output Class
#' @description S3 class for the output of \code{fetwfe()}.
#' @name fetwfe-class
NULL

#—--------------------------------------------------------------------
# coef() method (unchanged)
#—--------------------------------------------------------------------
#' @export
coef.fetwfe <- function(object, ...) object$beta_hat

#—--------------------------------------------------------------------
# print() method for fetwfe objects
#—--------------------------------------------------------------------
#' @export
print.fetwfe <- function(
	x,
	max_cohorts = getOption("fetwfe.max_cohorts", 10),
	order_by = c("cohort", "estimate", "abs_estimate", "pvalue", "none"),
	show_internal = FALSE,
	...
) {
	order_by <- match.arg(order_by)

	cat(
		"Fused Extended Two-Way Fixed Effects Results\n",
		"===========================================\n\n",
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
	if (!is.null(x$att_selected)) {
		cat(sprintf("  Selected:   %s\n", x$att_selected))
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
	cat(sprintf("  Selected size       : %d\n", x$lambda_star_model_size))
	cat(sprintf("  Lambda*             : %.4f\n", x$lambda_star))

	if (show_internal) {
		cat("\nInternal Details:\n")
		cat(
			"  X dims     :",
			paste(dim(x$internal$X_ints), collapse = " x "),
			"\n"
		)
		cat("  y length   :", length(x$internal$y), "\n")
		cat("  SEs computed:", x$internal$calc_ses, "\n")
	}

	invisible(x)
}

#—--------------------------------------------------------------------
# summary()
#—--------------------------------------------------------------------
#' @export
summary.fetwfe <- function(object, full_catt = FALSE, ...) {
	list(
		att = c(
			estimate = object$att_hat,
			se = object$att_se,
			p_value = object$att_p_value
		),
		att_selected = object$att_selected,
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
			lambda_star = object$lambda_star,
			model_size = object$lambda_star_model_size,
			sig_eps_sq = object$sig_eps_sq,
			sig_eps_c_sq = object$sig_eps_c_sq
		),
		alpha = object$alpha,
		se_type = object$se_type
	) |>
		structure(class = "summary.fetwfe")
}

#' @export
print.summary.fetwfe <- function(x, ...) {
	cat(
		"Summary of Fused Extended Two-Way Fixed Effects\n",
		"================================================\n\n",
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
		"Overall ATT: %.4f  (%s = %.4f, p = %s, %.0f%% CI = [%.4f, %.4f])\n",
		x$att["estimate"],
		se_label,
		x$att["se"],
		p_str,
		ci_pct,
		ci_low,
		ci_high
	))
	if (!is.null(x$att_selected)) {
		cat(sprintf("Selected: %s\n\n", x$att_selected))
	} else {
		cat("\n")
	}

	cat("CATT (preview):\n")
	.print_catt_tbl(x$catt)
	if (isTRUE(attr(x$catt, "truncated"))) {
		cat(sprintf("  ... + %d more cohorts.\n", attr(x$catt, "n_discarded")))
	}
	cat("\n")

	# cat("Model Info:\n")
	# print(x$model_info, row.names = FALSE, right = TRUE)

	## Model info
	cat("Model Details:\n")
	cat(sprintf("  Units (N)           : %d\n", x$model_info$N))
	cat(sprintf("  Time periods (T)    : %d\n", x$model_info$T))
	cat(sprintf("  Treated cohorts (R) : %d\n", x$model_info$R))
	cat(sprintf("  Covariates (d)      : %d\n", x$model_info$d))
	cat(sprintf("  Features (p)        : %d\n", x$model_info$p))
	cat(sprintf("  Selected size       : %d\n", x$model_info$model_size))
	cat(sprintf("  Lambda*             : %.4f\n", x$model_info$lambda_star))

	invisible(x)
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
	"N",
	"T",
	"R",
	"d",
	"p",
	"alpha",
	"indep_counts_used",
	"se_type",
	"y_mean",
	"response_col_name",
	"time_var",
	"unit_var",
	"treatment",
	"covs",
	"internal"
)

.EXPECTED_INTERNAL_SLOTS_FETWFE <- c(
	"X_ints",
	"y",
	"X_final",
	"y_final",
	"theta_hat",
	"calc_ses"
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
