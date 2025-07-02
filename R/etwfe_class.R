#' @title Extended Two-Way Fixed Effects Output Class
#' @description S3 class for the output of \code{etwfe()}.
#' @name etwfe-class
NULL

#—--------------------------------------------------------------------
# Utility: truncate / reorder cohort table
#—--------------------------------------------------------------------
.truncate_catt <- function(
	catt_df,
	max_cohorts = getOption("etwfe.max_cohorts", 10),
	order_by = c("cohort", "estimate", "abs_estimate", "pvalue", "none")
) {
	order_by <- match.arg(order_by)

	idx <- switch(
		order_by,
		cohort = order(catt_df$Cohort),
		estimate = order(catt_df$`Estimated TE`),
		abs_estimate = order(abs(catt_df$`Estimated TE`), decreasing = TRUE),
		pvalue = order(catt_df$P_value),
		none = seq_len(nrow(catt_df))
	)
	catt_df <- catt_df[idx, , drop = FALSE]

	if (nrow(catt_df) > max_cohorts) {
		attr(catt_df, "truncated") <- TRUE
		attr(catt_df, "n_discarded") <- nrow(catt_df) - max_cohorts
		catt_df <- catt_df[seq_len(max_cohorts), , drop = FALSE]
	} else {
		attr(catt_df, "truncated") <- FALSE
	}
	catt_df
}

#—--------------------------------------------------------------------
# Simple printer for cohort table (no colour)
#—--------------------------------------------------------------------
.print_catt_tbl <- function(df) {
	print(df, row.names = FALSE, right = TRUE)
	invisible()
}

#—--------------------------------------------------------------------
# coef() method
#—--------------------------------------------------------------------
#' @export
coef.etwfe <- function(object, ...) object$beta_hat

#—--------------------------------------------------------------------
# print() method for etwfe objects
#—--------------------------------------------------------------------
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
		"Overall Average Treatment Effect (ATT):\n  Estimate: %.4f\n",
		x$att_hat
	))
	cat(sprintf(
		"  Std. Error: %.4f\n  %.0f%% CI: [%.4f, %.4f]\n\n",
		x$att_se,
		ci_pct,
		ci_low,
		ci_high
	))

	## Cohort effects
	catt_df <- .truncate_catt(x$catt_df, max_cohorts, order_by)
	cat("Cohort Average Treatment Effects (CATT):\n")
	.print_catt_tbl(catt_df)
	if (isTRUE(attr(catt_df, "truncated")))
		cat(sprintf(
			"  ... and %d more cohorts.\n",
			attr(catt_df, "n_discarded")
		))
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

#—--------------------------------------------------------------------
# summary()
#—--------------------------------------------------------------------
#' @export
summary.etwfe <- function(object, full_catt = FALSE, ...) {
	list(
		att = c(estimate = object$att_hat, se = object$att_se),
		catt = if (full_catt) object$catt_df else
			.truncate_catt(object$catt_df, max_cohorts = 20),
		model_info = list(
			N = object$N,
			T = object$T,
			R = object$R,
			d = object$d,
			p = object$p,
			sig_eps_sq = object$sig_eps_sq,
			sig_eps_c_sq = object$sig_eps_c_sq
		),
		alpha = object$alpha
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
	cat(sprintf(
		"Overall ATT: %.4f  (SE = %.4f, %.0f%% CI = [%.4f, %.4f])\n\n",
		x$att["estimate"],
		x$att["se"],
		ci_pct,
		ci_low,
		ci_high
	))

	cat("CATT (preview):\n")
	.print_catt_tbl(x$catt)
	if (isTRUE(attr(x$catt, "truncated")))
		cat(sprintf("  ... + %d more cohorts.\n", attr(x$catt, "n_discarded")))
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
