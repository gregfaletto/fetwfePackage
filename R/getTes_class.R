#' @title Compute True Treatment Effects Output Class
#' @description S3 class for the output of \code{getTes()}.
#' @name FETWFE_tes-class
NULL

#' @export
print.FETWFE_tes <- function(x, ...) {
	cat(
		"True Treatment Effects (from FETWFE_coefs)\n",
		"==========================================\n",
		sep = ""
	)
	cat(sprintf("Overall true ATT: %.4f\n\n", x$att_true))
	cat("Cohort effects:\n")
	for (r in seq_along(x$actual_cohort_tes)) {
		cat(sprintf("  Cohort %d: %.4f\n", r, x$actual_cohort_tes[r]))
	}
	cat("\n")
	cat("Generated from:\n")
	cat(sprintf("  Cohorts (R)      : %d\n", x$R))
	cat(sprintf("  Time periods (T) : %d\n", x$T))
	cat(sprintf("  Covariates (d)   : %d\n", x$d))
	cat(sprintf(
		"  Seed             : %s\n",
		if (is.null(x$seed)) "<none>" else as.character(x$seed)
	))
	invisible(x)
}

#' @export
summary.FETWFE_tes <- function(object, ...) {
	tes <- object$actual_cohort_tes
	# Defensive: genCoefs() currently enforces R >= 2, so length(tes) >= 2 is
	# always true in practice. Guard is retained in case that validation ever
	# changes upstream (sd() of a length-1 vector is NA).
	out <- list(
		att_true = object$att_true,
		actual_cohort_tes = tes,
		R = object$R,
		T = object$T,
		d = object$d,
		seed = object$seed,
		cohort_te_stats = c(
			min = min(tes),
			max = max(tes),
			median = stats::median(tes),
			sd = if (length(tes) >= 2) stats::sd(tes) else NA_real_
		)
	)
	structure(out, class = "summary.FETWFE_tes")
}

#' @export
print.summary.FETWFE_tes <- function(x, ...) {
	cat(
		"Summary of True Treatment Effects (from FETWFE_coefs)\n",
		"=====================================================\n",
		sep = ""
	)
	cat(sprintf("Overall true ATT: %.4f\n\n", x$att_true))
	cat("Cohort effects:\n")
	for (r in seq_along(x$actual_cohort_tes)) {
		cat(sprintf("  Cohort %d: %.4f\n", r, x$actual_cohort_tes[r]))
	}
	cat("\n")
	cat("Cohort effect dispersion:\n")
	cat(sprintf("  min    : %.4f\n", x$cohort_te_stats["min"]))
	cat(sprintf("  max    : %.4f\n", x$cohort_te_stats["max"]))
	cat(sprintf("  median : %.4f\n", x$cohort_te_stats["median"]))
	cat(sprintf(
		"  sd     : %s\n",
		if (is.na(x$cohort_te_stats["sd"])) {
			"<NA>"
		} else {
			sprintf("%.4f", x$cohort_te_stats["sd"])
		}
	))
	cat("\n")
	cat("Generated from:\n")
	cat(sprintf("  Cohorts (R)      : %d\n", x$R))
	cat(sprintf("  Time periods (T) : %d\n", x$T))
	cat(sprintf("  Covariates (d)   : %d\n", x$d))
	cat(sprintf(
		"  Seed             : %s\n",
		if (is.null(x$seed)) "<none>" else as.character(x$seed)
	))
	invisible(x)
}
