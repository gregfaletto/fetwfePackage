#' @title Compute True Treatment Effects Output Class
#' @description S3 class for the output of \code{getTes()}.
#' @name FETWFE_tes-class
NULL

# Render the "Cohort assignment DGP" section shared by print.FETWFE_tes and
# print.summary.FETWFE_tes (#189). Gated on a non-marginal assignment_type, so
# marginal objects -- and pre-1.14.0 objects that predate the slot -- render
# nothing and their output is unchanged.
#' @keywords internal
#' @noRd
.cat_tes_assignment_dgp <- function(x) {
	if (is.null(x$assignment_type) || x$assignment_type == "marginal") {
		return(invisible(NULL))
	}
	cat("Cohort assignment DGP:\n")
	cat(sprintf("  %-16s: %s\n", "Type", x$assignment_type))
	if (!is.null(x$assignment_strength)) {
		cat(sprintf("  %-16s: %.4f\n", "Strength", x$assignment_strength))
	}
	cat(sprintf(
		"  %-16s: %s\n",
		"Cohort weights",
		paste(sprintf("%.4f", x$cohort_weights), collapse = ", ")
	))
	cat("\n")
	invisible(NULL)
}

# Resolve cohort labels for a FETWFE_tes (or summary.FETWFE_tes) object by
# calendar adoption time (cohort g adopts at time g + 1), matching
# tidy.FETWFE_tes() and the fitted estimators' catt_df$cohort (#261). Falls back
# to the simulator convention for pre-1.9.0 objects lacking the slot.
.resolve_tes_cohort_times <- function(x) {
	if (!is.null(x$cohort_times)) {
		x$cohort_times
	} else {
		as.integer(seq_along(x$actual_cohort_tes) + 1L)
	}
}

#' @export
print.FETWFE_tes <- function(x, ...) {
	cat(
		"True Treatment Effects (from FETWFE_coefs)\n",
		"==========================================\n",
		sep = ""
	)
	cat(sprintf("Overall true ATT: %.4f\n\n", x$att_true))
	cohort_times <- .resolve_tes_cohort_times(x)
	cat("Cohort effects:\n")
	for (g in seq_along(x$actual_cohort_tes)) {
		cat(sprintf(
			"  Cohort %d: %.4f\n",
			cohort_times[g],
			x$actual_cohort_tes[g]
		))
	}
	cat("\n")
	.cat_tes_assignment_dgp(x)
	cat("Generated from:\n")
	cat(sprintf("  Cohorts (G)      : %d\n", x$G))
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
	# genCoefs() supports a single treated cohort (G = 1), in which case
	# length(tes) == 1 and sd() of a length-1 vector is NA; the length(tes) >= 2
	# guard below keeps the spread statistics well defined in that case.
	out <- list(
		att_true = object$att_true,
		actual_cohort_tes = tes,
		G = object$G,
		R = object$G,
		T = object$T,
		d = object$d,
		seed = object$seed,
		cohort_weights = object$cohort_weights,
		cohort_times = object$cohort_times,
		assignment_type = object$assignment_type,
		assignment_strength = object$assignment_strength,
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
	cohort_times <- .resolve_tes_cohort_times(x)
	cat("Cohort effects:\n")
	for (g in seq_along(x$actual_cohort_tes)) {
		cat(sprintf(
			"  Cohort %d: %.4f\n",
			cohort_times[g],
			x$actual_cohort_tes[g]
		))
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
	.cat_tes_assignment_dgp(x)
	cat("Generated from:\n")
	cat(sprintf("  Cohorts (G)      : %d\n", x$G))
	cat(sprintf("  Time periods (T) : %d\n", x$T))
	cat(sprintf("  Covariates (d)   : %d\n", x$d))
	cat(sprintf(
		"  Seed             : %s\n",
		if (is.null(x$seed)) "<none>" else as.character(x$seed)
	))
	invisible(x)
}
