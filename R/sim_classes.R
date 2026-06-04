#' @title Simulated Panel-Data Class
#' @description S3 class for objects returned by \code{simulateData()}.
#'   Compact `print` method summarizes the panel's dimensions and cohort
#'   structure instead of dumping the full `N*T x p` design matrix
#'   (which the default \code{print.list} would do).
#' @name FETWFE_simulated-class
NULL

#' @export
print.FETWFE_simulated <- function(x, ...) {
	cat(sprintf(
		"Simulated FETWFE panel: N = %d units, T = %d periods\n",
		x$N,
		x$T
	))
	cat(sprintf(
		"  G = %d treated cohorts, d = %d covariates, p = %d design columns\n",
		x$G,
		x$d,
		x$p
	))
	cohort_sizes <- paste(
		sprintf(
			"cohort %d = %d",
			seq_len(x$G),
			x$assignments[seq_len(x$G) + 1L]
		),
		collapse = ", "
	)
	cat(sprintf(
		"  Cohort sizes: never-treated = %d, %s\n",
		x$N_UNTREATED,
		cohort_sizes
	))
	cat(sprintf(
		"  Noise variances: sig_eps_sq = %g, sig_eps_c_sq = %g\n",
		x$sig_eps_sq,
		x$sig_eps_c_sq
	))
	invisible(x)
}

#' @title FETWFE Coefficient-Vector Class
#' @description S3 class for objects returned by \code{genCoefs()}.
#'   Compact `print` method summarizes the coefficient vector and its
#'   sparsity pattern instead of dumping the full \code{beta} and
#'   \code{theta} vectors.
#' @name FETWFE_coefs-class
NULL

#' @export
print.FETWFE_coefs <- function(x, ...) {
	n_nonzero <- sum(x$theta != 0)
	p <- length(x$theta)
	cat(sprintf(
		"FETWFE coefficient vector: G = %d cohorts, T = %d periods, d = %d covariates\n",
		x$G,
		x$T,
		x$d
	))
	cat(sprintf(
		"  beta length: %d; theta nonzeros: %d of %d (%.1f%%)\n",
		length(x$beta),
		n_nonzero,
		p,
		100 * n_nonzero / p
	))
	cat(sprintf(
		"  seed: %s\n",
		if (is.null(x$seed)) "<none>" else as.character(x$seed)
	))
	cat(sprintf(
		"  fusion structure: %s\n",
		if (is.null(x$fusion_structure)) "<unknown>" else x$fusion_structure
	))
	invisible(x)
}
