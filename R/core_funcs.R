#' @title Prepare Transformed Design Matrix and Response for ETWFE Regression
#'
#' @description
#' Generates all matrix- and vector-level inputs required by the bridge/OLS
#' fitting step inside \code{\link{etwfe_core}}.
#' Its responsibilities include: estimating or assembling the covariance
#' matrix `\Omega`; performing the GLS whitening transformation;
#' (optionally) appending ridge-penalty rows; computing cohort probability
#' weights; and returning both raw and scaled versions of the final design
#' matrix.
#'
#' @param verbose Logical.  If \code{TRUE}, prints timing information for
#'   each major sub-task.  Default \code{FALSE}.
#' @param sig_eps_sq,sig_eps_c_sq Numeric or \code{NA}.  Observation-level
#'   and unit-level variance components.  If either is \code{NA},
#'   they are estimated via \code{estOmegaSqrtInv()}.
#' @param y Numeric vector of length \(N\times T\).  Centered response.
#' @param X_ints Numeric matrix.  Full design matrix in the \emph{original}
#'   parameterisation.
#' @param X_mod Numeric matrix.  Same dimensions as \code{X_ints}.  In the
#'   ETWFE path this is usually a transformed version, but for OLS it is
#'   often identical to \code{X_ints}.
#' @param N,T,R,d,p,num_treats Integers giving key problem dimensions:
#'   number of units, time periods, treated cohorts, covariates,
#'   total parameters, and base treatment-effect parameters, respectively.
#' @param add_ridge Logical.  Whether to append rows that implement a small
#'   L2 penalty on the \emph{untransformed} coefficients.  Default \code{FALSE}.
#' @param first_inds Integer vector (length \(R\)).  Column indices of the
#'   first base treatment-effect for each cohort within the block of
#'   \code{num_treats} columns.
#' @param in_sample_counts Integer vector of cohort sizes in the estimation
#'   sample (see \code{check_etwfe_core_inputs}).
#' @param indep_count_data_available Logical.  Indicates whether
#'   \code{indep_counts} contains valid cohort counts from an independent
#'   sample; affects the probability calculations below.
#' @param indep_counts Optional integer vector.  Cohort counts from an
#'   independent sample; only used when
#'   \code{indep_count_data_available = TRUE}.  Default \code{NA}.
#' @param is_fetwfe Logical.  If \code{TRUE}, a fusion transformation matrix
#'   has been applied upstream (this matters for how the optional ridge rows
#'   are constructed).  Required (no default): the prior default of
#'   \code{TRUE} produced GitHub #74 (BETWFE silently inherited the FETWFE
#'   augmentation basis).  All four caller sites
#'   (\code{fetwfe_core}, \code{etwfe_core}, \code{betwfe_core},
#'   \code{twfeCovs_core}) must pass this argument explicitly.
#' @param is_twfe_covs Logical.  If \code{TRUE}, columns will be removed and
#'   consolidated as required for `twfeCovs()`.  Default \code{FALSE}.
#'
#' @details
#' The routine carries out the following steps in order:
#' \enumerate{
#'   \item \strong{Variance-Component Handling}\newline
#'     Calls \code{estOmegaSqrtInv()} when either variance component is
#'     unknown; otherwise uses supplied values.
#'   \item \strong{GLS Transformation}\newline
#'     Forms \(\Omega = \sigma_\varepsilon^2 I_T + \sigma_{\varepsilon c}^2 J_T\)
#'     and multiplies both \code{y} and \code{X_mod} on the left by
#'     \(\sqrt{\sigma_\varepsilon^2}\,\Omega^{-1/2} \otimes I_N\).
#'   \item \strong{Optional Ridge Augmentation}\newline
#'     If \code{add_ridge = TRUE}, appends \(p\) rows to the scaled design
#'     matrix and response, thereby imposing a tiny L2 penalty equal to
#'     \code{lambda_ridge}.  When \code{is_fetwfe = TRUE} the rows are first
#'     premultiplied by the inverse fusion transformation so the penalty
#'     applies on the original coefficient scale.
#'   \item \strong{Cohort Probability Estimation}\newline
#'     Computes \eqn{\hat\pi_r = n_r / \sum_{s=1}^R n_s} from the in-sample
#'     counts; when an independent split is available the same is done for
#'     \code{indep_counts}.
#' }
#'
#' @return A list with everything \code{etwfe_core} needs next:
#'   \describe{
#'     \item{\code{X_final_scaled}}{Design matrix after GLS transform
#'       and column-wise scaling (and ridge rows if requested).}
#'     \item{\code{X_final}}{Design matrix after GLS transform but before scaling.}
#'     \item{\code{y_final}}{Response after GLS transform (and augmentation).}
#'     \item{\code{scale_center}, \code{scale_scale}}{Vectors used to undo
#'       column scaling.}
#'     \item{\code{cohort_probs}}{Vector of \(\hat\pi_r \mid \text{treated}\)
#'       from the estimation sample.}
#'     \item{\code{cohort_probs_overall}}{Vector of unconditional
#'       probabilities \(P(W=r)\).}
#'     \item{\code{indep_cohort_probs}, \code{indep_cohort_probs_overall}}{Same
#'       probabilities if an independent split was provided; otherwise \code{NA}.}
#'     \item{\code{sig_eps_sq}, \code{sig_eps_c_sq}}{Possibly estimated
#'       variance components carried forward.}
#'     \item{\code{lambda_ridge}}{Numeric scalar.  Value of the ridge penalty
#'       used (or \code{NA} if none).}
#'   }
#' @keywords internal
#' @noRd
prep_for_etwfe_regression <- function(
	verbose,
	sig_eps_sq,
	sig_eps_c_sq,
	y,
	X_ints,
	X_mod,
	N,
	T,
	R,
	d,
	p,
	num_treats,
	add_ridge,
	first_inds,
	in_sample_counts,
	indep_count_data_available,
	indep_counts = NA,
	is_fetwfe,
	is_twfe_covs = FALSE
) {
	gls <- .estimate_variance_and_gls(
		y = y,
		X_ints = X_ints,
		X_mod = X_mod,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		N = N,
		T = T,
		p = p,
		verbose = verbose
	)
	y_final <- gls$y_gls
	X_final <- gls$X_gls
	sig_eps_sq <- gls$sig_eps_sq
	sig_eps_c_sq <- gls$sig_eps_c_sq

	if (is_twfe_covs) {
		coll <- .collapse_design_for_twfe_covs(
			X_gls = X_final,
			N = N,
			T = T,
			R = R,
			d = d,
			num_treats = num_treats,
			first_inds = first_inds
		)
		X_final <- coll$X_collapsed
		p <- coll$p_short
	}

	X_final_scaled <- my_scale(X_final)
	stopifnot(ncol(X_final_scaled) == p)
	scale_center <- attr(X_final_scaled, "scaled:center")
	scale_scale <- attr(X_final_scaled, "scaled:scale")

	ridge <- .append_ridge_rows(
		X_scaled = X_final_scaled,
		y_gls = y_final,
		p = p,
		add_ridge = add_ridge,
		is_fetwfe = is_fetwfe,
		first_inds = first_inds,
		T = T,
		R = R,
		d = d,
		num_treats = num_treats,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		N = N
	)
	X_final_scaled <- ridge$X_scaled
	y_final <- ridge$y_final
	lambda_ridge <- ridge$lambda_ridge

	probs <- .compute_cohort_probs(
		in_sample_counts = in_sample_counts,
		indep_counts = indep_counts,
		N = N,
		R = R,
		indep_count_data_available = indep_count_data_available
	)

	list(
		X_final_scaled = X_final_scaled,
		X_final = X_final,
		y_final = y_final,
		scale_center = scale_center,
		scale_scale = scale_scale,
		cohort_probs = probs$cohort_probs,
		cohort_probs_overall = probs$cohort_probs_overall,
		indep_cohort_probs = probs$indep_cohort_probs,
		indep_cohort_probs_overall = probs$indep_cohort_probs_overall,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		lambda_ridge = lambda_ridge
	)
}

#' @title Estimate variance components and apply GLS whitening
#' @description If `sig_eps_sq` or `sig_eps_c_sq` is NA, estimate both via
#'   `estOmegaSqrtInv()` (REML on `y ~ X + (1 | unit)`). Then build
#'   `Omega = sig_eps_sq * I_T + sig_eps_c_sq * J_T`, take its
#'   matrix square root inverse, and apply the kronecker-product GLS
#'   transform to `y` and `X_mod`.
#' @param y Numeric vector of length `N*T`; raw response.
#' @param X_ints Numeric matrix used by REML when variance components are
#'   NA.
#' @param X_mod Numeric matrix to GLS-transform.
#' @param sig_eps_sq,sig_eps_c_sq Numeric or NA; row-level / unit-level
#'   variance components.
#' @param N,T,p Integers; units, time periods, design columns.
#' @param verbose Logical; if TRUE, print timing messages.
#' @return List with `y_gls`, `X_gls`, `sig_eps_sq`, `sig_eps_c_sq`.
#' @keywords internal
#' @noRd
.estimate_variance_and_gls <- function(
	y,
	X_ints,
	X_mod,
	sig_eps_sq,
	sig_eps_c_sq,
	N,
	T,
	p,
	verbose
) {
	if (verbose) {
		message("Getting omega sqrt inverse estimate...")
		t0 <- Sys.time()
	}

	if (is.na(sig_eps_sq) | is.na(sig_eps_c_sq)) {
		omega_res <- estOmegaSqrtInv(
			y,
			X_ints,
			N = N,
			T = T,
			p = p
		)

		sig_eps_sq <- omega_res$sig_eps_sq
		sig_eps_c_sq <- omega_res$sig_eps_c_sq

		rm(omega_res)

		if (verbose) {
			message("Done! Time to estimate noise variances:")
			message(Sys.time() - t0)
			t0 <- Sys.time()
		}
	}

	stopifnot(!is.na(sig_eps_sq) & !is.na(sig_eps_c_sq))

	# Closed-form Omega^(-1/2) for the equicorrelated-noise covariance
	# Omega = sig_eps_sq * I_T + sig_eps_c_sq * J_T (where J_T is the all-ones
	# matrix). Omega has eigenvalue (sig_eps_sq + T * sig_eps_c_sq) on the
	# 1-dimensional subspace span(1_T) and eigenvalue sig_eps_sq (multiplicity
	# T - 1) on its orthogonal complement. Hence
	#   Omega^(-1/2) =
	#     (1 / sqrt(sig_eps_sq + T * sig_eps_c_sq)) * (1_T 1_T^T / T)
	#     + (1 / sqrt(sig_eps_sq)) * (I_T - 1_T 1_T^T / T).
	# The closed form is exact (up to floating point), reduces the per-fit
	# cost from an `expm::sqrtm()` Schur-decomposition + solve() to a few
	# scalar sqrt()s plus a constant-matrix subtract, and drops the runtime
	# dependency on `expm` (now Suggests-only for the equivalence test).
	# Edge case `sig_eps_c_sq = 0`: the second eigenvalue collapses to
	# sig_eps_sq, the two scaled projectors recombine into
	# (1 / sqrt(sig_eps_sq)) * I_T, and the formula evaluates correctly.
	J_over_T <- matrix(1 / T, nrow = T, ncol = T)
	Omega_sqrt_inv <- (1 / sqrt(sig_eps_sq)) *
		(diag(T) - J_over_T) +
		(1 / sqrt(sig_eps_sq + T * sig_eps_c_sq)) * J_over_T

	if (verbose) {
		message("Time to get sqrt inverse matrix:")
		message(Sys.time() - t0)
	}

	y_gls <- kronecker(diag(N), sqrt(sig_eps_sq) * Omega_sqrt_inv) %*% y
	X_gls <- kronecker(diag(N), sqrt(sig_eps_sq) * Omega_sqrt_inv) %*% X_mod

	stopifnot(ncol(X_gls) == p)

	list(
		y_gls = y_gls,
		X_gls = X_gls,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq
	)
}

#' @title Collapse the twfeCovs design matrix
#' @description Drop treatment-effect interaction columns and the
#'   covariate-by-time / covariate-by-cohort interaction columns, then
#'   collapse the per-cohort treatment dummies into one column per
#'   cohort. Called from the orchestrator only when
#'   `is_twfe_covs = TRUE`.
#' @param X_gls Numeric matrix; GLS-transformed design.
#' @param N,T,R,d Integers; units, time periods, treated cohorts,
#'   covariates.
#' @param num_treats Integer; treatment-column count in the pre-collapse
#'   design.
#' @param first_inds Integer vector; first-index-within-cohort offsets.
#' @return List with `X_collapsed` (post-collapse design) and `p_short`
#'   (column count `= R + T - 1 + d + R`).
#' @keywords internal
#' @noRd
.collapse_design_for_twfe_covs <- function(
	X_gls,
	N,
	T,
	R,
	d,
	num_treats,
	first_inds
) {
	stopifnot(nrow(X_gls) == N * T)

	X <- X_gls[, 1:(R + T - 1 + d * (1 + R + T - 1) + num_treats)]

	first_treat_ind <- R + T - 1 + d * (1 + R + T - 1)
	treat_inds <- first_treat_ind:(first_treat_ind + num_treats - 1)
	stopifnot(length(treat_inds) == num_treats)

	X <- X[, c(1:(R + T - 1 + d), treat_inds)]
	stopifnot(nrow(X) == N * T)

	treat_inds_mat <- matrix(as.numeric(NA), nrow = N * T, ncol = R)
	for (r in 1:R) {
		inds_r <- .cohort_block_inds(r, R, first_inds, num_treats)
		cols_r <- R + T - 1 + d + inds_r
		treat_inds_mat[, r] <- rowSums(X[, cols_r, drop = FALSE])
	}
	stopifnot(all(!is.na(treat_inds_mat)))

	X_collapsed <- cbind(X[, 1:(R + T - 1 + d)], treat_inds_mat)

	p_short <- R + T - 1 + d + R
	stopifnot(ncol(X_collapsed) == p_short)

	list(X_collapsed = X_collapsed, p_short = p_short)
}

#' @title Append ridge-augmentation rows to the scaled design
#' @description When `add_ridge = TRUE`, build the augmentation matrix
#'   (`D_inverse` from the fusion transform for FETWFE; `diag(p)` for
#'   ETWFE / BETWFE / twfeCovs), then append
#'   `sqrt(lambda_ridge) * mat` rows to `X_scaled` and zeros to `y_gls`.
#'   When `add_ridge = FALSE`, returns inputs unchanged with
#'   `lambda_ridge = NA`. The long arg list reflects two semantic
#'   groups: `(N, T, p, sig_eps_sq, sig_eps_c_sq)` feed the
#'   `lambda_ridge` formula; `(is_fetwfe, first_inds, T, R, d,
#'   num_treats)` feed the `D_inverse` construction.
#' @param X_scaled Numeric matrix; the scaled, GLS-transformed design.
#' @param y_gls Numeric vector; the GLS-transformed response.
#' @param p Integer; column count of `X_scaled`.
#' @param add_ridge Logical.
#' @param is_fetwfe Logical; selects fusion-inverse vs identity.
#' @param first_inds,T,R,d,num_treats Args for
#'   `genFullInvFusionTransformMat()` (only used when
#'   `is_fetwfe = TRUE`).
#' @param sig_eps_sq,sig_eps_c_sq,N Numeric / integer; inputs to the
#'   `lambda_ridge = 1e-5 * (sig_eps_sq + sig_eps_c_sq) * sqrt(p/(N*T))`
#'   formula.
#' @return List with the (possibly-augmented) `X_scaled`, `y_final`,
#'   and `lambda_ridge`.
#' @keywords internal
#' @noRd
.append_ridge_rows <- function(
	X_scaled,
	y_gls,
	p,
	add_ridge,
	is_fetwfe,
	first_inds,
	T,
	R,
	d,
	num_treats,
	sig_eps_sq,
	sig_eps_c_sq,
	N
) {
	if (!add_ridge) {
		return(list(
			X_scaled = X_scaled,
			y_final = y_gls,
			lambda_ridge = as.numeric(NA)
		))
	}

	mat_to_multiply <- diag(p)

	if (is_fetwfe) {
		D_inverse <- genFullInvFusionTransformMat(
			first_inds = first_inds,
			T = T,
			R = R,
			d = d,
			num_treats = num_treats
		)

		stopifnot(ncol(D_inverse) == p)
		stopifnot(nrow(D_inverse) == p)

		mat_to_multiply <- D_inverse
	}

	lambda_ridge <- 0.00001 *
		(sig_eps_sq + sig_eps_c_sq) *
		sqrt(p / (N * T))

	X_scaled <- rbind(
		X_scaled,
		sqrt(lambda_ridge) * mat_to_multiply
	)
	y_final <- c(y_gls, rep(0, p))

	stopifnot(length(y_final) == N * T + p)
	stopifnot(nrow(X_scaled) == N * T + p)

	list(
		X_scaled = X_scaled,
		y_final = y_final,
		lambda_ridge = lambda_ridge
	)
}

#' @title Compute in-sample (and optionally independent) cohort probabilities
#' @description Convert raw cohort counts into two normalizations:
#'   within-treated (`cohort_probs`, sums to 1) and overall
#'   (`cohort_probs_overall`, equals `count_r / N`). If
#'   `indep_count_data_available = TRUE`, do the same for `indep_counts`;
#'   otherwise the indep-side outputs are NA.
#' @param in_sample_counts Integer vector of length `R + 1`; counts of
#'   never-treated + treated cohorts in the working panel.
#' @param indep_counts Integer vector of length `R + 1` or NA;
#'   independent-sample counts.
#' @param N Integer; total units.
#' @param R Integer; treated cohorts.
#' @param indep_count_data_available Logical.
#' @return List with `cohort_probs`, `cohort_probs_overall`,
#'   `indep_cohort_probs`, `indep_cohort_probs_overall`.
#' @keywords internal
#' @noRd
.compute_cohort_probs <- function(
	in_sample_counts,
	indep_counts,
	N,
	R,
	indep_count_data_available
) {
	cohort_probs <- in_sample_counts[2:(R + 1)] /
		sum(in_sample_counts[2:(R + 1)])

	stopifnot(all(!is.na(cohort_probs)))
	stopifnot(all(cohort_probs >= 0))
	stopifnot(all(cohort_probs <= 1))
	stopifnot(length(cohort_probs) == R)
	stopifnot(abs(sum(cohort_probs) - 1) < 1e-6)

	cohort_probs_overall <- in_sample_counts[2:(R + 1)] / N

	stopifnot(
		abs(1 - sum(cohort_probs_overall) - in_sample_counts[1] / N) < 1e-6
	)

	if (indep_count_data_available) {
		indep_cohort_probs <- indep_counts[2:(R + 1)] /
			sum(indep_counts[2:(R + 1)])

		stopifnot(all(!is.na(indep_cohort_probs)))
		stopifnot(all(indep_cohort_probs >= 0))
		stopifnot(all(indep_cohort_probs <= 1))
		stopifnot(length(indep_cohort_probs) == R)
		stopifnot(abs(sum(indep_cohort_probs) - 1) < 1e-6)

		indep_cohort_probs_overall <- indep_counts[2:(R + 1)] / N

		stopifnot(
			abs(
				1 -
					sum(indep_cohort_probs_overall) -
					indep_counts[1] / N
			) <
				1e-6
		)
	} else {
		indep_cohort_probs <- NA
		indep_cohort_probs_overall <- NA
	}

	list(
		cohort_probs = cohort_probs,
		cohort_probs_overall = cohort_probs_overall,
		indep_cohort_probs = indep_cohort_probs,
		indep_cohort_probs_overall = indep_cohort_probs_overall
	)
}


# getGramInv
#' @title Compute Inverse of Gram Matrix for Selected Features
#' @description Calculates the inverse of the Gram matrix formed by the selected
#'   features from the (potentially transformed) design matrix `X_final`. This
#'   is a key component in calculating standard errors for the estimated
#'   coefficients.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param X_final Numeric matrix; the final design matrix (e.g., after GLS and
#'   fusion transformations). Dimensions: `N*T` x `p_model`, where `p_model` is
#'   the total number of columns in this matrix.
#' @param sel_feat_inds Integer vector; indices of the features selected by the
#'   penalized regression, corresponding to columns in `X_final`.
#' @param treat_inds Integer vector; original indices (before selection) of the
#'   base treatment effect parameters. Used to subset `sel_feat_inds` to get
#'   only selected *treatment* features for the final Gram matrix.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param sel_treat_inds_shifted Integer vector; indices of the selected
#'   treatment effects within the `num_treats` block, shifted to start from 1.
#'   Used to check dimensions.
#' @param calc_ses Logical; if `FALSE`, the function may return `NA` for
#'   `gram_inv`.
#' @return A list containing:
#'   \item{gram_inv}{The inverse of the Gram matrix corresponding to the
#'     *selected treatment effect features*. Returns `NA` if `calc_ses` is
#'     `FALSE` or if the Gram matrix is found to be singular.}
#'   \item{calc_ses}{Logical, potentially modified to `FALSE` if the Gram
#'     matrix is singular.}
#' @details
#'   1. Subsets `X_final` to include only columns specified by `sel_feat_inds`.
#'   2. Centers these selected columns.
#'   3. Computes the Gram matrix: `(1/(N*T)) * t(X_sel_centered) %*% X_sel_centered`.
#'   4. Checks if the minimum eigenvalue of the Gram matrix is too small (close to zero).
#'      If so, it issues a warning and sets `calc_ses` to `FALSE`, returning `NA`
#'      for `gram_inv`.
#'   5. Otherwise, it computes the inverse of the Gram matrix.
#'   6. It then subsets this inverse Gram matrix to retain only the rows/columns
#'      that correspond to the *selected treatment effects* (identified via
#'      `sel_feat_inds` and `treat_inds`).
#' @keywords internal
#' @noRd
getGramInv <- function(
	N,
	T,
	X_final,
	treat_inds,
	num_treats,
	calc_ses,
	sel_feat_inds = NA,
	sel_treat_inds_shifted = NA
) {
	stopifnot(nrow(X_final) == N * T)
	if (any(!is.na(sel_feat_inds))) {
		X_sel <- X_final[, sel_feat_inds, drop = FALSE]
		p_sel <- length(sel_feat_inds)
	} else {
		X_sel <- X_final
		p_sel <- ncol(X_final)
	}

	stopifnot(length(treat_inds) == num_treats)

	# Centering X_sel even when estimating via OLS is harmless
	X_sel_centered <- scale(X_sel, center = TRUE, scale = FALSE)

	gram <- 1 / (N * T) * (t(X_sel_centered) %*% X_sel_centered)

	stopifnot(nrow(gram) == p_sel)
	stopifnot(ncol(gram) == p_sel)

	min_gram_eigen <- min(
		eigen(gram, symmetric = TRUE, only.values = TRUE)$values
	)

	if (min_gram_eigen < 10^(-16)) {
		warning(
			"Gram matrix corresponding to selected features is not invertible. Assumptions needed for inference are not satisfied. Standard errors will not be calculated."
		)
		return(list(gram_inv = NA, calc_ses = FALSE))
	}

	gram_inv <- solve(gram)

	if (any(!is.na(sel_feat_inds))) {
		# Get only the parts of gram_inv that have to do with treatment effects
		sel_treat_inds <- sel_feat_inds %in% treat_inds
		stopifnot(sum(sel_treat_inds) <= length(sel_feat_inds))
		stopifnot(length(sel_treat_inds) == length(sel_feat_inds))
	} else {
		sel_treat_inds <- rep(FALSE, p_sel)
		sel_treat_inds[treat_inds] <- TRUE
	}

	stopifnot(is.logical(sel_treat_inds))
	stopifnot(sum(sel_treat_inds) <= num_treats)
	stopifnot(length(sel_treat_inds) == p_sel)
	stopifnot(nrow(gram_inv) == ncol(gram_inv))

	stopifnot(all(!is.na(gram_inv)))

	gram_inv <- gram_inv[sel_treat_inds, sel_treat_inds]

	stopifnot(nrow(gram_inv) <= num_treats)
	stopifnot(nrow(gram_inv) == ncol(gram_inv))

	if (any(!is.na(sel_treat_inds_shifted))) {
		stopifnot(nrow(gram_inv) == length(sel_treat_inds_shifted))
		if (any(!is.na(sel_feat_inds))) {
			stopifnot(nrow(gram_inv) <= length(sel_feat_inds))
		}
	}

	return(list(gram_inv = gram_inv, calc_ses = calc_ses))
}


# estOmegaSqrtInv
#' @title Estimate Noise Variance Components by REML
#' @description Estimates the idiosyncratic error variance (`sig_eps_sq`) and
#'   the unit-level random effect variance (`sig_eps_c_sq`) via REML on the
#'   linear mixed-effects model `y ~ X + (1 | unit)`, using `lme4::lmer`.
#'   This matches the random-effects covariance structure (`Omega = sig_eps_sq * I_T
#'   + sig_eps_c_sq * 1_T 1_T'`) the package's GLS framework assumes. Called
#'   when these variances are not provided by the user.
#' @param y Numeric vector; the observed response variable, length `N*T`,
#'   in canonical `(unit, time)` order (rows `1..T` are unit 1, etc.).
#' @param X_ints Numeric matrix; the design matrix including all fixed effects,
#'   covariates, treatment dummies, and interactions. `N*T` rows.
#' @param N Integer; the total number of unique units.
#' @param T Integer; the total number of time periods.
#' @param p Integer; the number of columns in `X_ints`.
#' @return A list containing two named elements:
#'   \item{sig_eps_sq}{REML estimate of the idiosyncratic error variance.}
#'   \item{sig_eps_c_sq}{REML estimate of the unit-level random-effect variance.}
#' @details The function builds an `lme4`-friendly data frame with the
#'   response, a unit-membership factor, and `X_ints` columns, then calls
#'   `lme4::lmer(..., REML = TRUE)` to estimate the variance components via
#'   restricted maximum likelihood. Requires `lme4` to be installed
#'   (`Suggests:` dependency, gated by `requireNamespace`); errors with a
#'   clear message if not. lme4's "different scales" warning and
#'   "rank deficient" message — expected on FETWFE-scale designs where
#'   cohort dummies are 0/1 and covariates are arbitrary — are suppressed.
#'
#'   Prior versions of this package implemented the within-estimator
#'   procedure of Pesaran (2015, Section 26.5.1). That implementation
#'   contained two coupled bugs that caused `sig_eps_c_sq` to be returned
#'   as numerically zero, degenerating the downstream GLS step. The REML
#'   implementation here fixes both the bug and a methodological mismatch:
#'   the within-estimator drops time-invariant columns from `X_ints` (cohort
#'   dummies, time-invariant covariates, cohort-by-covariate interactions),
#'   which biases `sig_eps_c_sq` upward when the true coefficients on those
#'   columns are nonzero. REML handles them natively.
#' @references Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015).
#'   Fitting Linear Mixed-Effects Models Using lme4. *Journal of Statistical
#'   Software*, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.
#'
#'   Patterson, H. D., & Thompson, R. (1971). Recovery of inter-block
#'   information when block sizes are unequal. *Biometrika*, 58(3),
#'   545-554.
#'
#'   Pinheiro, J. C., & Bates, D. M. (2000). *Mixed-Effects Models in S
#'   and S-PLUS*. Springer.
#' @keywords internal
#' @noRd
estOmegaSqrtInv <- function(y, X_ints, N, T, p) {
	if (N * (T - 1) - p <= 0) {
		stop("Not enough units available to estimate the noise variance.")
	}
	stopifnot(N > 1)

	if (!requireNamespace("lme4", quietly = TRUE)) {
		stop(
			"Variance-component estimation requires the 'lme4' package. ",
			"Install it via install.packages('lme4'), or supply ",
			"`sig_eps_sq` and `sig_eps_c_sq` directly to the estimator.",
			call. = FALSE
		)
	}

	# The package's prepXints() places rows in canonical (unit, time) order,
	# so unit IDs are rep(1..N, each = T).
	unit_id <- factor(rep(seq_len(N), each = T))

	# Build a data frame for lme4. Rename X_ints columns to generic names so
	# the formula parser doesn't choke on special characters (`:`, `*`, `(`)
	# in the original interaction column names.
	X_named <- X_ints
	colnames(X_named) <- paste0("Xcol_", seq_len(p))
	df <- data.frame(
		y_resp = y,
		unit = unit_id,
		X_named,
		check.names = FALSE
	)
	rhs <- paste0("Xcol_", seq_len(p), collapse = " + ")
	fmla <- stats::as.formula(paste0("y_resp ~ ", rhs, " + (1 | unit)"))

	# Fit y ~ X + (1 | unit) by REML. lme4 emits a "different scales"
	# warning and a "rank deficient" message on most FETWFE-scale designs;
	# both are informational (cohort dummies are 0/1, covariates are
	# arbitrary; aliased columns are dropped automatically). Suppress both.
	fit <- tryCatch(
		suppressWarnings(suppressMessages(
			lme4::lmer(fmla, data = df, REML = TRUE)
		)),
		error = function(e) {
			stop(
				"lme4::lmer failed to fit the variance-component model: ",
				conditionMessage(e),
				".\nConsider supplying `sig_eps_sq` and `sig_eps_c_sq` ",
				"directly to the estimator.",
				call. = FALSE
			)
		}
	)

	# Extract variance components. VarCorr(fit)$unit is the 1x1 covariance
	# matrix of the random intercept; attr(., "sc") is the residual std dev.
	vc <- lme4::VarCorr(fit)
	sigma_c_sq_hat <- as.numeric(vc$unit[1, 1])
	sigma_hat_sq <- as.numeric(attr(vc, "sc"))^2

	return(list(
		sig_eps_sq = sigma_hat_sq,
		sig_eps_c_sq = sigma_c_sq_hat
	))
}

#' @title Validate Inputs for the ETWFE Core Estimator
#'
#' @description
#' Performs a complete set of argument, dimension, and consistency checks
#' for \code{\link{etwfe_core}}.
#' The function halts with informative error messages whenever a necessary
#' condition is violated (e.g., negative variances, inconsistent cohort
#' counts, or ill-posed significance levels).
#' When all checks pass it returns a compact list of derived quantities
#' used downstream by the core estimator.
#'
#' @param in_sample_counts Integer named vector. Length `R+1`.
#'   The first element must be the number of never-treated units; the
#'   remaining `R` elements give the number of units in each treated
#'   cohort.  Names must be unique and correspond to cohort identifiers.
#' @param N Integer. Total number of unique units in the (filtered) data.
#' @param T Integer. Total number of time periods.
#' @param sig_eps_sq Numeric scalar or \code{NA}.  Prespecified observation-level
#'   variance component (see Section 2 of Faletto 2025).  If \code{NA},
#'   it will later be estimated.
#' @param sig_eps_c_sq Numeric scalar or \code{NA}.  Prespecified unit-level
#'   random-effect variance component.  If \code{NA}, it will later be estimated.
#' @param indep_counts Optional integer vector.  Cohort counts from an
#'   independent sample used to obtain asymptotically exact SEs for the ATT.
#'   Must have the same length and ordering as \code{in_sample_counts}.
#'   Default is \code{NA}.
#' @param verbose Logical.  Forwarded verbosity flag.  Default \code{FALSE}.
#' @param alpha Numeric in \eqn{(0,1)}.  Significance level requested for
#'   confidence intervals.  Default \code{0.05}.
#' @param add_ridge Logical.  Whether a small ridge penalty will later be
#'   applied inside \code{etwfe_core}.  Used only to customize warning
#'   messages.  Default \code{FALSE}.
#'
#' @details
#' Key validations include:
#' \itemize{
#'   \item Ensuring the total of \code{in_sample_counts} equals \code{N}.
#'   \item Checking that at least one never-treated unit is present.
#'   \item Verifying that cohort names are unique and counts non-negative.
#'   \item Confirming that \eqn{1 \le R \le T-1}.
#'   \item Basic sanity checks for all numeric scalars (\code{alpha} inside \eqn{(0,1)}, non-negative variances, etc.).
#'   \item All structural requirements on \code{indep_counts} when supplied.
#' }
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{R}}{Integer. The number of treated cohorts `=length(in_sample_counts)-1`).}
#'     \item{\code{c_names}}{Character vector of cohort names (length `R`).}
#'     \item{\code{indep_count_data_available}}{Logical. \code{TRUE} if
#'       valid \code{indep_counts} were supplied, \code{FALSE} otherwise.}
#'   }
#' @keywords internal
#' @noRd
check_etwfe_core_inputs <- function(
	in_sample_counts,
	N,
	T,
	sig_eps_sq,
	sig_eps_c_sq,
	indep_counts,
	verbose,
	alpha,
	add_ridge
) {
	R <- length(in_sample_counts) - 1

	c_names <- names(in_sample_counts)[2:(R + 1)]

	stopifnot(N >= 2) # bare minimum, 2 units at 2 times

	stopifnot(T >= 2) # bare minimum, 2 units at 2 times

	if (any(!is.na(sig_eps_sq))) {
		stopifnot(is.numeric(sig_eps_sq) | is.integer(sig_eps_sq))
		stopifnot(length(sig_eps_sq) == 1)
		stopifnot(sig_eps_sq >= 0)
	}

	if (any(!is.na(sig_eps_c_sq))) {
		stopifnot(is.numeric(sig_eps_c_sq) | is.integer(sig_eps_c_sq))
		stopifnot(length(sig_eps_c_sq) == 1)
		stopifnot(sig_eps_c_sq >= 0)
	}

	stopifnot(sum(in_sample_counts) == N)
	stopifnot(all(in_sample_counts >= 0))
	if (in_sample_counts[1] == 0) {
		stop(
			"No never-treated units detected in data to fit model; estimating treatment effects is not possible"
		)
	}
	if (length(names(in_sample_counts)) != length(in_sample_counts)) {
		stop(
			"in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)"
		)
	}

	if (
		length(names(in_sample_counts)) !=
			length(unique(names(in_sample_counts)))
	) {
		stop(
			"in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)"
		)
	}

	stopifnot(R >= 1)
	stopifnot(R <= T - 1)

	indep_count_data_available <- FALSE
	if (any(!is.na(indep_counts))) {
		if (sum(indep_counts) != N) {
			stop(
				"Number of units in independent cohort count data does not equal number of units in data to be used to fit model."
			)
		}
		if (length(indep_counts) != length(in_sample_counts)) {
			stop(
				"Number of counts in independent counts does not match number of cohorts in data to be used to fit model."
			)
		}
		if (any(indep_counts <= 0)) {
			stop(
				"At least one cohort in the independent count data has 0 members"
			)
		}
		indep_count_data_available <- TRUE
	}

	stopifnot(is.logical(verbose))
	stopifnot(length(verbose) == 1)

	stopifnot(is.numeric(alpha))
	stopifnot(length(alpha) == 1)
	stopifnot(alpha > 0)
	stopifnot(alpha < 1)

	stopifnot(is.logical(add_ridge))
	stopifnot(length(add_ridge) == 1)

	return(list(
		R = R,
		c_names = c_names,
		indep_count_data_available = indep_count_data_available
	))
}

#' Build the **entire** inverse-fusion matrix \(D_N^{-1}\)
#'
#' Constructs the \eqn{p\times p} block–diagonal matrix that sends the sparse
#' coefficient vector \(\theta\) used in the bridge loss back to the original
#' coefficient scale \(\beta\), and is also used to multiply \(Z\) on the right
#' to transform the features into a transformed feature space:
#' \deqn{\beta \;=\;D_N^{-1}\,\theta.}
#' The block layout follows the factorisation proved in the
#' paper – repeated here using the helper generators already available in the
#' package.
#'
#' \preformatted{
#'       D_N^{-1} = diag(
#'         (D^{(1)}(R))^{-1},                          # 1. cohort FEs
#'         (D^{(1)}(T-1))^{-1},                        # 2. time  FEs
#'         I_d,                                        # 3. X main effects
#'         I_d ⊗ (D^{(1)}(R))^{-1},                   # 4. cohort × X
#'         I_d ⊗ (D^{(1)}(T-1))^{-1},                 # 5. time   × X
#'         (D^{(2)}(𝓡))^{-1},                         # 6. base τ_{r,t}
#'         I_d ⊗ (D^{(2)}(𝓡))^{-1} )                  # 7. τ_{r,t} × X
#' }
#'
#' @section Block dimensions:
#' \itemize{
#'   \item Cohort FEs: \(R\times R\)
#'   \item Time-period FEs: \((T-1)\times(T-1)\)
#'   \item Identitites: \(d\), \(dR\), \(d(T-1)\)
#'   \item Treatment blocks: \( \mathfrak W \times \mathfrak W\) with
#'         \(\mathfrak W = \texttt{num_treats}\)
#' }
#' All non–identity blocks contain only 0/1 entries, so the determinant of
#' the whole matrix is 1 (volume–preserving transform).
#'
#' @param first_inds Integer vector (length \code{R}).
#'   `first_inds[r]` is the 1-based column index of the first base
#'   treatment-effect parameter \(\tau_{r,0}\) for cohort \code{r}.
#' @param T Integer. Total number of time periods \(\ge 3\).
#' @param R Integer. Number of treated cohorts (\(\ge 1\)).
#'   The function stops if you accidentally pass \code{R = 0}.
#' @param d Integer. Number of time-invariant covariates (can be 0).
#' @param num_treats Integer.  Total number of base treatment-effect
#'   coefficients \(\mathfrak W = T R - R(R+1)/2\).
#'
#' @return A dense base-R matrix of size
#'   \eqn{p \times p} with
#'   \eqn{p = R + (T-1) + d + dR + d(T-1) + \mathfrak W + d\mathfrak W}.
#'
#' @examples
#' R  <- 3; T <- 6; d <- 2
#' nt <- getNumTreats(R, T)
#' Dinv <- genFullInvFusionTransformMat(getFirstInds(R,T), T, R, d, nt)
#' dim(Dinv)   # should be p × p
#'
#' @keywords internal
#' @noRd
genFullInvFusionTransformMat <- function(first_inds, T, R, d, num_treats) {
	##———— Safety checks ————————————————————————————————————————————
	stopifnot(is.numeric(R), length(R) == 1L, R >= 1L)
	stopifnot(is.numeric(T), length(T) == 1L, T >= 3L, R <= T - 1)
	stopifnot(is.numeric(d), length(d) == 1L, d >= 0L)
	stopifnot(length(first_inds) == R)

	##———— 1. Cohort fixed-effects block:  (D^{(1)}(R))^{-1} ————————
	block1 <- genBackwardsInvFusionTransformMat(R)

	##———— 2. Time fixed-effects block:    (D^{(1)}(T-1))^{-1} ———————
	block2 <- genBackwardsInvFusionTransformMat(T - 1)

	##———— 3. Covariate main effects:      I_d ————————————————
	block3 <- if (d > 0) diag(d) else NULL

	##———— 4. Cohort × X interactions:     I_d ⊗ (D^{(1)}(R))^{-1} ——
	block4 <- if (d > 0) {
		kronecker(diag(d), genBackwardsInvFusionTransformMat(R))
	} else {
		NULL
	}

	##———— 5. Time × X interactions:       I_d ⊗ (D^{(1)}(T-1))^{-1} —
	block5 <- if (d > 0) {
		kronecker(diag(d), genBackwardsInvFusionTransformMat(T - 1))
	} else {
		NULL
	}

	##———— 6. Base treatment effects:      (D^{(2)}(𝓡))^{-1} ————————
	block6 <- genInvTwoWayFusionTransformMat(
		n_vars = num_treats,
		first_inds = first_inds,
		R = R
	)

	##———— 7. Treatment × X interactions:  I_d ⊗ (D^{(2)}(𝓡))^{-1} ——
	block7 <- if (d > 0) {
		kronecker(
			diag(d),
			genInvTwoWayFusionTransformMat(
				n_vars = num_treats,
				first_inds = first_inds,
				R = R
			)
		)
	} else {
		NULL
	}

	## Gather present blocks in the same order as the theoretical expression
	blocks <- list(block1, block2, block3, block4, block5, block6, block7)
	blocks <- Filter(Negate(is.null), blocks) # drop NULLs for d = 0

	##———— Assemble block-diagonal matrix ————————————————
	## Matrix::bdiag() returns a sparse dgCMatrix.  We convert to base-R matrix
	## here because downstream (ridge-row augmentation) works with dense objects
	## via rbind().
	full_D_inv <- as.matrix(Matrix::bdiag(blocks))

	##———— Dimension cross-check ————————————————————————————
	p <- getP(R = R, T = T, d = d, num_treats = num_treats)
	stopifnot(nrow(full_D_inv) == p, ncol(full_D_inv) == p)

	return(full_D_inv)
}

#' @title Build the "no features selected" early-exit return list shared by
#'   `fetwfe_core()` and `betwfe_core()`
#'
#' @description
#' Internal helper that consolidates four near-duplicate early-exit blocks
#' that previously lived inline in `R/betwfe_core.R` and `R/fetwfe_core.R`
#' (BETWFE intercept-only, BETWFE no-treatment-selected, FETWFE
#' intercept-only, FETWFE no-treatment-selected). All four blocks construct
#' the same 33-/34-field return list with all-zero CATTs / ATTs and the same
#' `data.frame` of CATT placeholders; they differ only in the `beta_hat`
#' value, whether `theta_hat` is in the return, and the verbose-message
#' text. Caller computes the final `beta_hat` (including any
#' `add_ridge` adjustment and, for FETWFE, the
#' `untransformCoefImproved()` transformation back from `theta`) and
#' passes it in; helper builds the return list with the contractually-fixed
#' field ordering. Resolves drift surface that produced 9 stale `q < 1`
#' references in v1.9.5 (#56).
#'
#' @param message_text Character scalar; emitted via `message()` when
#'   `verbose` is `TRUE`. Per-block string (see callers).
#' @param verbose Logical; whether to emit `message_text`.
#' @param R,N,T,d,p Integers; problem dimensions.
#' @param c_names Character vector of length `R`; cohort labels for the
#'   `catt_df_to_ret` rows.
#' @param q Numeric; bridge regression exponent. Drives the
#'   `ret_se <- if (q < 1) 0 else NA` SE-shape gate.
#' @param beta_hat Numeric vector of length `p`; the final slope estimates
#'   (caller-computed). For BETWFE intercept-only this is
#'   `beta_hat[2:(p+1)]` (all zero); for BETWFE no-treatment this is
#'   `beta_hat_slopes` after caller-applied `add_ridge` adjustment; for
#'   FETWFE intercept-only this is `rep(0, p)`; for FETWFE no-treatment
#'   this is `untransformCoefImproved(theta_hat_slopes)` after
#'   caller-applied `add_ridge` adjustment.
#' @param treat_inds,treat_int_inds Integer vectors; treatment-effect column
#'   indices.
#' @param cohort_probs,indep_cohort_probs Numeric vectors; cohort-probability
#'   weights computed by `prep_for_etwfe_regression()`.
#' @param cohort_probs_overall,indep_cohort_probs_overall Numeric scalars;
#'   overall cohort-probability weights.
#' @param sig_eps_sq,sig_eps_c_sq Numeric scalars; variance components.
#' @param lambda.max,lambda.min,lambda_star Numeric scalars; bridge penalty
#'   path endpoints and the BIC-selected lambda.
#' @param lambda.max_model_size,lambda.min_model_size,lambda_star_model_size
#'   Integer scalars; model sizes at each lambda.
#' @param X_ints Numeric matrix; original-basis design.
#' @param y Numeric vector; centered response.
#' @param X_final,y_final Numeric matrix/vector; GLS-whitened design and
#'   response actually fed to the bridge regression.
#' @param theta_hat Optional numeric vector of length `p + 1` (includes
#'   intercept). For FETWFE blocks only; BETWFE blocks omit (default
#'   `NULL`).
#' @param include_theta Logical; if `TRUE`, inserts `theta_hat` into the
#'   return list between `catt_df` and `beta_hat` (FETWFE field order).
#'   Default `FALSE` (BETWFE field order: no `theta_hat`).
#'
#' @return A 33-field (BETWFE) or 34-field (FETWFE) list with the same
#'   field ordering as the pre-refactor inline blocks. The ordering is
#'   load-bearing: downstream code in `betwfe()` / `fetwfe()` and the S3
#'   class machinery in `R/*_class.R` indexes the result by name, but the
#'   snapshot tests in `tests/testthat/_snaps/print-method-snapshot.md` and
#'   the constructors in `R/class_helpers.R` are sensitive to the *names()*
#'   ordering.
#'
#' @keywords internal
#' @noRd
.build_selected_out_result <- function(
	message_text,
	verbose,
	R,
	c_names,
	q,
	beta_hat,
	treat_inds,
	treat_int_inds,
	cohort_probs,
	indep_cohort_probs,
	cohort_probs_overall,
	indep_cohort_probs_overall,
	sig_eps_sq,
	sig_eps_c_sq,
	lambda.max,
	lambda.max_model_size,
	lambda.min,
	lambda.min_model_size,
	lambda_star,
	lambda_star_model_size,
	X_ints,
	y,
	X_final,
	y_final,
	N,
	T,
	d,
	p,
	theta_hat = NULL,
	include_theta = FALSE
) {
	if (verbose) {
		message(message_text)
	}

	if (q < 1) {
		ret_se <- 0
	} else {
		ret_se <- NA
	}

	catt_df_to_ret <- data.frame(
		Cohort = c_names,
		`Estimated TE` = rep(0, R),
		SE = rep(ret_se, R),
		ConfIntLow = rep(ret_se, R),
		ConfIntHigh = rep(ret_se, R),
		P_value = rep(NA_real_, R),
		selected = rep(FALSE, R),
		check.names = FALSE
	)

	# Build the FULL union list with every possible field, then drop
	# the fields that don't belong to this block via name-vector
	# selection (per PR #92's `.summary_estimator_output()` pattern).
	# The `keep` vector encodes the contractual field ordering:
	# `theta_hat` slots between `catt_df` and `beta_hat` for FETWFE
	# blocks (3 + 4); BETWFE blocks (1 + 2) omit it.
	out <- list(
		in_sample_att_hat = 0,
		in_sample_att_se = ret_se,
		in_sample_att_se_no_prob = ret_se,
		indep_att_hat = 0,
		indep_att_se = ret_se,
		catt_hats = setNames(rep(0, R), c_names),
		catt_ses = setNames(rep(ret_se, R), c_names),
		catt_df = catt_df_to_ret,
		theta_hat = theta_hat,
		beta_hat = beta_hat,
		treat_inds = treat_inds,
		treat_int_inds = treat_int_inds,
		cohort_probs = cohort_probs,
		indep_cohort_probs = indep_cohort_probs,
		cohort_probs_overall = cohort_probs_overall,
		indep_cohort_probs_overall = indep_cohort_probs_overall,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		lambda.max = lambda.max,
		lambda.max_model_size = lambda.max_model_size,
		lambda.min = lambda.min,
		lambda.min_model_size = lambda.min_model_size,
		lambda_star = lambda_star,
		lambda_star_model_size = lambda_star_model_size,
		X_ints = X_ints,
		y = y,
		X_final = X_final,
		y_final = y_final,
		N = N,
		T = T,
		R = R,
		d = d,
		p = p,
		calc_ses = q < 1
	)

	keep <- c(
		"in_sample_att_hat",
		"in_sample_att_se",
		"in_sample_att_se_no_prob",
		"indep_att_hat",
		"indep_att_se",
		"catt_hats",
		"catt_ses",
		"catt_df",
		if (include_theta) "theta_hat",
		"beta_hat",
		"treat_inds",
		"treat_int_inds",
		"cohort_probs",
		"indep_cohort_probs",
		"cohort_probs_overall",
		"indep_cohort_probs_overall",
		"sig_eps_sq",
		"sig_eps_c_sq",
		"lambda.max",
		"lambda.max_model_size",
		"lambda.min",
		"lambda.min_model_size",
		"lambda_star",
		"lambda_star_model_size",
		"X_ints",
		"y",
		"X_final",
		"y_final",
		"N",
		"T",
		"R",
		"d",
		"p",
		"calc_ses"
	)
	out[keep]
}
