# GLS / linear-algebra machinery for the ETWFE / bridge fitting core. Split out
# of R/core_funcs.R for issue #186 (no behavior change). Variance-component
# estimation, GLS whitening, ridge-row augmentation, Gram inverse, and
# Omega^(-1/2) routines.

#' @title Estimate variance components and apply GLS whitening
#' @description If `sig_eps_sq` or `sig_eps_c_sq` is NA, estimate both via
#'   `estOmegaSqrtInv()` (REML on `y ~ X + (1 | unit)`) -- warning when exactly
#'   one was supplied, since REML estimates them jointly so the supplied value is
#'   discarded (#266). Then build
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
		# Both noise variances are estimated jointly by REML, so supplying only
		# ONE of them silently discards it. Warn when exactly one is given (#266).
		if (xor(is.na(sig_eps_sq), is.na(sig_eps_c_sq))) {
			supplied <- if (is.na(sig_eps_sq)) {
				sprintf("sig_eps_c_sq = %g", sig_eps_c_sq)
			} else {
				sprintf("sig_eps_sq = %g", sig_eps_sq)
			}
			warning(
				"Only one of `sig_eps_sq` / `sig_eps_c_sq` was supplied (",
				supplied,
				"); the two noise variances are estimated jointly (REML), so the ",
				"supplied value is ignored and both are re-estimated. Supply both ",
				"(or neither) to control the noise variances.",
				call. = FALSE
			)
		}
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

	# Block-apply form of (I_N kron (sqrt(sig_eps_sq) * Omega_sqrt_inv)).
	# The transform we're applying is the standard GLS whitening
	# Sigma^{-1/2} = I_N kron A, where A is the T-by-T per-unit factor
	# below. y and X_mod are laid out in (unit, time) order with T rows
	# per unit (the invariant `idCohorts()` enforces and `processCovs()`
	# re-asserts), so `matrix(y, nrow = T)` puts unit k in column k.
	# The identity (I_N kron A) %*% vec_{T,N}(M) = vec_{T,N}(A %*% M)
	# then says: don't form the Kronecker product; apply A on the left
	# to each unit's T-vector independently. The pre-fix code allocates
	# a (N*T)x(N*T) Kronecker matrix on every call (800 MB at N=2000),
	# then does a single dense multiply; the post-fix code allocates
	# nothing larger than the inputs. Verified bit-exact on the test
	# fixtures (#165). See R/utility.R::idCohorts and
	# R/design_matrix.R::processCovs for the row-ordering invariant.
	stopifnot(length(y) == N * T)
	stopifnot(nrow(X_mod) == N * T)
	A <- sqrt(sig_eps_sq) * Omega_sqrt_inv
	# Preserve the prior return shapes: y_gls is an (N*T) x 1 matrix
	# (kronecker(...) %*% y produces one), X_gls is an (N*T) x p matrix.
	# Downstream `grpreg::*()` accepts y as a vector OR a 1-column matrix,
	# but keeping the matrix form here avoids invalidating reference
	# fixtures that built expectations against the explicit-Kronecker
	# return shape (e.g. tests/testthat/test-est-omega-sqrt-inv.R:205).
	y_gls <- matrix(A %*% matrix(y, nrow = T), nrow = N * T, ncol = 1L)
	X_gls <- matrix(
		A %*% matrix(X_mod, nrow = T),
		nrow = N * T,
		ncol = ncol(X_mod)
	)

	stopifnot(ncol(X_gls) == p)

	list(
		y_gls = y_gls,
		X_gls = X_gls,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq
	)
}

#' @title Append ridge-augmentation rows to the scaled design
#' @description When `add_ridge = TRUE`, build the augmentation matrix
#'   (`D_inverse` from the fusion transform for FETWFE; `diag(p)` for
#'   ETWFE / BETWFE / twfeCovs), then append
#'   `sqrt(lambda_ridge) * mat` rows to `X_scaled` and zeros to `y_gls`.
#'   When `add_ridge = FALSE`, returns inputs unchanged with
#'   `lambda_ridge = NA`. The long arg list reflects two semantic
#'   groups: `(N, T, p, sig_eps_sq, sig_eps_c_sq)` feed the
#'   `lambda_ridge` formula; `(is_fetwfe, first_inds, T, G, d,
#'   num_treats)` feed the `D_inverse` construction.
#' @param X_scaled Numeric matrix; the scaled, GLS-transformed design.
#' @param y_gls Numeric vector; the GLS-transformed response.
#' @param p Integer; column count of `X_scaled`.
#' @param add_ridge Logical.
#' @param is_fetwfe Logical; selects fusion-inverse vs identity.
#' @param first_inds,T,G,d,num_treats Args for
#'   `genFullInvFusionTransformMat()` (only used when
#'   `is_fetwfe = TRUE`).
#' @param sig_eps_sq,sig_eps_c_sq,N Numeric / integer; inputs to the
#'   `lambda_ridge = 1e-5 * (sig_eps_sq + sig_eps_c_sq) * sqrt(p/(N*T))`
#'   formula.
#' @param fusion_structure Character; one of `"cohort"` (default) or
#'   `"event_study"`. Forwarded to `genFullInvFusionTransformMat()`
#'   (only used when `add_ridge = TRUE && is_fetwfe`).
#' @param d_inv_treat Optional `num_treats x num_treats` numeric matrix; the
#'   user-supplied already-inverted treatment-effect fusion block
#'   (`solve(fusion_matrix)`, #236), forwarded to
#'   `genFullInvFusionTransformMat()` so the custom block reaches the ridge
#'   penalty rows. Default `NULL`.
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
	G,
	d,
	num_treats,
	sig_eps_sq,
	sig_eps_c_sq,
	N,
	fusion_structure = "cohort",
	d_inv_treat = NULL
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
			G = G,
			d = d,
			num_treats = num_treats,
			fusion_structure = fusion_structure,
			d_inv_treat = d_inv_treat
		)

		stopifnot(ncol(D_inverse) == p)
		stopifnot(nrow(D_inverse) == p)

		mat_to_multiply <- D_inverse
	}

	# #185 SB5: sig_eps_sq > 0 here (validated for user-supplied input per SB2;
	# REML-estimated otherwise) and sig_eps_c_sq >= 0, so `lambda_ridge` is
	# strictly positive -- the former silent no-op (lambda_ridge = 0 when both
	# variance components were 0) is unreachable. Assert it so a future change
	# to the variance handling cannot silently reintroduce a no-op ridge under
	# add_ridge = TRUE.
	lambda_ridge <- 0.00001 *
		(sig_eps_sq + sig_eps_c_sq) *
		sqrt(p / (N * T))
	stopifnot(is.finite(lambda_ridge), lambda_ridge > 0)

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

	gram_eigvals <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
	min_gram_eigen <- min(gram_eigvals)
	max_gram_eigen <- max(gram_eigvals)

	gram_singular_msg <- "Gram matrix corresponding to selected features is not invertible. Assumptions needed for inference are not satisfied. Standard errors will not be calculated."

	# Guard Gram inversion with an rcond-aware relative tolerance (#205):
	# solve() aborts on reciprocal condition number, not absolute min
	# eigenvalue, and the old absolute floor (1e-16) sat below
	# .Machine$double.eps, letting near-singular Grams reach an unguarded
	# solve(). The tryCatch below backstops any residual rcond case so the
	# fit degrades to calc_ses = FALSE rather than aborting.
	if (
		min_gram_eigen < max(dim(gram)) * .Machine$double.eps * max_gram_eigen
	) {
		warning(gram_singular_msg)
		return(list(gram_inv = NA, calc_ses = FALSE))
	}

	gram_inv <- tryCatch(solve(gram), error = function(e) NULL)
	if (is.null(gram_inv)) {
		warning(gram_singular_msg)
		return(list(gram_inv = NA, calc_ses = FALSE))
	}

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
#'   "rank deficient" message --- expected on FETWFE-scale designs where
#'   cohort dummies are 0/1 and covariates are arbitrary --- are suppressed.
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
		stop(
			"Not enough units available to estimate the noise variance (REML needs ",
			"p < N(T - 1)). For a high-dimensional (p >= NT) fit, either supply ",
			"`sig_eps_sq` / `sig_eps_c_sq`, or re-fit with `gls = FALSE` to skip ",
			"variance-component estimation entirely (the debiasedATT() ",
			"cluster-robust standard error needs no whitening).",
			call. = FALSE
		)
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
