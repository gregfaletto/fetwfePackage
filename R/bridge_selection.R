# Bridge-regression coefficient selection for fetwfe_core() (#188).
#
# getBetaBIC() and getBetaCV() select the bridge-penalty lambda by BIC and by
# cross-validation, respectively; .untransform_scaled_theta() back-transforms
# the selected scaled coefficients to the original-data scale. All three are
# called only from .dispatch_bridge_selection() in R/utility.R. Relocated
# verbatim from R/fetwfe_core.R (which now holds checkFetwfeInputs() +
# fetwfe_core()); see issue #188.

# .untransform_scaled_theta
#' @title Back-transform scaled bridge coefficients to original-data scale
#' @description Shared rescaling step used by `getBetaBIC()` and `getBetaCV()`.
#'   Given a length-`(p + 1)` coefficient vector `theta_hat_scaled` (intercept
#'   at position 1, slopes at 2..p+1) on the my_scale()-centered/scaled
#'   design, returns the same shape on the original-data scale.
#'
#'   The back-transform: `beta_j = beta_scaled_j / scale_scale_j` for
#'   `j = 1..p`, and `intercept = intercept_scaled - sum(scale_center *
#'   (beta_scaled / scale_scale))`.
#' @keywords internal
#' @noRd
.untransform_scaled_theta <- function(
	theta_hat_scaled,
	p,
	scale_center,
	scale_scale
) {
	if (length(theta_hat_scaled) != p + 1L) {
		stop(
			"theta_hat_scaled length (",
			length(theta_hat_scaled),
			") does not match p + 1 (",
			p + 1L,
			")."
		)
	}
	if (length(scale_scale) != p) {
		stop("Length of scale_scale does not match number of predictors (p).")
	}
	if (length(scale_center) != p) {
		stop("Length of scale_center does not match number of predictors (p).")
	}
	adjusted <- theta_hat_scaled
	slopes_scaled <- theta_hat_scaled[2:(p + 1L)]
	adjusted[2:(p + 1L)] <- slopes_scaled / scale_scale
	adjusted[1L] <- theta_hat_scaled[1L] -
		sum(scale_center * (slopes_scaled / scale_scale))
	adjusted
}

# getBetaBIC
#' @title Select Optimal Coefficients using BIC from gBridge Fit
#' @description From a `gBridge` fit object (which contains solutions for a
#'   path of lambda penalties), this function selects the optimal set of
#'   coefficients based on the Bayesian Information Criterion (BIC). It also
#'   returns the chosen lambda index and the size of the selected model.
#'   Coefficients are returned on their original scale.
#' @param fit A `gBridge` fit object, typically the output from `grpreg::gBridge()`.
#' @param N Integer; the total number of unique units.
#' @param T Integer; the total number of time periods.
#' @param p Integer; the total number of predictor variables (excluding intercept)
#'   in the model matrix `X_mod`.
#' @param X_mod Numeric matrix; the design matrix (potentially transformed for
#'   FETWFE, and **not** yet GLStransformed or scaled/centered by `my_scale`) that was used to generate `y`.
#'   It's used here to calculate SSE on the original scale of `y`.
#' @param y Numeric vector; the original response variable (before GLS transform and centering)
#'   used to fit the model. Length `N*T`.
#' @param scale_center Numeric vector; the centering values used to scale `X_mod`
#'   before fitting `gBridge`. Length `p`.
#' @param scale_scale Numeric vector; the scaling values used to scale `X_mod`
#'   before fitting `gBridge`. Length `p`.
#' @return A list containing:
#'   \item{theta_hat}{Numeric vector of length `p+1`. The selected coefficients
#'     (including intercept at `theta_hat[1]`) on their original data scale.}
#'   \item{lambda_star_ind}{Integer; the index of the lambda value in `fit$lambda`
#'     that resulted in the best BIC.}
#'   \item{lambda_star_model_size}{Integer; the number of non-zero coefficients
#'     (excluding intercept) in the selected model.}
#' @details For every lambda in `fit$lambda` (all evaluated together rather than
#'   in a per-lambda loop):
#'   1. It extracts the intercepts (`eta_s`) and slopes (`beta_s`) on the scaled data.
#'   2. It converts these coefficients back to the original data scale using
#'      `scale_center` and `scale_scale`.
#'   3. It calculates the mean squared error using `sse_bridge()` with the
#'      original-scale coefficients, original `y`, and `X_mod` -- a single BLAS-3
#'      matrix multiply across all lambdas.
#'   4. It computes the BIC value: `N*T*log(SSE/(N*T)) + s*log(N*T)`, where `s`
#'      is the number of non-zero coefficients (excluding the always-present
#'      intercept; a constant offset that does not affect which lambda minimizes BIC).
#'      `SSE/(N*T)` is floored at `.Machine$double.eps * var(y)` (a
#'      response-scale-relative floor) so that an interpolating lambda gives a
#'      finite BIC instead of `-Inf` (#403).
#'   The set of coefficients corresponding to the minimum BIC is chosen. If multiple
#'   lambdas yield the same minimum BIC, the one resulting in the smallest model
#'   size (fewest non-zero coefficients) is selected.
#'   The final returned `theta_hat` also has its slopes and intercept adjusted back to the original scale.
#' @keywords internal
#' @noRd
getBetaBIC <- function(fit, N, T, p, X_mod, y, scale_center, scale_scale) {
	stopifnot(length(y) == N * T)
	stopifnot(nrow(fit$beta) == p + 1)

	## --- extract coefficients on the scaled data (all lambdas) -------
	eta_s <- fit$beta[1, ] # intercepts (scaled space), one per lambda
	beta_s <- fit$beta[2:(p + 1), , drop = FALSE] # slopes (scaled), p x n_lambda

	## --- convert to original scale -----------------------------------
	beta_hat <- beta_s / scale_scale # row i divided by scale_scale[i]
	eta_hat <- eta_s - colSums(scale_center * beta_hat)

	## --- residual MSE for every lambda in a single matmul ------------
	mse_hat <- sse_bridge(
		eta_hat,
		beta_hat,
		y = y,
		X_mod = X_mod,
		N = N,
		T = T
	)

	## --- model sizes and BIC for every lambda ------------------------
	# Number of selected *features* per lambda, excluding the always-present
	# (unpenalized) intercept in row 1 (#269). This feeds the BIC penalty below;
	# dropping the constant intercept shifts every BIC by the same amount, so the
	# argmin (the selected lambda) and the smallest-model tie-break are unchanged.
	model_sizes <- as.integer(colSums(fit$beta[-1, , drop = FALSE] != 0))

	# Coerce N * T to numeric before multiplying. `N * T` is integer
	# arithmetic and overflows at .Machine$integer.max ≈ 2.15e9
	# (panels with N * T > ~2.15e9, e.g. 50,000 x 50,000). On
	# overflow, R returns NA_integer_ and downstream
	# `which(BICs == min(BICs))` is empty, tripping the
	# `stopifnot(length(lambda_star_final_ind) == 1)` below with an
	# unhelpful message. The CV path already handles this safely (clip
	# to integer.max with a warning); mirror the same as.numeric()
	# coercion here. No regression test: the bug fires only at panel
	# sizes too large to construct in CI (#178).
	nt_double <- as.numeric(N) * as.numeric(T)
	# Floor mse_hat (and warn): an interpolating lambda (mse_hat == 0, which
	# sse_bridge()'s `>= 0` allows in the p >= NT / gls = FALSE corner) would give
	# log(0) = -Inf and win the argmin unconditionally, ignoring model size. The
	# floor keeps BICs finite/comparable -- so a sufficiently sparse
	# non-interpolating lambda can now win where -Inf always picked the
	# interpolator -- but a near-zero MSE may still not be overcome by the size
	# penalty, so the warning flags that the selected lambda may over-fit. (#403)
	#
	# The floor is RELATIVE to the response's own scale. `mse_hat` is in the raw
	# (unstandardized) units of `y`, so a fixed absolute cut is not
	# scale-equivariant: rescaling `y` by c rescales every mse_hat by c^2, and an
	# absolute .Machine$double.eps threshold would fire on a merely small-scale
	# response (a change of units) while leaving the fit itself untouched --
	# warning about "interpolation" on, and re-ranking the lambdas of, a fit that
	# is nowhere near interpolating. Scaling by var(y) makes the criterion
	# "residual MSE is zero to double precision relative to the variation being
	# explained" (R^2 >= 1 - eps), which is invariant to the units of `y`.
	# `.Machine$double.xmin` keeps the floor strictly positive (hence log() finite)
	# for a constant or vanishingly-dispersed response.
	var_y <- mean((y - mean(y))^2)
	mse_floor <- max(.Machine$double.eps * var_y, .Machine$double.xmin)
	# `var_y` overflows to Inf once max|y - mean(y)| exceeds ~1.3e154. An
	# infinite floor would floor EVERY lambda to the same value, collapsing the
	# whole BIC ranking onto the smallest-model tie-break below -- strictly worse
	# than the absolute floor this replaced. Fall back to the smallest positive
	# double: it floors nothing, and log() stays finite. (#403)
	if (!is.finite(mse_floor)) {
		mse_floor <- .Machine$double.xmin
	}
	if (any(mse_hat <= mse_floor)) {
		warning(
			"getBetaBIC(): at least one lambda interpolates the response ",
			"(residual MSE is zero to double precision relative to var(y)), ",
			"which would make its BIC -Inf; flooring MSE at ",
			".Machine$double.eps * var(y). This most often arises in the ",
			"p >= NT / gls = FALSE corner; the BIC-selected lambda may over-fit.",
			call. = FALSE
		)
	}
	BICs <- nt_double *
		log(pmax(mse_hat, mse_floor)) +
		model_sizes * log(nt_double)

	lambda_star_ind <- which(BICs == min(BICs))
	if (length(lambda_star_ind) == 1) {
		lambda_star_final_ind <- lambda_star_ind
		theta_hat <- fit$beta[, lambda_star_final_ind]
	} else {
		# Choose smallest model size among models with equal BIC
		model_sizes_star <- model_sizes[lambda_star_ind]
		min_model_size_ind <- which(model_sizes_star == min(model_sizes_star))
		lambda_star_final_ind <- lambda_star_ind[min_model_size_ind][1]
		stopifnot(length(lambda_star_final_ind) == 1)
		theta_hat <- fit$beta[, lambda_star_final_ind]
	}
	stopifnot(length(lambda_star_final_ind) == 1)
	stopifnot(length(theta_hat) == p + 1)
	stopifnot(all(!is.na(theta_hat)))

	adjusted_theta_hat <- .untransform_scaled_theta(
		theta_hat_scaled = theta_hat,
		p = p,
		scale_center = scale_center,
		scale_scale = scale_scale
	)

	return(list(
		theta_hat = adjusted_theta_hat,
		lambda_star_ind = lambda_star_final_ind,
		lambda_star_model_size = model_sizes[lambda_star_final_ind]
	))
}

# getBetaCV
#' @title Select Optimal Coefficients via 10-fold CV on `cv.grpreg`
#' @description Fits the bridge-penalized regression with `grpreg::cv.grpreg`
#'   (which performs k-fold CV over the same `grpreg`-derived lambda grid that
#'   `getBetaBIC()` would walk) and returns the coefficients at the CV-selected
#'   `lambda.min`, back-transformed to the original-data scale to match
#'   `getBetaBIC()`'s return contract.
#'
#'   The CV path is the v1.13.0 default for `fetwfe()` and `betwfe()`,
#'   replacing the BIC selection that previously was the only option. Phase B
#'   simulation studies (see issue #164 and the `lambda-selection-investigation-164`
#'   planning artifacts) showed that BIC selection produced systematically
#'   biased overall-ATT estimates with CI coverage collapsing to 0.00 at
#'   N = 2000 in the paper's second-simulation regime; CV restores
#'   near-nominal coverage in every tested regime.
#'
#'   Seed handling: if `cv_seed` is `NULL`, it defaults internally to
#'   `as.integer(N * T)` so consecutive calls on the same dataset are
#'   reproducible without user intervention. The seed is set via `set.seed()`
#'   immediately before the `cv.grpreg()` call so the fold assignment is
#'   deterministic.
#' @param X_final_scaled Numeric matrix; the design matrix after fusion +
#'   GLS transformations and `my_scale()` centering/scaling. Same matrix
#'   handed to `gBridge()` in the BIC path.
#' @param N Integer; number of unique units.
#' @param T Integer; number of time periods.
#' @param p Integer; number of predictors in `X_final_scaled` (excluding intercept).
#' @param scale_center Numeric vector of length `p`; the centering values
#'   `my_scale()` produced.
#' @param scale_scale Numeric vector of length `p`; the scaling values.
#' @param y_final Numeric vector; the GLS-transformed response.
#' @param gamma Numeric scalar; the bridge-penalty exponent `q`. Passed
#'   through to `cv.grpreg(penalty = "gBridge", gamma = gamma)`.
#' @param cv_folds Integer; number of folds for `cv.grpreg`.
#' @param cv_seed Integer or NULL; if NULL, defaults to `as.integer(N * T)`.
#' @return A list with the same shape as `getBetaBIC()`:
#'   \item{theta_hat}{Numeric vector of length `p + 1` with the selected
#'     coefficients (including intercept at position 1) on the original
#'     data scale.}
#'   \item{lambda_star_ind}{Integer; the index of the selected lambda in
#'     `cv_fit$fit$lambda`.}
#'   \item{lambda_star_model_size}{Integer; the number of non-zero
#'     coefficients (excluding the intercept) at the selected lambda.}
#' @keywords internal
#' @noRd
getBetaCV <- function(
	X_final_scaled,
	y_final,
	N,
	T,
	p,
	scale_center,
	scale_scale,
	gamma,
	cv_folds,
	cv_seed
) {
	if (is.null(cv_seed)) {
		# Default seed = N * T. Guard the (implausible) case where N * T
		# exceeds .Machine$integer.max ≈ 2.15e9: `as.integer()` would
		# silently coerce to `NA_integer_` and `set.seed(NA)` would then
		# error opaquely from the user's perspective ("supplied seed is
		# not a valid integer") because they never passed a seed. Clip to
		# the integer ceiling and warn so the user knows to override.
		nt_double <- as.numeric(N) * as.numeric(T)
		if (nt_double > .Machine$integer.max) {
			warning(
				"Default cv_seed (N * T = ",
				format(nt_double, scientific = FALSE),
				") exceeds .Machine$integer.max; clipping. ",
				"Pass cv_seed = <integer> to silence."
			)
			cv_seed <- .Machine$integer.max
		} else {
			cv_seed <- as.integer(nt_double)
		}
	}
	# `set.seed()` here drives the fold assignment for cv.grpreg(). We use
	# this rather than cv.grpreg's own `seed` argument because the latter
	# changed semantics across grpreg versions; the global RNG approach is
	# version-stable. `.with_preserved_rng()` saves and restores the caller's
	# .Random.seed so the call leaves the caller's RNG state untouched
	# (issue #177 — without this, every default-path fetwfe() / betwfe()
	# call would silently mutate the user's seed, a v1.13.0 regression
	# vs the v1.12.x BIC default).
	#
	# NOTE (#403): when add_ridge = TRUE, X_final_scaled already carries the p
	# appended sqrt(lambda_ridge)-scale ridge pseudo-rows (response 0), so they
	# participate in cv.grpreg()'s fold assignment and CV error. Accepted as-is:
	# at the default lambda_ridge (~1e-5 scale) the contamination is numerically
	# negligible. Excluding them (CV on the first N*T rows, then refit at the
	# selected lambda on the augmented design) would change the selected lambda,
	# so it is a separate behavior-changing enhancement, not a behavior-neutral
	# hardening fix.
	cv_fit <- .with_preserved_rng(cv_seed, {
		grpreg::cv.grpreg(
			X = X_final_scaled,
			y = y_final,
			penalty = "gBridge",
			gamma = gamma,
			nfolds = cv_folds
		)
	})

	lam_idx <- .cv_grpreg_at_min(cv_fit)
	stopifnot(nrow(cv_fit$fit$beta) == p + 1L)

	theta_hat_full <- cv_fit$fit$beta[, lam_idx]
	stopifnot(length(theta_hat_full) == p + 1L)
	stopifnot(all(!is.na(theta_hat_full)))

	adjusted_theta_hat <- .untransform_scaled_theta(
		theta_hat_scaled = theta_hat_full,
		p = p,
		scale_center = scale_center,
		scale_scale = scale_scale
	)

	list(
		theta_hat = adjusted_theta_hat,
		lambda_star_ind = lam_idx,
		# Exclude the intercept (row 1) to report selected-feature count (#269).
		lambda_star_model_size = sum(theta_hat_full[-1] != 0),
		# Carry the cv.grpreg fit object out — fetwfe_core() / betwfe_core()
		# need access to `fit$lambda`, `fit$beta`, and the four
		# lambda.max/min diagnostics that .fit_bridge_with_lambda_path()
		# would otherwise have computed.
		fit = cv_fit$fit,
		cv_seed = cv_seed
	)
}

# .fit_q1_nuisance
#' @title Fit the high-dimensional debiased nuisance as a `q = 1` fused lasso
#' @description The high-dimensional (`p >= NT`) debiased construction --- shared
#'   by `debiasedATT()` and the `simultaneousCIs(method = "bootstrap")` band
#'   center --- needs a nuisance `theta_hat` whose `l1` *rate* controls the
#'   Neyman-orthogonal remainder (the high-dimensional FETWFE theory). The paper
#'   establishes that rate for a **`q = 1` fused lasso** (the convex case,
#'   standard under the restricted-eigenvalue condition); the `q < 1` FETWFE
#'   bridge is super-efficient / non-uniform --- the wrong nuisance for a
#'   uniformly-valid CI (#303). This fits that nuisance on the SAME whitened,
#'   fusion-transformed design the bridge used, via `grpreg::cv.grpreg(penalty =
#'   "gBridge", gamma = 1)` (the group-bridge penalty at `gamma = 1` over
#'   singleton groups is the lasso), CV-selecting `lambda` --- the package's "CV
#'   on the raw lambda" convention (cf. `getBetaCV()`).
#'
#'   **Determinism contract.** `debiasedATT()` and `.simultaneous_cis_bootstrap()`
#'   must obtain the *identical* nuisance so the high-dim band center equals
#'   `debiasedATT()$att` to machine precision. The seed is derived from the data
#'   (`as.integer(N * T)`, `getBetaCV()`'s default), NOT the fit's `cv_seed`
#'   (which is `NA` for BIC-selected fits); both call sites pass the same
#'   `X_final` / `y_final` / `N` / `T` from the same fit, and
#'   `.with_preserved_rng()` leaves the caller's RNG untouched (#177).
#' @param X_final Numeric matrix; the `NT x p` unscaled, fusion + GLS-transformed
#'   design (`fit$internal$X_final`).
#' @param y_final Numeric vector of length `NT`; the GLS-transformed response
#'   (callers pass the `NT`-truncated `y_final`, matching the `add_ridge` design).
#' @param N,T Integer; the number of units and time periods.
#' @param cv_folds Integer; CV folds (fixed by the call sites for determinism).
#' @return Numeric vector of length `p + 1`: the `q = 1` nuisance coefficients
#'   (intercept at position 1) on the original-data scale.
#' @keywords internal
#' @noRd
.fit_q1_nuisance <- function(X_final, y_final, N, T, cv_folds = 10L) {
	p <- ncol(X_final)
	X_scaled <- my_scale(X_final)
	scale_center <- attr(X_scaled, "scaled:center")
	scale_scale <- attr(X_scaled, "scaled:scale")
	# Deterministic seed from the data only (NOT fit$cv_seed, which is NA for
	# BIC-selected fits). as.integer(N * T) is getBetaCV()'s default-seed
	# convention; identical at both call sites since N, T come from the same fit.
	nt_double <- as.numeric(N) * as.numeric(T)
	cv_seed <- if (nt_double > .Machine$integer.max) {
		.Machine$integer.max
	} else {
		as.integer(nt_double)
	}
	cv_fit <- .with_preserved_rng(cv_seed, {
		grpreg::cv.grpreg(
			X = X_scaled,
			y = y_final,
			penalty = "gBridge",
			gamma = 1,
			nfolds = cv_folds
		)
	})
	lam_idx <- .cv_grpreg_at_min(cv_fit)
	stopifnot(nrow(cv_fit$fit$beta) == p + 1L)
	theta_scaled <- cv_fit$fit$beta[, lam_idx]
	stopifnot(length(theta_scaled) == p + 1L, all(!is.na(theta_scaled)))
	.untransform_scaled_theta(
		theta_hat_scaled = theta_scaled,
		p = p,
		scale_center = scale_center,
		scale_scale = scale_scale
	)
}

#' @title Column index of cv.grpreg's lambda.min (exact match, guarded fallback)
#' @description The lockstep-critical step shared by the two cross-validation
#'   call sites (`getBetaCV()` and `.fit_q1_nuisance()`): map `cv.grpreg()`'s
#'   `lambda.min` to its column index in `cv_fit$fit$lambda`. Single-sourced
#'   (#401) because a fix applied to one copy would silently select the wrong
#'   lambda in the other.
#'
#'   **What the fallback is actually for (#431 item 5).** It is *not* for
#'   floating-point drift, which the earlier version of this block claimed.
#'   Read `grpreg::cv.grpreg` (measured against 3.5.0) at its tail: it sets
#'   `lambda.min = lambda[min]` where `lambda` is a *subset of `fit$lambda`* and
#'   `min` is a `which.min()` index, with no arithmetic applied anywhere. So
#'   `lambda.min` is bit-identical to a grid element by construction and cannot
#'   drift. `which(fit$lambda == lambda.min)` can therefore return a length
#'   other than 1 only if (a) the grid contains duplicate values, or (b)
#'   `grpreg` has broken that contract and `lambda.min` is absent from the grid
#'   altogether.
#'
#'   Case (a) selects the **earliest column of `fit$beta`** — duplicates are
#'   equal lambdas, so "the smallest lambda" would say nothing — and is
#'   harmless. Case (b) is a broken upstream contract and must not be papered
#'   over with the nearest grid point, which is what this function used to do
#'   silently; the `all.equal()` assertion below now catches it.
#'
#'   Neither case is reachable with `grpreg` 3.5.0's `gBridge` grid **as
#'   `fetwfe` calls it**. `grpreg:::setupLambda.gBridge`'s `lambda.min > 0`
#'   branch — the only one reached, since neither call site forwards
#'   `lambda.min` or `nlambda` — is a strictly monotone exponential sequence
#'   with `distinct == length` at every resolution tested. (Its other branch,
#'   taken when `lambda.min == 0`, *prepends* an explicit `0`; and the one
#'   degenerate case that would collapse the grid, `lambda.max == 0`, errors
#'   inside `seq()` before a grid exists at all.)
#'
#'   **Why the assertion is not a tautology.** The grid is geometric with ratio
#'   `lambda.min^(-1/(nlambda - 1))`, so adjacent lambdas differ by several
#'   percent — roughly 7% at `gBridge`'s `n > p` default (`lambda.min = 0.001`)
#'   and roughly 3% at its `n <= p` default (`0.05`), this package's usual
#'   regime. Either is six or seven orders of magnitude above `all.equal()`'s
#'   `1.5e-8` relative tolerance, so a wrong neighbour is rejected with enormous
#'   margin. Measured directly: the assertion accepts a `1e-8` relative
#'   perturbation of `lambda_star` and rejects `1e-7`.
#' @param cv_fit A `cv.grpreg` fit (with `$lambda.min` and `$fit$lambda`).
#' @return The integer column index (guaranteed length 1).
#' @keywords internal
#' @noRd
.cv_grpreg_at_min <- function(cv_fit) {
	lambda_star <- cv_fit$lambda.min
	lam_idx <- which(cv_fit$fit$lambda == lambda_star)
	if (length(lam_idx) != 1L) {
		# No exact match, or several. Take the nearest grid point, but assert
		# that it really IS lambda.min (#431 item 5) rather than returning the
		# nearest one silently. Deliberately inside the `if`: on the normal path
		# `lam_idx` came from an exact `==` match, so re-asserting it there
		# would be tautological.
		lam_idx <- which.min(abs(cv_fit$fit$lambda - lambda_star))
		stopifnot(isTRUE(all.equal(cv_fit$fit$lambda[lam_idx], lambda_star)))
	}
	# Residual role of this guard after the assertion above: on an EMPTY grid
	# the new assertion pre-empts it, but a zero-length `lambda.min` against a
	# non-empty grid still lands here, because `all.equal(numeric(0),
	# numeric(0))` is TRUE.
	stopifnot(length(lam_idx) == 1L)
	lam_idx
}
