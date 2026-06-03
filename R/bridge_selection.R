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
#'      is the number of non-zero coefficients (including intercept).
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
	# Number of fitted coefficients (including intercept), one per lambda.
	model_sizes <- as.integer(colSums(fit$beta != 0))

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
	BICs <- nt_double * log(mse_hat) + model_sizes * log(nt_double)

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
#'     coefficients (including intercept) at the selected lambda.}
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
	# version-stable. We save and restore the caller's .Random.seed via
	# on.exit() so the call leaves the caller's RNG state untouched
	# (issue #177 — without this, every default-path fetwfe() / betwfe()
	# call would silently mutate the user's seed, a v1.13.0 regression
	# vs the v1.12.x BIC default).
	old_rng <- if (
		exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	) {
		.GlobalEnv$.Random.seed
	} else {
		NULL
	}
	on.exit(
		{
			if (is.null(old_rng)) {
				if (
					exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				) {
					rm(".Random.seed", envir = .GlobalEnv)
				}
			} else {
				assign(".Random.seed", old_rng, envir = .GlobalEnv)
			}
		},
		add = TRUE
	)
	set.seed(cv_seed)

	cv_fit <- grpreg::cv.grpreg(
		X = X_final_scaled,
		y = y_final,
		penalty = "gBridge",
		gamma = gamma,
		nfolds = cv_folds
	)

	lambda_star <- cv_fit$lambda.min
	lam_idx <- which(cv_fit$fit$lambda == lambda_star)
	if (length(lam_idx) != 1L) {
		# Fall back to nearest match in case of floating-point comparison drift.
		lam_idx <- which.min(abs(cv_fit$fit$lambda - lambda_star))
	}
	stopifnot(length(lam_idx) == 1L)
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
		lambda_star_model_size = sum(theta_hat_full != 0),
		# Carry the cv.grpreg fit object out — fetwfe_core() / betwfe_core()
		# need access to `fit$lambda`, `fit$beta`, and the four
		# lambda.max/min diagnostics that .fit_bridge_with_lambda_path()
		# would otherwise have computed.
		fit = cv_fit$fit,
		cv_seed = cv_seed
	)
}
