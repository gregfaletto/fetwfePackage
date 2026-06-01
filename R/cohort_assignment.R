# Covariate-dependent cohort-assignment internals for the FETWFE
# simulator (#162). Four helpers:
#
#   .gen_assignment_coefs()    draws gamma (and cutpoints) given type/strength.
#   .compute_cohort_prob_matrix()    N x (R+1) per-unit propensity matrix.
#   .sample_cohort_assignments() per-unit categorical draws of W_i in 0..R.
#   .expected_cohort_probs()   length-(R+1) MC estimate of E[pi_r(X)].
#
# All four are @keywords internal @noRd. The factoring of
# .compute_cohort_prob_matrix() avoids duplicating the multinomial-softmax
# and ordered-cumprobs logic in both the sampler and the truth-
# derivation MC integrator (round-1 reviewer / sentinel convergence).

#' Validate a non-negative numeric scalar argument (with optional NULL allowance)
#'
#' Shared validator for the `assignment_strength` and
#' `assignment_interaction_strength` arguments on `genCoefs()` and the
#' `strength` argument on `.gen_assignment_coefs()`. Centralises the
#' "non-negative numeric scalar" check that existed at three sites
#' before #191 added `assignment_interaction_strength`. Errors with a
#' user-readable message naming the offending argument.
#'
#' @param x The value to validate.
#' @param name Character; the user-facing name of the argument (used in
#'   the error message).
#' @param allow_null Logical; if `TRUE`, `x = NULL` is accepted as a
#'   valid default. Used for `assignment_interaction_strength`'s
#'   `NULL` fall-through to `assignment_strength`. Default `FALSE`.
#'
#' @return Invisibly `NULL`. Called for its side-effect (`stop()` on
#'   invalid input).
#'
#' @keywords internal
#' @noRd
.validate_strength_arg <- function(x, name, allow_null = FALSE) {
	if (is.null(x)) {
		if (allow_null) {
			return(invisible(NULL))
		}
		stop(sprintf("%s must be a non-negative numeric scalar", name))
	}
	if (!is.numeric(x) || length(x) != 1 || x < 0) {
		if (allow_null) {
			stop(sprintf(
				"%s must be a non-negative numeric scalar (or NULL)",
				name
			))
		}
		stop(sprintf("%s must be a non-negative numeric scalar", name))
	}
	invisible(NULL)
}


#' Generate assignment-coefficient object for covariate-dependent cohorts
#'
#' Returns the (gamma, cutpoints) object that parameterizes the
#' propensity-score model used by `.compute_cohort_prob_matrix()`.
#'
#' Seed-offset convention: callers pass `seed` (typically the main
#' `genCoefs()` seed plus an offset). The package uses three deterministic
#' offsets so that the three random streams stay reproducible and
#' non-colliding:
#'
#'   - `seed`      -> main coefficient vector / theta draw (`genCoefsCore`).
#'   - `seed + 1L` -> assignment coefficients (this function).
#'   - `seed + 2L` -> Monte Carlo integration in `getTes()` for
#'                    `E[pi_r(X)]`.
#'
#' @param R Integer >= 2. Number of treated cohorts.
#' @param d Integer >= 1. Number of covariates. Only called when
#'   `type != "marginal"`, which requires `d >= 1`.
#' @param type Character. One of `"multinomial"` or `"ordered"`. The
#'   `"marginal"` case is handled by the caller and never reaches this
#'   function.
#' @param strength Non-negative numeric. Scales the (otherwise standard
#'   normal) gamma draw. `strength = 0` yields gamma = 0 so the
#'   cohort-probability distribution reduces to uniform 1/(R+1) by
#'   construction.
#' @param interactions Optional. A list of canonicalized length-2
#'   integer-pair vectors `list(c(j, k), ...)` whose elementwise products
#'   `X[, j] * X[, k]` augment the propensity linear predictor. The list
#'   is assumed to be canonicalized + deduplicated by the caller
#'   (`genCoefs()`). Defaults to `NULL` (no interactions; v1.14.0
#'   behavior is preserved byte-identically, including RNG ordering).
#'   New in 1.14.1.
#' @param interaction_strength Optional non-negative numeric scalar.
#'   Scales the Gaussian draws of the interaction coefficients
#'   independently of `strength`. Defaults to `NULL`, which falls
#'   through to `strength`. New in 1.14.1.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list with elements `type`, `strength`, `coefs`, `interactions`,
#'   `delta`, `interaction_strength`, and (for `"ordered"`) `cutpoints`.
#'   For multinomial, `coefs` is a `d x R` matrix; column `r` is
#'   `gamma_r` (the never-treated reference cohort has `gamma_0 = 0`
#'   implicitly), and `delta` is a `K x R` matrix of per-cohort
#'   interaction coefficients (or `NULL` when no interactions). For
#'   ordered, `coefs` is a length-`d` vector (one shared `gamma`),
#'   `delta` is a length-`K` vector (or `NULL`), and `cutpoints` is a
#'   length-`R` vector of cutpoints chosen so the marginal cohort
#'   probabilities are uniform 1/(R+1). `interaction_strength` is
#'   `NULL` when no interactions; otherwise the effective scaling
#'   factor used for the `delta` draws.
#'
#' @details
#' The ordered cutpoints are chosen by root-finding on the
#' marginal-uniform condition: for each `r` in `1..R`, solve
#' `E_X[plogis(alpha_r - gamma^T X - delta^T X_int)] = r / (R + 1)` via
#' `uniroot()` over a Monte Carlo reference sample of the augmented
#' linear predictor (standard McCullagh proportional-odds parameterization
#' with the subtraction convention: high augmented linpred shifts mass UP
#' the ordinal scale toward never-treated; low augmented linpred shifts
#' mass toward the earliest-adopting cohort). At `strength = 0` (gamma = 0)
#' AND no interactions the cutpoints collapse to `qlogis((1:R) / (R + 1))`
#' and the marginal probabilities are exactly uniform. At positive
#' strength the marginal probabilities are approximately uniform up to
#' MC noise.
#'
#' The `uniroot()` search interval is `c(-200, 200)` to accommodate
#' stress-test usage at very high strength (linpred ranges scale with
#' strength * sqrt(d); for stress-test usage at strength = 50 the
#' cutpoints can land near +/- 50, which would exceed a narrower
#' interval).
#'
#' The augmented linear predictor extends the standard multinomial-logit
#' / proportional-odds parameterization (see Faletto 2025 line 1016 for
#' the linear-in-X form). Interaction terms are user-specified pairwise
#' products of covariate columns and remain within the paper's framework
#' because the multinomial-logit MLE with interaction features still
#' satisfies Assumption (Psi-IF) (line 2018) via M-estimator theory.
#'
#' @keywords internal
#' @noRd
.gen_assignment_coefs <- function(
	R,
	d,
	type,
	strength,
	interactions = NULL,
	interaction_strength = NULL,
	seed = NULL
) {
	if (!is.null(seed)) {
		set.seed(seed)
	}
	.validate_strength_arg(strength, "assignment_strength")
	if (!(type %in% c("multinomial", "ordered"))) {
		stop(sprintf(
			"Internal error: .gen_assignment_coefs() called with type = '%s'",
			type
		))
	}
	if (d < 1) {
		stop(
			"Covariate-dependent cohort assignment requires d >= 1 (at least one covariate)"
		)
	}
	K <- if (is.null(interactions)) 0L else length(interactions)
	effective_int_strength <- if (is.null(interaction_strength)) {
		strength
	} else {
		interaction_strength
	}
	if (type == "multinomial") {
		coefs_mat <- matrix(
			rnorm(d * R) * strength,
			nrow = d,
			ncol = R
		)
		delta_mat <- if (K > 0L) {
			matrix(
				rnorm(K * R) * effective_int_strength,
				nrow = K,
				ncol = R
			)
		} else {
			NULL
		}
		return(list(
			type = "multinomial",
			strength = strength,
			coefs = coefs_mat,
			interactions = interactions,
			delta = delta_mat,
			interaction_strength = if (K > 0L) effective_int_strength else NULL
		))
	}
	# type == "ordered"
	gamma <- rnorm(d) * strength
	delta <- if (K > 0L) {
		rnorm(K) * effective_int_strength
	} else {
		NULL
	}
	# Root-find on the marginal-uniform condition: solve for alpha_r such
	# that mean(plogis(alpha_r - linpred_mc)) = r / (R + 1). At strength = 0
	# and no interactions, linpred_mc is identically 0 and the cutpoints
	# reduce to qlogis((1:R) / (R + 1)) by construction. With interactions,
	# the MC sample integrates over the augmented covariate distribution so
	# the marginal-uniform invariant is preserved up to MC noise. See
	# Faletto 2025 line 1016 for the linear-in-X form; this extends the
	# parameterization to include user-specified pairwise interactions.
	M_cutpoints <- 5000L
	X_mc <- matrix(
		rnorm(M_cutpoints * d),
		nrow = M_cutpoints,
		ncol = d
	)
	linpred_mc <- as.numeric(X_mc %*% gamma)
	if (K > 0L) {
		X_int_mc <- .build_x_int(X_mc, interactions)
		linpred_mc <- linpred_mc + as.numeric(X_int_mc %*% delta)
	}
	solve_cutpoint <- function(target_prob) {
		stats::uniroot(
			function(alpha) {
				mean(stats::plogis(alpha - linpred_mc)) - target_prob
			},
			interval = c(-200, 200)
		)$root
	}
	cutpoints <- vapply(
		seq_len(R),
		function(r) solve_cutpoint(r / (R + 1)),
		numeric(1)
	)
	list(
		type = "ordered",
		strength = strength,
		coefs = gamma,
		cutpoints = cutpoints,
		interactions = interactions,
		delta = delta,
		interaction_strength = if (K > 0L) effective_int_strength else NULL
	)
}


#' Build the per-pair pairwise-product augmentation matrix `X_int`
#'
#' Returns the `N x K` matrix whose column `k` is the elementwise product
#' `X[, j_k] * X[, k_k]` for the `k`-th pair `c(j_k, k_k)` in
#' `interactions`. Used only inside the propensity model — the outcome
#' design matrix never sees `X_int`.
#'
#' Centralized to one helper so the live propensity computation in
#' `.compute_cohort_prob_matrix()` and the ordered cutpoint root-finder's
#' Monte Carlo sample in `.gen_assignment_coefs()` use byte-identical
#' construction logic. See `.workflow/WORKFLOW_LESSONS.md` §14 (Class A
#' drift discipline).
#'
#' @param X Numeric matrix, `N x d`.
#' @param interactions A list of canonicalized length-2 integer-pair
#'   vectors. `NULL` or empty list returns `NULL`.
#'
#' @return Numeric matrix `N x K`, or `NULL` when no interactions are
#'   specified.
#'
#' @keywords internal
#' @noRd
.build_x_int <- function(X, interactions) {
	if (is.null(interactions) || length(interactions) == 0L) {
		return(NULL)
	}
	K <- length(interactions)
	X_int <- matrix(NA_real_, nrow = nrow(X), ncol = K)
	for (k in seq_len(K)) {
		pair <- interactions[[k]]
		X_int[, k] <- X[, pair[1L]] * X[, pair[2L]]
	}
	X_int
}


#' Compute the N x (R + 1) per-unit cohort-probability matrix
#'
#' Given per-unit covariates `X` (an `N x d` matrix) and the
#' `assignment_coefs` object produced by `.gen_assignment_coefs()`,
#' returns an `N x (R + 1)` matrix of cohort probabilities. Column 1 is
#' the never-treated reference cohort; columns `2..R+1` are the
#' treated cohorts `1..R`.
#'
#' Shared by both the sampler (`.sample_cohort_assignments()`) and the
#' MC integrator (`.expected_cohort_probs()`) so the multinomial-softmax
#' and ordered-cumprobs logic lives in exactly one place.
#'
#' @param X Numeric matrix, `N x d`.
#' @param assignment_coefs The list returned by `.gen_assignment_coefs()`.
#'   Must be non-`NULL`; the marginal path is dispatched separately and
#'   never reaches this function (cf. the `simulateDataCore()` dispatch
#'   on `assignment_type`).
#'
#' @return Numeric matrix, `N x (R + 1)`. Each row sums to 1.
#'
#' @keywords internal
#' @noRd
.compute_cohort_prob_matrix <- function(X, assignment_coefs) {
	stopifnot(is.matrix(X))
	stopifnot(!is.null(assignment_coefs))
	type <- assignment_coefs$type
	# Interactions (NULL or empty list -> no augmentation; v1.14.0 byte path).
	# See Faletto 2025 line 1016 for the linear-in-X form; this extends the
	# parameterization to include user-specified pairwise interactions.
	interactions <- assignment_coefs$interactions
	delta <- assignment_coefs$delta
	has_int <- !is.null(interactions) && length(interactions) > 0L
	if (type == "multinomial") {
		# linpred is N x R; the never-treated reference column is identically 0.
		linpred <- X %*% assignment_coefs$coefs
		if (has_int) {
			X_int <- .build_x_int(X, interactions)
			linpred <- linpred + X_int %*% delta
		}
		log_probs_unnorm <- cbind(0, linpred)
		# Subtract row maxima for numerical stability (invariant under softmax).
		row_max <- apply(log_probs_unnorm, 1, max)
		probs <- exp(log_probs_unnorm - row_max)
		probs <- probs / rowSums(probs)
		return(probs)
	}
	if (type == "ordered") {
		linpred <- as.numeric(X %*% assignment_coefs$coefs)
		if (has_int) {
			X_int <- .build_x_int(X, interactions)
			linpred <- linpred + as.numeric(X_int %*% delta)
		}
		cutpoints <- assignment_coefs$cutpoints
		# Ordinal scale W' = 1..R+1 in temporal-adoption order with never-
		# treated at the top:
		#
		#   slot 1 = cohort 1 (earliest treated)
		#   slot 2 = cohort 2
		#   ...
		#   slot R = cohort R (latest treated)
		#   slot R+1 = never-treated ("treatment timing -> infinity")
		#
		# Standard McCullagh proportional-odds parameterization with the
		# subtraction convention:
		#
		#   P(W' <= r | X) = plogis(alpha_r - gamma^T X)   for r = 1..R.
		#
		# Low gamma^T x makes the cumulative prob LARGE -> mass concentrates
		# on LOW W' (earliest-adopting cohorts). High gamma^T x makes the
		# cumulative prob SMALL -> mass concentrates on HIGH W' (later
		# cohorts and, at the extreme, never-treated). As gamma^T X
		# increases monotonically from -infinity to +infinity, the modal
		# cohort traverses the ordinal scale in order: cohort 1 -> cohort 2
		# -> ... -> cohort R -> never-treated. gamma is the propensity for
		# being LATER on the treatment-timing scale.
		cumprobs <- vapply(
			cutpoints,
			function(a) stats::plogis(a - linpred),
			numeric(length(linpred))
		)
		# vapply() returns a vector when length(linpred) == 1; force matrix shape.
		if (!is.matrix(cumprobs)) {
			cumprobs <- matrix(
				cumprobs,
				nrow = length(linpred),
				ncol = length(cutpoints)
			)
		}
		# Probabilities in ORDINAL order [cohort 1, cohort 2, ..., cohort R,
		# never-treated]:
		probs_ord <- cbind(
			cumprobs[, 1, drop = FALSE],
			cumprobs[, -1, drop = FALSE] -
				cumprobs[, -ncol(cumprobs), drop = FALSE],
			1 - cumprobs[, ncol(cumprobs)]
		)
		# Rearrange to the package's COHORT-LABEL order [cohort 0 (never-
		# treated), cohort 1, ..., cohort R] so the rest of the simulator
		# (which indexes cohort 0 = never as the first column) is unaffected.
		R_local <- ncol(probs_ord) - 1L
		probs <- cbind(
			probs_ord[, R_local + 1L, drop = FALSE],
			probs_ord[, seq_len(R_local), drop = FALSE]
		)
		return(probs)
	}
	stop(sprintf(
		"Internal error: .compute_cohort_prob_matrix() called with assignment_coefs$type = '%s'",
		type
	))
}


#' Sample per-unit cohort assignments from the propensity-score model
#'
#' For each row of `X`, draw a categorical cohort label `W_i` from
#' `Categorical(p_i)` where `p_i` is the corresponding row of
#' `.compute_cohort_prob_matrix(X, assignment_coefs)`. Vectorized via
#' `rowSums(u >= cumprobs_mat)` so the inner loop is implicit (no
#' per-row `apply` / `sample.int` overhead).
#'
#' @param X Numeric matrix, `N x d`.
#' @param assignment_coefs Non-`NULL` assignment-coefs object.
#'
#' @return Integer vector of length `N`, values in `0..R`.
#'
#' @keywords internal
#' @noRd
.sample_cohort_assignments <- function(X, assignment_coefs) {
	probs <- .compute_cohort_prob_matrix(X, assignment_coefs)
	N <- nrow(probs)
	u <- stats::runif(N)
	# Row-cumulative probs: rowSums(u >= cumprobs_mat) returns the number
	# of bins each u_i exceeds, equivalently the cohort label W_i in 0..R.
	# The last cumulative prob is exactly 1, so u_i strictly less than 1
	# always falls in some bin (u_i == 1 occurs with probability 0).
	cumprobs_mat <- t(apply(probs, 1, cumsum))
	assignments <- rowSums(u >= cumprobs_mat)
	# Cap at R defensively (in case of floating-point edge cases where
	# u_i == 1.0). Without this cap, an exact u_i == 1 would map to R + 1.
	assignments <- pmin(assignments, ncol(probs) - 1L)
	as.integer(assignments)
}


#' Monte Carlo estimate of `E[pi_r(X)]` for `r` in `0..R`
#'
#' Used by `getTes()` to compute the propensity-weighted population ATT
#' truth under non-marginal DGPs. Returns a length-`(R + 1)` numeric
#' vector summing to 1 (the marginal cohort probabilities under the
#' specified DGP and X distribution).
#'
#' @param assignment_coefs Non-`NULL` assignment-coefs object.
#' @param d Integer >= 1. Number of covariates (must match the
#'   dimension of `assignment_coefs$coefs`).
#' @param distribution Character. One of `"gaussian"` or `"uniform"`,
#'   matching `simulateData()`'s `distribution` argument.
#' @param M Integer. Number of Monte Carlo draws. Default 10000 gives
#'   ~0.5% standard error on each component.
#' @param seed Optional integer for reproducibility. Per the
#'   seed-offset convention in `.gen_assignment_coefs()`, this is
#'   typically `coefs_obj$seed + 2L`.
#'
#' @return Numeric vector of length `R + 1`, sums to 1.
#'
#' @keywords internal
#' @noRd
.expected_cohort_probs <- function(
	assignment_coefs,
	d,
	distribution = "gaussian",
	M = 10000L,
	seed = NULL
) {
	if (!is.null(seed)) {
		set.seed(seed)
	}
	X_mc <- if (distribution == "gaussian") {
		matrix(rnorm(M * d), nrow = M, ncol = d)
	} else if (distribution == "uniform") {
		a <- sqrt(3)
		matrix(stats::runif(M * d, min = -a, max = a), nrow = M, ncol = d)
	} else {
		stop(sprintf(
			"Unsupported distribution '%s'. Choose 'gaussian' or 'uniform'.",
			distribution
		))
	}
	probs <- .compute_cohort_prob_matrix(X_mc, assignment_coefs)
	colMeans(probs)
}


#' Draw auxiliary cohort counts under the covariate-dependent DGP
#'
#' Under non-marginal `assignment_type`, the `indep_counts` slot of the
#' simulated data must come from the SAME DGP as the main sample (BL2
#' fix per round-1 review). Drawing it from the uniform multinomial,
#' as the marginal path does, would silently break the asymptotic-
#' exactness guarantee for `indep_counts`-using downstream tests.
#'
#' Re-uses the same joint-(X, W) draw + rank-condition retry loop as
#' `.generateBaseEffectsCovariate()` but only returns the count vector
#' (not the design-matrix bits).
#'
#' @param N Integer. Number of units.
#' @param d Integer >= 1. Number of covariates.
#' @param T Integer. Number of time periods (kept for signature parity
#'   with `.drawCovariates()`; the count vector doesn't depend on T).
#' @param R Integer. Number of treated cohorts.
#' @param distribution Character. Covariate distribution.
#' @param guarantee_rank_condition Logical. If TRUE, retry until each
#'   cohort has at least `d + 1` units.
#' @param assignment_coefs Non-NULL assignment-coefs object.
#' @param max_iters Integer. Maximum retry attempts. Default 100.
#'
#' @return Integer vector of length `R + 1` of cohort counts, summing to N.
#' @keywords internal
#' @noRd
.draw_indep_assignments_covariate <- function(
	N,
	d,
	T,
	R,
	distribution,
	guarantee_rank_condition,
	assignment_coefs,
	max_iters = 100L
) {
	stopifnot(!is.null(assignment_coefs))
	res <- .retry_until_rank_ok(
		N = N,
		d = d,
		T = T,
		R = R,
		distribution = distribution,
		guarantee_rank_condition = guarantee_rank_condition,
		assignment_coefs = assignment_coefs,
		max_iters = max_iters,
		label = "auxiliary "
	)
	as.integer(res$counts)
}


#' Retry per-unit (X, W) draws until the cohort-count rank condition holds
#'
#' Shared helper between `.generateBaseEffectsCovariate()` (main-sample path)
#' and `.draw_indep_assignments_covariate()` (auxiliary `indep_counts` path).
#' Loops until either (a) the drawn cohort counts satisfy the rank condition
#' (each cohort has >= 1 unit, and >= d + 1 if `guarantee_rank_condition`),
#' or (b) the retry budget `max_iters` is exhausted (then stop with a
#' user-actionable error).
#'
#' Each iteration draws a fresh (X, W) pair via `.drawCovariates()` +
#' `.sample_cohort_assignments()`. The RNG state advances per iteration —
#' retries are deterministic given the initial seed but use fresh random
#' draws each time.
#'
#' @param N Integer. Number of units.
#' @param d Integer >= 1. Number of covariates.
#' @param T Integer. Number of time periods (for X_long expansion).
#' @param R Integer. Number of treated cohorts.
#' @param distribution Character. Covariate distribution.
#' @param guarantee_rank_condition Logical.
#' @param assignment_coefs Non-NULL assignment-coefs object.
#' @param max_iters Integer. Maximum retry attempts.
#' @param label Character. Inserted into the rank-condition-exhausted
#'   error message (empty string for the main sample, `"auxiliary "` for
#'   the indep_assignments path). Lets the same loop emit a contextually
#'   accurate message at each call site.
#'
#' @return List with `X` (N × d matrix), `X_long` (N*T × d matrix),
#'   `W` (length-N integer vector, values 0..R), and `counts` (length-(R+1)
#'   integer vector).
#' @keywords internal
#' @noRd
.retry_until_rank_ok <- function(
	N,
	d,
	T,
	R,
	distribution,
	guarantee_rank_condition,
	assignment_coefs,
	max_iters,
	label
) {
	attempt <- 0L
	repeat {
		attempt <- attempt + 1L
		cov <- .drawCovariates(
			N = N,
			d = d,
			T = T,
			distribution = distribution
		)
		W <- .sample_cohort_assignments(cov$X, assignment_coefs)
		counts <- tabulate(W + 1L, nbins = R + 1L)
		min_count <- min(counts)
		ok <- min_count >= 1L
		if (guarantee_rank_condition) {
			ok <- ok && (min_count >= d + 1L)
		}
		if (ok) {
			break
		}
		if (attempt >= max_iters) {
			msg <- sprintf(
				paste0(
					"Failed to draw %scohort assignments satisfying the rank condition ",
					"(min count >= d + 1) after %d attempts at assignment_type = '%s', ",
					"assignment_strength = %g, N = %d, R = %d, d = %d. ",
					"Under covariate-dependent assignment with strong covariate-cohort ",
					"coupling, the smallest treated cohorts can have arbitrarily small ",
					"counts. Try a larger N, a smaller assignment_strength, or set ",
					"guarantee_rank_condition = FALSE."
				),
				label,
				max_iters,
				assignment_coefs$type,
				assignment_coefs$strength,
				N,
				R,
				d
			)
			# When interactions are present (K > 0L), .gen_assignment_coefs()
			# stores the effective interaction strength on assignment_coefs.
			# Surface it so the user knows which knob to lower; the value is
			# the effective (post-fallthrough) value, which equals
			# assignment_strength when assignment_interaction_strength was
			# left at its NULL default (per round-1 Q2 + RC1).
			if (!is.null(assignment_coefs$interaction_strength)) {
				msg <- paste0(
					msg,
					sprintf(
						" (effective assignment_interaction_strength = %g; consider lowering it or assignment_strength.)",
						assignment_coefs$interaction_strength
					)
				)
			}
			stop(msg)
		}
	}
	list(
		X = cov$X,
		X_long = cov$X_long,
		W = W,
		counts = counts
	)
}
