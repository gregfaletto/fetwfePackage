# Internal helpers shared by the simulation pipeline:
# `generateBaseEffects()`, `genCohortTimeFE()`, `genAssignments()`,
# `genTreatVarsSim()`, `rfunc()`, and `testGenRandomDataInputs()`.
# Moved from R/gen_funcs.R in 1.9.25.

#' Draw covariates and their long-format expansion
#'
#' Helper shared between \code{generateBaseEffects()} (marginal cohort
#' assignment path) and \code{.generateBaseEffectsCovariate()}
#' (covariate-dependent cohort assignment path). Adding a new
#' distribution branch only needs to happen here.
#'
#' @param N Integer. Number of units.
#' @param d Integer. Number of covariates.
#' @param T Integer. Number of time periods.
#' @param distribution Character. \code{"gaussian"} or \code{"uniform"}.
#'
#' @return A list with \code{X} (the \code{N x d} unit-level matrix) and
#'   \code{X_long} (the row-repeated \code{(N * T) x d} long-format
#'   expansion).
#' @keywords internal
#' @noRd
.drawCovariates <- function(N, d, T, distribution = "gaussian") {
	if (distribution == "gaussian") {
		X <- matrix(rnorm(N * d), nrow = N, ncol = d)
	} else if (distribution == "uniform") {
		# U(-sqrt(3), sqrt(3)) so variance = 1 (matches standard normal).
		a <- sqrt(3)
		X <- matrix(runif(N * d, min = -a, max = a), nrow = N, ncol = d)
	} else {
		stop("Unsupported distribution. Please choose 'gaussian' or 'uniform'.")
	}
	X_long <- X[rep(seq_len(N), each = T), , drop = FALSE]
	list(X = X, X_long = X_long)
}


#' Generate Base Fixed Effects and Covariates for Simulation
#'
#' Creates cohort fixed effects, time fixed effects, and a long-format covariate matrix
#' for simulating panel data. Covariates are drawn based on the specified distribution.
#' This is the marginal-cohort-assignment path (cohort assignment independent
#' of covariates); the covariate-dependent path is in
#' \code{.generateBaseEffectsCovariate()}.
#'
#' @param N Integer. Number of units in the panel.
#' @param d Integer. Number of time-invariant covariates.
#' @param T Integer. Number of time periods.
#' @param G Integer. Number of treated cohorts.
#' @param distribution Character. Distribution to generate covariates.
#'   Defaults to \code{"gaussian"}. If set to \code{"uniform"}, covariates are drawn uniformly
#'   from \eqn{[-\sqrt{3}, \sqrt{3}]}.
#' @param guarantee_rank_condition (Optional). Logical. If TRUE, the returned
#' data set is guaranteed to have at least `d + 1` units per cohort, which is
#' necessary for the final design matrix to have full column rank. Default is
#' FALSE, in which case no such condition is enforced.
#'
#' @return A list containing:
#'   \item{cohort_fe}{A matrix of cohort fixed effects (dummy variables).}
#'   \item{time_fe}{A matrix of time fixed effects (dummy variables for periods 2 to T).}
#'   \item{X_long}{A long-format matrix of covariates, repeated for each time period.}
#'   \item{assignments}{A vector of counts indicating how many units fall into
#'         the never-treated group and each of the G treated cohorts.}
#'   \item{cohort_inds}{A list where each element contains the row indices in the
#'         long-format matrices corresponding to the units in a specific treated cohort.}
#' @keywords internal
#' @noRd
generateBaseEffects <- function(
	N,
	d,
	T,
	G,
	distribution = "gaussian",
	guarantee_rank_condition = FALSE
) {
	ret <- genCohortTimeFE(
		N = N,
		T = T,
		G = G,
		d = d,
		guarantee_rank_condition = guarantee_rank_condition
	)

	cov <- .drawCovariates(N = N, d = d, T = T, distribution = distribution)
	X_long <- cov$X_long

	stopifnot(ncol(ret$cohort_fe) == G)

	return(list(
		cohort_fe = ret$cohort_fe,
		time_fe = ret$time_fe,
		X_long = X_long,
		assignments = ret$assignments,
		cohort_inds = ret$inds
	))
}


#' Generate Base Fixed Effects and Covariates under Covariate-Dependent Cohort Assignment
#'
#' Sibling of \code{generateBaseEffects()} for the non-marginal
#' \code{assignment_type}s ("multinomial" / "ordered"). Differs in that
#' cohort membership is drawn per-unit from
#' \code{.sample_cohort_assignments(X, assignment_coefs)} rather than
#' from a uniform multinomial. Units are then sorted by drawn cohort
#' label so the row-block structure expected by downstream consumers
#' (\code{prepXints()}, cohort-block indexing in \code{simulateDataCore()})
#' is preserved.
#'
#' @param N Integer. Number of units.
#' @param d Integer >= 1. Number of covariates.
#' @param T Integer. Number of time periods.
#' @param G Integer. Number of treated cohorts.
#' @param distribution Character. Covariate distribution
#'   (\code{"gaussian"} / \code{"uniform"}).
#' @param guarantee_rank_condition Logical. If TRUE, retry the joint
#'   X + W draw until each cohort has at least \code{d + 1} units (up to
#'   \code{max_iters}; errors with a user-actionable message if exhausted).
#' @param assignment_coefs Non-NULL assignment-coefs object from
#'   \code{.gen_assignment_coefs()}.
#' @param max_iters Integer. Maximum number of joint X + W re-draws when
#'   \code{guarantee_rank_condition = TRUE}. Default 100.
#'
#' @return A list with the same shape as \code{generateBaseEffects()}'s
#'   return value: \code{cohort_fe}, \code{time_fe}, \code{X_long},
#'   \code{assignments}, \code{cohort_inds}. The covariates are sorted
#'   to match the row-block order of \code{cohort_fe}.
#' @keywords internal
#' @noRd
.generateBaseEffectsCovariate <- function(
	N,
	d,
	T,
	G,
	distribution,
	guarantee_rank_condition,
	assignment_coefs,
	max_iters = 100L
) {
	stopifnot(!is.null(assignment_coefs))
	stopifnot(d >= 1)
	stopifnot(N >= (G + 1))
	if (guarantee_rank_condition) {
		stopifnot(N >= (G + 1) * (d + 1))
	}

	# Draw (X, W) jointly with rank-condition retries via the shared
	# `.retry_until_rank_ok()` helper (also used by the auxiliary
	# `.draw_indep_assignments_covariate()` path). RNG state advances per
	# iteration — retries are deterministic given the initial seed.
	res <- .retry_until_rank_ok(
		N = N,
		d = d,
		T = T,
		G = G,
		distribution = distribution,
		guarantee_rank_condition = guarantee_rank_condition,
		assignment_coefs = assignment_coefs,
		max_iters = max_iters,
		label = ""
	)
	X <- res$X
	W <- res$W
	counts <- res$counts

	# Sort units by W (units are exchangeable; sorting is a permutation
	# that doesn't change the joint distribution of (X_i, W_i) but
	# restores the row-block invariant expected by downstream code).
	ord <- order(W)
	W <- W[ord]
	X <- X[ord, , drop = FALSE]
	X_long <- X[rep(seq_len(N), each = T), , drop = FALSE]

	# Now W is sorted: first counts[1] entries are 0 (never-treated),
	# next counts[2] are cohort 1, etc. Build cohort_fe and inds the
	# same way genCohortTimeFE() does.
	assignments <- counts
	stopifnot(sum(assignments) == N)
	stopifnot(length(assignments) == G + 1L)

	cohort_fe <- matrix(0, N * T, G)
	inds <- list()
	first_ind_g <- assignments[1] * T + 1
	for (g in seq_len(G)) {
		last_ind_g <- first_ind_g + assignments[g + 1] * T - 1
		cohort_fe[first_ind_g:last_ind_g, g] <- rep(1, assignments[g + 1] * T)
		inds[[g]] <- first_ind_g:last_ind_g
		first_ind_g <- last_ind_g + 1
	}
	stopifnot(last_ind_g == N * T)

	time_fe <- matrix(0, N * T, T - 1)
	for (t in seq_len(T - 1)) {
		rows_t <- (0:(N - 1)) * T + t + 1
		time_fe[rows_t, t] <- rep(1, N)
	}

	list(
		cohort_fe = cohort_fe,
		time_fe = time_fe,
		X_long = X_long,
		assignments = assignments,
		cohort_inds = inds
	)
}


#' Generate Cohort and Time Fixed Effects Matrices
#'
#' Internal function to create matrices of dummy variables for cohort and time fixed effects.
#' Assumes observations are arranged unit by unit, with T observations per unit.
#'
#' @param N Integer. Number of units.
#' @param T Integer. Number of time periods.
#' @param G Integer. Number of treated cohorts.
#' @param d Integer. Number of covariates (used to ensure minimum units per cohort).
#' @param guarantee_rank_condition (Optional). Logical. If TRUE, the returned
#' data set is guaranteed to have at least `d + 1` units per cohort, which is
#' necessary for the final design matrix to have full column rank. Default is
#' FALSE, in which case only `>= 1` unit per cohort is required (permitting small
#' cohorts and therefore rank-deficient, high-dimensional `p > NT` designs).
#'
#' @return A list containing:
#'   \item{cohort_fe}{An NT x G matrix of cohort dummy variables. The g-th column is 1
#'     if an observation belongs to the g-th treated cohort, 0 otherwise.}
#'   \item{time_fe}{An NT x (T-1) matrix of time dummy variables. The t-th column
#'     (for t from 1 to T-1) corresponds to time period (t+1), and is 1 if an
#'     observation is from that time period, 0 otherwise. Period 1 is the baseline.}
#'   \item{assignments}{An integer vector of length G+1. The first element is the
#'     count of never-treated units. Subsequent elements are counts of units in each
#'     of the G treated cohorts.}
#'   \item{inds}{A list of length G. Each element `inds[[g]]` contains the row indices
#'     in the NT-row matrices that correspond to units in the g-th treated cohort.}
#' @keywords internal
#' @noRd
genCohortTimeFE <- function(N, T, G, d, guarantee_rank_condition = FALSE) {
	# The observations will be arranged row-wise as blocks of T, one unit at
	# a time. So the first T rows correspond to all observations from the first
	# unit, and so on. Therefore the first N_UNTREATED*T rows contain all of
	# the observations from untreated units, then the next N_PER_COHORT*T
	# rows contain the observations from the first cohort, and so on.

	# Each of the G + 1 groups (G treated cohorts + the never-treated) needs at
	# least one unit so its fixed effect is defined; the remaining units are
	# allocated to the groups uniformly at random. Only when
	# guarantee_rank_condition = TRUE do we additionally require >= d + 1 units
	# per group (the OLS full-column-rank condition). Gating the stricter bound
	# behind the flag -- consistent with every other rank guard, e.g.
	# .generateBaseEffectsCovariate() and genAssignments() -- is what lets the
	# generator reach the regularized high-dimensional (p > NT) regime (#293).
	stopifnot(N >= (G + 1))
	if (guarantee_rank_condition) {
		stopifnot(N >= (G + 1) * (d + 1))
	}

	# Generate cohort assignments
	assignments <- genAssignments(
		N = N,
		G = G,
		guarantee_rank_condition = guarantee_rank_condition,
		d = d
	)

	# Cohort fixed effects
	cohort_fe <- matrix(0, N * T, G)
	first_ind_g <- assignments[1] * T + 1

	inds <- list()

	for (g in 1:G) {
		stopifnot(all(cohort_fe[, g] == 0))

		last_ind_g <- first_ind_g + assignments[g + 1] * T - 1

		stopifnot(last_ind_g <= N * T)
		stopifnot(length(first_ind_g:last_ind_g) == assignments[g + 1] * T)
		stopifnot(all(cohort_fe[first_ind_g:last_ind_g, ] == 0))

		# Now add cohort fixed effects in the rth column for these rows
		cohort_fe[first_ind_g:last_ind_g, g] <- rep(1, assignments[g + 1] * T)
		inds[[g]] <- first_ind_g:last_ind_g
		first_ind_g <- last_ind_g + 1

		stopifnot(length(inds[[g]]) == assignments[g + 1] * T)
	}

	stopifnot(last_ind_g == N * T)

	# Time fixed effects: only do 2 through T.
	time_fe <- matrix(0, N * T, T - 1)
	# No fixed effect for first time
	for (t in 1:(T - 1)) {
		stopifnot(all(time_fe[, t] == 0))
		# For each time (except the first), we have N rows of observations which
		# need to be assigned a time fixed effect. Identify those rows:
		# The (t + 1)st row of every block should have a 1 in the (t + 1)st
		# column
		rows_t <- (0:(N - 1)) * T + t + 1
		stopifnot(length(rows_t) == N)
		stopifnot(max(rows_t) <= N * T)
		time_fe[rows_t, t] <- rep(1, N)
	}

	return(list(
		cohort_fe = cohort_fe,
		time_fe = time_fe,
		assignments = assignments,
		inds = inds
	))
}


#' Generate Random Cohort Assignments
#'
#' Assigns N units to G+1 groups (G treated cohorts + 1 never-treated group)
#' ensuring each group has at least one unit.
#'
#' @param N Integer. Total number of units to assign.
#' @param G Integer. Number of treated cohorts.
#' @param guarantee_rank_condition (Optional). Logical. If TRUE, `d` must be
#' provided, and this function will ensure that the returned data set has at
#' least `d + 1` units per cohort, which is necessary for the final design
#' matrix to have full column rank. Default is FALSE, in which case no such
#' condition is enforced.
#' @param d (Optional). Integer. The total number of covariates in the data.
#' Only needs to be provied if guarantee_rank_condition is TRUE.
#'
#' @return An integer vector of length G+1, where the first element is the count
#'   of never-treated units, and subsequent elements are counts for each treated cohort.
#'   The sum of elements equals N.
#' @importFrom stats rmultinom
#' @keywords internal
#' @noRd
genAssignments <- function(N, G, guarantee_rank_condition = FALSE, d = NA) {
	# Make sure at least one observation in each cohort
	stopifnot(N >= G + 1)
	pass_condition <- FALSE
	while (!pass_condition) {
		assignments <- rmultinom(
			n = 1,
			size = N,
			prob = rep(1 / (G + 1), G + 1)
		)[,
			1
		]
		pass_condition <- all(assignments >= 1)
		if (guarantee_rank_condition) {
			stopifnot(all(!is.na(d)))
			stopifnot(d >= 0)
			pass_condition <- all(assignments >= d + 1)
		}
	}

	stopifnot(sum(assignments) == N)
	stopifnot(all(assignments >= 1))
	if (guarantee_rank_condition) {
		stopifnot(all(assignments >= d + 1))
	}
	stopifnot(length(assignments) == G + 1)
	return(assignments)
}


#' Generate Treatment Variable Matrix for Simulations
#'
#' Creates a matrix of treatment dummy variables for simulated panel data.
#' Treatment starts at period g+1 for cohort g.
#'
#' @param n_treats Integer. Total number of unique treatment (cohort x time) effects.
#'   Calculated as \eqn{T \times G - G(G+1)/2}.
#' @param N Integer. Number of units.
#' @param T Integer. Number of time periods.
#' @param G Integer. Number of treated cohorts.
#' @param assignments Integer vector from \code{genAssignments}, indicating unit counts
#'   per cohort (including never-treated).
#' @param cohort_inds List from \code{genCohortTimeFE}, containing row indices for each cohort.
#' @param N_UNTREATED Integer. Number of never-treated units.
#' @param first_inds_test Integer vector. Pre-calculated indices of the first treatment
#'   effect for each cohort (used for assertion).
#'
#' @return A list containing:
#'   \item{treat_mat_long}{An NT x n_treats matrix of treatment dummy variables.
#'     Each column corresponds to a specific cohort-time treatment indicator.}
#'   \item{first_inds}{An integer vector of length G, where `first_inds[g]` is the
#'     column index in `treat_mat_long` corresponding to the first treatment period
#'     for cohort g.}
#' @keywords internal
#' @noRd
genTreatVarsSim <- function(
	n_treats,
	N,
	T,
	G,
	assignments,
	cohort_inds,
	N_UNTREATED,
	first_inds_test
) {
	# Treatment indicators
	treat_mat_long <- matrix(0, N * T, n_treats)
	treat_ind <- 0
	total_feats_added <- 0
	first_inds <- rep(as.integer(NA), G)

	stopifnot(length(cohort_inds) == G)
	stopifnot(length(assignments) == G + 1)
	stopifnot(min(cohort_inds[[1]]) == N_UNTREATED * T + 1)

	all_inds_so_far <- integer()

	for (g in 1:G) {
		n_treats_g <- T - g

		stopifnot(length(cohort_inds[[g]]) == assignments[g + 1] * T)

		counter <- 0L

		for (t in (g + 1):T) {
			treat_ind <- treat_ind + 1

			stopifnot(treat_ind <= n_treats)

			if (t == g + 1) {
				# #185 SB6: keep `first_inds` an integer vector (its documented
				# type and `as.integer(NA)` allocation). `treat_ind` is promoted
				# to double by `treat_ind + 1`; coerce so the `identical()`
				# assertion against `first_inds_test = getFirstInds()` (now
				# integer) holds.
				first_inds[g] <- as.integer(treat_ind)

				stopifnot(treat_ind == first_inds_test[g])
			}

			first_ind_g <- min(cohort_inds[[g]]) + t - (g + 1) + g

			stopifnot(first_ind_g >= min(cohort_inds[[g]]))

			last_ind_g <- first_ind_g + (assignments[g + 1] - 1) * T

			stopifnot(last_ind_g <= max(cohort_inds[[g]]))
			stopifnot(
				(last_ind_g - first_ind_g) / T ==
					round((last_ind_g - first_ind_g) / T)
			)
			stopifnot((last_ind_g - first_ind_g) / T == assignments[g + 1] - 1)

			g_inds <- seq(first_ind_g, last_ind_g, by = T)

			all_inds_so_far <- c(all_inds_so_far, g_inds)

			stopifnot(length(g_inds) == assignments[g + 1])
			stopifnot(all(treat_mat_long[, treat_ind] == 0))

			treat_mat_long[g_inds, treat_ind] <- 1
			total_feats_added <- total_feats_added + 1
			counter <- counter + 1
		}

		stopifnot(counter == n_treats_g)

		treat_inds_g <- first_inds[g]:(first_inds[g] + length((g + 1):T) - 1)

		stopifnot(all(treat_inds_g <= n_treats))
		stopifnot(length(treat_inds_g) == n_treats_g)
		stopifnot(all(
			colSums(treat_mat_long[, treat_inds_g, drop = FALSE]) ==
				assignments[g + 1]
		))

		if (g < G) {
			stopifnot(last_ind_g == min(cohort_inds[[g + 1]]) - 1)
			stopifnot(last_ind_g == max(cohort_inds[[g]]))
		}
	}

	stopifnot(max(g_inds) == N * T)

	stopifnot(all(colSums(treat_mat_long) >= 1))
	stopifnot(total_feats_added == n_treats)

	stopifnot(all(!is.na(first_inds)))
	stopifnot(length(first_inds) == length(unique(first_inds)))
	stopifnot(first_inds[1] == 1)

	stopifnot(identical(first_inds, first_inds_test))

	return(list(treat_mat_long = treat_mat_long, first_inds = first_inds))
}


#' Generate Random Signs
#'
#' Generates a vector of -1s and 1s, with a specified probability for 1.
#'
#' @param n Integer. Number of signs to generate.
#' @param prob Numeric. Probability of generating a 1 (otherwise -1).
#'
#' @return An integer vector of length n containing -1s and 1s.
#' @importFrom stats rbinom
#' @keywords internal
#' @noRd
rfunc <- function(n, prob) {
	# sample(c(-1, 1), size=n, replace=TRUE)
	vec <- rbinom(n = n, size = 1, prob = prob)
	vec[vec == 0] <- -1
	stopifnot(all(vec %in% c(-1, 1)))
	return(vec)
}


#' Test Inputs for Random Data Generation
#'
#' Validates parameters for the `simulateDataCore` function, ensuring consistency
#' in dimensions and expected length of the beta coefficient vector.
#' TODO: This function assumes `gen_ints = TRUE` (i.e., full design matrix with all interactions).
#'
#' @param beta Numeric vector. The true coefficient vector used for data generation.
#' @param G Integer. Number of treated cohorts.
#' @param T Integer. Number of time periods.
#' @param d Integer. Number of time-invariant covariates.
#' @param N Integer. Number of units.
#' @param sig_eps_sq Numeric. Variance of idiosyncratic error.
#' @param sig_eps_c_sq Numeric. Variance of unit-level random effect.
#'
#' @return A list containing:
#'   \item{num_treats}{Integer. The calculated number of base treatment effects.}
#'   \item{p_expected}{Integer. The expected length of the `beta` vector given G, T, d,
#'     assuming a full design matrix with all interactions.}
#' @seealso \code{\link{getNumTreats}}, \code{\link{getP}}
#' @keywords internal
#' @noRd
testGenRandomDataInputs <- function(
	beta,
	G,
	T,
	d,
	N,
	sig_eps_sq,
	sig_eps_c_sq
) {
	stopifnot(G <= T - 1)
	stopifnot(T >= 2)
	stopifnot(G >= 1)
	stopifnot(N >= G + 1)
	stopifnot(sig_eps_sq > 0)
	stopifnot(sig_eps_c_sq >= 0)

	# Compute the number of treatment effects (common to both cases)
	num_treats <- getNumTreats(G = G, T = T)

	# --- Full design matrix with interactions ---
	# Expected number of columns:
	# If d > 0: p = G + (T - 1) + d + d*G + d*(T - 1) + num_treats + num_treats*d
	# If d == 0: p = G + (T - 1) + num_treats
	p_expected <- getP(G = G, T = T, d = d, num_treats = num_treats)

	if (length(beta) != p_expected) {
		stop(sprintf(
			"For gen_ints = TRUE, length(beta) must be %d",
			p_expected
		))
	}

	return(list(num_treats = num_treats, p_expected = p_expected))
}
