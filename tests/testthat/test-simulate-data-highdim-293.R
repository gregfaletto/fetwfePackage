# Tests for high-dimensional (p > NT) data generation (#293).
#
# Before #293, genCohortTimeFE() carried an UNCONDITIONAL guard
# `stopifnot(N >= (G + 1) * (d + 1))` (R/sim_helpers.R), so every group needed
# >= d + 1 units and the generator structurally lived in the p < NT regime --
# the high-dimensional debiasedATT() path (#31) could not be exercised on
# synthetic p > NT data with a known data-generating coefficient vector. The fix gates that stricter bound
# behind `guarantee_rank_condition` (consistent with every other rank guard), so
# `guarantee_rank_condition = FALSE` (the default) only needs >= 1 unit per
# group and small cohorts -> p > NT become reachable.

# A small-N, large-d config that lands the design well into p > NT.
.HD <- list(G = 3L, T = 5L, d = 22L, N = 8L) # N >= G+1 = 4

.make_hd_coefs <- function() {
	genCoefs(
		G = .HD$G,
		T = .HD$T,
		d = .HD$d,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
}

test_that("simulateData() reaches p > NT at small N with the default (rank condition off)", {
	coefs <- .make_hd_coefs()
	# Previously errored at the unconditional (G+1)(d+1) guard; now succeeds.
	dat <- simulateData(
		coefs,
		N = .HD$N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	expect_s3_class(dat, "FETWFE_simulated")
	pd <- dat$pdata
	expect_equal(length(unique(pd[[dat$unit_var]])), .HD$N)
	# The generated design is genuinely high-dimensional.
	expect_gt(dat$p, .HD$N * .HD$T)
	# Ground truth is carried on the object (the point of #293: synthetic data
	# with known coefficients for validating the p > NT method).
	expect_false(is.null(dat$coefs))
})

test_that("the p > NT simulated panel fits and debiasedATT() runs its high-dim path", {
	coefs <- .make_hd_coefs()
	dat <- simulateData(
		coefs,
		N = .HD$N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	# REML cannot estimate the variances at p > NT, so supply them.
	fit <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	X <- fit$internal$X_final
	expect_gt(ncol(X), nrow(X)) # p > NT
	db <- debiasedATT(fit)
	expect_true(db$converged)
	expect_lte(db$feasibility, db$lambda_node * (1 + 1e-9))
	expect_true(is.finite(db$att))
	expect_gt(db$se, 0)
})

test_that("guarantee_rank_condition = TRUE still enforces N >= (G+1)(d+1)", {
	coefs <- .make_hd_coefs()
	# The strict (full-rank) path is preserved: small N still errors.
	expect_error(
		simulateData(
			coefs,
			N = .HD$N,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = 1,
			guarantee_rank_condition = TRUE
		),
		"d \\+ 1"
	)
})

test_that("genCohortTimeFE() gates the rank bound on guarantee_rank_condition", {
	G <- .HD$G
	T <- .HD$T
	d <- .HD$d
	# Minimum feasible N (one unit per group) succeeds with the rank condition off
	# and yields a valid, fully-allocated set of small cohorts.
	set.seed(1)
	res <- genCohortTimeFE(N = G + 1L, T = T, G = G, d = d)
	expect_equal(sum(res$assignments), G + 1L)
	expect_true(all(res$assignments >= 1L)) # >= 1 unit per group, < d + 1
	expect_true(any(res$assignments < d + 1L))
	expect_equal(ncol(res$cohort_fe), G)
	expect_equal(nrow(res$cohort_fe), (G + 1L) * T)
	# With the rank condition on, the same small N errors.
	expect_error(
		genCohortTimeFE(
			N = G + 1L,
			T = T,
			G = G,
			d = d,
			guarantee_rank_condition = TRUE
		),
		"d \\+ 1"
	)
})

test_that("an all-singleton (N = G+1) high-dim panel generates and fits without erroring", {
	coefs <- .make_hd_coefs()
	# The extreme corner: exactly one unit per group (every treated cohort AND
	# the never-treated group is a singleton). Within-cohort covariate centering
	# then makes the treatment x covariate interaction columns exactly zero (no
	# NaN -- my_scale() guards zero-variance columns), so generation and a
	# regularized fit must still run to completion. The debiased SE may collapse
	# to 0 at this near-zero residual-d.o.f. corner; that is honest, not an error.
	# Reaching these lines without an error/NaN IS the assertion.
	dat <- simulateData(
		coefs,
		N = .HD$G + 1L,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	expect_equal(length(unique(dat$pdata[[dat$unit_var]])), .HD$G + 1L)
	fit <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	expect_gt(ncol(fit$internal$X_final), nrow(fit$internal$X_final)) # p > NT
	db <- debiasedATT(fit)
	expect_true(is.finite(db$att))
	expect_false(is.na(db$se))
	expect_gte(db$se, 0) # may be 0 at this degenerate corner; must not be NaN
})
