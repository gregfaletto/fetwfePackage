library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Tests for covariate-dependent cohort-assignment DGPs introduced in 1.14.0
# (issue #162).
#
# Six tests:
#   1. Marginal DGP byte-identity (lock-in for backward compatibility).
#   2. Multinomial DGP: empirical cohort proportions match the MC expectation.
#   3. Ordered DGP: empirical cohort proportions are approximately uniform
#      (the root-find cutpoint scheme guarantees this at strength = 0 and
#      keeps it approximately uniform at higher strength).
#   4. FETWFE smoke test: fits run cleanly on multinomial and ordered DGPs
#      (well-formedness, not load-bearing for numerical correctness).
#   5. Propensity-weighted att_true under multinomial DGP (LOAD-BEARING for
#      the propensity-weighting code path; uses an INDEPENDENT MC reference
#      computed inside the test with a different seed).
#   6. High-strength sanity: covariates strongly predict cohort at
#      strength = 10 (max-row P(W|X) is well above the uniform baseline).
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Test 1: marginal DGP byte-identity
# ------------------------------------------------------------------------------
test_that("marginal DGP: byte-identical to v1.13.x at default assignment_type", {
	# Reference values captured from the marginal-path output at branch
	# `feat-cohort-dgp-162` against the same call on `origin/main` SHA
	# ae515b6. Locks in the pre-1.14.0 behavior at the strictest possible
	# assertion. Any drift in the marginal-path RNG ordering fires here.
	coefs <- genCoefs(
		R = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 42
	)
	expect_identical(coefs$assignment_type, "marginal")
	expect_null(coefs$assignment_coefs)

	sim <- simulateData(coefs, N = 100, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	expect_equal(nrow(sim$pdata), 100 * 5)
	expect_identical(sim$assignments, c(31L, 29L, 18L, 22L))
	expect_identical(sim$indep_counts, c(24L, 30L, 19L, 27L))
	expect_equal(
		sim$pdata$y[1:3],
		c(1.249801325615, -3.167719220623, -3.455029448642),
		tolerance = 1e-9
	)
})

# ------------------------------------------------------------------------------
# Test 2: multinomial DGP — empirical cohort proportions match MC expectation
# ------------------------------------------------------------------------------
test_that("multinomial DGP: empirical cohort proportions converge to MC expectation", {
	skip_on_cran()
	# N = 10000 gives sqrt(p(1-p)/N) SE < 0.005 per cohort; tolerance 0.03
	# is comfortable.
	coefs <- genCoefs(
		R = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 1.0,
		seed = 42
	)
	sim <- simulateData(coefs, N = 10000, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	empirical_probs <- sim$assignments / sum(sim$assignments)
	expected_probs <- fetwfe:::.expected_cohort_probs(
		coefs$assignment_coefs,
		d = 2,
		distribution = "gaussian",
		M = 50000L,
		seed = 99L
	)
	expect_equal(
		as.numeric(empirical_probs),
		as.numeric(expected_probs),
		tolerance = 0.03
	)
})

# ------------------------------------------------------------------------------
# Test 3: ordered DGP — empirical cohort proportions ~ uniform 1/(R+1)
# ------------------------------------------------------------------------------
test_that("ordered DGP: empirical cohort proportions match cutpoint structure", {
	skip_on_cran()
	# The root-find cutpoint scheme produces approximately uniform 1/(R+1)
	# marginal probabilities at any strength (exactly uniform at
	# strength = 0; ~0.002 deviation at strength = 1 per round-2 plan
	# verification).
	coefs <- genCoefs(
		R = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "ordered",
		assignment_strength = 1.0,
		seed = 42
	)
	sim <- simulateData(coefs, N = 10000, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	empirical_probs <- sim$assignments / sum(sim$assignments)
	expect_equal(
		as.numeric(empirical_probs),
		rep(1 / 4, 4),
		tolerance = 0.03
	)
})

# ------------------------------------------------------------------------------
# Test 4: FETWFE smoke test on multinomial + ordered DGPs
# ------------------------------------------------------------------------------
test_that("FETWFE smoke test: fits run cleanly on multinomial and ordered DGPs", {
	skip_on_cran()
	# Well-formedness only (per round-1 sentinel NOTE: not load-bearing
	# for numerical correctness; that's covered by test #5 below).
	for (type in c("multinomial", "ordered")) {
		coefs <- genCoefs(
			R = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			assignment_type = type,
			assignment_strength = 1.0,
			seed = 42
		)
		sim <- simulateData(
			coefs,
			N = 200,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5
		)
		fit <- suppressWarnings(fetwfe(
			pdata = sim$pdata,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			response = "y",
			covs = c("cov1", "cov2")
		))
		expect_true(
			is.finite(fit$att_hat),
			info = paste0("att_hat not finite for type = ", type)
		)
		expect_true(
			is.finite(fit$att_se),
			info = paste0("att_se not finite for type = ", type)
		)
		expect_gt(fit$att_se, 0)
	}
})

# ------------------------------------------------------------------------------
# Test 5: propensity-weighted att_true under multinomial DGP (LOAD-BEARING)
# ------------------------------------------------------------------------------
test_that("getTes returns propensity-weighted att_true under multinomial DGP", {
	# Compares getTes()'s production output against an INDEPENDENT MC
	# integration computed inside the test with a different seed
	# (99L vs the production-side coefs$seed + 2L = 44L). This is the
	# load-bearing check for the propensity-weighting code path; it
	# verifies the production code does compute
	#   att_true = sum(cohort_weights * actual_cohort_tes)
	# with cohort_weights = E[pi_r(X)] / sum E[pi_r(X)] (treated only).
	coefs <- genCoefs(
		R = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 1.5,
		seed = 42
	)
	tes <- getTes(coefs)
	# Independent reference computed with a fresh MC sample under a
	# different seed.
	set.seed(99L)
	X_ref <- matrix(rnorm(50000 * 2), nrow = 50000, ncol = 2)
	probs_ref <- fetwfe:::.compute_cohort_prob_matrix(
		X_ref,
		coefs$assignment_coefs
	)
	expected_ref <- colMeans(probs_ref)
	treated_probs_ref <- expected_ref[2:4]
	cohort_weights_ref <- as.numeric(
		treated_probs_ref / sum(treated_probs_ref)
	)
	att_ref <- sum(cohort_weights_ref * tes$actual_cohort_tes)
	expect_equal(tes$att_true, att_ref, tolerance = 0.02)
	expect_equal(
		tes$cohort_weights,
		cohort_weights_ref,
		tolerance = 0.02
	)
	# Under non-marginal DGP, cohort_weights should NOT be uniform
	# (otherwise we'd not be exercising the propensity-weighting path).
	expect_false(
		isTRUE(all.equal(tes$cohort_weights, rep(1 / 3, 3), tolerance = 0.01))
	)
})

# ------------------------------------------------------------------------------
# Test 6: high-strength sanity — covariates strongly predict cohort
# ------------------------------------------------------------------------------
test_that("high-strength assignment: covariates strongly predict cohort", {
	# At assignment_strength = 10, cohort assignment becomes nearly
	# deterministic given X. Verify: max-row-probability of P(W|X) is
	# well above the uniform-baseline 1/(R+1) = 0.25 at most units.
	coefs <- genCoefs(
		R = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 10.0,
		seed = 42
	)
	X_test <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
	probs <- fetwfe:::.compute_cohort_prob_matrix(
		X_test,
		coefs$assignment_coefs
	)
	max_row_probs <- apply(probs, 1, max)
	expect_gt(mean(max_row_probs), 0.85)
})
