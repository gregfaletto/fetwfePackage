# Tests for #398: getTes() integrates the propensity-weighted truth over the
# covariate distribution given by its new `distribution` argument, instead of a
# hardcoded Gaussian. Only bites covariate-dependent assignment (uniform vs
# gaussian X give different E[pi_g(X)]); marginal assignment is unaffected.

# A covariate-dependent (multinomial) DGP whose gaussian and uniform propensity
# integrals differ; and a marginal DGP where they cannot.
.cf_cov <- genCoefs(
	G = 3,
	T = 5,
	d = 2,
	density = 0.5,
	eff_size = 2,
	assignment_type = "multinomial",
	assignment_strength = 2,
	seed = 11
)
.cf_marg <- genCoefs(
	G = 3,
	T = 5,
	d = 2,
	density = 0.5,
	eff_size = 2,
	seed = 11
)

test_that("getTes() `distribution` arg contract: default back-compat, effect, marginal no-op, validation (#398)", {
	# Default is byte-identical to the pre-#398 behavior (hardcoded gaussian).
	expect_identical(getTes(.cf_cov), getTes(.cf_cov, "gaussian"))

	# The argument has an effect under covariate-dependent assignment.
	g <- getTes(.cf_cov, "gaussian")
	u <- getTes(.cf_cov, "uniform")
	expect_false(isTRUE(all.equal(u$cohort_weights, g$cohort_weights)))
	expect_false(isTRUE(all.equal(u$att_true, g$att_true)))

	# Under marginal assignment the truth does not integrate over X, so the
	# argument is ignored (cohort weights are uniform 1/G either way).
	expect_identical(getTes(.cf_marg, "uniform"), getTes(.cf_marg, "gaussian"))

	# Unsupported values error at the front door (match.arg).
	expect_error(getTes(.cf_cov, "poisson"))
})

test_that("getTes(distribution = 'uniform') matches the uniform-X sampling it claims to describe (#398)", {
	# End-to-end: the propensity-weighted truth getTes() reports for a uniform
	# panel must match the cohort frequencies a large uniform-X draw actually
	# produces -- which the pre-#398 gaussian-integrated truth does NOT.
	sim <- simulateData(
		.cf_cov,
		N = 30000,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		distribution = "uniform",
		seed = 7
	)
	# sim$assignments is the length-(G+1) per-cohort count in the package's
	# never-treated-first order: index 1 is the never-treated count (guarded
	# below), indices 2..(G+1) are the treated cohorts. Use assignments[-1] --
	# tie-proof, unlike matching on the count value.
	expect_identical(sim$assignments[1], sim$N_UNTREATED)
	emp_treated <- sim$assignments[-1] / sum(sim$assignments[-1])

	wu <- getTes(.cf_cov, "uniform")$cohort_weights
	wg <- getTes(.cf_cov, "gaussian")$cohort_weights
	l1_u <- sum(abs(emp_treated - wu))
	l1_g <- sum(abs(emp_treated - wg))

	# The uniform truth is much closer to the uniform sampling than the gaussian
	# truth is (~0.006 vs ~0.033 here), and matches within MC tolerance.
	expect_lt(l1_u, l1_g)
	expect_lt(l1_u, 0.02)
})
