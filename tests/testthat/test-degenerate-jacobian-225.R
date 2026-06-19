library(testthat)
library(fetwfe)

# Issue #225 (remaining items, deferred from #227): two regression tests for the
# recently-added inference surface that previously had no isolated coverage.
#
#   Item 1 -- the K=0 "degenerate support" branch of .simultaneous_cis_impl()
#     (R/simultaneous_cis.R:405-432): when bridge selection zeroes EVERY treatment
#     effect, simultaneousCIs() returns collapsed point bands at zero, with the
#     pointwise critical value and NA adjusted p-values. No test reached it.
#
#   Item 2 -- the K=1 second-variance ("Sigma_2") term computed two ways: the
#     joint .assemble_joint_cov_var2() (R/variance_machinery.R:1259, used by the
#     simultaneous-CI worker) vs the scalar getSecondVarTermOLS()
#     (R/variance_machinery.R:251, used by eventStudy()). They are the same
#     quantity in different parameterizations; nothing cross-checked them.
#
# Both are test-only drift sentinels; the shipped behavior is correct. Each is
# verified to FAIL on a deliberately-broken version of the guarded path.

# ---------------------------------------------------------------------------
# Item 1: K=0 degenerate-support branch.
#
# Trigger: a real fit with NO true treatment effect (eff_size = 0) + high noise,
# so bridge selection (q = 0.5) zeroes every treatment effect while still
# computing standard errors (calc_ses = TRUE) -- reaching the empty-support
# branch rather than the earlier calc_ses = FALSE stop().
#
# NOTE: q = 0.5 is load-bearing. With q >= 1 this fixture singularizes the Gram
# on the (degenerate) support, giving calc_ses = FALSE and the earlier stop() --
# NOT the K=0 branch. Do not "simplify" q here.
# ---------------------------------------------------------------------------
.deg_sim_225 <- local({
	sc <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 0, seed = 1)
	simulateData(sc, N = 60, sig_eps_sq = 5, sig_eps_c_sq = 1, seed = 1)
})

.deg_fit_225 <- suppressMessages(fetwfe(
	pdata = .deg_sim_225$pdata,
	time_var = .deg_sim_225$time_var,
	unit_var = .deg_sim_225$unit_var,
	treatment = .deg_sim_225$treatment,
	covs = .deg_sim_225$covs,
	response = .deg_sim_225$response,
	sig_eps_sq = .deg_sim_225$sig_eps_sq,
	sig_eps_c_sq = .deg_sim_225$sig_eps_c_sq,
	q = 0.5,
	verbose = FALSE
))

test_that("simultaneousCIs() K=0 degenerate branch returns collapsed point bands (#225)", {
	skip_if_not_installed("mvtnorm")

	# Fixture invariant: selection really did zero every treatment effect, and SEs
	# were still requested. If a future penalty-default change keeps a nonzero
	# effect, these turn silent rot into a loud failure.
	expect_true(all(.deg_fit_225$beta_hat[.deg_fit_225$treat_inds] == 0))
	expect_true(.deg_fit_225$calc_ses)

	G <- .deg_fit_225$G
	T_ <- .deg_fit_225$T
	num_treats <- length(.deg_fit_225$treat_inds)
	expected_K <- c(
		event_study = T_ - 1L,
		cohort = G,
		all_post_treatment = num_treats
	)

	for (fam in c("event_study", "cohort", "all_post_treatment")) {
		expect_message(
			sci <- simultaneousCIs(.deg_fit_225, family = fam),
			"zeroed out every treatment effect"
		)
		expect_s3_class(sci, "simultaneous_cis")
		expect_identical(sci$K, as.integer(expected_K[[fam]]))

		# Every estimate is zero (the support is empty).
		expect_true(all(sci$ci$estimate == 0))

		# The critical value falls back to the pointwise value (NOT Bonferroni).
		expect_equal(sci$critical_value, sci$pointwise_critical_value)
		expect_equal(sci$critical_value, qnorm(1 - sci$alpha / 2))

		# All four bands collapse to a point at the (zero) estimate.
		expect_equal(sci$ci$simultaneous_ci_low, sci$ci$estimate)
		expect_equal(sci$ci$simultaneous_ci_high, sci$ci$estimate)
		expect_equal(sci$ci$pointwise_ci_low, sci$ci$estimate)
		expect_equal(sci$ci$pointwise_ci_high, sci$ci$estimate)

		# Degenerate => no adjusted (max-T) p-values.
		expect_true(all(is.na(sci$adjusted_p_values)))
	}
})

# ---------------------------------------------------------------------------
# Item 2: K=1 Sigma_2 cross-parameterization identity.
#
# For a single event time e (valid cohort set V_e, |V_e| > 1), the second
# variance term is computed two ways that must agree exactly:
#   PATH A getSecondVarTermOLS(): scalar, via a psi_e_mat selector + V_e-masked
#          cohort probabilities (the eventStudy() construction).
#   PATH B .assemble_joint_cov_var2(): the [1,1] entry of the K=1 joint Term-2,
#          via a per_effect_masked Jacobian (the simultaneous-CI worker's path).
# Use the OLS family (etwfe): theta_sel = tes and d_inv = identity, so the two
# parameterizations map transparently. The identity is exact (not asymptotic).
# ---------------------------------------------------------------------------
.var2_sim_225 <- local({
	sc <- genCoefs(G = 3, T = 5, d = 2, density = 0.6, eff_size = 2, seed = 11)
	simulateData(sc, N = 150, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 11)
})

.var2_fit_225 <- suppressMessages(etwfe(
	pdata = .var2_sim_225$pdata,
	time_var = .var2_sim_225$time_var,
	unit_var = .var2_sim_225$unit_var,
	treatment = .var2_sim_225$treatment,
	covs = .var2_sim_225$covs,
	response = .var2_sim_225$response,
	sig_eps_sq = .var2_sim_225$sig_eps_sq,
	sig_eps_c_sq = .var2_sim_225$sig_eps_c_sq,
	verbose = FALSE
))

test_that(".assemble_joint_cov_var2 (K=1) matches getSecondVarTermOLS on an event time (#225)", {
	fit <- .var2_fit_225
	G <- fit$G
	T_ <- fit$T
	N <- fit$N
	num_treats <- length(fit$treat_inds)
	tes <- fit$beta_hat[fit$treat_inds]
	cpo <- fit$cohort_probs_overall

	offs <- fetwfe:::.resolve_cohort_offsets_and_first_inds(
		fit,
		G = G,
		T = T_
	)
	first_inds <- offs$first_inds
	coi <- offs$cohort_offsets_int
	Sigma_pi_hat <- fetwfe:::.multinomial_cov(cpo[1:G])

	# Event time e = 0 pools all G cohorts here; assert |V_e| > 1 so the cross
	# check can't reduce to a vacuous 0 == 0 (single-cohort event times give 0).
	e <- 0L
	V_e <- which(coi <= T_ - e)
	expect_gt(length(V_e), 1L)

	# PATH A: scalar Term-2 via getSecondVarTermOLS (the eventStudy() construction).
	psi_e_mat <- matrix(0, nrow = num_treats, ncol = G)
	for (g in V_e) {
		psi_e_mat[first_inds[g] + e, g] <- 1
	}
	masked_probs <- numeric(G)
	masked_probs[V_e] <- cpo[V_e]
	v2_scalar <- fetwfe:::getSecondVarTermOLS(
		psi_mat = psi_e_mat,
		tes = tes,
		cohort_probs_overall = masked_probs,
		num_treats = num_treats,
		N = N,
		T = T_,
		G = G
	)

	# PATH B: K=1 diagonal of the joint Term-2 (the simultaneous-CI worker's path).
	J_e <- fetwfe:::.build_jacobian(
		cohort_probs_overall = cpo,
		G = G,
		d_inv_treat_sel = diag(num_treats),
		mode = "per_effect_masked",
		V_e = V_e,
		first_inds = first_inds,
		e = e
	)
	v2_joint <- fetwfe:::.assemble_joint_cov_var2(
		J_list = list(J_e),
		theta_sel = tes,
		Sigma_pi_hat = Sigma_pi_hat,
		N = N,
		T = T_
	)

	expect_equal(dim(v2_joint), c(1L, 1L))
	# Exact structural identity (default tolerance; not bit-identity across BLAS).
	expect_equal(v2_joint[1, 1], v2_scalar)
	# Guard: this event time genuinely has a positive second-variance term.
	expect_gt(v2_scalar, 0)
})
