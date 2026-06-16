# Tests for the high-dimensional (full p >= NT) regime of the simultaneousCIs()
# multiplier bootstrap (#142 Phase 2). There the analytic Gram inverse need not
# exist, so the bootstrap uses the FULL-design desparsified construction
# (debiasedATT()'s #31 nodewise directions generalized to K effects -- uniformly
# valid, NOT post-selection). EXPERIMENTAL. The simulator can make full p > NT
# panels since #293 (small cohorts).

.make_highdim_sim_fit <- function() {
	# G=3, T=5, d=22 => full p = 390 > NT = 200 at N = 40.
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 22,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	dat <- simulateData(
		coefs,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	dat$indep_counts <- NA
	fetwfe(
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
}

hd_fit <- .make_highdim_sim_fit()

test_that("the high-dim fixture really is full p >= NT with valid SEs", {
	expect_gte(ncol(hd_fit$internal$X_final), nrow(hd_fit$internal$X_final))
	expect_true(isTRUE(hd_fit$internal$calc_ses))
})

test_that("the bootstrap runs in the high-dim regime and returns nodewise diagnostics", {
	for (fam in c("cohort", "all_post_treatment")) {
		bo <- simultaneousCIs(
			hd_fit,
			family = fam,
			method = "bootstrap",
			B = 1000,
			seed = 1
		)
		expect_s3_class(bo, "simultaneous_cis")
		expect_identical(bo$regime, "high-dimensional")
		expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
		expect_true(all(
			c("feasibility", "converged", "lambda_node") %in% names(bo)
		))
		# valid studentized sup-t: critical value in [pointwise, Bonferroni]
		expect_gte(bo$critical_value, bo$pointwise_critical_value - 1e-8)
		expect_lte(bo$critical_value, bo$bonferroni_critical_value + 1e-8)
	}
})

test_that("at the default lambda_c the nodewise directions converge and are feasible", {
	bo <- simultaneousCIs(
		hd_fit,
		family = "cohort",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_true(all(bo$converged))
	expect_true(all(bo$feasibility <= bo$lambda_node * (1 + 1e-9)))
})

test_that("too-small lambda_c yields infeasible directions and warns (experimental)", {
	expect_warning(
		simultaneousCIs(
			hd_fit,
			family = "cohort",
			method = "bootstrap",
			B = 500,
			seed = 1,
			lambda_c = 0.05
		),
		"feasibility|experimental"
	)
})

test_that("high-dim bootstrap is deterministic given a seed", {
	a <- simultaneousCIs(
		hd_fit,
		family = "cohort",
		method = "bootstrap",
		B = 800,
		seed = 7
	)
	b <- simultaneousCIs(
		hd_fit,
		family = "cohort",
		method = "bootstrap",
		B = 800,
		seed = 7
	)
	expect_identical(a$critical_value, b$critical_value)
})

test_that("event_study bootstrap uses the desparsified path in the high-dim regime (#299)", {
	# #299: a full p >= NT event_study fit routes through the desparsified
	# `targets` path (the regression-only guard was lifted), so it carries BOTH the
	# full-design desparsified regression channel and the propensity channel
	# `F_pi`, with `regime == "high-dimensional"` + nodewise diagnostics -- not the
	# fixed-p selected-support fallback it used before #299.
	bo <- simultaneousCIs(
		hd_fit,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(bo, "simultaneous_cis")
	expect_identical(bo$regime, "high-dimensional")
	# nodewise diagnostics attached (the desparsified-path signature)
	expect_true(all(
		c("feasibility", "converged", "lambda_node") %in% names(bo)
	))
	expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
})

test_that("high-dim event_study bootstrap is deterministic given a seed", {
	a <- simultaneousCIs(
		hd_fit,
		family = "event_study",
		method = "bootstrap",
		B = 800,
		seed = 11
	)
	b <- simultaneousCIs(
		hd_fit,
		family = "event_study",
		method = "bootstrap",
		B = 800,
		seed = 11
	)
	expect_identical(a$critical_value, b$critical_value)
	expect_identical(a$ci$estimate, b$ci$estimate)
	expect_identical(a$ci$simultaneous_ci_low, b$ci$simultaneous_ci_low)
})

test_that("too-small lambda_c warns for high-dim event_study (experimental)", {
	expect_warning(
		simultaneousCIs(
			hd_fit,
			family = "event_study",
			method = "bootstrap",
			B = 500,
			seed = 1,
			lambda_c = 0.05
		),
		"feasibility|experimental"
	)
})

# ------------------------------------------------------------------------------
# Debiasing anchor (#299, mutation-checkable cross-check vs the validated
# debiasedATT()). The high-dim band center is the DEBIASED estimate
# `estimates + colSums(F_reg)/(N*T)` (Theorem 6.6). For a 1 x num_treats `custom`
# contrast that matches debiasedATT()'s overall-ATT direction (cohort g loads its
# treat_inds cells with cohort_probs[g] / (#cells in cohort g)), that center must
# equal debiasedATT(fit)$att. RHS is the independently-validated single-effect
# debiased ATT; LHS is the simultaneousCIs high-dim center -- an independent path,
# not a round-trip. Mutation-checkable against the `colSums(F_mat)/(N*T)` line.
# ------------------------------------------------------------------------------
test_that("high-dim debiased band center equals debiasedATT()$att for the overall-ATT contrast", {
	G <- hd_fit$G
	T_ <- hd_fit$T
	treat_inds <- hd_fit$treat_inds
	num_treats <- length(treat_inds)
	cohort_probs <- hd_fit$cohort_probs
	# cohort of each treatment cell: T-1 cells for cohort 1, ..., T-G for cohort G
	cohort_of_treat <- rep(seq_len(G), times = (T_ - 1):(T_ - G))
	contrast <- numeric(num_treats)
	for (g in seq_len(G)) {
		idx <- which(cohort_of_treat == g)
		contrast[idx] <- cohort_probs[g] / length(idx)
	}
	# overall-ATT direction: weights sum to sum(cohort_probs[1:G])
	expect_equal(sum(contrast), sum(cohort_probs[seq_len(G)]))
	contrast_mat <- matrix(contrast, nrow = 1L)

	sc <- simultaneousCIs(
		hd_fit,
		family = "custom",
		contrasts = contrast_mat,
		method = "bootstrap",
		B = 200,
		seed = 1
	)
	expect_identical(sc$regime, "high-dimensional")
	db <- debiasedATT(hd_fit)
	# the load-bearing cross-check: equal to ~1e-9 (observed |diff| ~4e-16).
	expect_equal(sc$ci$estimate, db$att, tolerance = 1e-9)
})

# ------------------------------------------------------------------------------
# Band / adjusted-p duality under debiasing (#299, plan-review item 1). The
# adjusted-p `t_stat` MUST use the debiased center `estimates_used`, not the bridge
# `estimates`, or the band and its adjusted p-values disagree on a bridge-zeroed
# effect. Fixture: seed = 2 (same G=3,T=5,d=22,N=40 high-dim DGP), where the
# bridge zeroes SOME-but-not-all event-study effects (bridge e3 == 0 exactly, its
# debiased center is non-zero). The mutation-checkable signal: e3's adjusted
# p-value is finite and < 1 here (debiased |center|/se ~1.5); reverting the
# `t_stat` to abs(estimates) would make e3's t_stat abs(0)/se = 0 and its adjusted
# p ~1. NOTE: seed = 6 zeroes ALL effects and trips the upstream all-zero
# early-exit (R/simultaneous_cis.R:551), so it never reaches this path.
# ------------------------------------------------------------------------------
test_that("high-dim event_study band/adjusted-p duality holds with the debiased center (#299)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 22,
		density = 0.5,
		eff_size = 2,
		seed = 2
	)
	dat <- simulateData(
		coefs,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 2
	)
	dat$indep_counts <- NA
	fit2 <- fetwfe(
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
	expect_gte(ncol(fit2$internal$X_final), nrow(fit2$internal$X_final))
	# the fixture genuinely zeroes SOME (not all) treatment cells in the bridge,
	# so the path reaches the bootstrap AND has a bridge-zero effect.
	n_sel <- sum(fit2$internal$theta_hat[-1][fit2$treat_inds] != 0)
	expect_gt(n_sel, 0L)
	expect_lt(n_sel, length(fit2$treat_inds))

	bo <- simultaneousCIs(
		fit2,
		family = "event_study",
		method = "bootstrap",
		B = 4000,
		seed = 1
	)
	expect_identical(bo$regime, "high-dimensional")
	# bridge ES estimate for the last event time is exactly 0; the debiased center
	# is not -- the configuration that makes the t_stat source decisive.
	an <- simultaneousCIs(fit2, family = "event_study", method = "analytic")
	bridge_zero <- which(abs(an$ci$estimate) < 1e-10)
	expect_gt(length(bridge_zero), 0L)
	expect_true(any(abs(bo$ci$estimate[bridge_zero]) > 1e-6))

	# Exact dual: an effect's simultaneous band excludes 0 iff its adjusted p-value
	# is < alpha (uses the debiased center on BOTH sides).
	excl0 <- (bo$ci$simultaneous_ci_low > 0) | (bo$ci$simultaneous_ci_high < 0)
	sig <- bo$adjusted_p_values < bo$alpha
	expect_identical(excl0, sig)

	# Mutation anchor: a bridge-zero effect's adjusted p-value is driven by its
	# DEBIASED center, so it is < 1 (here ~0.28); under the bug (abs(estimates))
	# it would be ~1 (mean(boot_max >= 0)).
	expect_true(all(bo$adjusted_p_values[bridge_zero] < 0.95))
})
