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

test_that("event_study is still rejected for the bootstrap in the high-dim regime", {
	expect_error(
		simultaneousCIs(hd_fit, family = "event_study", method = "bootstrap"),
		"event_study"
	)
})
