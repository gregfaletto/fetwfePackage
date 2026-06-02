# Deprecation-contract tests for the user-facing cohort-count rename
# R -> G (#41). Locks: (a) the deprecated `R` argument warns AND maps to
# `G` (identical result); (b) the canonical `G` path is warning-free;
# (c) both-supplied-conflict errors, both-equal warns-not-errors;
# (d) the dual-populated `$G` field equals `$R` on every classed object;
# (e) the package's own pipeline never self-triggers the deprecation
# warning. See PLAN.md Edit 10.

# ------------------------------------------------------------------------------
# Argument deprecation: genCoefs / genCoefsCore / simulateDataCore
# ------------------------------------------------------------------------------

test_that("genCoefs(R = ) warns about deprecation", {
	expect_warning(
		genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 1, seed = 1),
		"deprecated"
	)
})

test_that("genCoefs(G = ) emits no deprecation warning", {
	expect_no_warning(
		genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 1, seed = 1)
	)
})

test_that("genCoefs(R = k) equals genCoefs(G = k)", {
	g <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 1, seed = 1)
	r <- suppressWarnings(
		genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 1, seed = 1)
	)
	expect_identical(r$beta, g$beta)
	expect_identical(r$theta, g$theta)
	expect_identical(r$G, g$G)
	expect_identical(r$R, g$R)
})

test_that("genCoefs positional first arg binds to G (no warning)", {
	g <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 1, seed = 1)
	expect_no_warning(
		p <- genCoefs(3, 6, 2, 0.5, 1, seed = 1)
	)
	expect_identical(p$beta, g$beta)
})

test_that("genCoefs errors when G and R disagree", {
	expect_error(
		genCoefs(G = 3, R = 4, T = 6, d = 2, density = 0.5, eff_size = 1),
		"only one"
	)
})

test_that("genCoefs accepts both G and R when they agree (still warns)", {
	expect_warning(
		coefs <- genCoefs(
			G = 3,
			R = 3,
			T = 6,
			d = 2,
			density = 0.5,
			eff_size = 1,
			seed = 1
		),
		"deprecated"
	)
	expect_equal(coefs$G, 3)
})

test_that("genCoefs errors when the cohort count is missing entirely", {
	expect_error(
		genCoefs(T = 6, d = 2, density = 0.5, eff_size = 1),
		"requires the cohort count"
	)
})

test_that("genCoefsCore(R = ) warns and equals genCoefsCore(G = )", {
	g <- genCoefsCore(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 1,
		seed = 1
	)
	expect_warning(
		r <- genCoefsCore(
			R = 3,
			T = 6,
			d = 2,
			density = 0.5,
			eff_size = 1,
			seed = 1
		),
		"deprecated"
	)
	expect_identical(r$beta, g$beta)
})

test_that("simulateDataCore(R = ) warns and equals simulateDataCore(G = )", {
	cf <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 1, seed = 1)
	g <- simulateDataCore(
		N = 50,
		T = 6,
		G = 3,
		d = 2,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		beta = cf$beta,
		seed = 7
	)
	expect_warning(
		r <- simulateDataCore(
			N = 50,
			T = 6,
			R = 3,
			d = 2,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			beta = cf$beta,
			seed = 7
		),
		"deprecated"
	)
	expect_identical(r$pdata, g$pdata)
	expect_identical(r$G, g$G)
})

# ------------------------------------------------------------------------------
# Dual-populated $G field on every classed return ($G == $R, $R still works)
# ------------------------------------------------------------------------------

test_that("$G equals $R on genCoefs / simulateData / getTes objects", {
	cf <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 1, seed = 1)
	expect_equal(cf$G, cf$R)
	expect_equal(cf$G, 3)

	sim <- simulateData(cf, N = 80, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	expect_equal(sim$G, sim$R)
	expect_equal(sim$G, 3)

	tes <- getTes(cf)
	expect_equal(tes$G, tes$R)
	expect_equal(tes$G, 3)
})

test_that("$G equals $R on a fit from each of the four estimators", {
	cf <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 1, seed = 1)
	sim <- simulateData(cf, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)

	fit_fetwfe <- fetwfeWithSimulatedData(sim)
	expect_equal(fit_fetwfe$G, fit_fetwfe$R)

	fit_etwfe <- etwfeWithSimulatedData(sim)
	expect_equal(fit_etwfe$G, fit_etwfe$R)

	fit_betwfe <- betwfeWithSimulatedData(sim)
	expect_equal(fit_betwfe$G, fit_betwfe$R)

	fit_twfeCovs <- twfeCovsWithSimulatedData(sim)
	expect_equal(fit_twfeCovs$G, fit_twfeCovs$R)
})

test_that("ci_type = 'simultaneous' fit exposes $G == $R (cross-feature with #199)", {
	cf <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 1, seed = 1)
	sim <- simulateData(cf, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	fit <- fetwfeWithSimulatedData(sim, ci_type = "simultaneous")
	expect_equal(fit$G, fit$R)
	expect_equal(fit$ci_type, "simultaneous")
})

# ------------------------------------------------------------------------------
# The package never self-triggers its own deprecation warning
# ------------------------------------------------------------------------------

test_that("a representative pipeline emits no deprecation warning", {
	withr::local_options(warn = 1)
	expect_no_condition(
		{
			cf <- genCoefs(
				G = 3,
				T = 6,
				d = 2,
				density = 0.5,
				eff_size = 1,
				seed = 1
			)
			sim <- simulateData(cf, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
			fit <- fetwfeWithSimulatedData(sim)
			invisible(capture.output({
				print(fit)
				print(summary(fit))
			}))
			getTes(cf)
			broom::tidy(fit)
			eventStudy(fit)
		},
		class = "warning"
	)
})
