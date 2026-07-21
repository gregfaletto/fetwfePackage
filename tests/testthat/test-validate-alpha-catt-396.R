# Tests for #396: (1) every CI accessor rejects an invalid `alpha` via the shared
# .validate_alpha_arg() guard -- eventStudy()/plot() formerly accepted it and
# silently returned inverted intervals; (2) plot(type = "catt") orders cohorts
# numerically, not alphabetically; (3) simultaneousCIs() validates its nodewise
# (riesz) controls like debiasedATT() does.

.fit_396 <- local({
	cf <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 1)
	sim <- simulateData(
		cf,
		N = 200,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	fetwfeWithSimulatedData(sim)
})

test_that("every CI accessor rejects an invalid alpha with the shared (0, 1) message (#396)", {
	msg <- "must be a single number in \\(0, 1\\)"
	# The two formerly-unguarded bug sites: eventStudy() and plot() (both views).
	expect_error(eventStudy(.fit_396, alpha = 1.5), msg)
	expect_error(eventStudy(.fit_396, alpha = NA_real_), msg)
	expect_error(plot(.fit_396, type = "catt", alpha = 1.5), msg)
	expect_error(plot(.fit_396, type = "event_study", alpha = 1.5), msg)
	# The three consolidated sites keep rejecting it (now via the shared helper).
	expect_error(cohortTimeATTs(.fit_396, alpha = 1.5), msg)
	expect_error(debiasedATT(.fit_396, alpha = 1.5), msg)
	expect_error(simultaneousCIs(.fit_396, alpha = 1.5), msg)
	# Boundary / type violations all rejected.
	expect_error(eventStudy(.fit_396, alpha = 0), msg)
	expect_error(eventStudy(.fit_396, alpha = 1), msg)
	expect_error(eventStudy(.fit_396, alpha = c(0.05, 0.1)), msg)
	# Valid alpha (default and explicit) is unaffected: the guards are no-ops.
	expect_s3_class(eventStudy(.fit_396), "eventStudy")
	expect_s3_class(eventStudy(.fit_396, alpha = 0.1), "eventStudy")
	expect_s3_class(cohortTimeATTs(.fit_396, alpha = 0.1), "cohortTimeATTs")
})

test_that("plot(type = 'catt') orders the cohort axis numerically, not alphabetically (#396)", {
	# >= 10 cohorts so the string order ("10" before "2") diverges from numeric.
	cf <- genCoefs(G = 10, T = 12, d = 1, density = 0.5, eff_size = 2, seed = 2)
	sim <- simulateData(
		cf,
		N = 600,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 2
	)
	fit <- fetwfeWithSimulatedData(sim)
	p <- plot(fit, type = "catt")
	xl <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()
	# Numeric order: "10"/"11" come AFTER "9", not after "1".
	expect_identical(xl, as.character(sort(as.numeric(xl))))
	expect_gt(match("10", xl), match("9", xl))
})

test_that("simultaneousCIs() validates its nodewise (riesz) controls (#396)", {
	expect_error(
		simultaneousCIs(.fit_396, riesz_max_iter = 0),
		"riesz_max_iter"
	)
	expect_error(
		simultaneousCIs(.fit_396, riesz_max_iter = c(1L, 2L)),
		"riesz_max_iter"
	)
	expect_error(simultaneousCIs(.fit_396, riesz_tol = -1), "riesz_tol")
	expect_error(simultaneousCIs(.fit_396, riesz_tol = NA_real_), "riesz_tol")
	# Valid controls unaffected.
	expect_s3_class(
		simultaneousCIs(.fit_396, riesz_max_iter = 5000L, riesz_tol = 1e-9),
		"simultaneous_cis"
	)
})
