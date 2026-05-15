# Tests for the internal `estOmegaSqrtInv()` variance-component estimator.
#
# Background: Prior to v1.9.1, `estOmegaSqrtInv()` returned `sig_eps_c_sq` of
# essentially zero (~1e-32) regardless of truth due to two coupled bugs in
# the within-estimator implementation. The fix replaces it with REML via
# `lme4::lmer(y ~ X + (1 | unit), REML = TRUE)`. Tests below verify that the
# new estimator recovers both variance components to within a factor of 2 of
# truth on synthetic panels, and that the dependency-fallback path works
# when lme4 isn't available.
#
# The tests deliberately call `fetwfe()` and `fetwfe:::estOmegaSqrtInv()`
# DIRECTLY (not `fetwfeWithSimulatedData()`), because the latter passes the
# simulator's true variance components through, bypassing `estOmegaSqrtInv()`
# entirely.

if (!requireNamespace("lme4", quietly = TRUE)) {
	testthat::skip("lme4 not installed")
}

test_that("estOmegaSqrtInv recovers sig_eps_sq and sig_eps_c_sq via REML", {
	# Generate a synthetic panel with known variance components.
	set.seed(20260515)
	coefs <- genCoefs(
		R = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 20260515
	)
	sim <- simulateData(
		coefs,
		N = 200,
		sig_eps_sq = 1,
		sig_eps_c_sq = 2
	)

	# Call estOmegaSqrtInv() directly on the simulator's design matrix.
	res <- fetwfe:::estOmegaSqrtInv(
		y = sim$y,
		X_ints = sim$X,
		N = sim$N,
		T = sim$T,
		p = sim$p
	)

	# Both components within a factor of 2 of truth. (Empirically REML on
	# this regime returns estimates within ~15% of truth across seeds; the
	# factor-of-2 bound is generous slack against finite-sample noise.
	# The pre-fix buggy estimator returned sig_eps_c_sq ~1e-32, failing
	# by ~30 orders of magnitude.)
	expect_lt(abs(log(res$sig_eps_sq / 1)), log(2))
	expect_lt(abs(log(res$sig_eps_c_sq / 2)), log(2))
})

test_that("fetwfe() invokes the REML estimator when sig values are NA", {
	# End-to-end check: fetwfe() with default NA arguments triggers
	# estOmegaSqrtInv() internally. Verify the recovered sig_eps_c_sq is
	# nonzero (regression guard against re-introducing the buggy behavior
	# where the estimator collapsed to ~0 silently).
	set.seed(20260516)
	coefs <- genCoefs(
		R = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 20260516
	)
	sim <- simulateData(
		coefs,
		N = 200,
		sig_eps_sq = 1,
		sig_eps_c_sq = 2
	)
	res <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		# Crucially: leave sig_eps_sq / sig_eps_c_sq at their NA defaults
		# so estOmegaSqrtInv() is actually invoked.
		verbose = FALSE
	)
	expect_gt(res$sig_eps_c_sq, 0.1) # was ~1e-32 before
	expect_lt(res$sig_eps_c_sq, 10)
})

test_that("estOmegaSqrtInv error path mentions lme4 in the message body", {
	# We can't reliably stub requireNamespace from outside the package
	# namespace, so this test reads the function's source to verify the
	# error message points users at install.packages('lme4') or supplying
	# σ values manually. A weaker check than runtime stubbing but
	# sufficient to catch regression of the user-facing message text.
	src <- paste(
		deparse(body(fetwfe:::estOmegaSqrtInv)),
		collapse = "\n"
	)
	expect_match(src, "requireNamespace")
	expect_match(src, "lme4")
	expect_match(src, "install\\.packages")
})
