# Regression test for #324: simultaneousCIs() inherits the fit's `alpha`.
#
# Before the fix, the accessor hard-coded `alpha = 0.05`, so a fit built at a
# non-default level (e.g. 0.10) showed 90% `catt_df` simultaneous bands, yet
# `simultaneousCIs(fit)` returned 95% bands -- breaking the documented identity
# ("the `catt_df` bounds ARE the `simultaneousCIs(fit, family=...)` bounds", the
# `ci_type` docs) and diverging from every other accessor (`eventStudy()` already
# resolves the level via `x$alpha`). The fix defaults `alpha = NULL` across the
# generic + 4 S3 methods and resolves `if (is.null(alpha)) alpha <- .alpha_of(x)`
# once in `.simultaneous_cis_impl()`, before validation.
#
# Mutation-checkable: reverting the signature default to `0.05` (or deleting the
# `if (is.null(alpha)) alpha <- .alpha_of(x)` resolution) makes the inherited
# call fall back to 0.05, failing the `sci$alpha == 0.10` and catt_df-identity
# assertions below. All four S3 methods are byte-identical wrappers over the same
# impl, so fetwfe + twfeCovs cover the fusion and OLS dispatch paths.

make_alpha_panel <- function(seed = 101, G = 4, T = 6, d = 2, N = 200) {
	coefs <- genCoefs(
		G = G,
		T = T,
		d = d,
		density = 0.5,
		eff_size = 2,
		seed = seed
	)
	simulateData(coefs, N = N, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = seed)
}

test_that("simultaneousCIs() inherits the fit's alpha when not supplied (#324)", {
	sim <- make_alpha_panel()

	check_inherits <- function(fit, label) {
		# fixture invariant: the fit was built at the non-default alpha = 0.10.
		expect_equal(fit$alpha, 0.10, info = label)

		# (1) the inherited call (no `alpha`) uses the fit's 0.10, NOT the old
		#     hard-coded 0.05.
		sci <- simultaneousCIs(fit, family = "cohort") # method = "analytic" default
		expect_equal(sci$alpha, 0.10, info = label)

		# (2) documented identity: the inherited bounds ARE the fit's `catt_df`
		#     simultaneous bounds (positional, cohort-block order; degenerate rows
		#     match by construction since catt_df's bounds come from this band).
		expect_equal(
			sci$ci$simultaneous_ci_low,
			fit$catt_df$ci_low,
			info = label
		)
		expect_equal(
			sci$ci$simultaneous_ci_high,
			fit$catt_df$ci_high,
			info = label
		)

		# (3) inherited == explicit-0.10 (the inheritance is exactly `x$alpha`).
		sci_explicit <- simultaneousCIs(fit, family = "cohort", alpha = 0.10)
		expect_equal(
			sci$ci$simultaneous_ci_low,
			sci_explicit$ci$simultaneous_ci_low,
			info = label
		)

		# (4) an explicit alpha still overrides -> the 90% band is strictly
		#     narrower than the 95% band on finite-SE cohorts.
		sci05 <- simultaneousCIs(fit, family = "cohort", alpha = 0.05)
		expect_equal(sci05$alpha, 0.05, info = label)
		fin <- is.finite(fit$catt_df$se) & fit$catt_df$se > 0
		skip_if_not(sum(fin) >= 1L)
		w10 <- (sci$ci$simultaneous_ci_high - sci$ci$simultaneous_ci_low)[fin]
		w05 <- (sci05$ci$simultaneous_ci_high - sci05$ci$simultaneous_ci_low)[
			fin
		]
		expect_true(all(w10 < w05), info = label)
	}

	check_inherits(fetwfeWithSimulatedData(sim, alpha = 0.10), "fetwfe")
	check_inherits(twfeCovsWithSimulatedData(sim, alpha = 0.10), "twfeCovs")
})
