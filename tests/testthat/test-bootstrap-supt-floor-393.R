# Tests for #393: the K >= 2 sup-t multiplier-bootstrap critical value AND its
# dual single-step max-T adjusted p-value are floored at the pointwise Gaussian
# value, so `simultaneousCIs(method = "bootstrap")` is never anti-conservative
# relative to its own pointwise band/p-values.
#
# Two floors, two failure modes (mutation-proven to bite separately):
#   (1) band crit floor  -> `.simultaneous_bootstrap_crit()` (crit >= z);
#   (2) adjusted-p floor -> the caller (adjusted_p >= pointwise two-sided p).
# The raw Monte-Carlo estimate of each dips below the pointwise value under
# strong cross-effect correlation (a sup over K >= 1 standardized effects
# dominates a single |Z| in the population, but the B-draw estimate need not).

test_that(".simultaneous_bootstrap_crit() floors the K>=2 sup-t crit at qnorm(1-alpha/2) (#393)", {
	alpha <- 0.05
	z <- stats::qnorm(1 - alpha / 2)
	N <- 200L
	# Three near-perfectly correlated IF columns (common base + tiny noise): the
	# sup-t statistic collapses to ~|Z|, so its 1-alpha MC quantile fluctuates
	# around z and dips below it on a sizeable fraction (~44%) of bootstrap draws.
	set.seed(393)
	v <- stats::rnorm(N)
	F_mat <- sapply(1:3, function(k) v + 0.01 * stats::rnorm(N))
	# Loop bootstrap seeds: which seeds dip is BLAS/RNG-dependent, so a single
	# pinned seed is fragile; looping guarantees the floor is exercised (>= 1
	# dipping seed) on any platform. Every returned crit must be >= z.
	crits <- vapply(
		1:8,
		function(s) {
			fetwfe:::.simultaneous_bootstrap_crit(
				F_mat,
				n = N,
				alpha = alpha,
				B = 1000,
				multiplier = "rademacher",
				seed = s
			)$crit
		},
		numeric(1)
	)
	expect_true(all(crits >= z))
})

test_that("simultaneousCIs(method='bootstrap') is never anti-conservative vs pointwise (#393)", {
	# A plain 2-cohort fit whose bootstrap adjusted p-values dip below the
	# pointwise p on ~75% of seeds (the reachable-in-normal-use half of #393; the
	# band crit itself rarely dips on real fits, so the crit assertion here is an
	# invariant guard, and the direct test above is the crit floor's biting test).
	coefs <- genCoefs(
		G = 2,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 5
	)
	sim <- simulateData(
		coefs,
		N = 200,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 5
	)
	fit <- fetwfeWithSimulatedData(sim)
	alpha <- fit$alpha
	for (bs in 1:8) {
		sc <- suppressWarnings(simultaneousCIs(
			fit,
			family = "cohort",
			method = "bootstrap",
			B = 1000,
			seed = bs
		))
		est <- sc$ci$estimate
		# Recover the per-effect SE exactly from the simultaneous half-width.
		se <- (sc$ci$simultaneous_ci_high - est) / sc$critical_value
		pw_p <- 2 * stats::pnorm(-abs(est / se))
		# (1) Band: the simultaneous crit / band never inside the pointwise one.
		expect_gte(sc$critical_value, sc$pointwise_critical_value)
		expect_true(all(
			sc$ci$simultaneous_ci_low <= sc$ci$pointwise_ci_low + 1e-10
		))
		expect_true(all(
			sc$ci$simultaneous_ci_high >= sc$ci$pointwise_ci_high - 1e-10
		))
		# (2) p-value: the adjusted (simultaneous, family-wise) p never below the
		# pointwise two-sided p. This is the assertion the p-value floor restores;
		# it fails pre-fix on the dipping seeds.
		ok <- is.na(sc$adjusted_p_values) |
			(sc$adjusted_p_values >= pw_p - 1e-10)
		expect_true(all(ok))
	}
})
