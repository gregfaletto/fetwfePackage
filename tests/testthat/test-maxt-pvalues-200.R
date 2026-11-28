# Tests for issue #200: single-step max-T multiplicity-adjusted p-values, the
# exact dual of the simultaneous CI band. Under ci_type = "simultaneous" (the
# default) the catt_df / eventStudy() `p_value` is the family-wise max-T
# adjusted p-value; under "pointwise" it stays the per-effect Wald p-value.

# ---------------------------------------------------------------------------
# 1. The core helper `.maxt_adjusted_p_nd` matches an independent pmvnorm
#    recomputation, and the adjustment strictly inflates the pointwise p.
# ---------------------------------------------------------------------------
test_that(".maxt_adjusted_p_nd equals an independent pmvnorm recomputation", {
	skip_if_not_installed("mvtnorm")
	K <- 5L
	rho <- matrix(0.45, K, K)
	diag(rho) <- 1
	estimates <- c(2.5, 1.8, 0.9, 3.1, 2.0)
	ses <- rep(1, K)

	set.seed(1)
	adj <- .maxt_adjusted_p_nd(estimates, ses, rep(TRUE, K), rho)

	set.seed(1)
	indep <- vapply(
		seq_len(K),
		function(k) {
			z <- abs(estimates[k] / ses[k])
			1 -
				as.numeric(mvtnorm::pmvnorm(
					lower = rep(-z, K),
					upper = rep(z, K),
					corr = rho,
					algorithm = mvtnorm::GenzBretz()
				))
		},
		numeric(1)
	)

	expect_equal(adj, indep, tolerance = 1e-12)

	# Each value lies between its pointwise p-value and its Bonferroni bound
	# (#551), and the adjustment is strictly visible for at least one effect.
	pw <- 2 * stats::pnorm(-abs(estimates / ses))
	expect_true(all(adj >= pw))
	expect_true(all(adj <= pmin(1, K * pw)))
	expect_true(any(adj > pw + 1e-4))
})

# ---------------------------------------------------------------------------
# 2. The adjusted p-value is the exact dual of the qmvnorm band:
#    z_k > crit  <=>  adjusted_p_k < alpha (away from the MC-noise boundary).
# ---------------------------------------------------------------------------
test_that(".maxt_adjusted_p_nd is the exact dual of the qmvnorm critical value", {
	skip_if_not_installed("mvtnorm")
	K <- 5L
	rho <- matrix(0.45, K, K)
	diag(rho) <- 1
	estimates <- c(2.5, 1.8, 0.9, 3.1, 2.0)
	ses <- rep(1, K)

	set.seed(1)
	adj <- .maxt_adjusted_p_nd(estimates, ses, rep(TRUE, K), rho)
	set.seed(1)
	crit <- mvtnorm::qmvnorm(
		1 - 0.05,
		corr = rho,
		tail = "both.tails",
		algorithm = mvtnorm::GenzBretz()
	)$quantile

	z <- abs(estimates / ses)
	# Exclude any effect within the GenzBretz noise floor of the boundary.
	near <- abs(z - crit) < 0.02
	expect_equal((z > crit)[!near], (adj < 0.05)[!near])
})

# ---------------------------------------------------------------------------
# 3. Degenerate effects scatter to NA; each non-degenerate slot is finite and
#    lies within its own effect's bounds (#551).
# ---------------------------------------------------------------------------
test_that(".maxt_adjusted_p_nd scatters NA into degenerate slots", {
	skip_if_not_installed("mvtnorm")
	K <- 5L
	rho_full <- matrix(0.4, K, K)
	diag(rho_full) <- 1
	estimates <- c(2.5, 1.8, 0.9, 3.1, 2.0)
	ses <- rep(1, K)
	nondeg <- c(TRUE, TRUE, FALSE, TRUE, FALSE)
	rho_nd <- rho_full[nondeg, nondeg]

	set.seed(1)
	adj <- .maxt_adjusted_p_nd(estimates, ses, nondeg, rho_nd)
	expect_true(all(is.na(adj[!nondeg])))
	expect_true(all(!is.na(adj[nondeg])))
	pw <- 2 * stats::pnorm(-abs(estimates / ses))
	expect_true(all(adj[nondeg] >= pw[nondeg]))
	expect_true(all(adj[nondeg] <= pmin(1, sum(nondeg) * pw[nondeg])))
})

# ---------------------------------------------------------------------------
# Shared fixture: a fit with multiple non-degenerate cohorts.
# ---------------------------------------------------------------------------
make_maxt_fit <- function(se_type = "default") {
	sc <- genCoefs(G = 4, T = 6, d = 2, density = 0.5, eff_size = 1.5, seed = 7)
	sim <- simulateData(sc, N = 120, sig_eps_sq = 2, sig_eps_c_sq = 1, seed = 7)
	fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		covs = sim$covs,
		response = sim$response,
		sig_eps_sq = sim$sig_eps_sq,
		sig_eps_c_sq = sim$sig_eps_c_sq,
		q = 0.5,
		se_type = se_type,
		verbose = FALSE
	)
}

# ---------------------------------------------------------------------------
# 4. Integration: catt_df p_value == simultaneousCIs()$adjusted_p_values, is
#    >= the pointwise p, and satisfies the band duality.
# ---------------------------------------------------------------------------
test_that("catt_df p_value is the adjusted dual under simultaneous (cohort)", {
	skip_if_not_installed("mvtnorm")
	fit <- make_maxt_fit()
	expect_identical(fit$ci_type, "simultaneous")
	cd <- fit$catt_df
	sci <- suppressMessages(simultaneousCIs(fit, family = "cohort"))

	# Same source: the fit-time finalizer copies the worker's adjusted p. The
	# values reach the far tail, so they are compared as a ratio (#551).
	fin <- is.finite(cd$p_value)
	expect_identical(fin, is.finite(sci$adjusted_p_values))
	expect_true(sum(fin) >= 2L) # fixture has >= 2 non-degenerate cohorts
	expect_equal(
		cd$p_value[fin] / sci$adjusted_p_values[fin],
		rep(1, sum(fin)),
		tolerance = 1e-6
	)

	# Between the pointwise p-value and its Bonferroni bound (#551).
	pw <- 2 * stats::pnorm(-abs(cd$estimate / cd$se))
	m <- sum(fin)
	expect_true(all(cd$p_value[fin] >= pw[fin] * (1 - 1e-6)))
	expect_true(all(cd$p_value[fin] <= pmin(1, m * pw[fin]) * (1 + 1e-6)))

	# Band duality: a cohort's band excludes 0 iff its adjusted p < alpha.
	crit <- sci$critical_value
	z <- abs(cd$estimate / cd$se)
	near <- !fin | abs(z - crit) < 0.02
	sig_band <- cd$ci_low > 0 | cd$ci_high < 0
	expect_equal(sig_band[!near], (cd$p_value < fit$alpha)[!near])
})

# ---------------------------------------------------------------------------
# 5. Event-study p_value follows ci_type + band duality.
# ---------------------------------------------------------------------------
test_that("eventStudy p_value is adjusted under simultaneous, Wald under pointwise", {
	skip_if_not_installed("mvtnorm")
	fit <- make_maxt_fit()
	es_s <- eventStudy(fit)
	es_p <- eventStudy(fit, ci_type = "pointwise")

	fin <- is.finite(es_s$p_value) & is.finite(es_p$p_value)
	expect_gte(sum(fin), 2L) # fixture has >= 2 non-degenerate event times
	# Pointwise is the Wald p-value; simultaneous is >= it.
	expect_equal(
		es_p$p_value[fin],
		(2 * stats::pnorm(-abs(es_p$estimate / es_p$se)))[fin]
	)
	expect_true(all(es_s$p_value[fin] >= es_p$p_value[fin] * (1 - 1e-6)))
})

# ---------------------------------------------------------------------------
# 6. Conservative se_type: the Bonferroni dual, min(1, K * pointwise_p) on the
#    band's Cauchy-Schwarz ses. NOTE: this pins the adjusted-p / band-se DUAL,
#    not the Cauchy-Schwarz SE formula itself (it recovers se_band from the
#    band's own CI, a round-trip, and uses family = "cohort" where Sigma_2 = 0).
#    The SE formula is pinned in test-conservative-band-se-225.R (#225).
# ---------------------------------------------------------------------------
test_that("conservative se_type yields the Bonferroni-dual adjusted p-value", {
	skip_if_not_installed("mvtnorm")
	fit_c <- make_maxt_fit("conservative")
	sci <- suppressMessages(simultaneousCIs(fit_c, family = "cohort"))
	K <- sci$K
	# Recover the band's per-effect se from its pointwise CI half-width.
	se_band <- (sci$ci$pointwise_ci_high - sci$ci$estimate) /
		sci$pointwise_critical_value
	pw_band <- 2 * stats::pnorm(-abs(sci$ci$estimate / se_band))
	expect_param <- pmin(1, K * pw_band)
	fin <- is.finite(sci$adjusted_p_values)
	expect_gte(sum(fin), 2L) # fixture has >= 2 non-degenerate cohorts
	expect_equal(
		sci$adjusted_p_values[fin],
		expect_param[fin],
		tolerance = 1e-9
	)
})

# ---------------------------------------------------------------------------
# 7. The overall ATT is K = 1, so its p-value stays the pointwise Wald value
#    (no multiplicity), regardless of ci_type.
# ---------------------------------------------------------------------------
test_that("overall-ATT p_value (K=1) stays the pointwise Wald value", {
	fit <- make_maxt_fit()
	expect_equal(
		fit$att_p_value,
		2 * stats::pnorm(-abs(fit$att_hat / fit$att_se))
	)
})

# ---------------------------------------------------------------------------
# 8. Determinism.
# ---------------------------------------------------------------------------
test_that("adjusted p-values are deterministic across calls", {
	skip_if_not_installed("mvtnorm")
	fit <- make_maxt_fit()
	es1 <- eventStudy(fit)
	es2 <- eventStudy(fit)
	expect_equal(es1$p_value, es2$p_value)
	sci1 <- suppressMessages(simultaneousCIs(fit, family = "cohort"))
	sci2 <- suppressMessages(simultaneousCIs(fit, family = "cohort"))
	expect_equal(sci1$adjusted_p_values, sci2$adjusted_p_values)
})

# ---------------------------------------------------------------------------
# 9. tidy() carries the (ci_type-dependent) p_value through unchanged.
# ---------------------------------------------------------------------------
test_that("tidy() p.value matches catt_df p_value (cohort rows)", {
	skip_if_not_installed("mvtnorm")
	skip_if_not_installed("broom")
	fit <- make_maxt_fit()
	td <- broom::tidy(fit)
	# .tidy_estimator_output builds p.value = c(att_pvalue, cohort_pvalues), so
	# the cohort p-values are the trailing nrow(catt_df) entries.
	cohort_p <- tail(td$p.value, nrow(fit$catt_df))
	expect_equal(cohort_p, fit$catt_df$p_value)
})

# ---------------------------------------------------------------------------
# 10. Across the tail, each non-degenerate adjusted p-value lies between its
#     pointwise p-value and its Bonferroni bound min(1, m * p), with
#     m = sum(nondeg) (#551). Both sides are the same scalar arithmetic on the
#     same inputs, so the bounds are asserted with no tolerance.
# ---------------------------------------------------------------------------
test_that(".maxt_adjusted_p_nd stays within [p, min(1, m p)] across the tail (#551)", {
	skip_if_not_installed("mvtnorm")
	equicorr <- function(r, k) {
		rho <- matrix(r, k, k)
		diag(rho) <- 1
		rho
	}
	others <- c(1.8, 0.9, 3.1, 2.0)
	ses <- rep(1, 5L)
	slices <- list(
		identity = list(nondeg = rep(TRUE, 5L), rho = equicorr(0, 5L)),
		rho_045 = list(nondeg = rep(TRUE, 5L), rho = equicorr(0.45, 5L)),
		rho_090 = list(nondeg = rep(TRUE, 5L), rho = equicorr(0.9, 5L)),
		# A degenerate slot between non-degenerate ones, so that m = 4 < K = 5.
		degenerate_slot = list(
			nondeg = c(TRUE, FALSE, TRUE, TRUE, TRUE),
			rho = equicorr(0.45, 4L)
		)
	)
	for (nm in names(slices)) {
		nondeg <- slices[[nm]]$nondeg
		rho <- slices[[nm]]$rho
		m <- sum(nondeg)
		below_floor <- character(0)
		above_cap <- character(0)
		for (z1 in seq(3, 12, by = 0.25)) {
			estimates <- c(z1, others)
			set.seed(1)
			adj <- .maxt_adjusted_p_nd(estimates, ses, nondeg, rho)
			pw <- 2 * stats::pnorm(-abs(estimates / ses))
			if (any(adj[nondeg] < pw[nondeg])) {
				below_floor <- c(below_floor, format(z1))
			}
			if (any(adj[nondeg] > pmin(1, m * pw[nondeg]))) {
				above_cap <- c(above_cap, format(z1))
			}
		}
		# The z values, if any, at which a bound is broken.
		expect_identical(below_floor, character(0), info = nm)
		expect_identical(above_cap, character(0), info = nm)
	}

	# #551's own case: block 1's fixture with its first effect at z = 9.8.
	K <- 5L
	set.seed(1)
	adj <- .maxt_adjusted_p_nd(
		c(9.8, others),
		ses,
		rep(TRUE, K),
		equicorr(0.45, K)
	)
	expect_identical(adj[1], 2 * stats::pnorm(-9.8))
})

# ---------------------------------------------------------------------------
# Shared fixture for blocks 11-13.
# ---------------------------------------------------------------------------
make_bounds_panel_551 <- function() {
	sc <- genCoefs(G = 3, T = 5, density = 0.5, eff_size = 2, d = 2, seed = 1)
	simulateData(sc, N = 120, sig_eps_sq = 0.5, sig_eps_c_sq = 0.5, seed = 1)
}

# Asserts that each finite entry of `p` lies within [pw, min(1, m * pw)] at a
# relative tolerance of 1e-6, where pw is the pointwise p-value built from
# `estimate` and `se`. Returns the finite mask.
expect_p_within_bounds_551 <- function(p, estimate, se, m, label) {
	pw <- 2 * stats::pnorm(-abs(estimate / se))
	fin <- is.finite(p)
	expect_true(all(p[fin] >= pw[fin] * (1 - 1e-6)), info = label)
	expect_true(
		all(p[fin] <= pmin(1, m * pw[fin]) * (1 + 1e-6)),
		info = label
	)
	fin
}

# ---------------------------------------------------------------------------
# 11. Integration (#551): every estimator's catt_df p-values, and the
#     eventStudy() p-values of the three estimators that have one, lie within
#     [p, min(1, m p)], with m the band's count of non-degenerate effects.
# ---------------------------------------------------------------------------
test_that("catt_df and eventStudy() p-values lie within [p, min(1, m p)] for every estimator (#551)", {
	skip_if_not_installed("mvtnorm")
	sim <- make_bounds_panel_551()
	fits <- list(
		fetwfe = fetwfeWithSimulatedData(sim),
		etwfe = etwfeWithSimulatedData(sim),
		betwfe = betwfeWithSimulatedData(sim),
		twfeCovs = twfeCovsWithSimulatedData(sim)
	)
	for (nm in names(fits)) {
		fit <- fits[[nm]]
		cd <- fit$catt_df
		sci <- suppressMessages(simultaneousCIs(fit, family = "cohort"))
		m <- sum(is.finite(sci$adjusted_p_values))
		expect_gte(m, 2L)
		fin <- expect_p_within_bounds_551(
			cd$p_value,
			cd$estimate,
			cd$se,
			m,
			paste(nm, "cohort")
		)
		# catt_df copies the band's p-values, so the counts match.
		expect_identical(sum(fin), m, info = nm)
		if (nm == "twfeCovs") {
			next
		}
		es <- eventStudy(fit)
		sci_es <- suppressMessages(simultaneousCIs(fit, family = "event_study"))
		m_es <- sum(is.finite(sci_es$adjusted_p_values))
		fin_es <- expect_p_within_bounds_551(
			es$p_value,
			es$estimate,
			es$se,
			m_es,
			paste(nm, "event study")
		)
		expect_gte(sum(fin_es), 2L)
	}

	# The etwfe fit: cohort 3 sits in the far tail, where the adjusted p-value
	# equals the pointwise one, and event time 2 is at most its Bonferroni bound.
	cd <- fits$etwfe$catt_df
	r3 <- as.character(cd$cohort) == "3"
	expect_identical(sum(r3), 1L)
	pw3 <- 2 * stats::pnorm(-abs(cd$estimate[r3] / cd$se[r3]))
	expect_equal(cd$p_value[r3] / pw3, 1, tolerance = 1e-6)

	es <- eventStudy(fits$etwfe)
	m_es <- sum(is.finite(
		suppressMessages(simultaneousCIs(
			fits$etwfe,
			family = "event_study"
		))$adjusted_p_values
	))
	r2 <- es$event_time == 2L
	expect_identical(sum(r2), 1L)
	pw2 <- 2 * stats::pnorm(-abs(es$estimate[r2] / es$se[r2]))
	expect_lte(es$p_value[r2], min(1, m_es * pw2) * (1 + 1e-6))
})

# ---------------------------------------------------------------------------
# 12. At a small alpha the analytic critical value is at most the Bonferroni
#     critical value over the non-degenerate effects, and strictly above the
#     pointwise one (#551).
# ---------------------------------------------------------------------------
test_that("the analytic critical value stays within its bounds at a small alpha (#551)", {
	skip_if_not_installed("mvtnorm")
	sim <- make_bounds_panel_551()
	for (nm in c("etwfe", "fetwfe")) {
		fit <- get(paste0(nm, "WithSimulatedData"))(sim)
		sci <- suppressMessages(simultaneousCIs(
			fit,
			family = "event_study",
			alpha = 0.001
		))
		m <- sum(is.finite(sci$adjusted_p_values))
		expect_gte(m, 2L)
		expect_lte(sci$critical_value, stats::qnorm(1 - 0.001 / (2 * m)))
		expect_lt(sci$pointwise_critical_value, sci$critical_value)
	}
})

# simultaneousCIs() with mvtnorm::qmvnorm()'s quantile overwritten by `q`. The
# real qmvnorm() still runs, so the RNG stream the adjusted p-values draw from
# is unchanged.
sci_with_forced_quantile_551 <- function(fit, family, q) {
	real_qmvnorm <- mvtnorm::qmvnorm
	testthat::with_mocked_bindings(
		suppressMessages(simultaneousCIs(fit, family = family)),
		qmvnorm = function(...) {
			out <- real_qmvnorm(...)
			out$quantile <- q
			out
		},
		.package = "mvtnorm"
	)
}

# ---------------------------------------------------------------------------
# 13. A quantile forced past either end comes out as that end, the pointwise
#     or the Bonferroni critical value over the m non-degenerate effects, and
#     the band and the adjusted p-values still agree on every non-degenerate
#     effect (#551). That agreement needs no such effect with its z in
#     (qnorm(1 - alpha/2), qnorm(1 - alpha/(2m))], which is asserted first.
# ---------------------------------------------------------------------------
test_that("a forced quantile is clamped to its bounds, keeping band duality (#551)", {
	skip_if_not_installed("mvtnorm")
	sim <- make_bounds_panel_551()
	fit <- etwfeWithSimulatedData(sim)
	alpha <- fit$alpha
	sci_hi <- sci_with_forced_quantile_551(fit, "event_study", 10)
	sci_lo <- sci_with_forced_quantile_551(fit, "event_study", 1.5)
	nd <- is.finite(sci_hi$adjusted_p_values)
	m <- sum(nd)
	lo <- stats::qnorm(1 - alpha / 2)
	hi <- stats::qnorm(1 - alpha / (2 * m))

	# Precondition, from eventStudy()'s own estimates and SEs.
	es <- eventStudy(fit)
	z <- abs(es$estimate / es$se)
	expect_gte(m, 2L)
	expect_identical(sum(is.finite(z)), m)
	expect_false(any(z > lo & z <= hi, na.rm = TRUE))

	expect_identical(sci_hi$critical_value, hi)
	expect_identical(sci_lo$critical_value, lo)
	for (sci in list(sci_hi, sci_lo)) {
		excludes_zero <- sci$ci$simultaneous_ci_low > 0 |
			sci$ci$simultaneous_ci_high < 0
		expect_identical(
			excludes_zero[nd],
			sci$adjusted_p_values[nd] < alpha,
			info = format(sci$critical_value)
		)
	}
})

# ---------------------------------------------------------------------------
# 14. The upper end of the critical value's clamp counts the non-degenerate
#     effects m, not the family size K (#551). On this panel the cohort family
#     has a degenerate effect.
# ---------------------------------------------------------------------------
test_that("the critical value's upper bound counts the non-degenerate effects (#551)", {
	skip_if_not_installed("mvtnorm")
	coefs <- genCoefs(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 2026
	)
	sim <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 2026
	)
	fit <- fetwfeWithSimulatedData(sim)
	sci <- sci_with_forced_quantile_551(fit, "cohort", 10)
	m <- sum(is.finite(sci$adjusted_p_values))
	expect_gte(m, 2L)
	expect_lt(m, sci$K)
	expect_identical(sci$critical_value, stats::qnorm(1 - fit$alpha / (2 * m)))
})
