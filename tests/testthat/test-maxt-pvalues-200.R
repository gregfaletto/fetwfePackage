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

	# Adjustment only widens, and is strictly visible for at least one effect.
	pw <- 2 * stats::pnorm(-abs(estimates / ses))
	expect_true(all(adj >= pw - 1e-9))
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
# 3. Degenerate effects scatter to NA; non-degenerate slots are finite.
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

	# Same source: the fit-time finalizer copies the worker's adjusted p.
	expect_equal(cd$p_value, sci$adjusted_p_values)

	fin <- is.finite(cd$p_value)
	expect_true(sum(fin) >= 2L) # fixture has >= 2 non-degenerate cohorts
	pw <- 2 * stats::pnorm(-abs(cd$estimate / cd$se))
	expect_true(all(cd$p_value[fin] >= pw[fin] - 1e-9))

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
	# Pointwise is the Wald p-value; simultaneous is >= it.
	expect_equal(
		es_p$p_value[fin],
		(2 * stats::pnorm(-abs(es_p$estimate / es_p$se)))[fin]
	)
	expect_true(all(es_s$p_value[fin] >= es_p$p_value[fin] - 1e-9))
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
# 8. Determinism: repeated calls give byte-identical adjusted p-values.
# ---------------------------------------------------------------------------
test_that("adjusted p-values are deterministic across calls", {
	skip_if_not_installed("mvtnorm")
	fit <- make_maxt_fit()
	es1 <- eventStudy(fit)
	es2 <- eventStudy(fit)
	expect_identical(es1$p_value, es2$p_value)
	sci1 <- suppressMessages(simultaneousCIs(fit, family = "cohort"))
	sci2 <- suppressMessages(simultaneousCIs(fit, family = "cohort"))
	expect_identical(sci1$adjusted_p_values, sci2$adjusted_p_values)
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
