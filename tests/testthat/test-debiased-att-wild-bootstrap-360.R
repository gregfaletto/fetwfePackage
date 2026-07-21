library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #360 / #363: the debiased overall ATT via a studentized score /
# influence-function WILD CLUSTER BOOTSTRAP (method = "bootstrap"),
# re-signing the per-unit influence summands (V1 regression `unit_scores` + V2
# cohort-weight, built from a_att_G = (catt - att_hat)/S) WITHOUT a refit per
# replicate. Point estimate and reported `se` are unchanged; only the critical
# value changes, and it is FLOORED at qnorm(1 - alpha/2) (#363) so the interval is
# never narrower than the analytic Wald interval -- this is NOT a few-clusters
# remedy (an unrestricted no-refit score bootstrap narrows below the Gaussian
# under heterogeneous cluster influence). Default method = "analytic" is
# byte-identical to the previous behavior.
# ------------------------------------------------------------------------------

# Single-sample fixtures (the wild bootstrap requires single-sample fits;
# fetwfeWithSimulatedData() uses indep_counts, which is unsupported -- see below).
.wb_fixedp <- local({
	cf <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 1)
	sim <- simulateData(
		cf,
		N = 60,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		q = 0.5,
		verbose = FALSE
	)
})

.wb_highdim <- local({
	cf <- genCoefs(G = 3, T = 5, d = 22, density = 0.5, eff_size = 2, seed = 1)
	sim <- simulateData(
		cf,
		N = 30,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		q = 0.5,
		gls = FALSE, # Omega-free high-dim path (REML can't run at p >= NT)
		verbose = FALSE
	)
})

test_that("analytic default is unchanged and adds no bootstrap fields (#360)", {
	for (fit in list(.wb_fixedp, .wb_highdim)) {
		a0 <- debiasedATT(fit)
		a1 <- debiasedATT(fit, method = "analytic")
		# method = "analytic" is the default and is identical.
		expect_identical(a0, a1)
		# No bootstrap fields on the analytic object.
		expect_false(any(
			c("method", "crit_value", "B", "multiplier") %in% names(a0)
		))
	}
	# fixed-p analytic keeps EXACTLY the documented names (the #326 contract).
	expect_identical(
		names(debiasedATT(.wb_fixedp)),
		c("att", "se", "ci_low", "ci_high", "var_reg", "var_weight")
	)
})

test_that("wild bootstrap leaves att and se unchanged, only the CI changes (#360)", {
	for (fit in list(.wb_fixedp, .wb_highdim)) {
		a <- debiasedATT(fit)
		w <- debiasedATT(fit, method = "bootstrap", seed = 1)
		expect_equal(w$att, a$att)
		expect_equal(w$se, a$se) # analytic sandwich SE is reused verbatim
		expect_identical(w$method, "bootstrap")
		expect_true(all(
			c("crit_value", "B", "multiplier", "alpha") %in% names(w)
		))
		# CI is centered on att with half-width crit * se.
		expect_equal(w$ci_high - w$att, w$crit_value * w$se)
		expect_equal(w$att - w$ci_low, w$crit_value * w$se)
		expect_gt(w$crit_value, 0)
	}
})

test_that("the V2 cohort-weight influence function reproduces var_weight (the anchor) (#360)", {
	# The correctness gate: the per-unit V2 summands psi2 built from the
	# closed-form gradient a_att_G = (catt - att_hat)/S must satisfy
	# sum(psi2^2)/n^2 == var_weight (== .plugin_v2). Mutating the gradient breaks
	# this, so the test is sensitive to a wrong V2 construction.
	for (fit in list(.wb_fixedp, .wb_highdim)) {
		db <- debiasedATT(fit)
		G <- fit$G
		Tt <- fit$T
		N <- fit$N
		n <- N * Tt
		S <- sum(fit$cohort_probs_overall[seq_len(G)])
		a_att_G <- (fit$catt_df$estimate - fit$att_hat) / S
		psi2 <- as.numeric(fetwfe:::.build_propensity_if(
			cohort_probs_overall = fit$cohort_probs_overall,
			G = G,
			N = N,
			T = Tt,
			A = matrix(a_att_G, nrow = G, ncol = 1L)
		))
		expect_equal(sum(psi2^2) / n^2, db$var_weight)
	}
})

test_that("wild bootstrap is reproducible with a seed and ambient without (#360)", {
	fit <- .wb_fixedp
	w1 <- debiasedATT(fit, method = "bootstrap", seed = 42)
	w2 <- debiasedATT(fit, method = "bootstrap", seed = 42)
	expect_identical(w1$crit_value, w2$crit_value)
	expect_identical(w1$ci_low, w2$ci_low)
	# A different seed (almost surely) gives a different critical value.
	w3 <- debiasedATT(fit, method = "bootstrap", seed = 7)
	expect_false(isTRUE(all.equal(w1$crit_value, w3$crit_value)))
	# seed = NULL draws from the ambient stream and still returns a finite crit.
	set.seed(1)
	w_null <- debiasedATT(fit, method = "bootstrap")
	expect_true(is.finite(w_null$crit_value) && w_null$crit_value > 0)
	# Seeded draw does not disturb the ambient RNG stream (.with_preserved_rng).
	set.seed(100)
	u_before <- runif(1)
	set.seed(100)
	invisible(debiasedATT(fit, method = "bootstrap", seed = 42))
	u_after <- runif(1)
	expect_identical(u_before, u_after)
})

test_that("all three multiplier types produce finite critical values (#360)", {
	fit <- .wb_fixedp
	for (m in c("rademacher", "mammen", "webb")) {
		w <- debiasedATT(
			fit,
			method = "bootstrap",
			multiplier = m,
			seed = 1
		)
		expect_identical(w$multiplier, m)
		expect_true(is.finite(w$crit_value) && w$crit_value > 0)
	}
})

test_that(".draw_multipliers webb has the right support and moments (#360)", {
	set.seed(1)
	w <- fetwfe:::.draw_multipliers(1L, 2e5L, "webb")
	expect_setequal(
		round(sort(unique(as.numeric(w))), 6),
		round(c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5)), 6)
	)
	expect_equal(mean(w), 0, tolerance = 0.02) # mean 0
	expect_equal(mean(w^2), 1, tolerance = 0.02) # unit variance
})

test_that("wild bootstrap now supports indep_counts (two-sample) fits (#402)", {
	# #360 originally REFUSED the bootstrap on two-sample fits; #402 lifted that
	# gate after verifying the per-unit two-sample propensity IF reproduces the
	# two-sample var_weight (the runtime V2 anchor enforces this to 1e-6). The
	# comprehensive proof (IF reconstruction vs att_var_2) lives in
	# test-align-bootstrap-accessors-402.R; here we pin that the call now succeeds.
	cf <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 3)
	sim <- simulateData(
		cf,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 3
	)
	fit_ic <- fetwfeWithSimulatedData(sim, q = 0.5, verbose = FALSE)
	expect_true(isTRUE(fit_ic$indep_counts_used))
	# Analytic still works on indep_counts fits.
	expect_s3_class(debiasedATT(fit_ic), "debiased_att")
	# The wild bootstrap now runs (previously errored with "indep_counts").
	res <- debiasedATT(fit_ic, method = "bootstrap", seed = 1)
	expect_s3_class(res, "debiased_att")
	expect_identical(res$method, "bootstrap")
})

test_that("wild bootstrap validates B and seed (#360)", {
	fit <- .wb_fixedp
	expect_error(
		debiasedATT(fit, method = "bootstrap", B = -1),
		"`B`"
	)
	expect_error(
		debiasedATT(fit, method = "bootstrap", B = 2.5),
		"`B`"
	)
	expect_error(
		debiasedATT(fit, method = "bootstrap", seed = c(1, 2)),
		"`seed`"
	)
	# A non-integer seed is rejected (it would be silently truncated by set.seed).
	expect_error(
		debiasedATT(fit, method = "bootstrap", seed = 1.5),
		"`seed`"
	)
	expect_error(
		debiasedATT(fit, method = "bootstrap", multiplier = "gaussian")
	)
})

test_that("the default wild-bootstrap multiplier is webb (#360, review 1)", {
	# Webb is kept as the default (a studentized statistic, vs rademacher's constant
	# denominator); given the #363 floor the multiplier choice only matters in the
	# near-homogeneous case where the bootstrap widens above the floor.
	w <- debiasedATT(.wb_fixedp, method = "bootstrap", seed = 1)
	expect_identical(w$multiplier, "webb")
})

test_that("wild bootstrap handles a single-cohort fit, var_weight = 0 (#360, review 6)", {
	# One treated cohort => no cohort-weight variance => psi2 == 0 => the bootstrap
	# degrades cleanly to a V1-only wild bootstrap (anchor 0 == 0 holds, no NaN).
	cf <- genCoefs(G = 1, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 4)
	sim <- simulateData(
		cf,
		N = 60,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 4
	)
	fit <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		q = 0.5,
		verbose = FALSE
	)
	a <- debiasedATT(fit)
	expect_equal(a$var_weight, 0)
	w <- debiasedATT(fit, method = "bootstrap", seed = 1)
	expect_equal(w$se, a$se)
	expect_true(is.finite(w$crit_value) && w$crit_value > 0)
	expect_equal(w$ci_high - w$att, w$crit_value * w$se)
})

test_that("webb is reachable from simultaneousCIs() too (#360, review 5)", {
	sc <- suppressWarnings(simultaneousCIs(
		.wb_fixedp,
		family = "cohort",
		method = "bootstrap",
		B = 200,
		seed = 1,
		multiplier = "webb"
	))
	expect_s3_class(sc, "simultaneous_cis")
	expect_true(is.finite(sc$critical_value) && sc$critical_value > 0)
})

test_that("print.debiased_att surfaces the bootstrap method and level (#360)", {
	fit <- .wb_fixedp
	w <- debiasedATT(
		fit,
		method = "bootstrap",
		multiplier = "webb",
		seed = 1
	)
	out <- capture.output(print(w))
	expect_true(any(grepl("wild cluster bootstrap", out)))
	expect_true(any(grepl("webb weights", out)))
	# Level read from stored alpha (0.05 -> 95%), NOT inferred from crit * se.
	expect_true(any(grepl("95% CI", out)))
	# The analytic print does NOT show the bootstrap line.
	out_a <- capture.output(print(debiasedATT(fit)))
	expect_false(any(grepl("wild cluster bootstrap", out_a)))
})

test_that("the wild-bootstrap crit is floored at the Gaussian quantile (#363)", {
	z <- stats::qnorm(0.975)
	# The (floored) crit is never below z on any fit.
	for (fit in list(.wb_fixedp, .wb_highdim)) {
		w <- debiasedATT(fit, method = "bootstrap", seed = 1)
		expect_gte(w$crit_value, z)
	}
	# The floor BITES: on a heterogeneous few-cluster fit whose pre-floor crit is
	# ~1.44 (< z -- the anti-conservative narrowing this fix removes, #363), the
	# returned crit equals z exactly, so the interval is the analytic Wald interval.
	# Remove the `max(crit, z)` floor and this returns ~1.44 != z.
	cf <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 1)
	sim <- simulateData(
		cf,
		N = 18,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	fit_het <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		q = 0.5,
		gls = FALSE,
		verbose = FALSE
	)
	w_het <- debiasedATT(
		fit_het,
		method = "bootstrap",
		multiplier = "webb",
		seed = 1,
		B = 2000
	)
	expect_equal(w_het$crit_value, z)
	a_het <- debiasedATT(fit_het)
	expect_equal(w_het$ci_low, a_het$att - z * a_het$se)
	expect_equal(w_het$ci_high, a_het$att + z * a_het$se)
	# The floor is alpha-dependent (qnorm(1 - alpha/2), not a hardcoded 1.96): at
	# alpha = 0.1 the same heterogeneous fit floors to qnorm(0.95). A mutant that
	# floored at a literal qnorm(0.975) would pass at the default alpha but fail
	# here.
	w_het_90 <- debiasedATT(
		fit_het,
		method = "bootstrap",
		multiplier = "webb",
		seed = 1,
		B = 2000,
		alpha = 0.1
	)
	expect_equal(w_het_90$crit_value, stats::qnorm(0.95))
})
