# Coverage backfill for the high-dimensional desparsified inference cluster
# (#328, from the 2026-06-19 review Lens 4). Test hardening only -- no production
# change. Item 1 (scattered x lambda_c="cv") is already covered by
# test-aatt-scattered-cv-323.R; gap "T=2 x G>=2" is infeasible (G <= T-1) and
# "no-never-treated x high-dim" is not a distinct path (allow_no_never_treated
# auto-truncates to low-dim), so neither is tested here.

# Cheap full-p>=NT fixture (mirrors test-simultaneous-bootstrap-highdim-142.R):
# G=3,T=5,d=22,N=40 -> p=390 >= NT=200. sig_eps_sq / sig_eps_c_sq must be supplied
# (REML cannot estimate them when p >= N(T-1)).
.hd_fit_328 <- function(q) {
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
	fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = q,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
}

# ------------------------------------------------------------------------------
# Gap 2: q in {1, 2} through the high-dim debiased / bootstrap paths. The high-dim
# regime accepts ANY q (the gate is regime-aware, R/debiased_att.R ~338) because
# it re-fits its OWN q=1 fused-lasso nuisance from the FULL design (#303), so the
# debiased POINT estimate and the regression variance channel are INVARIANT to the
# input fit's q. (The SE's cohort-weight channel `var_weight` is q-dependent -- the
# fit's stored `att_var_2` when present (q < 1) else the `.plugin_v2()` recompute
# (q >= 1, where `att_var_2` is NA); R/debiased_att.R:386-395 -- so the SE itself
# varies with q, pinned below.)
# ------------------------------------------------------------------------------
test_that("high-dim debiasedATT() accepts q in {0.5,1,2} and the point estimate is q-invariant (#328)", {
	skip_on_cran()
	f05 <- .hd_fit_328(0.5)
	expect_gte(ncol(f05$internal$X_final), nrow(f05$internal$X_final)) # p >= NT
	f1 <- .hd_fit_328(1)
	f2 <- .hd_fit_328(2)

	# the bridge fits genuinely differ across q (so invariance below is non-trivial)
	expect_false(isTRUE(all.equal(f05$att_hat, f1$att_hat)))
	expect_false(isTRUE(all.equal(f05$att_hat, f2$att_hat)))

	d05 <- debiasedATT(f05)
	d1 <- debiasedATT(f1)
	d2 <- debiasedATT(f2)

	# every q is accepted in the high-dim regime, with a finite positive SE.
	for (d in list(d05, d1, d2)) {
		expect_s3_class(d, "debiased_att")
		expect_true(is.finite(d$att) && is.finite(d$se) && d$se > 0)
		expect_false(is.null(d$lambda_node)) # high-dim diagnostics present
	}

	# LOAD-BEARING: the debiased point estimate AND the regression variance channel
	# are EXACTLY invariant to the input q (they depend only on the full design and
	# the internal q=1 nuisance, not the bridge `theta_hat`). Mutating the high-dim
	# nuisance to reuse the bridge `theta_hat` breaks this.
	expect_equal(d1$att, d05$att, tolerance = 1e-9)
	expect_equal(d2$att, d05$att, tolerance = 1e-9)
	expect_equal(d1$var_reg, d05$var_reg, tolerance = 1e-9)
	expect_equal(d2$var_reg, d05$var_reg, tolerance = 1e-9)
	# the SE's cohort-weight channel `var_weight` is q-dependent (att_var_2 for
	# q < 1, the .plugin_v2() recompute for q >= 1), so the SE legitimately varies
	# with q -- confirm it actually moves (not q-invariant).
	expect_false(isTRUE(all.equal(d05$se, d1$se)))
})

test_that("high-dim simultaneousCIs(bootstrap) accepts q in {0.5,1,2}; band center is q-invariant (#328)", {
	skip_on_cran()
	mkband <- function(q) {
		suppressWarnings(simultaneousCIs(
			.hd_fit_328(q),
			family = "cohort",
			method = "bootstrap",
			B = 200,
			seed = 1
		))
	}
	b05 <- mkband(0.5)
	b1 <- mkband(1)
	b2 <- mkband(2)
	for (b in list(b05, b1, b2)) {
		expect_identical(b$regime, "high-dimensional")
		expect_true(all(is.finite(b$ci$simultaneous_ci_low)))
		expect_true(all(is.finite(b$ci$simultaneous_ci_high)))
	}
	# the desparsified band center is the debiased estimate -> q-invariant.
	expect_equal(b1$ci$estimate, b05$ci$estimate, tolerance = 1e-9)
	expect_equal(b2$ci$estimate, b05$ci$estimate, tolerance = 1e-9)
})

# ------------------------------------------------------------------------------
# Gap 3: minimal panels (d=0 no covariates; G=1 single cohort) through the
# debiased / band paths. Both are low-dimensional (p < NT) but were exercised by
# no debiasedATT()/simultaneousCIs() test before.
# ------------------------------------------------------------------------------
test_that("debiasedATT() and simultaneousCIs() run on a d=0 (no-covariate) panel (#328)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 0,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	dat <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	fit <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE
	)
	expect_lt(ncol(fit$internal$X_final), nrow(fit$internal$X_final)) # fixed-p
	db <- debiasedATT(fit)
	expect_true(is.finite(db$att) && is.finite(db$se) && db$se > 0)
	expect_null(db$lambda_node) # fixed-p: no high-dim diagnostics

	sc <- simultaneousCIs(
		fit,
		family = "cohort",
		method = "bootstrap",
		seed = 1
	)
	expect_identical(sc$regime, "fixed-p")
	expect_identical(nrow(sc$ci), 3L) # G = 3 cohorts
	fin <- is.finite(sc$ci$simultaneous_ci_low)
	expect_true(all(
		sc$ci$simultaneous_ci_low[fin] <= sc$ci$estimate[fin] &
			sc$ci$estimate[fin] <= sc$ci$simultaneous_ci_high[fin]
	))
})

test_that("debiasedATT() and simultaneousCIs() run on a G=1 (single-cohort) panel (#328)", {
	coefs <- genCoefs(
		G = 1,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	dat <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	fit <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE
	)
	db <- debiasedATT(fit)
	expect_true(is.finite(db$att) && is.finite(db$se) && db$se > 0)

	sc <- simultaneousCIs(fit, family = "cohort")
	expect_identical(nrow(sc$ci), 1L) # single cohort -> K = 1
	# K = 1 bypasses the correlation machinery: the simultaneous critical value
	# collapses to the pointwise qnorm(1 - alpha/2).
	expect_equal(sc$critical_value, sc$pointwise_critical_value)
	expect_equal(sc$critical_value, stats::qnorm(1 - 0.05 / 2))
})

# ------------------------------------------------------------------------------
# Gap 5: schema pins for tidy.simultaneous_cis and debiasedATTWithSimulatedData.
# ------------------------------------------------------------------------------
test_that("tidy.simultaneous_cis returns the documented column schema (#328)", {
	coefs <- genCoefs(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 3
	)
	dat <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 3
	)
	fit <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE
	)
	expected_cols <- c(
		"term",
		"estimate",
		"conf.low",
		"conf.high",
		"pointwise.conf.low",
		"pointwise.conf.high"
	)
	for (fam in c("cohort", "event_study")) {
		for (meth in c("analytic", "bootstrap")) {
			sci <- suppressWarnings(simultaneousCIs(
				fit,
				family = fam,
				method = meth,
				seed = 1
			))
			td <- tidy(sci)
			info <- paste(fam, meth)
			expect_s3_class(td, "data.frame")
			expect_identical(names(td), expected_cols, info = info)
			expect_identical(nrow(td), nrow(sci$ci), info = info)
			expect_type(td$term, "character")
			expect_true(
				all(vapply(td[-1], is.numeric, NA)),
				info = info
			)
			# the tidy columns mirror the source object's ci frame.
			expect_equal(td$estimate, sci$ci$estimate, info = info)
			expect_equal(td$conf.low, sci$ci$simultaneous_ci_low, info = info)
		}
	}
})

test_that("debiasedATTWithSimulatedData() returns the debiased_att schema, fixed-p and high-dim (#328)", {
	skip_on_cran()
	fields <- c("att", "se", "ci_low", "ci_high", "var_reg", "var_weight")
	check_schema <- function(dws, label) {
		expect_s3_class(dws, "debiased_att")
		expect_true(all(fields %in% names(dws)), info = label)
		expect_true(
			is.finite(dws$att) && is.finite(dws$se) && dws$se > 0,
			info = label
		)
		expect_true(
			dws$ci_low <= dws$att && dws$att <= dws$ci_high,
			info = label
		)
	}

	# fixed-p: also byte-identical to the explicit two-step pipeline.
	coefs_lp <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 2
	)
	dat_lp <- simulateData(
		coefs_lp,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 2
	)
	dws_lp <- debiasedATTWithSimulatedData(dat_lp, q = 0.5)
	check_schema(dws_lp, "fixed-p")
	expect_equal(
		unclass(dws_lp),
		unclass(debiasedATT(fetwfeWithSimulatedData(dat_lp, q = 0.5))),
		tolerance = 1e-9
	)

	# high-dim: the wrapper passes the simulator's true sigmas through, so the
	# p >= NT fit is constructible; the result carries the high-dim diagnostics.
	coefs_hd <- genCoefs(
		G = 3,
		T = 5,
		d = 22,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	dat_hd <- simulateData(
		coefs_hd,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	dws_hd <- debiasedATTWithSimulatedData(dat_hd, q = 0.5)
	check_schema(dws_hd, "high-dim")
	expect_false(is.null(dws_hd$lambda_node)) # high-dim diagnostics present
})
