# Tests for the multiplier-bootstrap method of simultaneousCIs() (#142, Phase 1).
#
# The bootstrap perturbs the per-unit influence-function matrix F (N x K) instead
# of computing the analytic mvtnorm::qmvnorm() critical value. Phase 1 covers the
# regression-channel families (cohort / all_post_treatment / custom), where the
# cohort-probability variance term is zero. The correctness anchor: F is built so
# that crossprod(F)/(NT)^2 equals the cluster-robust Sigma, hence on a
# se_type = "cluster" fit the bootstrap per-effect SEs equal the analytic
# per-effect SEs to machine precision; and the bootstrap critical value matches
# the analytic qmvnorm one up to Monte-Carlo error.

.boot_fit <- function(
	se_type = "cluster",
	G = 3,
	T = 5,
	d = 2,
	N = 150,
	seed = 1
) {
	coefs <- genCoefs(
		G = G,
		T = T,
		d = d,
		density = 0.6,
		eff_size = 1.5,
		seed = seed
	)
	dat <- simulateData(
		coefs,
		N = N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = seed
	)
	dat$indep_counts <- NA
	fetwfeWithSimulatedData(dat, q = 0.5, se_type = se_type)
}

.se_from <- function(sc) {
	# per-effect SE backed out of the pointwise band
	(sc$ci$pointwise_ci_high - sc$ci$estimate) / sc$pointwise_critical_value
}

test_that("bootstrap per-effect SEs match the analytic cluster SEs to machine precision (the anchor)", {
	fit <- .boot_fit(se_type = "cluster")
	an <- simultaneousCIs(fit, family = "cohort", method = "analytic")
	bo <- simultaneousCIs(
		fit,
		family = "cohort",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	# crossprod(F)/n^2 == cluster Sigma / cadjust, and the bootstrap SE restores
	# cadjust, so the per-effect SEs coincide with the analytic cluster SEs. The
	# SEs do NOT depend on the random multipliers, so this is exact.
	expect_equal(.se_from(bo), .se_from(an), tolerance = 1e-9)
})

test_that("bootstrap critical value matches the analytic qmvnorm one up to Monte-Carlo error", {
	fit <- .boot_fit(se_type = "cluster")
	for (fam in c("cohort", "all_post_treatment")) {
		an <- simultaneousCIs(fit, family = fam, method = "analytic")
		bo <- simultaneousCIs(
			fit,
			family = fam,
			method = "bootstrap",
			B = 5000,
			seed = 1
		)
		expect_lt(abs(bo$critical_value - an$critical_value), 0.12)
	}
})

test_that("bootstrap critical value sits in [pointwise, Bonferroni]", {
	fit <- .boot_fit(se_type = "cluster")
	bo <- simultaneousCIs(
		fit,
		family = "all_post_treatment",
		method = "bootstrap",
		B = 4000,
		seed = 7
	)
	expect_gte(bo$critical_value, bo$pointwise_critical_value - 1e-8)
	expect_lte(bo$critical_value, bo$bonferroni_critical_value + 1e-8)
	# strictly tighter than Bonferroni when K > 1 and effects are correlated
	expect_lt(bo$critical_value, bo$bonferroni_critical_value)
})

test_that("a single non-degenerate effect uses the exact pointwise critical value", {
	fit <- .boot_fit(se_type = "cluster")
	ctr <- matrix(c(1, rep(0, length(fit$treat_inds) - 1L)), nrow = 1)
	bo <- simultaneousCIs(
		fit,
		family = "custom",
		contrasts = ctr,
		method = "bootstrap",
		B = 1000,
		seed = 1
	)
	expect_equal(bo$K, 1L)
	expect_equal(bo$critical_value, stats::qnorm(1 - 0.05 / 2))
})

test_that("bootstrap is deterministic given a seed, and varies without one", {
	fit <- .boot_fit(se_type = "cluster")
	a <- simultaneousCIs(
		fit,
		family = "all_post_treatment",
		method = "bootstrap",
		B = 1000,
		seed = 42
	)
	b <- simultaneousCIs(
		fit,
		family = "all_post_treatment",
		method = "bootstrap",
		B = 1000,
		seed = 42
	)
	expect_identical(a$critical_value, b$critical_value)
	c2 <- simultaneousCIs(
		fit,
		family = "all_post_treatment",
		method = "bootstrap",
		B = 1000,
		seed = 99
	)
	expect_false(isTRUE(all.equal(a$critical_value, c2$critical_value)))
})

test_that("seed is restored: a seeded bootstrap does not perturb the ambient RNG stream", {
	fit <- .boot_fit(se_type = "cluster")
	set.seed(123)
	before <- stats::runif(1)
	set.seed(123)
	invisible(simultaneousCIs(
		fit,
		family = "cohort",
		method = "bootstrap",
		B = 500,
		seed = 7
	))
	after <- stats::runif(1)
	expect_identical(before, after)
})

test_that("mammen multipliers also produce an in-range critical value", {
	fit <- .boot_fit(se_type = "cluster")
	bo <- simultaneousCIs(
		fit,
		family = "all_post_treatment",
		method = "bootstrap",
		B = 4000,
		seed = 1,
		multiplier = "mammen"
	)
	expect_gte(bo$critical_value, bo$pointwise_critical_value - 1e-8)
	expect_lte(bo$critical_value, bo$bonferroni_critical_value + 1e-8)
})

test_that("method = 'analytic' is the default and byte-identical to omitting it", {
	fit <- .boot_fit(se_type = "cluster")
	expect_equal(
		simultaneousCIs(fit, family = "cohort"),
		simultaneousCIs(fit, family = "cohort", method = "analytic")
	)
	# the bootstrap return adds method/B/seed/multiplier and is otherwise the
	# same S3 shape
	bo <- simultaneousCIs(
		fit,
		family = "cohort",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(bo, "simultaneous_cis")
	expect_true(all(c("method", "B", "seed", "multiplier") %in% names(bo)))
	expect_identical(bo$method, "bootstrap")
})

test_that("event_study is supported by the bootstrap method (#142 Phase 3)", {
	# Phase 3 added the per-unit propensity influence-function channel, so
	# event_study + method = "bootstrap" now returns a valid band (it used to
	# stop()). See test-simultaneous-bootstrap-eventstudy-142.R for the anchor /
	# SE-match / coverage tests; this just confirms it no longer errors and is
	# consistent with the analytic path.
	fit <- .boot_fit(se_type = "cluster")
	bo <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(bo, "simultaneous_cis")
	expect_identical(bo$family, "event_study")
	expect_identical(bo$regime, "fixed-p")
	expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
	# and still works analytically
	expect_s3_class(
		simultaneousCIs(fit, family = "event_study", method = "analytic"),
		"simultaneous_cis"
	)
})

test_that("bootstrap validates B and seed", {
	fit <- .boot_fit(se_type = "cluster")
	expect_error(
		simultaneousCIs(fit, family = "cohort", method = "bootstrap", B = 0),
		"B"
	)
	expect_error(
		simultaneousCIs(
			fit,
			family = "cohort",
			method = "bootstrap",
			seed = c(1, 2)
		),
		"seed"
	)
})

test_that("adjusted p-values are the dual of the simultaneous band (both sides)", {
	# Weak-/sparse-effect DGP so the all_post_treatment family contains BOTH
	# rejected and non-rejected cells, exercising both directions of the iff
	# (a strong-effect DGP rejects every cell and never bites the non-reject side).
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 0.4,
		seed = 2
	)
	dat <- simulateData(
		coefs,
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 2
	)
	dat$indep_counts <- NA
	fit <- fetwfeWithSimulatedData(dat, q = 0.5, se_type = "cluster")
	bo <- simultaneousCIs(
		fit,
		family = "all_post_treatment",
		method = "bootstrap",
		B = 4000,
		seed = 3
	)
	excludes_zero <- bo$ci$simultaneous_ci_low > 0 |
		bo$ci$simultaneous_ci_high < 0
	nd <- !is.na(bo$adjusted_p_values)
	# outside band iff adjusted p < alpha -- over all non-degenerate effects, so
	# both the reject and non-reject sides bite whenever the family has each.
	expect_equal(excludes_zero[nd], (bo$adjusted_p_values < bo$alpha)[nd])
})

test_that("the bootstrap method works for the non-FETWFE estimators (etwfe, betwfe)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.6,
		eff_size = 1.5,
		seed = 1
	)
	dat <- simulateData(
		coefs,
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	dat$indep_counts <- NA
	# etwfe: full design, no selection (sel_feat_inds = NA path);
	# betwfe: beta-space selected support.
	for (fit in list(
		etwfeWithSimulatedData(dat),
		betwfeWithSimulatedData(dat)
	)) {
		bo <- simultaneousCIs(
			fit,
			family = "cohort",
			method = "bootstrap",
			B = 1000,
			seed = 1
		)
		expect_s3_class(bo, "simultaneous_cis")
		expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
		expect_gte(bo$critical_value, bo$pointwise_critical_value - 1e-8)
		expect_lte(bo$critical_value, bo$bonferroni_critical_value + 1e-8)
	}
})

test_that(".simultaneous_bootstrap_crit handles degenerate / single-effect F without warnings", {
	# All-degenerate F (zero-variance columns): no draw, no max()-over-empty
	# warning, crit falls back to the exact pointwise quantile.
	F0 <- matrix(0, 60, 3)
	expect_no_warning(
		r0 <- .simultaneous_bootstrap_crit(
			F0,
			n = 240,
			alpha = 0.05,
			B = 100,
			multiplier = "rademacher",
			seed = 1
		)
	)
	expect_equal(r0$crit, stats::qnorm(1 - 0.05 / 2))
	expect_null(r0$boot_max)
	# Exactly one non-degenerate effect: same fallback, no warning.
	set.seed(1)
	F1 <- cbind(stats::rnorm(60), 0, 0)
	expect_no_warning(
		r1 <- .simultaneous_bootstrap_crit(F1, n = 240, alpha = 0.05, B = 100)
	)
	expect_equal(r1$crit, stats::qnorm(1 - 0.05 / 2))
})

test_that("bootstrap simultaneous bands attain near-nominal family-wise coverage", {
	skip_on_cran()
	# A small but honest family-wise coverage check for the cohort family: in what
	# fraction of simulated panels does the simultaneous band cover ALL true cohort
	# ATTs jointly? Nominal is 0.95; the lower bound is generous to absorb the
	# modest nsim and finite N (this is a sanity floor, not a precise estimate).
	G <- 2L
	T <- 4L
	nsim <- 40L
	coefs <- genCoefs(
		G = G,
		T = T,
		d = 1,
		density = 0.8,
		eff_size = 1.5,
		seed = 1
	)
	num_treats <- getNumTreats(G = G, T = T)
	true_catt <- getActualCohortTes(
		G = G,
		first_inds = getFirstInds(G = G, T = T),
		treat_inds = getTreatInds(
			G = G,
			T = T,
			d = 1L,
			num_treats = num_treats
		),
		coefs = coefs$beta,
		num_treats = num_treats
	)
	covered <- logical(nsim)
	for (s in seq_len(nsim)) {
		dat <- simulateData(
			coefs,
			N = 100,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = 100L + s
		)
		dat$indep_counts <- NA
		fit <- fetwfeWithSimulatedData(dat, q = 0.5, se_type = "cluster")
		bo <- simultaneousCIs(
			fit,
			family = "cohort",
			method = "bootstrap",
			B = 400,
			seed = 1
		)
		covered[s] <- all(
			true_catt >= bo$ci$simultaneous_ci_low &
				true_catt <= bo$ci$simultaneous_ci_high
		)
	}
	expect_gt(mean(covered), 0.80)
})
