# ==============================================================================
# Single-treated-cohort (G = 1) support (#112).
#
# Re-relaxes the historical `G >= 2` floor (re-tightened in #109) so that
# fetwfe / etwfe / betwfe / twfeCovs and genCoefs / simulateData accept panels
# with a SINGLE treated cohort (one cohort that adopts treatment, plus
# never-treated units).
#
# NOTE (#251): #112 originally also required at least two treatment effects
# (a `num_treats >= 2` guard in prep_for_etwfe_core()), rejecting num_treats = 1
# panels (T = 2 and late-adopting single cohorts). #251 removed that guard so
# those panels now FIT. The num_treats = 1 tests below (section 4) were
# converted from rejection assertions to fits accordingly; the T == 2 warning
# behavior is covered in detail in test-2x2-support-251.R.
#
# Coverage:
#   1. All four estimators recover a known ATT at G = 1 (numerical, not just
#      finite) via a simulateData(G = 1) fixture; twfeCovs is biased by design
#      (asserted finite + G == 1 only).
#   2. The three fusion structures (cohort / event_study / custom fusion_matrix)
#      all fit at G = 1, and cohort == event_study EXACTLY (no cross-cohort
#      dimension to fuse, so both reduce to within-cohort dynamic-effect fusion).
#   3. A hand-built SCATTERED-offset (offset 3, num_treats = 4) single-cohort
#      panel fits -- genCoefs / simulateData force offset 2, so the scattered
#      and late-adopting fixtures are built by hand.
#   4. num_treats = 1 now FITS (#251 removed the #112 guard): T = 2 (fetwfe
#      warns "no fusion") AND a late-adopting single cohort (offset = T) both
#      produce a valid fit. (Warning asymmetries are asserted in
#      test-2x2-support-251.R.)
#   5. No-never-treated single-cohort panels fail via TWO distinct guards:
#      adopting at t = 2 trips the time-period check (truncation leaves 1
#      period); adopting at t >= 3 trips the cohort check in
#      .truncate_if_no_never_treated() (truncation leaves 0 retained cohorts).
#   6. An R >= 2 no-never-treated panel that truncates to a SINGLE retained
#      cohort now SUCCEEDS as a valid G = 1 sub-panel (the real regression
#      target for the utility.R `length(retained_cohorts) < 1` relaxation).
#   7. The three accessors (eventStudy / cohortStudy / simultaneousCIs) and the
#      conservative / cluster SE paths all degrade correctly at G = 1.
#
# Mutation-checks for each relaxed/added guard (re-tighten the source guard ->
# the named test FAILs) were performed by the implementer and recorded in the
# PR; they are noted in comments next to the relevant tests below.
#
# All assertions are on the canonical `$G` (never the deprecated `$R` compat
# field, #41/#201).
# ==============================================================================

# ------------------------------------------------------------------------------
# Shared fixtures.
# ------------------------------------------------------------------------------

# A simulateData(G = 1) panel with a known overall ATT (= the single cohort's
# effect). Offset 2 (earliest adoption), T = 6 -> num_treats = 5.
.sc112_sim <- function(
	N = 200,
	T = 6,
	d = 2,
	density = 0.5,
	eff_size = 2,
	seed = 20260605
) {
	coefs <- genCoefs(
		G = 1,
		T = T,
		d = d,
		density = density,
		eff_size = eff_size,
		seed = seed
	)
	tes <- getTes(coefs)
	sim <- simulateData(coefs, N = N, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	list(sim = sim, att_true = tes$att_true)
}

.sc112_fit <- function(sim, estimator = fetwfe, ...) {
	estimator(
		pdata = sim$pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = sim$covs,
		response = "y",
		verbose = FALSE,
		...
	)
}

# A hand-built single-cohort panel: `frac` of the units adopt at calendar time
# `adopt` (offset = adopt), the rest are never-treated. When `tau` is non-zero
# the treated post-periods get an additive effect of size `tau`, so ATT
# recovery is non-tautological. The seed stream is fixed for reproducibility.
.sc112_make_panel <- function(
	N = 240,
	T = 6,
	adopt = 3,
	tau = 0,
	frac = 0.5,
	seed = 99
) {
	set.seed(seed)
	units <- sprintf("u%03d", seq_len(N))
	treated_units <- units[seq_len(round(N * frac))]
	rows <- do.call(
		rbind,
		lapply(units, function(u) {
			is_tr <- u %in% treated_units
			unit_effect <- rnorm(1)
			time_shocks <- rnorm(T)
			cov1 <- rnorm(T)
			cov2 <- stats::runif(T)
			treat <- as.integer(is_tr & (seq_len(T) >= adopt))
			y <- unit_effect +
				time_shocks +
				0.5 * cov1 +
				tau * treat +
				rnorm(T, sd = 0.3)
			data.frame(
				time = as.integer(seq_len(T)),
				unit = as.character(u),
				treatment = treat,
				cov1 = cov1,
				cov2 = cov2,
				y = y
			)
		})
	)
	rows[order(rows$unit, rows$time), , drop = FALSE]
}

# A no-never-treated panel: ALL units adopt (no never-treated). `adopts` may
# name one cohort (single-cohort cases 5a/5b) or two (case 6).
.sc112_make_no_nt <- function(N = 200, T = 6, adopts = 2, tau = 0, seed = 7) {
	set.seed(seed)
	units <- sprintf("u%03d", seq_len(N))
	n_cohorts <- length(adopts)
	rows <- do.call(
		rbind,
		lapply(seq_along(units), function(i) {
			u <- units[i]
			a <- adopts[((i - 1L) %% n_cohorts) + 1L]
			cov1 <- rnorm(T)
			cov2 <- stats::runif(T)
			treat <- as.integer(seq_len(T) >= a)
			y <- rnorm(1) +
				rnorm(T) +
				0.5 * cov1 +
				tau * treat +
				rnorm(T, sd = 0.3)
			data.frame(
				time = as.integer(seq_len(T)),
				unit = as.character(u),
				treatment = treat,
				cov1 = cov1,
				cov2 = cov2,
				y = y
			)
		})
	)
	rows[order(rows$unit, rows$time), , drop = FALSE]
}

# ------------------------------------------------------------------------------
# 1. All four estimators recover a known ATT at G = 1.
#
# Mutation target for the cohort-count relaxation (input_prep.R `G < 1`): if the
# guard is re-tightened to `G < 2`, all four fits stop with the no-treated-cohort
# error and this test FAILs.
# ------------------------------------------------------------------------------

test_that("fetwfe recovers the known ATT at G = 1", {
	skip_on_cran()
	fx <- .sc112_sim()
	res <- .sc112_fit(fx$sim, fetwfe)
	expect_s3_class(res, "fetwfe")
	expect_equal(res$G, 1)
	expect_true(is.finite(res$att_hat))
	expect_true(is.finite(res$att_se))
	expect_gt(res$att_se, 0)
	# Numerical recovery (true ATT = 2.4; FETWFE att_hat ~ 2.41 at this fixture).
	expect_equal(res$att_hat, fx$att_true, tolerance = 0.3)
})

test_that("etwfe recovers the known ATT at G = 1", {
	skip_on_cran()
	fx <- .sc112_sim()
	res <- .sc112_fit(fx$sim, etwfe)
	expect_s3_class(res, "etwfe")
	expect_equal(res$G, 1)
	expect_true(is.finite(res$att_hat))
	expect_equal(res$att_hat, fx$att_true, tolerance = 0.3)
})

test_that("betwfe recovers the known ATT at G = 1", {
	skip_on_cran()
	fx <- .sc112_sim()
	res <- .sc112_fit(fx$sim, betwfe)
	expect_s3_class(res, "betwfe")
	expect_equal(res$G, 1)
	expect_true(is.finite(res$att_hat))
	expect_equal(res$att_hat, fx$att_true, tolerance = 0.3)
})

test_that("twfeCovs fits at G = 1 (biased by design, so finiteness only)", {
	skip_on_cran()
	fx <- .sc112_sim()
	res <- .sc112_fit(fx$sim, twfeCovs)
	expect_s3_class(res, "twfeCovs")
	expect_equal(res$G, 1)
	# twfeCovs is the naive TWFE-with-covariates benchmark; it does NOT recover
	# the ATT under staggered adoption, so we assert only finiteness here.
	expect_true(is.finite(res$att_hat))
	expect_true(is.finite(res$att_se))
})

# ------------------------------------------------------------------------------
# 2. The three fusion structures at G = 1, plus the cohort == event_study
#    invariant. With a single cohort there is no cross-cohort dimension, so the
#    cohort and event_study penalties reduce to the same within-cohort
#    dynamic-effect fusion and must give byte-identical fits.
# ------------------------------------------------------------------------------

test_that("cohort and event_study fusion are identical at G = 1", {
	skip_on_cran()
	fx <- .sc112_sim()
	res_cohort <- .sc112_fit(fx$sim, fetwfe, fusion_structure = "cohort")
	res_es <- .sc112_fit(fx$sim, fetwfe, fusion_structure = "event_study")

	expect_equal(res_cohort$G, 1)
	expect_equal(res_es$G, 1)
	# The defining invariant at G = 1: no cross-cohort fusion direction exists,
	# so the two penalties coincide exactly.
	expect_equal(res_cohort$att_hat, res_es$att_hat)
	expect_equal(res_cohort$beta_hat, res_es$beta_hat)
	expect_equal(res_cohort$catt_df$estimate, res_es$catt_df$estimate)
	expect_equal(nrow(res_cohort$catt_df), 1L)
})

test_that("a custom fusion_matrix fits at G = 1", {
	skip_on_cran()
	fx <- .sc112_sim()
	# num_treats = T - 1 = 5 for the offset-2 single cohort; a num_treats x
	# num_treats fusion matrix (here the identity, i.e. a ridge-on-theta penalty)
	# is the documented #236 shape.
	num_treats <- 5L
	D <- diag(num_treats)
	res <- .sc112_fit(fx$sim, fetwfe, fusion_matrix = D)
	expect_s3_class(res, "fetwfe")
	expect_equal(res$G, 1)
	expect_true(is.finite(res$att_hat))
})

# ------------------------------------------------------------------------------
# 3. Scattered-offset single cohort (offset 3 -> num_treats = 4). This exercises
#    the getFirstIndsFromOffsets() path (#174) at G = 1, which genCoefs /
#    simulateData cannot produce (they force offset 2).
# ------------------------------------------------------------------------------

test_that("a scattered-offset (offset 3) single cohort fits and recovers ATT", {
	skip_on_cran()
	# T = 6, single cohort adopts at t = 3 -> offset 3, num_treats = 4, with a
	# true additive effect tau = 2.
	pdata <- .sc112_make_panel(N = 240, T = 6, adopt = 3, tau = 2, seed = 311)
	res <- fetwfe(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	)
	expect_s3_class(res, "fetwfe")
	expect_equal(res$G, 1)
	# Single cohort adopting at t = 3 -> 4 post-treatment periods (t = 3..6).
	expect_true(is.finite(res$att_hat))
	expect_equal(res$att_hat, 2, tolerance = 0.3)
})

# ------------------------------------------------------------------------------
# 4. num_treats = 1 now FITS (#251 removed the #112 `num_treats >= 2` guard).
#
# Both panels previously stopped with "at least two treatment effects"; they now
# produce a valid single-cohort fit. The T = 2 case is the unique no-fusion
# configuration, so fetwfe warns there (asserted in test-2x2-support-251.R; here
# we just suppress it and confirm the fit). The late-adopting case (offset = T,
# T >= 3) still fuses the time fixed effects, so it does NOT warn.
# ------------------------------------------------------------------------------

test_that("T = 2 (num_treats = 1) now fits at G = 1", {
	# T = 2, single cohort adopting at t = 2 -> exactly one post-treatment
	# period -> num_treats = 1.
	pdata <- .sc112_make_panel(N = 60, T = 2, adopt = 2, tau = 0, seed = 12)
	res <- suppressWarnings(fetwfe(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	))
	expect_s3_class(res, "fetwfe")
	expect_equal(res$G, 1)
	expect_equal(res$T, 2)
	# Single post-treatment period -> a single treatment effect (one catt row).
	expect_equal(nrow(res$catt_df), 1L)
	expect_true(is.finite(res$att_hat))
})

test_that("a late-adopting single cohort (offset = T, num_treats = 1) now fits", {
	# T = 5, single cohort adopts at the FINAL period t = 5 -> one
	# post-treatment period -> num_treats = 1, despite T >= 3. getNumTreats(G, T)
	# would report T - 1 = 4 here; the offset-aware count produces a single
	# treatment effect.
	pdata <- .sc112_make_panel(N = 80, T = 5, adopt = 5, tau = 0, seed = 13)
	res <- fetwfe(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	)
	expect_s3_class(res, "fetwfe")
	expect_equal(res$G, 1)
	expect_equal(res$T, 5)
	# Adopting in the final period -> exactly one post-treatment period.
	expect_equal(nrow(res$catt_df), 1L)
	expect_true(is.finite(res$att_hat))
})

# ------------------------------------------------------------------------------
# 5. No-never-treated single-cohort panels fail via TWO distinct guards.
# ------------------------------------------------------------------------------

test_that("no-never-treated single cohort adopting at t = 2 trips the time-period check", {
	# All units adopt at t = 2, no never-treated. Truncating to the largest
	# all-untreated sub-panel leaves only t = 1 -> 1 time period.
	pdata <- .sc112_make_no_nt(N = 40, T = 6, adopts = 2, seed = 71)
	expect_error(
		suppressWarnings(fetwfe(
			pdata = pdata,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE
		)),
		"would leave only"
	)
})

test_that("no-never-treated single cohort adopting at t = 3 trips the cohort check", {
	# All units adopt at t = 3, no never-treated. Truncation keeps t = 1, 2
	# (>= 2 periods, so the time-period check passes) but leaves 0 retained
	# cohorts -> the .truncate_if_no_never_treated() cohort check fires.
	pdata <- .sc112_make_no_nt(N = 40, T = 6, adopts = 3, seed = 72)
	expect_error(
		suppressWarnings(fetwfe(
			pdata = pdata,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE
		)),
		"treated cohort"
	)
})

# ------------------------------------------------------------------------------
# 6. R >= 2 no-never-treated panel truncating to a SINGLE retained cohort now
#    SUCCEEDS (the real regression target for the utility.R relaxation from
#    `< 2` to `< 1`).
#
# Mutation target for the utility.R `length(retained_cohorts) < 1` relaxation:
# re-tighten to `< 2` and this test FAILs (the 1-retained-cohort truncation is
# rejected instead of fitting).
# ------------------------------------------------------------------------------

test_that("two-cohort no-never-treated panel truncating to one cohort fits (G = 1)", {
	# Half the units adopt at t = 3, half at t = 5; no never-treated. r_max = 5,
	# truncation keeps t = 1..4 (>= 2 periods) and retains the t = 3 cohort only
	# -> a valid single-cohort sub-panel.
	pdata <- .sc112_make_no_nt(
		N = 200,
		T = 6,
		adopts = c(3, 5),
		tau = 2,
		seed = 81
	)
	res <- suppressWarnings(fetwfe(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	))
	expect_s3_class(res, "fetwfe")
	# Truncated down to a single retained cohort.
	expect_equal(res$G, 1)
	expect_true(res$T < 6)
	expect_true(is.finite(res$att_hat))
})

# ------------------------------------------------------------------------------
# 7. Accessors and SE paths degrade correctly at G = 1.
# ------------------------------------------------------------------------------

test_that("eventStudy, cohortStudy, and simultaneousCIs work at G = 1", {
	skip_on_cran()
	fx <- .sc112_sim()
	res <- .sc112_fit(fx$sim, fetwfe)

	es <- eventStudy(res)
	expect_s3_class(es, "eventStudy")
	# Single cohort -> every event time is covered by exactly one cohort.
	expect_true(all(es$n_cohorts == 1L))
	expect_true(all(is.finite(es$estimate)))

	cs <- cohortStudy(res)
	expect_s3_class(cs, "cohortStudy")
	# Exactly one cohort row at G = 1.
	expect_equal(nrow(cs), 1L)

	# event_study family (K > 1) and cohort family (K = 1, the critical-value
	# bypass) both produce valid simultaneous CIs at G = 1.
	sci_es <- simultaneousCIs(res, family = "event_study")
	expect_s3_class(sci_es, "simultaneous_cis")
	expect_equal(nrow(sci_es$ci), 5L)
	expect_true(all(is.finite(sci_es$ci$simultaneous_ci_low)))
	expect_true(all(is.finite(sci_es$ci$simultaneous_ci_high)))

	sci_cohort <- simultaneousCIs(res, family = "cohort")
	expect_s3_class(sci_cohort, "simultaneous_cis")
	expect_equal(sci_cohort$K, 1L)

	sci_all <- simultaneousCIs(res, family = "all_post_treatment")
	expect_s3_class(sci_all, "simultaneous_cis")
})

test_that("conservative and cluster SE paths work at G = 1", {
	skip_on_cran()
	fx <- .sc112_sim()

	res_cons <- .sc112_fit(fx$sim, fetwfe, se_type = "conservative")
	expect_equal(res_cons$G, 1)
	expect_true(is.finite(res_cons$att_se))
	expect_gt(res_cons$att_se, 0)

	res_clus <- .sc112_fit(fx$sim, fetwfe, se_type = "cluster")
	expect_equal(res_clus$G, 1)
	expect_true(is.finite(res_clus$att_se))
	expect_gt(res_clus$att_se, 0)
})

test_that("add_ridge and REML (sig_eps_sq = NA) paths work at G = 1", {
	skip_on_cran()
	fx <- .sc112_sim()

	# add_ridge at the standard multi-period single-cohort fixture (T = 6,
	# num_treats = 5). The 2x2 add_ridge path (where the old `T >= 3L` stopifnot
	# in genFullInvFusionTransformMat used to bite) is covered in
	# test-2x2-support-251.R.
	res_ridge <- .sc112_fit(fx$sim, fetwfe, add_ridge = TRUE)
	expect_equal(res_ridge$G, 1)
	expect_true(is.finite(res_ridge$att_hat))

	# REML noise-variance estimation (sig_eps_sq = NA) is cohort-count
	# independent.
	res_reml <- .sc112_fit(
		fx$sim,
		fetwfe,
		sig_eps_sq = NA,
		sig_eps_c_sq = NA
	)
	expect_equal(res_reml$G, 1)
	expect_true(is.finite(res_reml$att_hat))
})

# ------------------------------------------------------------------------------
# 8. genCoefs / simulateData round trip at G = 1 (the DGP side of #112).
# ------------------------------------------------------------------------------

test_that("genCoefs and simulateData round-trip at G = 1", {
	coefs <- genCoefs(
		G = 1,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 4242
	)
	expect_s3_class(coefs, "FETWFE_coefs")
	expect_equal(coefs$G, 1)

	tes <- getTes(coefs)
	expect_length(tes$actual_cohort_tes, 1L)
	# Overall ATT equals the lone cohort's effect.
	expect_equal(tes$att_true, tes$actual_cohort_tes[1])

	sim <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	expect_equal(sim$G, 1)
	expect_s3_class(sim$pdata, "data.frame")
	# One treated cohort -> exactly one unique post-adoption start time among
	# the treated units.
	treated <- sim$pdata[sim$pdata$treatment == 1L, , drop = FALSE]
	first_treat <- stats::aggregate(time ~ unit, data = treated, FUN = min)
	expect_equal(length(unique(first_treat$time)), 1L)
})
