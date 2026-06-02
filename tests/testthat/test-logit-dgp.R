library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Tests for covariate-dependent cohort-assignment DGPs introduced in 1.14.0
# (issue #162) and the assignment_interactions augmentation in 1.14.1
# (issue #191).
#
# Eleven tests:
#   1. Marginal DGP byte-identity (lock-in for backward compatibility).
#   2. Multinomial DGP: empirical cohort proportions match the MC expectation.
#   3. Ordered DGP: empirical cohort proportions are approximately uniform
#      (the root-find cutpoint scheme guarantees this at strength = 0 and
#      keeps it approximately uniform at higher strength).
#   4. FETWFE smoke test: fits run cleanly on multinomial and ordered DGPs
#      (well-formedness, not load-bearing for numerical correctness).
#   5. Propensity-weighted att_true under multinomial DGP (LOAD-BEARING for
#      the propensity-weighting code path; uses an INDEPENDENT MC reference
#      computed inside the test with a different seed).
#   6. High-strength sanity: covariates strongly predict cohort at
#      strength = 10 (max-row P(W|X) is well above the uniform baseline).
#   7. (#191) Interactions-augmented multinomial DGP yields propensity a
#      linear-in-X logit cannot fit (LOAD-BEARING; seed-pinned).
#   8. (#191) getTes() accounts for interactions in MC integration
#      (LOAD-BEARING; independent MC reference at fresh seed).
#   9. (#191) Marginal-uniform invariant holds under ordered DGP with
#      interactions + paired no-interactions baseline contrast
#      (LOAD-BEARING at strength = 2, tolerance = 0.015; both assertions
#      load-bearing at seed=42).
#  10. (#191) Canonicalization + deduplication of pair list (lightweight).
#  11. (#191) v1.14.0 backward-compat: hand-constructed FETWFE_coefs missing
#      the new $assignment_interaction_strength + $assignment_coefs$interactions
#      + $assignment_coefs$delta slots round-trips through simulateData()
#      and getTes() without error.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Test 1: marginal DGP byte-identity
# ------------------------------------------------------------------------------
test_that("marginal DGP: byte-identical to v1.13.x at default assignment_type", {
	# Reference values captured from the marginal-path output at branch
	# `feat-cohort-dgp-162` against the same call on `origin/main` SHA
	# ae515b6. Locks in the pre-1.14.0 behavior at the strictest possible
	# assertion. Any drift in the marginal-path RNG ordering fires here.
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 42
	)
	expect_identical(coefs$assignment_type, "marginal")
	expect_null(coefs$assignment_coefs)

	sim <- simulateData(coefs, N = 100, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	expect_equal(nrow(sim$pdata), 100 * 5)
	expect_identical(sim$assignments, c(31L, 29L, 18L, 22L))
	expect_identical(sim$indep_counts, c(24L, 30L, 19L, 27L))
	expect_equal(
		sim$pdata$y[1:3],
		c(1.249801325615, -3.167719220623, -3.455029448642),
		tolerance = 1e-9
	)
})

# ------------------------------------------------------------------------------
# Test 2: multinomial DGP — empirical cohort proportions match MC expectation
# ------------------------------------------------------------------------------
test_that("multinomial DGP: empirical cohort proportions converge to MC expectation", {
	skip_on_cran()
	# N = 10000 gives sqrt(p(1-p)/N) SE < 0.005 per cohort; tolerance 0.03
	# is comfortable.
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 1.0,
		seed = 42
	)
	sim <- simulateData(coefs, N = 10000, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	empirical_probs <- sim$assignments / sum(sim$assignments)
	expected_probs <- fetwfe:::.expected_cohort_probs(
		coefs$assignment_coefs,
		d = 2,
		distribution = "gaussian",
		M = 50000L,
		seed = 99L
	)
	expect_equal(
		as.numeric(empirical_probs),
		as.numeric(expected_probs),
		tolerance = 0.03
	)
})

# ------------------------------------------------------------------------------
# Test 3: ordered DGP — empirical cohort proportions ~ uniform 1/(R+1)
# ------------------------------------------------------------------------------
test_that("ordered DGP: empirical cohort proportions match cutpoint structure", {
	skip_on_cran()
	# The root-find cutpoint scheme produces approximately uniform 1/(R+1)
	# marginal probabilities at any strength (exactly uniform at
	# strength = 0; ~0.002 deviation at strength = 1 per round-2 plan
	# verification).
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "ordered",
		assignment_strength = 1.0,
		seed = 42
	)
	sim <- simulateData(coefs, N = 10000, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	empirical_probs <- sim$assignments / sum(sim$assignments)
	expect_equal(
		as.numeric(empirical_probs),
		rep(1 / 4, 4),
		tolerance = 0.03
	)
})

# ------------------------------------------------------------------------------
# Test 4: FETWFE smoke test on multinomial + ordered DGPs
# ------------------------------------------------------------------------------
test_that("FETWFE smoke test: fits run cleanly on multinomial and ordered DGPs", {
	skip_on_cran()
	# Well-formedness only (per round-1 sentinel NOTE: not load-bearing
	# for numerical correctness; that's covered by test #5 below).
	for (type in c("multinomial", "ordered")) {
		coefs <- genCoefs(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			assignment_type = type,
			assignment_strength = 1.0,
			seed = 42
		)
		sim <- simulateData(
			coefs,
			N = 200,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5
		)
		fit <- suppressWarnings(fetwfe(
			pdata = sim$pdata,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			response = "y",
			covs = c("cov1", "cov2")
		))
		expect_true(
			is.finite(fit$att_hat),
			info = paste0("att_hat not finite for type = ", type)
		)
		expect_true(
			is.finite(fit$att_se),
			info = paste0("att_se not finite for type = ", type)
		)
		expect_gt(fit$att_se, 0)
	}
})

# ------------------------------------------------------------------------------
# Test 5: propensity-weighted att_true under multinomial DGP (LOAD-BEARING)
# ------------------------------------------------------------------------------
test_that("getTes returns propensity-weighted att_true under multinomial DGP", {
	# Compares getTes()'s production output against an INDEPENDENT MC
	# integration computed inside the test with a different seed
	# (99L vs the production-side coefs$seed + 2L = 44L). This is the
	# load-bearing check for the propensity-weighting code path; it
	# verifies the production code does compute
	#   att_true = sum(cohort_weights * actual_cohort_tes)
	# with cohort_weights = E[pi_r(X)] / sum E[pi_r(X)] (treated only).
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 1.5,
		seed = 42
	)
	tes <- getTes(coefs)
	# Independent reference computed with a fresh MC sample under a
	# different seed.
	set.seed(99L)
	X_ref <- matrix(rnorm(50000 * 2), nrow = 50000, ncol = 2)
	probs_ref <- fetwfe:::.compute_cohort_prob_matrix(
		X_ref,
		coefs$assignment_coefs
	)
	expected_ref <- colMeans(probs_ref)
	treated_probs_ref <- expected_ref[2:4]
	cohort_weights_ref <- as.numeric(
		treated_probs_ref / sum(treated_probs_ref)
	)
	att_ref <- sum(cohort_weights_ref * tes$actual_cohort_tes)
	expect_equal(tes$att_true, att_ref, tolerance = 0.02)
	expect_equal(
		tes$cohort_weights,
		cohort_weights_ref,
		tolerance = 0.02
	)
	# Under non-marginal DGP, cohort_weights should NOT be uniform
	# (otherwise we'd not be exercising the propensity-weighting path).
	expect_false(
		isTRUE(all.equal(tes$cohort_weights, rep(1 / 3, 3), tolerance = 0.01))
	)
})

# ------------------------------------------------------------------------------
# Test 6: high-strength sanity — covariates strongly predict cohort
# ------------------------------------------------------------------------------
test_that("high-strength assignment: covariates strongly predict cohort", {
	# At assignment_strength = 10, cohort assignment becomes nearly
	# deterministic given X. Verify: max-row-probability of P(W|X) is
	# well above the uniform-baseline 1/(R+1) = 0.25 at most units.
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 10.0,
		seed = 42
	)
	X_test <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
	probs <- fetwfe:::.compute_cohort_prob_matrix(
		X_test,
		coefs$assignment_coefs
	)
	max_row_probs <- apply(probs, 1, max)
	expect_gt(mean(max_row_probs), 0.85)
})

# ------------------------------------------------------------------------------
# Test 7 (#191, LOAD-BEARING, seed-pinned): interactions-augmented multinomial
# DGP yields propensity a linear-in-X logit cannot fit.
#
# Fit per-cohort one-vs-rest binomial logits on the simulated data — once with
# `~ cov1 + cov2` (linear-in-X), once with `~ cov1 + cov2 + I(cov1 * cov2)`
# (augmented). At the pinned seed = 42 with assignment_strength =
# assignment_interaction_strength = 2.0 the linear-truth gap is empirically
# 67 and the augmented gap is empirically 2005; the threshold of 500 provides
# seed-specific separation (not cross-seed safety — the test is seed-pinned).
# The cross-seed linear-truth gap varies in [18, 593] (30-seed round-2 sweep)
# but the seed-pinned test pass/fail is unaffected.
# ------------------------------------------------------------------------------
test_that("interactions-augmented multinomial DGP injects propensity nonlinearity (#191)", {
	skip_on_cran()
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 2.0,
		assignment_interactions = list(c(1L, 2L)),
		assignment_interaction_strength = 2.0,
		seed = 42
	)
	sim <- suppressWarnings(simulateData(
		coefs,
		N = 5000,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	))
	# Per-unit cohort labels (first observation per unit; treatment is an
	# absorbing state, so the unit's cohort is encoded in the panel's
	# `treatment` column — we extract via first-period rows).
	pdata <- sim$pdata
	first_period <- min(pdata$time)
	first_rows <- pdata[pdata$time == first_period, ]
	# `cohort` is encoded as the unit's adoption-period label (or NA for
	# never-treated). simulateData() exposes per-unit assignment via the
	# panel structure; reconstruct W from the treatment column.
	cohort_by_unit <- tapply(
		pdata$treatment,
		pdata$unit,
		function(v) {
			treated_inds <- which(v == 1)
			if (length(treated_inds) == 0L) {
				0L
			} else {
				as.integer(min(treated_inds))
			}
		}
	)
	# Align X with cohort label.
	X1 <- first_rows$cov1
	X2 <- first_rows$cov2
	cohort_aligned <- cohort_by_unit[as.character(first_rows$unit)]

	# One-vs-rest binomial logits, summed deviance across cohorts.
	linear_dev <- 0
	augmented_dev <- 0
	cohort_levels <- unique(cohort_aligned)
	treated_levels <- sort(cohort_levels[cohort_levels > 0L])
	for (r in treated_levels) {
		y <- as.integer(cohort_aligned == r)
		df_r <- data.frame(y = y, cov1 = X1, cov2 = X2)
		fit_lin <- suppressWarnings(stats::glm(
			y ~ cov1 + cov2,
			data = df_r,
			family = stats::binomial()
		))
		fit_aug <- suppressWarnings(stats::glm(
			y ~ cov1 + cov2 + I(cov1 * cov2),
			data = df_r,
			family = stats::binomial()
		))
		linear_dev <- linear_dev + fit_lin$deviance
		augmented_dev <- augmented_dev + fit_aug$deviance
	}
	# A meaningful deviance reduction means the augmented model explains
	# nonlinearity that the linear-in-X logit cannot. Threshold 500 chosen
	# for seed-specific separation at seed = 42.
	expect_gt(linear_dev - augmented_dev, 500)
})

# ------------------------------------------------------------------------------
# Test 8 (#191, LOAD-BEARING): getTes() accounts for interactions in MC
# integration.
#
# Compare getTes()'s production output against an INDEPENDENT MC reference
# computed inside the test with a fresh seed (99L), matching the precedent
# set by Test 5. Pinned at strength = int_strength = 2.0, seed = 42 for
# consistency with Tests 7 + 9.
# ------------------------------------------------------------------------------
test_that("getTes() accounts for interactions in MC integration (#191)", {
	skip_on_cran()
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 2.0,
		assignment_interactions = list(c(1L, 2L)),
		assignment_interaction_strength = 2.0,
		seed = 42
	)
	tes <- getTes(coefs)
	# Independent reference computed with a fresh MC sample.
	set.seed(99L)
	X_ref <- matrix(rnorm(50000 * 2), nrow = 50000, ncol = 2)
	probs_ref <- fetwfe:::.compute_cohort_prob_matrix(
		X_ref,
		coefs$assignment_coefs
	)
	expected_ref <- colMeans(probs_ref)
	treated_probs_ref <- expected_ref[2:4]
	cohort_weights_ref <- as.numeric(
		treated_probs_ref / sum(treated_probs_ref)
	)
	att_ref <- sum(cohort_weights_ref * tes$actual_cohort_tes)
	expect_equal(tes$att_true, att_ref, tolerance = 0.02)
	expect_equal(
		tes$cohort_weights,
		cohort_weights_ref,
		tolerance = 0.02
	)
})

# ------------------------------------------------------------------------------
# Test 9 (#191, LOAD-BEARING at strength = 2 + tol = 0.015): marginal-uniform
# invariant holds under ordered DGP with interactions; paired with a
# no-interactions baseline contrast.
#
# Two assertions, both load-bearing at seed = 42:
# (a) Marginal-uniform: at strength = int_strength = 2.0, the buggy un-
#     augmented cutpoint solver gives max_dev = 0.0188 (fails tolerance =
#     0.015) while the correct path gives max_dev = 0.0097 (passes with ~55%
#     headroom). Catches Failure-mode A (root-finder forgets augmentation).
# (b) Baseline contrast: a no-interactions baseline panel at the same seed
#     differs from the with-interactions panel by mean|diff| ~ 10 at seed =
#     42 (correct case). Catches Failure-mode C (both root-finder and
#     runtime silently drop the augmentation -> WITH and WITHOUT are
#     byte-identical -> mean|diff| == 0).
# ------------------------------------------------------------------------------
test_that("ordered DGP with interactions preserves marginal uniformity + baseline contrast (#191)", {
	skip_on_cran()
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "ordered",
		assignment_strength = 2.0,
		assignment_interactions = list(c(1L, 2L)),
		assignment_interaction_strength = 2.0,
		seed = 42
	)
	sim_ord <- suppressWarnings(simulateData(
		coefs,
		N = 10000,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	))
	empirical_probs <- sim_ord$assignments / sum(sim_ord$assignments)
	max_dev <- max(abs(as.numeric(empirical_probs) - 1 / 4))
	# (a) Marginal-uniform invariant.
	expect_lt(max_dev, 0.015)

	# (b) No-interactions baseline contrast — identical inputs except
	# assignment_interactions = NULL. Same seed for reproducibility.
	coefs_baseline <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "ordered",
		assignment_strength = 2.0,
		seed = 42
	)
	sim_ord_baseline <- suppressWarnings(simulateData(
		coefs_baseline,
		N = 10000,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	))
	diff_means <- mean(abs(
		as.numeric(sim_ord$assignments) -
			as.numeric(sim_ord_baseline$assignments)
	))
	# Threshold of 5 has ~2x headroom against the empirical correct-case
	# value of ~10 at seed = 42; fires at exactly 0 under Failure-mode C.
	expect_gt(diff_means, 5)
})

# ------------------------------------------------------------------------------
# Test 10 (#191, lightweight): canonicalization + deduplication of pair list.
# ------------------------------------------------------------------------------
test_that("assignment_interactions canonicalizes + dedupes pairs (#191)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_interactions = list(c(2, 1), c(1, 2)),
		seed = 42
	)
	expect_identical(
		coefs$assignment_coefs$interactions,
		list(c(1L, 2L))
	)

	# Self-interactions (quadratic) are allowed.
	coefs_quad <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_interactions = list(c(1, 1)),
		seed = 42
	)
	expect_identical(
		coefs_quad$assignment_coefs$interactions,
		list(c(1L, 1L))
	)

	# Index outside [1, d] errors with the informative message.
	expect_error(
		genCoefs(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			assignment_type = "multinomial",
			assignment_interactions = list(c(0, 1))
		),
		regexp = "outside \\[1, d"
	)
	expect_error(
		genCoefs(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			assignment_type = "multinomial",
			assignment_interactions = list(c(3, 1))
		),
		regexp = "outside \\[1, d"
	)

	# verbose = TRUE emits a dedup message when canonicalization removes
	# pairs (per post-execution review item #6 maintainer-elevation:
	# verbose-gated message preserves silent default while giving users
	# an opt-in signal that dedup happened).
	expect_message(
		genCoefs(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			assignment_type = "multinomial",
			assignment_interactions = list(c(2, 1), c(1, 2)),
			verbose = TRUE,
			seed = 42
		),
		regexp = "Deduplicated 1 assignment_interactions pair"
	)

	# verbose = FALSE (default) stays silent even when dedup happens.
	expect_silent(
		genCoefs(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			assignment_type = "multinomial",
			assignment_interactions = list(c(2, 1), c(1, 2)),
			verbose = FALSE,
			seed = 42
		)
	)
})

# ------------------------------------------------------------------------------
# Test 11 (#191, lightweight): v1.14.0 backward-compat round-trip.
#
# Mirrors the v1.13.x backward-compat defensive pattern at
# R/gen_data.R:122-129 and R/gen_coefs.R:345-352. A hand-constructed
# FETWFE_coefs object missing the new 1.14.1 slots ($interaction_strength,
# $assignment_coefs$interactions, $assignment_coefs$delta) must round-trip
# through simulateData() and getTes() without error so the package's
# defensive is.null() handling does not regress.
# ------------------------------------------------------------------------------
test_that("v1.14.0 FETWFE_coefs round-trips through simulateData + getTes (#191)", {
	# Construct a minimal v1.14.0-shape multinomial coefs object. We use a
	# small R = 3, T = 5, d = 2 panel matching the test surface.
	R <- 3L
	T <- 5L
	d <- 2L

	# Pull beta + theta off a freshly-constructed v1.14.1 object, then
	# strip the new top-level + sub-slots to simulate a v1.14.0 object.
	template <- genCoefs(
		G = R,
		T = T,
		d = d,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 1.0,
		seed = 42
	)
	mock <- list(
		beta = template$beta,
		theta = template$theta,
		R = R,
		T = T,
		d = d,
		seed = 42L,
		assignment_type = "multinomial",
		assignment_strength = 1.0,
		# v1.14.0 NO $assignment_interaction_strength,
		# v1.14.0 NO $assignment_coefs$interactions / $delta / $interaction_strength
		assignment_coefs = list(
			type = "multinomial",
			strength = 1.0,
			coefs = template$assignment_coefs$coefs
		)
	)
	class(mock) <- "FETWFE_coefs"

	# Both calls must succeed (the defensive is.null() handling in the
	# downstream code catches missing slots).
	sim <- expect_no_error(
		suppressWarnings(simulateData(
			mock,
			N = 100,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5
		))
	)
	tes <- expect_no_error(getTes(mock))
	# Sanity: returned objects have the expected shapes.
	expect_true(is.finite(tes$att_true))
	expect_equal(length(tes$actual_cohort_tes), R)
})
