# Tests for issues #141 + #146 (v1.12.0): the same-data Gaussianity upgrade.
#
# Three things to lock in:
#
# 1. The default `att_se` value uses the TIGHT Gaussian variance
#    `sqrt(att_var_1 + att_var_2)` from paper Theorem (c').
#    Verify: byte-equal to the direct computation from the exposed
#    `att_var_1` / `att_var_2` slots.
#
# 2. `se_type = "conservative"` recovers the OLD default (Cauchy-Schwarz
#    upper bound). Verify: byte-equal to
#    `sqrt(att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2))`.
#    Also verify the conservative SE is STRICTLY WIDER than the default
#    (regression guard against either branch silently doing the wrong
#    thing).
#
# 3. The new `internal$variance_components` slot is populated correctly
#    with the paper-notation names (`V_1`, `V_2`, `tilde_v_N`, the four
#    catalogued `tilde_v_N_*` variants).
#
# Plus a coverage sanity check (skip_on_cran) -- the substantive
# validation that the tight Gaussian default lands near nominal 95%
# coverage on simulated data with the (Psi-IF)-satisfying default
# cohort-sample-proportions estimator.

test_that("default `att_se` is the tight Gaussian variance sqrt(V1 + V2)", {
	# Use a fixture without indep_counts so the same-data path is
	# exercised (otherwise indep_probs = TRUE forces the tight formula
	# regardless of se_type and the test wouldn't catch a default-flipped
	# regression).
	set.seed(141)
	sim_coefs <- genCoefs(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 141
	)
	sim <- simulateData(
		sim_coefs,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	res <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		covs = sim$covs,
		response = sim$response,
		sig_eps_sq = sim$sig_eps_sq,
		sig_eps_c_sq = sim$sig_eps_c_sq,
		q = 0.5
	)

	expect_false(res$indep_counts_used)
	expect_identical(res$se_type, "default")

	vc <- res$internal$variance_components
	expect_true(is.list(vc))
	expect_true(is.finite(vc$att_var_1))
	expect_true(is.finite(vc$att_var_2))
	expect_gt(vc$att_var_1, 0)
	expect_gt(vc$att_var_2, 0)

	expected_default <- sqrt(vc$att_var_1 + vc$att_var_2)
	expect_equal(res$att_se, expected_default, tolerance = 1e-12)
})

test_that("se_type = 'conservative' recovers the Cauchy-Schwarz bound (pre-v1.12.0 default)", {
	set.seed(141)
	sim_coefs <- genCoefs(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 141
	)
	sim <- simulateData(
		sim_coefs,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	common_args <- list(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		covs = sim$covs,
		response = sim$response,
		sig_eps_sq = sim$sig_eps_sq,
		sig_eps_c_sq = sim$sig_eps_c_sq,
		q = 0.5
	)

	res_def <- do.call(fetwfe, common_args)
	res_cons <- do.call(fetwfe, c(common_args, list(se_type = "conservative")))

	expect_identical(res_def$se_type, "default")
	expect_identical(res_cons$se_type, "conservative")

	# Point estimates should be identical -- only the SE combination
	# differs.
	expect_equal(res_def$att_hat, res_cons$att_hat)
	expect_equal(res_def$catt_hats, res_cons$catt_hats)

	# Verify the conservative formula numerically.
	vc <- res_cons$internal$variance_components
	expected_cons <- sqrt(
		vc$att_var_1 +
			vc$att_var_2 +
			2 * sqrt(vc$att_var_1 * vc$att_var_2)
	)
	expect_equal(res_cons$att_se, expected_cons, tolerance = 1e-12)

	# The conservative SE must be STRICTLY wider than the default
	# (regression guard).
	expect_gt(res_cons$att_se, res_def$att_se)
})

test_that("variance_components block exposes paper-notation V_1, V_2, and tilde_v_N family", {
	set.seed(141)
	sim_coefs <- genCoefs(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 141
	)
	sim <- simulateData(
		sim_coefs,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	res <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		covs = sim$covs,
		response = sim$response,
		sig_eps_sq = sim$sig_eps_sq,
		sig_eps_c_sq = sim$sig_eps_c_sq,
		q = 0.5
	)

	vc <- res$internal$variance_components
	expected_slots <- c(
		"att_var_1",
		"att_var_2",
		"V_1",
		"V_2",
		"tilde_v_N",
		"hat_v_N",
		"tilde_v_N_C",
		"tilde_v_N_C_pi_hat",
		"tilde_v_N_C_pi_hat_cons",
		"tilde_v_N_cons",
		"se_type",
		"indep_counts_used"
	)
	expect_true(setequal(names(vc), expected_slots))

	# V_1 = N * att_var_1 (per-unit paper scale).
	expect_equal(vc$V_1, res$N * vc$att_var_1)
	expect_equal(vc$V_2, res$N * vc$att_var_2)

	# tilde_v_N := hat_v_N / T (paper line 2004), hat_v_N = T * tilde_v_N.
	expect_equal(vc$hat_v_N, res$T * vc$tilde_v_N)

	# Headline `tilde_v_N` matches the variance underlying `att_se`:
	# att_se = sqrt(tilde_v_N / N) per paper Eq. conf.int.form.
	expect_equal(res$att_se, sqrt(vc$tilde_v_N / res$N), tolerance = 1e-12)

	# The six catalogued unit-scaled variance estimators.
	# `tilde_v_N_C` is the fixed-pi exact (no propensity-score noise):
	# equal to V_1 = N * att_var_1.
	expect_equal(vc$tilde_v_N_C, res$N * vc$att_var_1)
	# `tilde_v_N_C_pi_hat` is the random-pi exact (tight Gaussian under
	# (Psi-IF)): equal to V_1 + V_2.
	expect_equal(vc$tilde_v_N_C_pi_hat, res$N * (vc$att_var_1 + vc$att_var_2))
	# `tilde_v_N_C_pi_hat_cons` is the random-pi conservative Cauchy-
	# Schwarz upper bound: V_1 + V_2 + 2 sqrt(V_1 V_2) at the per-unit
	# scale.
	expect_equal(
		vc$tilde_v_N_C_pi_hat_cons,
		res$N *
			(vc$att_var_1 +
				vc$att_var_2 +
				2 * sqrt(vc$att_var_1 * vc$att_var_2))
	)
	# `tilde_v_N_cons` numerically matches `tilde_v_N_C_pi_hat_cons` at
	# the per-cohort-proportions weighting the package uses.
	expect_equal(vc$tilde_v_N_cons, vc$tilde_v_N_C_pi_hat_cons)

	# Echo'd metadata.
	expect_identical(vc$se_type, "default")
	expect_false(vc$indep_counts_used)
})

test_that("variance_components block: conservative se_type swaps the headline tilde_v_N", {
	set.seed(141)
	sim_coefs <- genCoefs(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 141
	)
	sim <- simulateData(
		sim_coefs,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	res_def <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		covs = sim$covs,
		response = sim$response,
		sig_eps_sq = sim$sig_eps_sq,
		sig_eps_c_sq = sim$sig_eps_c_sq,
		q = 0.5
	)
	res_cons <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		covs = sim$covs,
		response = sim$response,
		sig_eps_sq = sim$sig_eps_sq,
		sig_eps_c_sq = sim$sig_eps_c_sq,
		q = 0.5,
		se_type = "conservative"
	)

	vc_def <- res_def$internal$variance_components
	vc_cons <- res_cons$internal$variance_components

	# Under default, headline `tilde_v_N` = `tilde_v_N_C_pi_hat`.
	expect_equal(vc_def$tilde_v_N, vc_def$tilde_v_N_C_pi_hat)
	# Under conservative, headline `tilde_v_N` = `tilde_v_N_C_pi_hat_cons`.
	expect_equal(vc_cons$tilde_v_N, vc_cons$tilde_v_N_C_pi_hat_cons)

	# But the underlying `att_var_1`, `att_var_2`, `V_1`, `V_2`, and the
	# full catalogue of `tilde_v_N_*` slots are identical across the two
	# fits -- only the combination chosen as the headline changes.
	expect_equal(vc_def$att_var_1, vc_cons$att_var_1)
	expect_equal(vc_def$att_var_2, vc_cons$att_var_2)
	expect_equal(vc_def$V_1, vc_cons$V_1)
	expect_equal(vc_def$V_2, vc_cons$V_2)
	expect_equal(vc_def$tilde_v_N_C, vc_cons$tilde_v_N_C)
	expect_equal(vc_def$tilde_v_N_C_pi_hat, vc_cons$tilde_v_N_C_pi_hat)
	expect_equal(
		vc_def$tilde_v_N_C_pi_hat_cons,
		vc_cons$tilde_v_N_C_pi_hat_cons
	)
})

test_that("variance_components block: indep_counts path uses tight formula regardless of se_type", {
	# Two-sample path: both `default` and `conservative` should fall
	# through to the asymptotically-exact `sqrt(att_var_1 + att_var_2)`
	# formula.
	set.seed(141)
	sim_coefs <- genCoefs(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 141
	)
	sim <- simulateData(
		sim_coefs,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	# fetwfeWithSimulatedData passes indep_counts (from the simulator).
	res_def <- fetwfeWithSimulatedData(sim, q = 0.5)
	res_cons <- fetwfeWithSimulatedData(
		sim,
		q = 0.5,
		se_type = "conservative"
	)

	expect_true(res_def$indep_counts_used)
	expect_true(res_cons$indep_counts_used)

	# `att_se` is identical: the indep path uses the tight formula in
	# both branches.
	expect_equal(res_def$att_se, res_cons$att_se)

	# The metadata still echoes the user's se_type choice (even though
	# numerically inert).
	expect_identical(res_def$internal$variance_components$se_type, "default")
	expect_identical(
		res_cons$internal$variance_components$se_type,
		"conservative"
	)
})

test_that("variance_components block: q >= 1 path returns NA-filled slots", {
	# Under ridge (q = 2), `calc_ses = FALSE` and `att_se` is NA. The
	# variance-components block should be NA-filled in lockstep.
	set.seed(141)
	sim_coefs <- genCoefs(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 141
	)
	sim <- simulateData(
		sim_coefs,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	res <- fetwfeWithSimulatedData(sim, q = 2)
	expect_true(is.na(res$att_se))

	vc <- res$internal$variance_components
	expect_true(is.na(vc$att_var_1))
	expect_true(is.na(vc$att_var_2))
	expect_true(is.na(vc$V_1))
	expect_true(is.na(vc$V_2))
	expect_true(is.na(vc$tilde_v_N))
	expect_true(is.na(vc$hat_v_N))
	expect_true(is.na(vc$tilde_v_N_C))
	expect_true(is.na(vc$tilde_v_N_C_pi_hat))
	expect_true(is.na(vc$tilde_v_N_C_pi_hat_cons))
	expect_true(is.na(vc$tilde_v_N_cons))
})

test_that("etwfe, betwfe, twfeCovs all expose the variance_components block", {
	# Parity check across the four estimator classes (issue #144 made
	# `$internal` a canonical access path for all four; #141/#146
	# extends this to `$internal$variance_components`).
	set.seed(141)
	sim_coefs <- genCoefs(
		G = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 141
	)
	sim <- simulateData(
		sim_coefs,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	fits <- list(
		fetwfe = fetwfeWithSimulatedData(sim, q = 0.5),
		etwfe = etwfeWithSimulatedData(sim),
		betwfe = betwfeWithSimulatedData(sim, q = 0.5),
		twfeCovs = twfeCovsWithSimulatedData(sim)
	)

	for (cls in names(fits)) {
		vc <- fits[[cls]]$internal$variance_components
		expect_true(
			is.list(vc),
			info = paste0(cls, ": variance_components is not a list")
		)
		expected_slots <- c(
			"att_var_1",
			"att_var_2",
			"V_1",
			"V_2",
			"tilde_v_N",
			"hat_v_N",
			"tilde_v_N_C",
			"tilde_v_N_C_pi_hat",
			"tilde_v_N_C_pi_hat_cons",
			"tilde_v_N_cons",
			"se_type",
			"indep_counts_used"
		)
		expect_true(
			setequal(names(vc), expected_slots),
			info = paste0(cls, ": variance_components slot list mismatch")
		)
	}
})

# ------------------------------------------------------------------------------
# Coverage sanity check: the substantive validation of #141.
#
# Under the new tight Gaussian default, empirical 95% CI coverage on
# simulated data with the (Psi-IF)-satisfying default cohort-sample-
# proportions estimator should land near nominal (>= 0.92 over 60 reps
# at N = 500 -- consistent with the SPEC's target). Pre-v1.12.0 the
# conservative default over-covered (typically 1.00 over the same
# fixture); the new default lands near 0.95.
#
# We use etwfe() (OLS, no bridge-penalty shrinkage) so coverage isolates
# the variance-formula effect from the bridge-regression bias. The
# fetwfe()/betwfe() bridge-penalty fits have non-trivial finite-sample
# bias that competes with the SE for coverage (independent of the SE
# formula); see paper Theorem 6.1's asymptotic-unbiasedness statement
# vs the finite-sample behavior at N = 500.
#
# Skip on CRAN (slow; runs ~50 etwfe fits at N=500).
# ------------------------------------------------------------------------------

test_that("tight Gaussian default delivers near-nominal 95% CI coverage on simulated data", {
	skip_on_cran()

	n_reps <- 60
	N_per <- 500
	covered_default <- logical(n_reps)
	covered_conservative <- logical(n_reps)
	att_ses_default <- numeric(n_reps)
	att_ses_conservative <- numeric(n_reps)
	att_hats <- numeric(n_reps)
	true_atts <- numeric(n_reps)

	z <- qnorm(0.975)
	for (i in seq_len(n_reps)) {
		seed_i <- 1000L + i
		coefs <- genCoefs(
			G = 3,
			T = 6,
			d = 2,
			density = 0.5,
			eff_size = 0.5,
			seed = seed_i
		)
		sim <- simulateData(
			coefs,
			N = N_per,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5
		)

		# Population ATT is `cohort_tes %*% pi_cond`, with `pi_cond` the
		# conditional cohort-membership proportions in the realized
		# sample. This is the estimand `att_hat` targets (paper's
		# `att.def`).
		te_truth <- getTes(coefs)
		ic <- sim$indep_counts
		pi_cond <- ic[2:length(ic)] / sum(ic[2:length(ic)])
		true_att <- sum(te_truth$actual_cohort_tes * pi_cond)

		# etwfe() so the SE-formula effect is isolated from bridge
		# bias. Same-data path (no indep_counts passed) so we exercise
		# the (Psi-IF)-tight-vs-conservative branch.
		res_def <- etwfe(
			pdata = sim$pdata,
			time_var = sim$time_var,
			unit_var = sim$unit_var,
			treatment = sim$treatment,
			covs = sim$covs,
			response = sim$response,
			sig_eps_sq = sim$sig_eps_sq,
			sig_eps_c_sq = sim$sig_eps_c_sq
		)
		res_cons <- etwfe(
			pdata = sim$pdata,
			time_var = sim$time_var,
			unit_var = sim$unit_var,
			treatment = sim$treatment,
			covs = sim$covs,
			response = sim$response,
			sig_eps_sq = sim$sig_eps_sq,
			sig_eps_c_sq = sim$sig_eps_c_sq,
			se_type = "conservative"
		)

		covered_default[i] <- (true_att >=
			res_def$att_hat - z * res_def$att_se) &&
			(true_att <= res_def$att_hat + z * res_def$att_se)
		covered_conservative[i] <- (true_att >=
			res_cons$att_hat - z * res_cons$att_se) &&
			(true_att <= res_cons$att_hat + z * res_cons$att_se)
		att_ses_default[i] <- res_def$att_se
		att_ses_conservative[i] <- res_cons$att_se
		att_hats[i] <- res_def$att_hat
		true_atts[i] <- true_att
	}

	emp_coverage_default <- mean(covered_default)
	emp_coverage_conservative <- mean(covered_conservative)

	# Target: tight Gaussian default near 0.95. Lower bound 0.92 matches
	# the SPEC; with N_per = 500 and n_reps = 60 the binomial 1-sigma
	# variability around 0.95 is ~0.03, so 0.92 is a soft floor.
	expect_gte(
		emp_coverage_default,
		0.92,
		label = sprintf(
			"default-tight empirical coverage = %.3f (N_per=%d, n_reps=%d, mean att_se=%.4f)",
			emp_coverage_default,
			N_per,
			n_reps,
			mean(att_ses_default)
		)
	)
	# Tight Gaussian shouldn't over-cover catastrophically (would
	# suggest residual conservative-bound contamination somewhere).
	expect_lte(
		emp_coverage_default,
		0.995,
		label = sprintf(
			"default-tight empirical coverage = %.3f (suspiciously high)",
			emp_coverage_default
		)
	)
	# Conservative CIs should over-cover relative to the tight CIs:
	# they're STRICTLY WIDER on every rep (same point estimate, wider
	# SE), so the conservative coverage is the tight coverage plus the
	# marginal coverage from the wider intervals. >= the tight rate is
	# the structural guarantee.
	expect_gte(
		emp_coverage_conservative,
		emp_coverage_default,
		label = sprintf(
			"conservative coverage (%.3f) should be >= default coverage (%.3f)",
			emp_coverage_conservative,
			emp_coverage_default
		)
	)
})
