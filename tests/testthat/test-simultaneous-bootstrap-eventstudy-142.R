# Tests for the event_study family of the simultaneousCIs() multiplier bootstrap
# (#142, Phase 3). event_study is the one family whose simultaneous-band variance
# carries a non-zero cohort-probability (propensity) channel Sigma_2: each
# event-time effect pools across cohorts weighted by the estimated cohort
# probabilities pi_hat_g = N_g/N, so estimating those probabilities contributes
# variance. The bootstrap perturbs the per-unit propensity influence function
# F_pi[i,k] = T * a_k'(e_{W_i} - pi_hat) (a_k = J_k theta_sel) with its OWN
# independent multiplier stream, combined with the regression channel as a
# two-channel bootstrap matching the analytic Sigma = Sigma_1 + Sigma_2 (paper
# Thm C.1 Step 5 -- zero cross-term). The load-bearing correctness anchor:
# colSums(F_pi^2)/(NT)^2 == diag(Sigma_2) (an independent code path from the
# analytic delta-method sandwich, so not a round-trip); hence on a
# se_type = "cluster" fit the bootstrap event_study SEs equal the analytic ones
# to machine precision, and the bootstrap critical value matches qmvnorm up to
# Monte-Carlo error. In fixed-p the propensity channel composes with the
# selected-support regression channel; since #299 it also composes with the
# high-dim desparsified regression channel (high-dim event_study uses both).

.es_boot_fit <- function(
	se_type = "cluster",
	G = 3,
	T = 5,
	d = 2,
	N = 150,
	density = 0.6,
	eff_size = 1.5,
	seed = 1,
	add_ridge = FALSE,
	q = 0.5,
	keep_indep_counts = FALSE
) {
	coefs <- genCoefs(
		G = G,
		T = T,
		d = d,
		density = density,
		eff_size = eff_size,
		seed = seed
	)
	dat <- simulateData(
		coefs,
		N = N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = seed
	)
	if (!keep_indep_counts) {
		dat$indep_counts <- NA
	}
	fetwfeWithSimulatedData(
		dat,
		q = q,
		se_type = se_type,
		add_ridge = add_ridge
	)
}

.es_se_from <- function(sc) {
	# per-effect SE backed out of the pointwise band
	(sc$ci$pointwise_ci_high - sc$ci$estimate) / sc$pointwise_critical_value
}

# Reconstruct the dispatch state .simultaneous_cis_impl() builds for the
# event_study family, so the propensity-IF builder can be exercised directly
# (the mutation-checkable anchor below). Mirrors R/simultaneous_cis.R:421-707.
.es_dispatch_state <- function(fit) {
	N <- fit$N
	T_ <- fit$T
	G <- fit$G
	p <- fit$p
	treat_inds <- fit$treat_inds
	num_treats <- length(treat_inds)
	cohort_probs_overall <- fit$cohort_probs_overall
	theta_hat_full <- fit$internal$theta_hat
	offs <- fetwfe:::.resolve_event_study_offsets_and_first_inds(
		fit,
		G = G,
		T = T_
	)
	cohort_offsets_int <- offs$cohort_offsets_int
	first_inds <- offs$first_inds
	theta_hat_slopes <- theta_hat_full[2:(p + 1)]
	sel_treat_inds_shifted <- which(theta_hat_slopes[treat_inds] != 0)
	theta_sel <- theta_hat_slopes[treat_inds][sel_treat_inds_shifted]
	d_inv_treat <- fetwfe:::.gen_inv_treat_block(
		num_treats = num_treats,
		first_inds = first_inds,
		G = G,
		fusion_structure = fit$fusion_structure,
		d_inv_treat = fit$internal$d_inv_treat
	)
	d_inv_treat_sel <- d_inv_treat[, sel_treat_inds_shifted, drop = FALSE]
	psi_tes_mat <- fetwfe:::.build_psi_tes_for_family(
		family = "event_study",
		contrasts = NULL,
		G = G,
		T = T_,
		num_treats = num_treats,
		cohort_offsets_int = cohort_offsets_int,
		first_inds = first_inds,
		cohort_probs_overall = cohort_probs_overall
	)
	K <- nrow(psi_tes_mat)
	J_list <- fetwfe:::.build_j_list_for_family(
		family = "event_study",
		K = K,
		G = G,
		T = T_,
		num_treats = num_treats,
		cohort_offsets_int = cohort_offsets_int,
		first_inds = first_inds,
		cohort_probs_overall = cohort_probs_overall,
		d_inv_treat_sel = d_inv_treat_sel
	)
	list(
		J_list = J_list,
		theta_sel = theta_sel,
		cohort_probs_overall = cohort_probs_overall,
		G = G,
		N = N,
		T = T_,
		K = K
	)
}

# ------------------------------------------------------------------------------
# Correctness anchor (mutation-checkable, independent of the analytic builder).
# colSums(.build_propensity_if(...)^2)/(NT)^2 == diag(.assemble_joint_cov_var2(...)).
# LHS is the new centered-one-hot assembly; RHS is the analytic delta-method
# sandwich T*theta'J_k'Sigma_pi J_k theta/(NT). Neither is defined in terms of the
# other -- not a round-trip. Fails on pre-change source (no .build_propensity_if).
# ------------------------------------------------------------------------------
test_that(".build_propensity_if reproduces diag(Sigma_2) (the propensity anchor)", {
	fit <- .es_boot_fit(se_type = "cluster")
	st <- .es_dispatch_state(fit)
	F_pi <- fetwfe:::.build_propensity_if(
		J_list = st$J_list,
		theta_sel = st$theta_sel,
		cohort_probs_overall = st$cohort_probs_overall,
		G = st$G,
		N = st$N,
		T = st$T
	)
	expect_identical(dim(F_pi), c(st$N, st$K))
	Sigma_pi <- fetwfe:::.multinomial_cov(
		st$cohort_probs_overall[seq_len(st$G)]
	)
	Sigma_2 <- fetwfe:::.assemble_joint_cov_var2(
		J_list = st$J_list,
		theta_sel = st$theta_sel,
		Sigma_pi_hat = Sigma_pi,
		N = st$N,
		T = st$T
	)
	lhs <- colSums(F_pi^2) / (st$N * st$T)^2
	expect_equal(lhs, diag(Sigma_2), tolerance = 1e-9)
	# the channel is genuinely non-zero (the family DOES exercise Sigma_2)
	expect_gt(max(diag(Sigma_2)), 0)
})

# ------------------------------------------------------------------------------
# End-to-end SE match: combines the regression anchor (Phase 1) + the propensity
# anchor. The SEs do not depend on the random multipliers, so this is exact. The
# asymmetric cadjust (N/(N-1) on the regression channel only, none on the
# multinomial Sigma_2) is the unique scaling that matches the analytic SE.
# ------------------------------------------------------------------------------
test_that("bootstrap event_study SEs match the analytic cluster SEs to machine precision", {
	fit <- .es_boot_fit(se_type = "cluster")
	an <- simultaneousCIs(fit, family = "event_study", method = "analytic")
	bo <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_equal(.es_se_from(bo), .es_se_from(an), tolerance = 1e-9)
})

test_that("bootstrap event_study critical value matches the analytic qmvnorm one", {
	fit <- .es_boot_fit(se_type = "cluster")
	an <- simultaneousCIs(fit, family = "event_study", method = "analytic")
	bo <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 5000,
		seed = 1
	)
	expect_lt(abs(bo$critical_value - an$critical_value), 0.12)
})

test_that("bootstrap event_study critical value sits in [pointwise, Bonferroni]", {
	fit <- .es_boot_fit(se_type = "cluster")
	bo <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 4000,
		seed = 7
	)
	expect_s3_class(bo, "simultaneous_cis")
	expect_identical(bo$family, "event_study")
	expect_identical(bo$regime, "fixed-p")
	expect_gte(bo$critical_value, bo$pointwise_critical_value - 1e-8)
	expect_lte(bo$critical_value, bo$bonferroni_critical_value + 1e-8)
	# strictly tighter than Bonferroni when K > 1 and effects are correlated
	expect_lt(bo$critical_value, bo$bonferroni_critical_value)
})

# ------------------------------------------------------------------------------
# Determinism / RNG transparency, and independence-yet-reproducibility of the
# two multiplier streams through the single `seed`.
# ------------------------------------------------------------------------------
test_that("event_study bootstrap is deterministic given a seed and varies without one", {
	fit <- .es_boot_fit(se_type = "cluster")
	a <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 1000,
		seed = 42
	)
	b <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 1000,
		seed = 42
	)
	expect_identical(a$critical_value, b$critical_value)
	expect_identical(a$ci, b$ci)
	c2 <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 1000,
		seed = 99
	)
	expect_false(isTRUE(all.equal(a$critical_value, c2$critical_value)))
})

test_that("a seeded event_study bootstrap does not perturb the ambient RNG stream", {
	fit <- .es_boot_fit(se_type = "cluster")
	set.seed(123)
	before <- stats::runif(1)
	set.seed(123)
	invisible(simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 7
	))
	after <- stats::runif(1)
	expect_identical(before, after)
})

# ------------------------------------------------------------------------------
# New-option interaction matrix: family = "event_study" x method = "bootstrap"
# crossed with the bootstrap options and the three estimators.
# ------------------------------------------------------------------------------
test_that("event_study bootstrap works across multiplier / add_ridge / se_type / indep_counts", {
	# mammen multiplier
	fit <- .es_boot_fit(se_type = "cluster")
	bm <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 4000,
		seed = 1,
		multiplier = "mammen"
	)
	expect_gte(bm$critical_value, bm$pointwise_critical_value - 1e-8)
	expect_lte(bm$critical_value, bm$bonferroni_critical_value + 1e-8)

	# add_ridge = TRUE
	fit_r <- .es_boot_fit(se_type = "cluster", add_ridge = TRUE)
	br <- simultaneousCIs(
		fit_r,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(br, "simultaneous_cis")
	expect_true(all(is.finite(br$ci$simultaneous_ci_low)))

	# se_type = "default" (the bootstrap is se_type-agnostic: it always uses the
	# empirical per-unit IF, so it returns a band regardless of the fit's se_type)
	fit_d <- .es_boot_fit(se_type = "default")
	bd <- simultaneousCIs(
		fit_d,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(bd, "simultaneous_cis")
	expect_true(all(is.finite(bd$ci$simultaneous_ci_low)))

	# indep_counts supplied (indep_counts_used = TRUE path)
	fit_ic <- .es_boot_fit(se_type = "cluster", keep_indep_counts = TRUE)
	expect_true(isTRUE(fit_ic$indep_counts_used))
	bic <- simultaneousCIs(
		fit_ic,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(bic, "simultaneous_cis")
	expect_true(all(is.finite(bic$ci$simultaneous_ci_low)))
})

test_that("q >= 1 errors with the has_valid_ses stop before the event_study dispatch", {
	# q >= 1 => calc_ses = FALSE => no valid SEs => simultaneousCIs() stops before
	# reaching the bootstrap dispatch (so event_study inherits the stop, not a band).
	fit_q1 <- .es_boot_fit(se_type = "cluster", q = 1)
	expect_error(
		simultaneousCIs(
			fit_q1,
			family = "event_study",
			method = "bootstrap",
			B = 100,
			seed = 1
		),
		"calc_ses|standard errors"
	)
})

test_that("event_study bootstrap works for etwfe and betwfe", {
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
	for (fit in list(
		etwfeWithSimulatedData(dat),
		betwfeWithSimulatedData(dat)
	)) {
		bo <- simultaneousCIs(
			fit,
			family = "event_study",
			method = "bootstrap",
			B = 1000,
			seed = 1
		)
		expect_s3_class(bo, "simultaneous_cis")
		expect_identical(bo$regime, "fixed-p")
		expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
		expect_gte(bo$critical_value, bo$pointwise_critical_value - 1e-8)
		expect_lte(bo$critical_value, bo$bonferroni_critical_value + 1e-8)
	}
})

# ------------------------------------------------------------------------------
# Edge cases.
# ------------------------------------------------------------------------------
test_that("single treated cohort (G = 1): the propensity builder is exercised and the band is graceful", {
	# G = 1 => Sigma_2 collapses (Sigma_pi is the scalar pi(1-pi)). The
	# A = [a_1 ... a_K] matrix must stay G x K (the bare vapply collapses to a
	# vector when G == 1); this directly exercises the builder, not just the
	# analytic path.
	fit <- .es_boot_fit(se_type = "cluster", G = 1, seed = 3)
	expect_equal(fit$G, 1)
	st <- .es_dispatch_state(fit)
	F_pi <- fetwfe:::.build_propensity_if(
		J_list = st$J_list,
		theta_sel = st$theta_sel,
		cohort_probs_overall = st$cohort_probs_overall,
		G = st$G,
		N = st$N,
		T = st$T
	)
	# builder must return the full N x K shape (the matrix(., nrow = G) guard)
	expect_identical(dim(F_pi), c(st$N, st$K))
	Sigma_pi <- fetwfe:::.multinomial_cov(
		st$cohort_probs_overall[seq_len(st$G)]
	)
	Sigma_2 <- fetwfe:::.assemble_joint_cov_var2(
		J_list = st$J_list,
		theta_sel = st$theta_sel,
		Sigma_pi_hat = Sigma_pi,
		N = st$N,
		T = st$T
	)
	expect_equal(
		colSums(F_pi^2) / (st$N * st$T)^2,
		diag(Sigma_2),
		tolerance = 1e-9
	)
	# the full dispatch returns a graceful band too
	bo <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(bo, "simultaneous_cis")
	expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
})

test_that("event_study bootstrap works with a never-treated group and with no covariates (d = 0)", {
	# never-treated group present (the default DGP has one); d = 0 (no covariates)
	fit <- .es_boot_fit(se_type = "cluster", d = 0)
	bo <- simultaneousCIs(
		fit,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(bo, "simultaneous_cis")
	expect_identical(bo$regime, "fixed-p")
	expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
})

# ------------------------------------------------------------------------------
# High-dimensional (full p >= NT): since #299 event_study routes through the
# DESPARSIFIED `targets` path (both the full-design regression channel and the
# propensity channel `F_pi`), not the fixed-p selected-support fallback. No
# analytic ground truth at p >= NT, so assert finiteness + regime + diagnostics.
# ------------------------------------------------------------------------------
test_that("high-dim event_study bootstrap uses the desparsified path and returns a finite band (#299)", {
	# G=3, T=5, d=22 => full p = 390 > NT = 200 at N = 40; selected support ~42 < NT.
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
	dat$indep_counts <- NA
	hd_fit <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	# genuinely full p >= NT, but a low-dim selected support
	expect_gte(ncol(hd_fit$internal$X_final), nrow(hd_fit$internal$X_final))
	expect_lt(
		sum(hd_fit$internal$theta_hat[-1] != 0),
		nrow(hd_fit$internal$X_final)
	)
	bo <- simultaneousCIs(
		hd_fit,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(bo, "simultaneous_cis")
	# the proof it uses the desparsified path: regime is high-dimensional, with
	# the nodewise diagnostics attached.
	expect_identical(bo$regime, "high-dimensional")
	expect_true(all(
		c("feasibility", "converged", "lambda_node") %in% names(bo)
	))
	expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
	expect_gte(bo$critical_value, bo$pointwise_critical_value - 1e-8)
	expect_lte(bo$critical_value, bo$bonferroni_critical_value + 1e-8)
})

# ------------------------------------------------------------------------------
# Calibrated family-wise coverage (skip_on_cran): on a non-degenerate event_study
# DGP (Sigma_2 > 0), in what fraction of simulated panels does the simultaneous
# band cover ALL true event-time effects jointly? Nominal 0.95; floor 0.80. A
# 50%-miscalibrated F_pi (here: a 50%-miscalibrated propensity channel via a
# scaled critical value) must fail this floor -- so the propensity channel is
# tested for calibration, not merely finiteness.
# ------------------------------------------------------------------------------
test_that("event_study bootstrap bands attain near-nominal family-wise coverage", {
	skip_on_cran()
	G <- 3L
	T_ <- 5L
	d_ <- 1L
	N_ <- 200L
	nsim <- 40L
	coefs <- genCoefs(
		G = G,
		T = T_,
		d = d_,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	num_treats <- getNumTreats(G = G, T = T_)
	treat_inds_truth <- getTreatInds(
		G = G,
		T = T_,
		d = d_,
		num_treats = num_treats
	)
	tes_truth <- coefs$beta[treat_inds_truth]
	first_inds_truth <- getFirstInds(G = G, T = T_)
	event_times <- 0:(T_ - 2L)

	covered <- logical(nsim)
	for (s in seq_len(nsim)) {
		dat <- simulateData(
			coefs,
			N = N_,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = 100L + s
		)
		dat$indep_counts <- NA
		fit <- fetwfeWithSimulatedData(dat, q = 0.5, se_type = "cluster")
		bo <- simultaneousCIs(
			fit,
			family = "event_study",
			method = "bootstrap",
			B = 400,
			seed = 1
		)
		# true event-time effects: probability-weighted pooling of the true
		# per-cell effects across the cohorts treated by event time e (matches
		# the truth construction in test-simultaneous-cis.R Test 3).
		cpo <- fit$cohort_probs_overall
		offs <- fetwfe:::.resolve_event_study_offsets_and_first_inds(
			fit,
			G = G,
			T = T_
		)
		coh_off <- offs$cohort_offsets_int
		truth_e <- numeric(length(event_times))
		for (kk in seq_along(event_times)) {
			e <- event_times[kk]
			V_e <- which(coh_off <= T_ - e)
			weights_Ve <- cpo[V_e] / sum(cpo[V_e])
			truth_e[kk] <- sum(
				weights_Ve * tes_truth[first_inds_truth[V_e] + e]
			)
		}
		covered[s] <- all(
			truth_e >= bo$ci$simultaneous_ci_low &
				truth_e <= bo$ci$simultaneous_ci_high
		)
	}
	expect_gt(mean(covered), 0.80)
})
