library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Helpers for event_study tests.
# ------------------------------------------------------------------------------
make_es_panel <- function(
	seed = 7,
	R = 2,
	T = 4,
	d = 0,
	N = 500,
	eff_size = 2
) {
	set.seed(seed)
	coefs <- genCoefs(
		R = R,
		T = T,
		d = d,
		density = 0.5,
		eff_size = eff_size,
		seed = seed
	)
	simulateData(coefs, N = N, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
}

# ------------------------------------------------------------------------------
# Test 1: structure of the returned data frame.
# ------------------------------------------------------------------------------
test_that("event_study returns expected structure", {
	sim <- make_es_panel(R = 3, T = 6, d = 0, N = 200)
	res <- etwfeWithSimulatedData(sim)
	es <- event_study(res)

	expect_s3_class(es, "fetwfe_event_study")
	expect_s3_class(es, "data.frame")
	expect_equal(
		colnames(es),
		c(
			"event_time",
			"n_cohorts",
			"estimate",
			"se",
			"ci_low",
			"ci_high",
			"p_value"
		)
	)
	# T = 6, so max_event = 4, length = 5
	expect_equal(nrow(es), 5L)
	expect_equal(es$event_time, 0:4)
	# Estimate and SE columns are numeric, finite, non-negative SE
	expect_true(all(is.finite(es$estimate)))
	expect_true(all(is.finite(es$se)))
	expect_true(all(es$se >= 0))
})

# ------------------------------------------------------------------------------
# Test 2: ETWFE point estimate math against a hand-computed reference.
# R = 2, T = 4, d = 0; cohorts adopt at calendar times 2 and 3.
# Event time 0 cells: cohort 1 at t=2, cohort 2 at t=3.
# Pooled estimate = (n_1 / (n_1 + n_2)) * tes[idx(1, 0)] +
#                   (n_2 / (n_1 + n_2)) * tes[idx(2, 0)]
# where idx(r_idx, e) = first_inds[r_idx] + e.
# ------------------------------------------------------------------------------
test_that("event_study ETWFE point estimates match hand-computed weights x tes", {
	sim <- make_es_panel(R = 2, T = 4, d = 0, N = 500)
	res <- etwfeWithSimulatedData(sim)
	es <- event_study(res)

	tes <- unname(res$beta_hat[res$treat_inds])
	first_inds <- c(1L, 4L) # T*R - R(R+1)/2 = 8 - 3 = 5 cells; cohort 1: 1:3, cohort 2: 4:5
	cohort_probs_overall <- unname(res$cohort_probs_overall)

	# Event time 0: V_0 = {1, 2}; cells (1,0) and (2,0)
	S0 <- sum(cohort_probs_overall[c(1, 2)])
	w1 <- cohort_probs_overall[1] / S0
	w2 <- cohort_probs_overall[2] / S0
	expected_e0 <- w1 * tes[first_inds[1] + 0] + w2 * tes[first_inds[2] + 0]
	expect_equal(es$estimate[1], expected_e0, tolerance = 1e-10)

	# Event time 1: V_1 = {1, 2}; cells (1,1) and (2,1)
	expected_e1 <- w1 * tes[first_inds[1] + 1] + w2 * tes[first_inds[2] + 1]
	expect_equal(es$estimate[2], expected_e1, tolerance = 1e-10)

	# Event time 2: V_2 = {1} only (cohort 2 needs r+1+e = 3+2 = 5 > T = 4)
	expected_e2 <- tes[first_inds[1] + 2]
	expect_equal(es$estimate[3], expected_e2, tolerance = 1e-10)
	expect_equal(es$n_cohorts, c(2L, 2L, 1L))
})

# ------------------------------------------------------------------------------
# Test 3: var_2 = 0 when |V_e| = 1 (single-cohort degeneracy).
# At e = T - 2, only cohort 1 is valid. The cohort-probability variance term
# must vanish; the combined SE equals sqrt(var_1) alone.
# ------------------------------------------------------------------------------
test_that("event_study var_2 vanishes at the last event time (single-cohort)", {
	sim <- make_es_panel(R = 2, T = 4, d = 0, N = 500)
	res <- etwfeWithSimulatedData(sim)
	es <- event_study(res)
	# At e = 2 (last event time), only cohort 1 contributes.
	# var_1(2) = sig_eps_sq * gram_inv[idx(1,2), idx(1,2)] / (N*T)
	# We don't recompute it independently here; instead, check that the SE
	# is exactly sqrt(var_1) by re-deriving var_1 and seeing it matches se^2.
	N <- res$N
	T <- res$T
	R <- res$R
	num_treats <- length(res$treat_inds)
	sig_eps_sq <- res$sig_eps_sq
	res_gram <- fetwfe:::getGramInv(
		N = N,
		T = T,
		X_final = res$X_final,
		sel_feat_inds = NA,
		treat_inds = res$treat_inds,
		num_treats = num_treats,
		sel_treat_inds_shifted = seq_len(num_treats),
		calc_ses = TRUE
	)
	gram_inv <- res_gram$gram_inv
	first_inds <- c(1L, 4L)
	idx_1_2 <- first_inds[1] + 2 # = 3
	var_1_e2 <- sig_eps_sq * gram_inv[idx_1_2, idx_1_2] / (N * T)
	expect_equal(es$se[3]^2, var_1_e2, tolerance = 1e-10)
})

# ------------------------------------------------------------------------------
# Test 4: BETWFE selected-out cells contribute 0 to the estimate.
# When the bridge zeros tes[idx] in beta-space, the weighted sum at each event
# time is zero when all valid cells at that event time are selected out.
# SE handling is implementation-dependent (0 when the support is non-empty,
# NA when the bridge dropped everything), so we focus on the estimate, which
# is the well-defined invariant.
# ------------------------------------------------------------------------------
test_that("event_study BETWFE: estimate is 0 when all cells at e are selected out", {
	sim <- make_es_panel(
		R = 3,
		T = 6,
		d = 0,
		N = 200,
		eff_size = 0.05
	)
	res <- betwfeWithSimulatedData(sim, q = 0.5)
	es <- event_study(res)

	tes <- res$beta_hat[res$treat_inds]
	first_inds <- getFirstInds(R = res$R, T = res$T)
	for (k in seq_len(nrow(es))) {
		e <- es$event_time[k]
		V_e <- which(seq_len(res$R) <= res$T - 1L - e)
		cells <- first_inds[V_e] + e
		if (all(tes[cells] == 0)) {
			expect_equal(es$estimate[k], 0, tolerance = 1e-12)
			# SE is 0 (when some support exists) or NA (when bridge dropped
			# everything). Either is consistent with "no information from this
			# event time"; only assert it's not a spurious non-zero number.
			expect_true(es$se[k] == 0 || is.na(es$se[k]))
		}
	}
})

# ------------------------------------------------------------------------------
# Test 5: FETWFE all-zero-theta-treatment-block early-exit.
# When the bridge zeros out the entire theta-treatment block, att_hat = 0 and
# all event-time estimates should be 0.
# ------------------------------------------------------------------------------
test_that("event_study FETWFE: all-zero-theta-block produces all-zero estimates", {
	# Use very small eff_size to push the bridge toward all-zero
	sim <- make_es_panel(
		R = 2,
		T = 4,
		d = 0,
		N = 200,
		eff_size = 0.01
	)
	res <- fetwfeWithSimulatedData(sim, q = 0.5)
	# Construct the all-zero case deterministically by post-patching BOTH
	# theta_hat (the bridge-space coefficient) AND beta_hat (the
	# back-transformed coefficient that event_study.fetwfe() actually
	# reads at R/event_study.R:274). The fit-time all-zero early-exit at
	# R/fetwfe_core.R:904-939 zeroes both representations together; the
	# test mirrors that. Previously this block conditionally skipped if
	# the bridge happened not to zero the treatment block on the test
	# seed — fragile under a seed bump or genCoefs change. See issue #57.
	res$internal$theta_hat[2:(res$p + 1)][res$treat_inds] <- 0
	res$beta_hat[res$treat_inds] <- 0
	stopifnot(
		all(res$internal$theta_hat[2:(res$p + 1)][res$treat_inds] == 0),
		all(res$beta_hat[res$treat_inds] == 0)
	)
	es <- event_study(res)
	expect_true(all(es$estimate == 0))
	expect_true(all(es$se == 0 | is.na(es$se)))
})

# ------------------------------------------------------------------------------
# Test 6: se_type = "cluster" returns finite SEs and substitutes the sandwich.
# Quick correctness check: under cluster, the SE differs from the default SE
# (typically larger under AR(1)-style noise, but on F1-conforming sim the two
# can be close). We just check finiteness and that the cluster path runs.
# ------------------------------------------------------------------------------
test_that("event_study se_type = 'cluster' returns finite SEs", {
	sim <- make_es_panel(R = 3, T = 6, d = 0, N = 200)
	res_def <- fetwfeWithSimulatedData(sim, q = 0.5, se_type = "default")
	res_cls <- fetwfeWithSimulatedData(sim, q = 0.5, se_type = "cluster")
	es_def <- event_study(res_def)
	es_cls <- event_study(res_cls)
	# Both produce finite SEs at event times where at least one cell is selected
	finite_def <- is.finite(es_def$se) & es_def$se > 0
	finite_cls <- is.finite(es_cls$se) & es_cls$se > 0
	expect_true(any(finite_def))
	expect_true(any(finite_cls))
	# Point estimates are the same (only variance differs by se_type)
	expect_equal(es_def$estimate, es_cls$estimate, tolerance = 1e-10)
})

test_that(".compute_cluster_robust_sandwich matches hand-rolled Liang-Zeger on a small fixture", {
	# Bug-discrimination test for the cluster-robust SE machinery.
	# The integration test above only verifies finiteness; this one
	# verifies the actual numerical formula. A 50% miscalibration in
	# .compute_cluster_robust_sandwich would pass the integration check
	# but fail here. Mutation test: change `cadjust <- N / (N - 1)` in
	# R/ols_calcs.R:623 to `cadjust <- 1` and re-run; this test fails
	# (diff ~6e-4 vs 1e-10 tolerance). See issue #57.
	set.seed(31)
	N <- 25L
	T <- 4L
	p_S <- 3L
	X_S <- matrix(rnorm(N * T * p_S), nrow = N * T, ncol = p_S)
	residuals <- rnorm(N * T)

	V_pkg <- fetwfe:::.compute_cluster_robust_sandwich(
		X_S,
		residuals,
		N = N,
		T = T
	)

	# Hand-rolled Liang-Zeger sandwich (Liang & Zeger 1986):
	#   bread:    (X'X)^{-1} on column-demeaned X
	#   meat:     sum_i (Xi' eps_i)(Xi' eps_i)'
	#   sandwich: cadjust * bread * meat * bread, cadjust = N / (N - 1)
	Xc <- scale(X_S, center = TRUE, scale = FALSE)
	bread <- solve(crossprod(Xc))
	meat <- matrix(0, p_S, p_S)
	for (i in seq_len(N)) {
		rows_i <- ((i - 1L) * T + 1L):(i * T)
		Xi <- Xc[rows_i, , drop = FALSE]
		eps_i <- residuals[rows_i]
		score_i <- crossprod(Xi, eps_i)
		meat <- meat + tcrossprod(score_i)
	}
	V_hand <- (N / (N - 1)) * bread %*% meat %*% bread

	expect_equal(V_pkg, V_hand, tolerance = 1e-10)
})

# ------------------------------------------------------------------------------
# Test 7: plot methods return ggplot objects with expected layer structure.
# ------------------------------------------------------------------------------
test_that("plot.fetwfe / plot.etwfe / plot.betwfe return ggplot objects", {
	skip_if_not_installed("ggplot2")
	sim <- make_es_panel(R = 3, T = 6, d = 0, N = 200)
	res_f <- fetwfeWithSimulatedData(sim, q = 0.5)
	res_e <- etwfeWithSimulatedData(sim)
	res_b <- betwfeWithSimulatedData(sim, q = 0.5)
	p_f <- plot(res_f)
	p_e <- plot(res_e)
	p_b <- plot(res_b)
	expect_s3_class(p_f, "ggplot")
	expect_s3_class(p_e, "ggplot")
	expect_s3_class(p_b, "ggplot")
	# Each plot has the three layers we added (hline, errorbar, point)
	expect_equal(length(p_f$layers), 3L)
	expect_equal(length(p_e$layers), 3L)
	expect_equal(length(p_b$layers), 3L)
})

# ------------------------------------------------------------------------------
# Test 8: dispatch errors on unsupported classes.
# ------------------------------------------------------------------------------
test_that("event_study errors on objects of unsupported class", {
	expect_error(event_study(list(foo = 1)), "fetwfe.*etwfe.*betwfe")
	expect_error(event_study("not a fit"), "fetwfe.*etwfe.*betwfe")
})
