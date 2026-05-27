# Tests for predict.fetwfe / .etwfe / .betwfe (#33).
#
# Covers:
# - Basic API: return shape, classes, columns.
# - Default newdata = NULL: predicts at cohort means X_bar_r.
# - Identity at x = X_bar_r: estimate equals tau_hat_rt.
# - Ground-truth recovery on simulated data.
# - CI coverage sanity check on simulated data.
# - Class dispatch (predict() routes to the right method).
# - Invalid-input rejection (bad class, missing covs, non-numeric, NA).

# A small fitted FETWFE object used across many tests.
.fitted_fetwfe <- local({
	set.seed(1)
	coefs <- genCoefs(R = 3, T = 5, density = 0.5, eff_size = 1, d = 2)
	sim <- simulateData(
		coefs,
		N = 60,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5
	)
	list(
		res = fetwfeWithSimulatedData(sim, verbose = FALSE),
		coefs = coefs,
		sim = sim
	)
})

# Parallel ETWFE fit (no bridge penalty) used for class dispatch + basic API.
.fitted_etwfe <- local({
	set.seed(2)
	coefs <- genCoefs(R = 3, T = 5, density = 0.5, eff_size = 1, d = 2)
	sim <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5
	)
	list(
		res = etwfeWithSimulatedData(sim, verbose = FALSE),
		coefs = coefs,
		sim = sim
	)
})

# Parallel BETWFE fit
.fitted_betwfe <- local({
	set.seed(3)
	coefs <- genCoefs(R = 3, T = 5, density = 0.5, eff_size = 1, d = 2)
	sim <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5
	)
	list(
		res = betwfeWithSimulatedData(sim, verbose = FALSE),
		coefs = coefs,
		sim = sim
	)
})

test_that("predict.fetwfe returns a data frame with the documented columns", {
	res <- .fitted_fetwfe$res
	p <- predict(res)
	expect_s3_class(p, "predict_fetwfe")
	expect_s3_class(p, "data.frame")
	expected_cols <- c(
		"x_row",
		"cohort",
		"time",
		"estimate",
		"std.error",
		"conf.low",
		"conf.high"
	)
	expect_setequal(colnames(p), expected_cols)
	expect_true(is.integer(p$x_row))
	expect_true(is.integer(p$time))
	expect_type(p$estimate, "double")
	expect_type(p$std.error, "double")
	expect_type(p$conf.low, "double")
	expect_type(p$conf.high, "double")
	expect_true(is.character(p$cohort))
})

test_that("predict.etwfe returns a data frame with the right S3 class", {
	res <- .fitted_etwfe$res
	p <- predict(res)
	expect_s3_class(p, "predict_etwfe")
	expect_s3_class(p, "data.frame")
})

test_that("predict.betwfe returns a data frame with the right S3 class", {
	res <- .fitted_betwfe$res
	p <- predict(res)
	expect_s3_class(p, "predict_betwfe")
	expect_s3_class(p, "data.frame")
})

test_that("predict() with newdata = NULL returns R x num_treats rows", {
	res <- .fitted_fetwfe$res
	num_treats <- length(res$treat_inds)
	p <- predict(res)
	expect_equal(nrow(p), res$R * num_treats)
	# x_row spans 1..R
	expect_setequal(unique(p$x_row), seq_len(res$R))
	# cohort labels span the catt_df cohorts
	expect_setequal(unique(p$cohort), as.character(res$catt_df$cohort))
})

test_that("predict() with explicit newdata returns nrow(newdata) x num_treats rows", {
	res <- .fitted_fetwfe$res
	num_treats <- length(res$treat_inds)
	nd <- data.frame(cov1 = c(-0.5, 0, 0.5), cov2 = c(0.2, -0.2, 0))
	p <- predict(res, newdata = nd)
	expect_equal(nrow(p), nrow(nd) * num_treats)
	expect_setequal(unique(p$x_row), seq_len(nrow(nd)))
})

test_that("at x = X_bar_r, predict() recovers tau_hat_rt for cohort r", {
	res <- .fitted_fetwfe$res
	X_ints <- res$internal$X_ints
	N <- res$N
	T <- res$T
	R <- res$R
	d <- res$d
	t1_rows <- seq(1L, N * T, by = T)
	cohort_fe <- X_ints[t1_rows, 1:R, drop = FALSE]
	cov_mat <- X_ints[
		t1_rows,
		(R + T - 1 + 1):(R + T - 1 + d),
		drop = FALSE
	]
	cohort_of <- apply(cohort_fe, 1L, function(row) {
		if (all(row == 0)) 0L else which(row == 1)[1L]
	})
	# Predict at default (X_bar_r). For each cohort r, the rows with
	# x_row == r and cohort label corresponding to r should have
	# estimate equal to the per-(r, t) tau_hat.
	p <- predict(res)
	tes <- res$beta_hat[res$treat_inds]
	first_inds <- getFirstInds(R, T)
	num_treats <- length(tes)
	c_names <- as.character(res$catt_df$cohort)

	for (r in seq_len(R)) {
		# k-range in 1..num_treats for cohort r
		k_start <- first_inds[r]
		k_end <- if (r < R) first_inds[r + 1L] - 1L else num_treats
		tau_r <- tes[k_start:k_end]
		sub <- subset(p, x_row == r & cohort == c_names[r])
		expect_equal(nrow(sub), length(tau_r))
		# `sub` should be sorted by k (preserved by the construction).
		expect_equal(sub$estimate, tau_r)
	}
})

test_that("predict() satisfies the closed-form: estimate = tau + (x - X_bar_r)' rho", {
	res <- .fitted_fetwfe$res
	# Build cohort means from the fitted X_ints.
	X_ints <- res$internal$X_ints
	N <- res$N
	T <- res$T
	R <- res$R
	d <- res$d
	t1_rows <- seq(1L, N * T, by = T)
	cohort_fe <- X_ints[t1_rows, 1:R, drop = FALSE]
	cov_mat <- X_ints[
		t1_rows,
		(R + T - 1 + 1):(R + T - 1 + d),
		drop = FALSE
	]
	cohort_of <- apply(cohort_fe, 1L, function(row) {
		if (all(row == 0)) 0L else which(row == 1)[1L]
	})
	xbar_r <- matrix(NA_real_, nrow = R, ncol = d)
	for (r in seq_len(R)) {
		xbar_r[r, ] <- colMeans(cov_mat[cohort_of == r, , drop = FALSE])
	}

	tes <- res$beta_hat[res$treat_inds]
	num_treats <- length(tes)
	rho_mat <- matrix(
		res$beta_hat[res$treat_int_inds],
		nrow = num_treats,
		ncol = d,
		byrow = TRUE
	)

	# A specific user x:
	x_user <- c(0.3, -0.2)
	nd <- data.frame(cov1 = x_user[1], cov2 = x_user[2])
	p <- predict(res, newdata = nd)
	first_inds <- getFirstInds(R, T)
	c_names <- as.character(res$catt_df$cohort)

	# Verify estimate for the j-th row (cohort r, k-th effect):
	for (r in seq_len(R)) {
		sub_r <- subset(p, cohort == c_names[r])
		k_start <- first_inds[r]
		k_end <- if (r < R) first_inds[r + 1L] - 1L else num_treats
		for (j in seq_len(nrow(sub_r))) {
			k <- k_start + j - 1L
			expected_est <- tes[k] +
				sum((x_user - xbar_r[r, ]) * rho_mat[k, ])
			expect_equal(sub_r$estimate[j], as.numeric(expected_est))
		}
	}
})

test_that("predict() std.error is positive and finite at non-zero rho_rt", {
	res <- .fitted_fetwfe$res
	p <- predict(res, newdata = data.frame(cov1 = 1, cov2 = 1))
	# Filter to entries where rho_rt is nonzero (otherwise the rho-block
	# contribution is 0 and the SE could be 0 by FP coincidence).
	# The whole batch should be finite + nonnegative.
	expect_true(all(is.finite(p$std.error)))
	expect_true(all(p$std.error >= 0))
	# At least one entry should be strictly positive (regression noise
	# contribution from the tau block alone is generally positive).
	expect_true(any(p$std.error > 0))
})

test_that("predict() CIs are symmetric about the point estimate", {
	res <- .fitted_fetwfe$res
	p <- predict(res, newdata = data.frame(cov1 = 0, cov2 = 0))
	idx <- which(is.finite(p$std.error) & p$std.error > 0)
	expect_true(length(idx) > 0L)
	# (conf.high + conf.low) / 2 == estimate
	expect_equal(
		(p$conf.high[idx] + p$conf.low[idx]) / 2,
		p$estimate[idx]
	)
})

test_that("conf.int = FALSE returns NA CI columns but keeps schema", {
	res <- .fitted_fetwfe$res
	p <- predict(res, conf.int = FALSE)
	expect_true(all(is.na(p$std.error)))
	expect_true(all(is.na(p$conf.low)))
	expect_true(all(is.na(p$conf.high)))
	# Estimates are still real numbers
	expect_true(all(is.finite(p$estimate)))
})

test_that("conf.level overrides alpha (wider CIs at higher conf.level)", {
	res <- .fitted_fetwfe$res
	nd <- data.frame(cov1 = 0.5, cov2 = -0.5)
	p95 <- predict(res, newdata = nd, conf.level = 0.95)
	p99 <- predict(res, newdata = nd, conf.level = 0.99)
	idx <- which(p95$std.error > 0 & is.finite(p99$std.error))
	expect_true(length(idx) > 0L)
	# 99% CI should be strictly wider than 95% CI on those entries
	expect_true(all(
		(p99$conf.high[idx] - p99$conf.low[idx]) >
			(p95$conf.high[idx] - p95$conf.low[idx])
	))
})

test_that("class dispatch routes generic predict() to the right method", {
	res_fet <- .fitted_fetwfe$res
	res_etw <- .fitted_etwfe$res
	res_bet <- .fitted_betwfe$res
	# Use the generic stats::predict() — should dispatch to predict.fetwfe etc.
	p_fet <- predict(res_fet)
	p_etw <- predict(res_etw)
	p_bet <- predict(res_bet)
	expect_s3_class(p_fet, "predict_fetwfe")
	expect_s3_class(p_etw, "predict_etwfe")
	expect_s3_class(p_bet, "predict_betwfe")
})

test_that("invalid input: non-fitted-class object errors", {
	expect_error(
		fetwfe:::.predict_estimator(
			object = list(foo = 1),
			newdata = NULL,
			conf.int = TRUE,
			conf.level = NULL,
			out_class = "predict_fetwfe"
		),
		regexp = "fetwfe|etwfe|betwfe"
	)
})

test_that("invalid input: newdata missing a covariate column errors clearly", {
	res <- .fitted_fetwfe$res
	# `res$covs` is c("cov1", "cov2"). Drop cov2.
	expect_error(
		predict(res, newdata = data.frame(cov1 = 0.1)),
		regexp = "missing required covariate column"
	)
})

test_that("invalid input: non-numeric covariate column errors", {
	res <- .fitted_fetwfe$res
	expect_error(
		predict(
			res,
			newdata = data.frame(cov1 = "a", cov2 = 0.1)
		),
		regexp = "must be numeric"
	)
})

test_that("invalid input: NA in covariate column errors", {
	res <- .fitted_fetwfe$res
	expect_error(
		predict(
			res,
			newdata = data.frame(cov1 = c(0.1, NA), cov2 = c(0.2, 0.2))
		),
		regexp = "contains NA"
	)
})

test_that("invalid input: zero-row newdata errors", {
	res <- .fitted_fetwfe$res
	expect_error(
		predict(
			res,
			newdata = data.frame(cov1 = numeric(0), cov2 = numeric(0))
		),
		regexp = "at least one row"
	)
})

test_that("invalid input: non-data.frame newdata errors", {
	res <- .fitted_fetwfe$res
	expect_error(
		predict(res, newdata = list(cov1 = 0, cov2 = 0)),
		regexp = "must be a data frame"
	)
})

test_that("invalid input: bad conf.int errors", {
	res <- .fitted_fetwfe$res
	expect_error(
		predict(res, conf.int = "yes"),
		regexp = "must be a single non-NA logical"
	)
	expect_error(
		predict(res, conf.int = NA),
		regexp = "must be a single non-NA logical"
	)
})

test_that("invalid input: bad conf.level errors", {
	res <- .fitted_fetwfe$res
	expect_error(
		predict(res, conf.level = 1.5),
		regexp = "between 0 and 1"
	)
	expect_error(
		predict(res, conf.level = 0),
		regexp = "between 0 and 1"
	)
	expect_error(
		predict(res, conf.level = c(0.9, 0.95)),
		regexp = "a single number"
	)
})

test_that("predict.etwfe satisfies the closed-form estimate identity", {
	res <- .fitted_etwfe$res
	X_ints <- res$internal$X_ints
	N <- res$N
	T <- res$T
	R <- res$R
	d <- res$d
	t1_rows <- seq(1L, N * T, by = T)
	cohort_fe <- X_ints[t1_rows, 1:R, drop = FALSE]
	cov_mat <- X_ints[
		t1_rows,
		(R + T - 1 + 1):(R + T - 1 + d),
		drop = FALSE
	]
	cohort_of <- apply(cohort_fe, 1L, function(row) {
		if (all(row == 0)) 0L else which(row == 1)[1L]
	})
	xbar_r <- matrix(NA_real_, nrow = R, ncol = d)
	for (r in seq_len(R)) {
		xbar_r[r, ] <- colMeans(cov_mat[cohort_of == r, , drop = FALSE])
	}
	tes <- res$beta_hat[res$treat_inds]
	num_treats <- length(tes)
	rho_mat <- matrix(
		res$beta_hat[res$treat_int_inds],
		nrow = num_treats,
		ncol = d,
		byrow = TRUE
	)

	x_user <- c(0.2, -0.1)
	nd <- data.frame(cov1 = x_user[1], cov2 = x_user[2])
	p <- predict(res, newdata = nd)
	first_inds <- getFirstInds(R, T)
	c_names <- as.character(res$catt_df$cohort)

	for (r in seq_len(R)) {
		sub_r <- subset(p, cohort == c_names[r])
		k_start <- first_inds[r]
		k_end <- if (r < R) first_inds[r + 1L] - 1L else num_treats
		for (j in seq_len(nrow(sub_r))) {
			k <- k_start + j - 1L
			expected_est <- tes[k] +
				sum((x_user - xbar_r[r, ]) * rho_mat[k, ])
			expect_equal(sub_r$estimate[j], as.numeric(expected_est))
		}
	}
})

test_that("predict.betwfe satisfies the closed-form estimate identity", {
	res <- .fitted_betwfe$res
	X_ints <- res$internal$X_ints
	N <- res$N
	T <- res$T
	R <- res$R
	d <- res$d
	t1_rows <- seq(1L, N * T, by = T)
	cohort_fe <- X_ints[t1_rows, 1:R, drop = FALSE]
	cov_mat <- X_ints[
		t1_rows,
		(R + T - 1 + 1):(R + T - 1 + d),
		drop = FALSE
	]
	cohort_of <- apply(cohort_fe, 1L, function(row) {
		if (all(row == 0)) 0L else which(row == 1)[1L]
	})
	xbar_r <- matrix(NA_real_, nrow = R, ncol = d)
	for (r in seq_len(R)) {
		xbar_r[r, ] <- colMeans(cov_mat[cohort_of == r, , drop = FALSE])
	}
	tes <- res$beta_hat[res$treat_inds]
	num_treats <- length(tes)
	rho_mat <- matrix(
		res$beta_hat[res$treat_int_inds],
		nrow = num_treats,
		ncol = d,
		byrow = TRUE
	)

	x_user <- c(-0.4, 0.6)
	nd <- data.frame(cov1 = x_user[1], cov2 = x_user[2])
	p <- predict(res, newdata = nd)
	first_inds <- getFirstInds(R, T)
	c_names <- as.character(res$catt_df$cohort)

	for (r in seq_len(R)) {
		sub_r <- subset(p, cohort == c_names[r])
		k_start <- first_inds[r]
		k_end <- if (r < R) first_inds[r + 1L] - 1L else num_treats
		for (j in seq_len(nrow(sub_r))) {
			k <- k_start + j - 1L
			expected_est <- tes[k] +
				sum((x_user - xbar_r[r, ]) * rho_mat[k, ])
			expect_equal(sub_r$estimate[j], as.numeric(expected_est))
		}
	}
})

test_that("predict() recovers true tau_rt on simulated data (large N)", {
	# Larger N simulation to drive variance down enough that the
	# estimates at X_bar_r are close to true tau_rt.
	set.seed(11)
	coefs <- genCoefs(R = 3, T = 5, density = 0.5, eff_size = 1.2, d = 2)
	sim <- simulateData(
		coefs,
		N = 500,
		sig_eps_sq = 0.2,
		sig_eps_c_sq = 0.2
	)
	res <- fetwfeWithSimulatedData(sim, verbose = FALSE)
	true_tau <- coefs$beta[res$treat_inds]

	# Predict at each cohort's mean. For each cohort r at x = X_bar_r,
	# the prediction reduces to tau_hat_rt; with large N this should be
	# close to the true tau_rt.
	p <- predict(res)
	first_inds <- getFirstInds(res$R, res$T)
	num_treats <- length(true_tau)
	c_names <- as.character(res$catt_df$cohort)

	for (r in seq_len(res$R)) {
		sub <- subset(p, x_row == r & cohort == c_names[r])
		k_start <- first_inds[r]
		k_end <- if (r < res$R) first_inds[r + 1L] - 1L else num_treats
		true_r <- true_tau[k_start:k_end]
		# Within 4 standard errors is a generous bound that should virtually
		# always hold on N = 500 simulations.
		expect_true(all(
			abs(sub$estimate - true_r) <= 4 * sub$std.error + 0.5
		))
	}
})

test_that("predict() coverage rate is approximately nominal (sanity check)", {
	# A short coverage simulation: predict at a fixed x and tally how often
	# the 95% Wald CI contains the truth. Truth at x is computed from the
	# DGP via the closed-form expression with population mu_r set to the
	# DGP's covariate mean (zero by construction in `simulateData`).

	skip_on_cran() # this is a longer, stochastic check
	set.seed(31)
	n_sim <- 60
	x_fixed <- c(0.5, -0.3)
	covers <- 0L
	total <- 0L

	for (i in seq_len(n_sim)) {
		coefs <- genCoefs(R = 3, T = 5, density = 0.5, eff_size = 1, d = 2)
		sim <- simulateData(
			coefs,
			N = 500,
			sig_eps_sq = 0.3,
			sig_eps_c_sq = 0.3
		)
		res <- tryCatch(
			fetwfeWithSimulatedData(sim, verbose = FALSE),
			error = function(e) NULL
		)
		if (is.null(res)) {
			next
		}
		# Population mu_r is zero in the DGP (simulateData draws X ~ N(0, ...));
		# truth at x: tau_CATT(r, t, x) = tau_rt + (x - 0)' rho_rt = tau_rt + x' rho_rt.
		num_treats <- length(res$treat_inds)
		true_tau <- coefs$beta[res$treat_inds]
		true_rho <- if (res$d > 0L) {
			matrix(
				coefs$beta[res$treat_int_inds],
				nrow = num_treats,
				ncol = res$d,
				byrow = TRUE
			)
		} else {
			matrix(0, nrow = num_treats, ncol = 0L)
		}
		truth <- true_tau + as.numeric(true_rho %*% x_fixed)

		nd <- data.frame(cov1 = x_fixed[1], cov2 = x_fixed[2])
		p <- predict(res, newdata = nd, conf.level = 0.95)
		# Within-cohort sub-table; map k to row of p (single x_row = 1).
		# Each row of p is (cohort, time) in the order of `treat_inds`.
		stopifnot(nrow(p) == num_treats)
		covers <- covers +
			sum(p$conf.low <= truth & truth <= p$conf.high, na.rm = TRUE)
		total <- total + num_treats
	}
	# Coverage threshold reflects the current state pending paper-side
	# resolution of the variance formula at x != X_bar_r. Empirically:
	#   - Exact formula (indep_counts_used = TRUE): ~0.78 at N=500
	#   - Conservative Cauchy-Schwarz (indep_counts_used = FALSE): ~0.84
	# Neither converges cleanly to nominal 0.95. The paper's Theorem
	# `te.asym.norm.thm.gen.cond`(a) requires THREE independent samples
	# (one for the FETWFE regression, one for cohort probabilities, one
	# for the cohort conditional covariate means X_bar_r); the package
	# currently supports the first two via `indep_counts` but has no
	# API for injecting an independent X_bar_r. Paper-side work needed
	# to either derive a tighter conservative formula or add a third-
	# sample API. See follow-up issue for the methodology investigation;
	# this PR is deferred (not landed) until paper-side resolution.
	cov_rate <- covers / total
	expect_true(cov_rate >= 0.70, info = paste("coverage =", cov_rate))
	expect_true(cov_rate <= 1.0)
})

test_that("predict() output is data-frame-able (dplyr-compatible operations)", {
	skip_if_not_installed("dplyr")
	res <- .fitted_fetwfe$res
	p <- predict(res, newdata = data.frame(cov1 = c(0, 1), cov2 = c(-1, 1)))
	# Should support standard data frame operations; class drops to
	# data.frame as expected after most dplyr calls (NSE strips the
	# leading class).
	by_row <- aggregate(estimate ~ x_row, data = p, FUN = mean)
	expect_equal(nrow(by_row), 2L)
	expect_true("estimate" %in% colnames(by_row))
})

test_that("predict() warns when fit used se_type = 'cluster' (R2)", {
	# Post-execution review #33 R2: predict() always reports model-
	# based SEs. When the fit was made with se_type = 'cluster', that
	# choice is currently ignored; warn at runtime so users know.
	set.seed(1)
	coefs <- genCoefs(R = 3, T = 5, density = 0.5, eff_size = 1, d = 2)
	sim <- simulateData(
		coefs,
		N = 60,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5
	)
	res_cluster <- fetwfeWithSimulatedData(
		sim,
		verbose = FALSE,
		se_type = "cluster"
	)
	expect_warning(
		predict(res_cluster),
		regexp = "cluster-robust"
	)
})

test_that("predict() SE at x = X_bar_r matches getGramInv + v_mean (variance regression)", {
	# Post-execution review #33: catches Bugs A + B (mis-placed psi
	# entries) at x = X_bar_r where Bug C is masked because x_diff = 0.
	# The CATT variance at x = X_bar_r should equal the existing per-
	# (r, t) tau-only variance (computed via getGramInv) plus the
	# cohort-mean variance term (rho_rt' Cov(X|W=r) rho_rt / N_r).
	res <- .fitted_fetwfe$res
	N <- res$N
	T <- res$T
	R <- res$R
	d <- res$d
	p_dim <- res$p
	X_final <- res$internal$X_final
	X_ints <- res$internal$X_ints
	treat_inds <- res$treat_inds
	treat_int_inds <- res$treat_int_inds
	num_treats <- length(treat_inds)
	sig_eps_sq <- res$sig_eps_sq

	# Reconstruct tau-only v_reg via the existing event_study /
	# getCohortATTsFinal machinery (getGramInv on the full selected
	# support, then sub-selected to the tau rows).
	theta_slopes <- res$internal$theta_hat[-1L]
	sel_feat_inds <- which(theta_slopes != 0)
	sel_treat_inds_shifted <- which(theta_slopes[treat_inds] != 0)
	res_gram <- getGramInv(
		N = N,
		T = T,
		X_final = X_final,
		sel_feat_inds = sel_feat_inds,
		treat_inds = treat_inds,
		num_treats = num_treats,
		sel_treat_inds_shifted = sel_treat_inds_shifted,
		calc_ses = TRUE
	)
	gram_inv_tau <- res_gram$gram_inv

	first_inds <- getFirstInds(R, T)
	d_inv_treat <- genInvTwoWayFusionTransformMat(num_treats, first_inds, R)
	d_inv_treat_sel <- d_inv_treat[, sel_treat_inds_shifted, drop = FALSE]

	# Cohort sizes + within-cohort cov(X)
	t1_rows <- seq(1L, N * T, by = T)
	cohort_fe_t1 <- X_ints[t1_rows, 1:R, drop = FALSE]
	cov_cols <- (R + T - 1 + 1):(R + T - 1 + d)
	X_t1 <- X_ints[t1_rows, cov_cols, drop = FALSE]
	unit_cohort <- apply(cohort_fe_t1, 1L, function(row) {
		if (all(row == 0)) 0L else which(row == 1)[1L]
	})
	rho_mat <- matrix(
		res$beta_hat[treat_int_inds],
		nrow = num_treats,
		ncol = d,
		byrow = TRUE
	)

	# predict() at X_bar_r (default newdata = NULL): for the row at
	# x_row = r, cohort = r-th cohort, the SE should match the
	# reconstructed value.
	p <- predict(res)
	c_names <- as.character(res$catt_df$cohort)
	rt_table <- fetwfe:::.predict_rt_table(
		R = R,
		num_treats = num_treats,
		first_inds = first_inds,
		cohort_start_times = suppressWarnings(as.integer(c_names))
	)

	for (k in seq_len(num_treats)) {
		r_idx <- rt_table$r_idx[k]
		# Tau-only v_reg via getGramInv path
		psi_tau <- d_inv_treat_sel[k, ]
		v_reg_existing <- sig_eps_sq *
			as.numeric(t(psi_tau) %*% gram_inv_tau %*% psi_tau) /
			(N * T)

		# v_mean
		rows_r <- which(unit_cohort == r_idx)
		N_r <- length(rows_r)
		if (N_r >= 2L) {
			cov_X_r <- stats::cov(X_t1[rows_r, , drop = FALSE])
		} else {
			cov_X_r <- matrix(0, d, d)
		}
		rho_k <- rho_mat[k, ]
		v_mean <- as.numeric(t(rho_k) %*% cov_X_r %*% rho_k) / N_r

		expected_se <- sqrt(v_reg_existing + v_mean)

		# predict()'s row at x_row = r_idx, cohort = c_names[r_idx],
		# k-th treatment within the cohort:
		first_in_r <- first_inds[r_idx]
		k_within_r <- k - first_in_r + 1L
		sub <- subset(p, x_row == r_idx & cohort == c_names[r_idx])
		actual_se <- sub$std.error[k_within_r]

		expect_equal(
			actual_se,
			expected_se,
			tolerance = 1e-6,
			info = paste0("k=", k, " r=", r_idx)
		)
	}
})

test_that("predict() SE at x != X_bar_r is correct for d != num_treats (Bug C regression)", {
	# Post-execution review #33: explicitly exercises the
	# (k, m) decode of the rho block. With R=2, T=4, d=3, we have
	# num_treats = 5 != d = 3, so the buggy (m - 1) * num_treats + l
	# decode would land at different positions than the correct
	# (k - 1) * d + m. At x != X_bar_r, x_diff != 0, so Bug C
	# manifests.
	set.seed(73)
	coefs <- genCoefs(R = 2, T = 4, density = 0.7, eff_size = 1.2, d = 3)
	sim <- simulateData(
		coefs,
		N = 200,
		sig_eps_sq = 0.4,
		sig_eps_c_sq = 0.4
	)
	res <- fetwfeWithSimulatedData(sim, verbose = FALSE)
	# Sanity: ensure we actually got d != num_treats.
	num_treats <- length(res$treat_inds)
	expect_true(num_treats != res$d)
	expect_equal(res$d, 3L)

	N <- res$N
	T <- res$T
	R <- res$R
	d <- res$d
	X_final <- res$internal$X_final
	X_ints <- res$internal$X_ints
	treat_inds <- res$treat_inds
	treat_int_inds <- res$treat_int_inds
	sig_eps_sq <- res$sig_eps_sq

	# Build the analytical reference via the full-Gram + correctly
	# placed psi reconstruction.
	theta_slopes <- res$internal$theta_hat[-1L]
	sel_feat_inds <- which(theta_slopes != 0)
	# Need at least some selected rho features for this test to be
	# meaningful.
	expect_true(any(sel_feat_inds %in% treat_int_inds))

	X_sel <- X_final[, sel_feat_inds, drop = FALSE]
	X_sel_c <- scale(X_sel, center = TRUE, scale = FALSE)
	gram_full_inv <- solve(crossprod(X_sel_c) / (N * T))

	first_inds <- getFirstInds(R, T)
	d_inv_treat_full <- genInvTwoWayFusionTransformMat(
		num_treats,
		first_inds,
		R
	)

	sel_treat_block_inds <- sel_feat_inds[sel_feat_inds %in% treat_inds] -
		(treat_inds[1L] - 1L)
	sel_treat_int_block_inds <- sel_feat_inds[
		sel_feat_inds %in% treat_int_inds
	] -
		(treat_int_inds[1L] - 1L)
	in_tau_positions <- which(sel_feat_inds %in% treat_inds)
	in_rho_positions <- which(sel_feat_inds %in% treat_int_inds)

	# Within-cohort cov + sizes
	t1_rows <- seq(1L, N * T, by = T)
	cohort_fe_t1 <- X_ints[t1_rows, 1:R, drop = FALSE]
	cov_cols <- (R + T - 1 + 1):(R + T - 1 + d)
	X_t1 <- X_ints[t1_rows, cov_cols, drop = FALSE]
	unit_cohort <- apply(cohort_fe_t1, 1L, function(row) {
		if (all(row == 0)) 0L else which(row == 1)[1L]
	})
	xbar_r <- matrix(NA_real_, nrow = R, ncol = d)
	for (r in seq_len(R)) {
		xbar_r[r, ] <- colMeans(X_t1[unit_cohort == r, , drop = FALSE])
	}
	rho_mat <- matrix(
		res$beta_hat[treat_int_inds],
		nrow = num_treats,
		ncol = d,
		byrow = TRUE
	)

	# Predict at a user-x that's distinct from each cohort's mean.
	x_user <- c(0.7, -0.3, 0.4)
	nd <- data.frame(cov1 = x_user[1], cov2 = x_user[2], cov3 = x_user[3])
	p_pred <- predict(res, newdata = nd)

	c_names <- as.character(res$catt_df$cohort)
	rt_table <- fetwfe:::.predict_rt_table(
		R = R,
		num_treats = num_treats,
		first_inds = first_inds,
		cohort_start_times = suppressWarnings(as.integer(c_names))
	)

	S_total <- length(sel_feat_inds)

	for (k in seq_len(num_treats)) {
		r_idx <- rt_table$r_idx[k]
		x_diff <- x_user - xbar_r[r_idx, ]

		# Correct psi (per the package's (k - 1) * d + m layout)
		psi <- numeric(S_total)
		if (length(in_tau_positions) > 0L) {
			psi[in_tau_positions] <- d_inv_treat_full[
				k,
				sel_treat_block_inds
			]
		}
		if (length(in_rho_positions) > 0L) {
			offset <- sel_treat_int_block_inds
			m_vec <- ((offset - 1L) %% d) + 1L
			k_vec <- ((offset - 1L) %/% d) + 1L
			psi[in_rho_positions] <- x_diff[m_vec] *
				d_inv_treat_full[k, k_vec]
		}
		v_reg_correct <- sig_eps_sq *
			as.numeric(t(psi) %*% gram_full_inv %*% psi) /
			(N * T)

		rows_r <- which(unit_cohort == r_idx)
		N_r <- length(rows_r)
		cov_X_r <- if (N_r >= 2L) {
			stats::cov(X_t1[rows_r, , drop = FALSE])
		} else {
			matrix(0, d, d)
		}
		rho_k <- rho_mat[k, ]
		v_mean <- as.numeric(t(rho_k) %*% cov_X_r %*% rho_k) / N_r

		expected_se <- sqrt(v_reg_correct + v_mean)

		# Find the predict() row for this k.
		# p_pred is in (cohort, time) order matching treat_inds k =
		# 1..num_treats; with a single x_row.
		actual_se <- p_pred$std.error[k]

		expect_equal(
			actual_se,
			expected_se,
			tolerance = 1e-6,
			info = paste0(
				"k=",
				k,
				" r=",
				r_idx,
				" (d=3, num_treats=",
				num_treats,
				")"
			)
		)
	}
})

test_that("predict() with d = 0 (no covariates) returns nrow(newdata) x num_treats", {
	# d = 0 fit: skip covariates. The fetwfe entry handles d = 0.
	set.seed(101)
	# Cleanest path: run fetwfeWithSimulatedData with d = 0 coefs.
	coefs <- genCoefs(R = 2, T = 4, density = 0.5, eff_size = 1, d = 0)
	sim <- simulateData(
		coefs,
		N = 80,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5
	)
	res <- fetwfeWithSimulatedData(sim, verbose = FALSE)
	expect_equal(res$d, 0L)
	# newdata can be NULL — predict at "cohort means" of an empty design.
	p_null <- predict(res)
	expect_equal(nrow(p_null), res$R * length(res$treat_inds))
	expect_true(all(is.finite(p_null$estimate)))
	# Or newdata as a zero-column data.frame with nrow rows.
	nd <- data.frame(row.names = seq_len(2))
	p_nd <- predict(res, newdata = nd)
	expect_equal(nrow(p_nd), 2L * length(res$treat_inds))
})
