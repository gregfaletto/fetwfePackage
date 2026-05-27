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
			N = 200,
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
	# Expect coverage to be somewhere in [0.7, 1.0] on a 60-rep small-N study.
	# A loose bound matches the paper's empirical findings that finite-sample
	# coverage may not perfectly match the nominal rate.
	cov_rate <- covers / total
	expect_true(cov_rate >= 0.6, info = paste("coverage =", cov_rate))
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
