# Direct return-shape tests for the helpers extracted from
# prep_for_etwfe_regression() in #81. Each helper is also exercised
# end-to-end by every fit-side test; these tests lock the return-list
# names/shapes so future return-shape drift in an individual helper is
# caught at the helper boundary, not 20 stack frames later.

# ------------------------------------------------------------------------------
# .estimate_variance_and_gls(): when variance components are supplied,
# the helper should pass them through unchanged and return a transformed
# (y_gls, X_gls) of the right shape.
# ------------------------------------------------------------------------------
test_that(".estimate_variance_and_gls returns expected list shape (#81)", {
	set.seed(1)
	N <- 5
	T <- 4
	p <- 3
	y <- rnorm(N * T)
	X_mod <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)
	X_ints <- X_mod # unused when variance components are supplied

	out <- fetwfe:::.estimate_variance_and_gls(
		y = y,
		X_ints = X_ints,
		X_mod = X_mod,
		sig_eps_sq = 1.0,
		sig_eps_c_sq = 0.5,
		N = N,
		T = T,
		p = p,
		verbose = FALSE
	)

	expect_named(out, c("y_gls", "X_gls", "sig_eps_sq", "sig_eps_c_sq"))
	expect_equal(length(out$y_gls), N * T)
	expect_equal(dim(out$X_gls), c(N * T, p))
	expect_equal(out$sig_eps_sq, 1.0)
	expect_equal(out$sig_eps_c_sq, 0.5)
})

# ------------------------------------------------------------------------------
# .collapse_design_for_twfe_covs(): given a synthetic design matrix sized
# to satisfy the column-slicing arithmetic, the helper should return a
# (N*T) x p_short matrix with p_short = R + T - 1 + d + R.
# ------------------------------------------------------------------------------
test_that(".collapse_design_for_twfe_covs returns expected shape (#81)", {
	R <- 2L
	T <- 4L
	N <- 3L
	d <- 1L
	first_inds <- fetwfe:::getFirstInds(R, T)
	num_treats <- fetwfe:::getNumTreats(R, T)

	n_cols <- R + T - 1 + d * (1 + R + T - 1) + num_treats
	X_gls <- matrix(seq_len(N * T * n_cols), nrow = N * T, ncol = n_cols)

	out <- fetwfe:::.collapse_design_for_twfe_covs(
		X_gls = X_gls,
		N = N,
		T = T,
		G = R,
		d = d,
		num_treats = num_treats,
		first_inds = first_inds
	)

	expect_named(out, c("X_collapsed", "p_short", "treat_inds"))
	expect_equal(out$p_short, R + T - 1 + d + R)
	expect_equal(ncol(out$X_collapsed), out$p_short)
	expect_equal(nrow(out$X_collapsed), N * T)
	# treat_inds is the trailing R (= G) treatment columns -- the single source
	# of truth for the collapsed-design layout consumed by the core (#337).
	expect_equal(out$treat_inds, (R + T - 1 + d + 1):out$p_short)
	expect_equal(length(out$treat_inds), R)
})

# ------------------------------------------------------------------------------
# .append_ridge_rows(): with add_ridge = FALSE, the helper should pass
# inputs through unchanged and return lambda_ridge = NA.
# ------------------------------------------------------------------------------
test_that(".append_ridge_rows is a no-op when add_ridge = FALSE (#81)", {
	N <- 4
	T <- 3
	p <- 5
	X_scaled <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)
	y_gls <- rnorm(N * T)

	out <- fetwfe:::.append_ridge_rows(
		X_scaled = X_scaled,
		y_gls = y_gls,
		p = p,
		add_ridge = FALSE,
		is_fetwfe = FALSE,
		first_inds = 1L, # unused
		T = T,
		G = 1L,
		d = 0L,
		num_treats = 1L,
		sig_eps_sq = 1.0,
		sig_eps_c_sq = 0.5,
		N = N
	)

	expect_named(out, c("X_scaled", "y_final", "lambda_ridge"))
	expect_identical(out$X_scaled, X_scaled)
	expect_identical(out$y_final, y_gls)
	expect_true(is.na(out$lambda_ridge))
})

# ------------------------------------------------------------------------------
# .compute_cohort_probs(): the two normalizations (within-treated and
# overall) should match the documented formulae on a simple count vector.
# ------------------------------------------------------------------------------
test_that(".compute_cohort_probs returns expected shape and values (#81)", {
	# 100 total: 50 never-treated, 30 cohort 1, 20 cohort 2.
	in_sample_counts <- c(50L, 30L, 20L)
	out <- fetwfe:::.compute_cohort_probs(
		in_sample_counts = in_sample_counts,
		indep_counts = NA,
		N = 100L,
		G = 2L,
		indep_count_data_available = FALSE
	)

	expect_named(
		out,
		c(
			"cohort_probs",
			"cohort_probs_overall",
			"indep_cohort_probs",
			"indep_cohort_probs_overall"
		)
	)
	expect_equal(out$cohort_probs, c(30 / 50, 20 / 50))
	expect_equal(out$cohort_probs_overall, c(30 / 100, 20 / 100))
	expect_true(is.na(out$indep_cohort_probs))
	expect_true(is.na(out$indep_cohort_probs_overall))
})
