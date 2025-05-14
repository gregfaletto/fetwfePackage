library(testthat)
library(fetwfe)

# Helper functions to compute expected dimensions:
compute_num_treats <- function(T, R) {
	T * R - (R * (R + 1)) / 2
}

# For full interactions:
compute_p_int <- function(T, R, d) {
	num_treats <- compute_num_treats(T, R)
	R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d
}

# For no interactions:
compute_p_no_int <- function(T, R, d) {
	num_treats <- compute_num_treats(T, R)
	R + (T - 1) + d + num_treats
}

# ------------------------------------------------------------------------------
# Test 1: Output structure and dimensions when gen_ints = TRUE
# ------------------------------------------------------------------------------
test_that("simulateDataCore (with interactions) returns expected output", {
	N <- 120
	T_val <- 30
	R_val <- 5
	d_val <- 12
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	p_int <- compute_p_int(T_val, R_val, d_val)
	beta_int <- rnorm(p_int)

	res <- simulateDataCore(
		N = N,
		T = T_val,
		R = R_val,
		d = d_val,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		beta = beta_int,
		seed = 123,
		gen_ints = TRUE
	)

	# Check list names (X now instead of X_int)
	expected_names <- c(
		"X",
		"y",
		"coefs",
		"covs",
		"first_inds",
		"N_UNTREATED",
		"assignments",
		"indep_counts",
		"p",
		"N",
		"T",
		"R",
		"d",
		"sig_eps_sq",
		"sig_eps_c_sq"
	)
	for (nm in expected_names) {
		expect_true(nm %in% names(res))
	}

	# Check dimensions of X and y
	expect_equal(dim(res$X), c(N * T_val, p_int))
	expect_equal(length(res$y), N * T_val)
})

# ------------------------------------------------------------------------------
# Test 2: Output structure and dimensions when gen_ints = FALSE
# ------------------------------------------------------------------------------
test_that("simulateDataCore (without interactions) returns expected output", {
	N <- 120
	T_val <- 30
	R_val <- 5
	d_val <- 12
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	p_no_int <- compute_p_no_int(T_val, R_val, d_val)

	p_int <- compute_p_int(T_val, R_val, d_val)
	beta_int <- rnorm(p_int)

	res <- simulateDataCore(
		N,
		T_val,
		R_val,
		d_val,
		sig_eps_sq,
		sig_eps_c_sq,
		beta_int,
		seed = 123,
		gen_ints = FALSE
	)

	expect_equal(dim(res$X), c(N * T_val, p_no_int))
	expect_equal(length(res$y), N * T_val)
	expect_equal(res$p, p_no_int)
})

# ------------------------------------------------------------------------------
# Test 3: Reproducibility test (seed)
# ------------------------------------------------------------------------------
test_that("simulateDataCore is reproducible with same seed", {
	N <- 120
	T_val <- 30
	R_val <- 5
	d_val <- 12
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	p_no_int <- compute_p_no_int(T_val, R_val, d_val)

	p_int <- compute_p_int(T_val, R_val, d_val)
	beta_int <- rnorm(p_int)

	res1 <- simulateDataCore(
		N,
		T_val,
		R_val,
		d_val,
		sig_eps_sq,
		sig_eps_c_sq,
		beta_int,
		seed = 456,
		gen_ints = FALSE
	)
	res2 <- simulateDataCore(
		N,
		T_val,
		R_val,
		d_val,
		sig_eps_sq,
		sig_eps_c_sq,
		beta_int,
		seed = 456,
		gen_ints = FALSE
	)

	expect_equal(res1$X, res2$X)
	expect_equal(res1$y, res2$y)
	expect_equal(res1$first_inds, res2$first_inds)
	expect_equal(res1$assignments, res2$assignments)
	expect_equal(res1$indep_counts, res2$indep_counts)
	expect_equal(res1$actual_cohort_tes, res2$actual_cohort_tes)
	expect_equal(res1$att_true, res2$att_true)
})

# ------------------------------------------------------------------------------
# Test 4: Error when beta has incorrect length for gen_ints = TRUE
# ------------------------------------------------------------------------------
test_that("simulateDataCore errors when beta has wrong length (with interactions)", {
	N <- 120
	T_val <- 30
	R_val <- 5
	d_val <- 12
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	p_int <- compute_p_int(T_val, R_val, d_val)
	beta_wrong <- rnorm(p_int - 1)

	expect_error(
		simulateDataCore(
			N,
			T_val,
			R_val,
			d_val,
			sig_eps_sq,
			sig_eps_c_sq,
			beta_wrong,
			seed = 123,
			gen_ints = TRUE
		),
		"length\\(beta\\) must be"
	)
})

# ------------------------------------------------------------------------------
# Test 5: Error when beta has incorrect length for gen_ints = FALSE
# ------------------------------------------------------------------------------
test_that("simulateDataCore errors when beta has wrong length (without interactions)", {
	N <- 120
	T_val <- 30
	R_val <- 5
	d_val <- 12
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	p_no_int <- compute_p_no_int(T_val, R_val, d_val)
	beta_wrong <- rnorm(p_no_int + 2)

	expect_error(
		simulateDataCore(
			N,
			T_val,
			R_val,
			d_val,
			sig_eps_sq,
			sig_eps_c_sq,
			beta_wrong,
			seed = 123,
			gen_ints = FALSE
		),
		"length\\(beta\\) must be"
	)
})

# ------------------------------------------------------------------------------
# Test 6: Check that outputs (X and y) are non-missing.
# ------------------------------------------------------------------------------
test_that("simulateDataCore returns non-missing X and y", {
	N <- 120
	T_val <- 30
	R_val <- 5
	d_val <- 12
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	p_no_int <- compute_p_no_int(T_val, R_val, d_val)

	p_int <- compute_p_int(T_val, R_val, d_val)
	beta_int <- rnorm(p_int)

	res <- simulateDataCore(
		N,
		T_val,
		R_val,
		d_val,
		sig_eps_sq,
		sig_eps_c_sq,
		beta_int,
		seed = 789,
		gen_ints = FALSE
	)
	expect_true(all(!is.na(res$X)))
	expect_true(all(!is.na(res$y)))
})


# ------------------------------------------------------------------------------
# Test 7: When distribution is "uniform", covariates are bounded between -sqrt(3) and sqrt(3)
# ------------------------------------------------------------------------------
test_that("Covariates are uniformly bounded when distribution = 'uniform'", {
	N <- 120
	T_val <- 30
	R_val <- 5
	d_val <- 12
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	# For gen_ints = FALSE, expected p = R + (T-1) + d + num_treats.
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	sim_data <- simulateDataCore(
		N,
		T_val,
		R_val,
		d_val,
		sig_eps_sq,
		sig_eps_c_sq,
		beta,
		seed = 123,
		gen_ints = FALSE,
		distribution = "uniform"
	)

	# In the no-interactions case, the design matrix X is built as:
	# X = [cohort_fe, time_fe, X_long, treat_mat_long]
	# The covariate part is in columns (R + (T-1) + 1) to (R + (T-1) + d)
	base_cols <- R_val + (T_val - 1)
	cov_cols <- seq(from = base_cols + 1, to = base_cols + d_val)
	covariates <- sim_data$X[, cov_cols, drop = FALSE]

	a <- sqrt(3)
	expect_true(all(covariates >= -a - 1e-6))
	expect_true(all(covariates <= a + 1e-6))
})

# ------------------------------------------------------------------------------
# Test 8: Covariates are constant over time for each unit.
# ------------------------------------------------------------------------------
test_that("Covariates are constant over time within each unit", {
	N <- 50
	T_val <- 20
	R_val <- 3
	d_val <- 5
	sig_eps_sq <- 2
	sig_eps_c_sq <- 2
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	sim_data <- simulateDataCore(
		N,
		T_val,
		R_val,
		d_val,
		sig_eps_sq,
		sig_eps_c_sq,
		beta,
		seed = 456,
		gen_ints = FALSE,
		distribution = "gaussian"
	)

	# The design matrix X = [cohort_fe, time_fe, X_long, treat_mat_long].
	# The covariate part is in columns (R + (T-1) + 1):(R + (T-1) + d)
	base_cols <- R_val + (T_val - 1)
	cov_cols <- seq(from = base_cols + 1, to = base_cols + d_val)
	covariates <- sim_data$X[, cov_cols, drop = FALSE]

	# Because covariates were generated as X_long = X (N x d) repeated each T times,
	# for each unit (each block of T rows) the covariate values should be identical.
	X_long_mat <- matrix(NA, nrow = N, ncol = d_val)
	for (i in 1:N) {
		rows <- ((i - 1) * T_val + 1):(i * T_val)
		# For each covariate, check that all values in the block are equal.
		X_long_mat[i, ] <- apply(
			covariates[rows, , drop = FALSE],
			2,
			function(x) {
				x[1]
			}
		)
		for (j in 1:d_val) {
			expect_equal(covariates[rows, j], rep(X_long_mat[i, j], T_val))
		}
	}
})

# ------------------------------------------------------------------------------
# Test 9: The output can be used to create a panel data frame
#         that is accepted by fetwfe().
# ------------------------------------------------------------------------------
test_that("Output from simulateDataCore can be passed to fetwfe()", {
	N <- 30
	T_val <- 10
	R_val <- 2
	d_val <- 2
	sig_eps_sq <- 1
	sig_eps_c_sq <- 1
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	sim_data <- simulateDataCore(
		N,
		T_val,
		R_val,
		d_val,
		sig_eps_sq,
		sig_eps_c_sq,
		beta,
		seed = 789,
		gen_ints = FALSE,
		distribution = "gaussian"
	)

	# Now call fetwfe() with this panel data frame.
	# Note: Since the simulated data come from a simplified process,
	# we can use an empty list for indep_counts and leave noise variance as provided.
	result <- fetwfe(
		pdata = sim_data$pdata,
		time_var = sim_data$time_var,
		unit_var = sim_data$unit_var,
		treatment = sim_data$treatment,
		covs = sim_data$covs,
		response = sim_data$response,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		verbose = FALSE
	)

	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
})

# ------------------------------------------------------------------------------
# Test 10: The output can be used to create a panel data frame
#         that is accepted by fetwfe() even when d = 0.
# ------------------------------------------------------------------------------
test_that("Output from simulateDataCore can be passed to fetwfe()", {
	N <- 30
	T_val <- 10
	R_val <- 2
	d_val <- 0
	sig_eps_sq <- 1
	sig_eps_c_sq <- 1
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	sim_data <- simulateDataCore(
		N,
		T_val,
		R_val,
		d_val,
		sig_eps_sq,
		sig_eps_c_sq,
		beta,
		seed = 789,
		gen_ints = TRUE,
		distribution = "uniform"
	)

	# Now call fetwfe() with this panel data frame.
	# Note: Since the simulated data come from a simplified process,
	# we can use an empty list for indep_counts and leave noise variance as provided.
	result <- fetwfe(
		pdata = sim_data$pdata,
		time_var = sim_data$time_var,
		unit_var = sim_data$unit_var,
		treatment = sim_data$treatment,
		covs = sim_data$covs,
		response = sim_data$response,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		verbose = FALSE
	)

	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
})

# ------------------------------------------------------------------------------
# Test 11: Error when R >= T
# ------------------------------------------------------------------------------
test_that("simulateDataCore errors when R >= T", {
	# For example, set T = 5 and R = 5 (since R should be <= T-1).
	N <- 30
	T_val <- 5
	R_val <- 5
	d_val <- 2
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	expect_error(
		simulateDataCore(
			N,
			T_val,
			R_val,
			d_val,
			1,
			1,
			beta,
			seed = 123,
			gen_ints = FALSE,
			distribution = "gaussian"
		),
		regexp = "R <= T - 1" # Expect an error message indicating R must be <= T-1
	)
})

# ------------------------------------------------------------------------------
# Test 12: Error when T < 3
# ------------------------------------------------------------------------------
test_that("simulateDataCore errors when T < 3", {
	# For instance, T = 2 should trigger an error (we require at least T = 3).
	N <- 30
	T_val <- 2
	R_val <- 1
	d_val <- 2
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	expect_error(
		simulateDataCore(
			N,
			T_val,
			R_val,
			d_val,
			1,
			1,
			beta,
			seed = 123,
			gen_ints = FALSE,
			distribution = "gaussian"
		),
		regexp = "T >= 3" # Expect an error message indicating T must be at least 3
	)
})

# ------------------------------------------------------------------------------
# Test 13: Error when N < R
# ------------------------------------------------------------------------------
test_that("simulateDataCore errors when N < R", {
	# For example, if there are 3 units and 4 treated cohorts, that's impossible.
	N <- 3
	T_val <- 5
	R_val <- 4
	d_val <- 2
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	expect_error(
		simulateDataCore(
			N,
			T_val,
			R_val,
			d_val,
			1,
			1,
			beta,
			seed = 123,
			gen_ints = FALSE,
			distribution = "gaussian"
		),
		regexp = "N >= R" # Expect an error message indicating N must be at least R
	)
})

# ------------------------------------------------------------------------------
# Test 14: Input from genCoefsCore() works.
# ------------------------------------------------------------------------------
test_that("Input from genCoefsCore() works", {
	N <- 30
	T_val <- 3
	R_val <- 2
	d_val <- 2
	sig_eps_sq <- 1
	sig_eps_c_sq <- 1
	density = 0.2
	eff_size <- 1

	res <- genCoefsCore(
		R = R_val,
		T = T_val,
		d = d_val,
		density = density,
		eff_size = eff_size
	)

	p_int <- compute_p_int(T = T_val, R = R_val, d = d_val)

	expect_equal(p_int, length(res$beta))

	sim_data <- simulateDataCore(
		N = N,
		T = T_val,
		R = R_val,
		d = d_val,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		beta = res$beta,
		gen_ints = TRUE,
		seed = 789
	)

	# Check list names (X now instead of X_int)
	expected_names <- c(
		"X",
		"y",
		"coefs",
		"covs",
		"first_inds",
		"N_UNTREATED",
		"assignments",
		"indep_counts",
		"p",
		"N",
		"T",
		"R",
		"d",
		"sig_eps_sq",
		"sig_eps_c_sq"
	)
	for (nm in expected_names) {
		expect_true(nm %in% names(sim_data))
	}

	# Check dimensions of X and y
	expect_equal(nrow(sim_data$X), N * T_val)
	expect_equal(ncol(sim_data$X), p_int)
	expect_equal(length(sim_data$y), N * T_val)
})
