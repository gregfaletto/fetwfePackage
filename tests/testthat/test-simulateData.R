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
# Test 1: Output structure and dimensions
# ------------------------------------------------------------------------------
test_that("simulateData (with interactions) returns expected output", {
	N <- 30
	T_val <- 6
	R_val <- 5
	d_val <- 4
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	p_int <- compute_p_int(T_val, R_val, d_val)
	beta_int <- rnorm(p_int)

	obj <- list(beta = beta_int, R = R_val, T = T_val, d = d_val, seed = 843)
	class(obj) <- "FETWFE_coefs"

	res <- obj |>
		simulateData(
			N = N,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq
		)

	# Check list names (X now instead of X_int)
	expected_names <- c(
		"pdata",
		"X",
		"y",
		"covs",
		"time_var",
		"unit_var",
		"treatment",
		"response",
		"coefs",
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
# Test 2: Reproducibility test (seed)
# ------------------------------------------------------------------------------
test_that("simulateData is reproducible with same seed", {
	N <- 30
  T_val <- 6
	R_val <- 5
  d_val <- 4
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5

	coefs_obj <- genCoefs(
		R = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 456
	)

	res1 <- coefs_obj |>
		simulateData(
			N = N,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq
		)
	res2 <- coefs_obj |>
		simulateData(
			N = N,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq
		)

	expect_equal(res1$X, res2$X)
	expect_equal(res1$y, res2$y)
	expect_equal(res1$first_inds, res2$first_inds)
	expect_equal(res1$assignments, res2$assignments)
	expect_equal(res1$indep_counts, res2$indep_counts)
})

# ------------------------------------------------------------------------------
# Test 3: Error when beta has incorrect length
# ------------------------------------------------------------------------------
test_that("simulateData errors when beta has wrong length", {
	N <- 30
  T_val <- 6
	R_val <- 5
  d_val <- 4
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	p_int <- compute_p_int(T_val, R_val, d_val)
	# Create an invalid beta vector (one element short)
	beta_wrong <- rnorm(p_int - 1)

	# Manually construct a bad coefs object.
	bad_coefs <- list(
		beta = beta_wrong,
		R = R_val,
		T = T_val,
		d = d_val,
		seed = 123
	)
	class(bad_coefs) <- "FETWFE_coefs"

	expect_error(
		simulateData(
			coefs_obj = bad_coefs,
			N = N,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq
		),
		"length\\(beta\\)"
	)
})


# ------------------------------------------------------------------------------
# Test 4: Pipe test - Verify that the output of simulateData() can be directly piped into fetwfeWithSimulatedData().
# ------------------------------------------------------------------------------
test_that("Output from simulateData can be piped into fetwfeWithSimulatedData", {
	coefs_obj <- genCoefs(
		R = 5,
		T = 6,
		d = 4,
		density = 0.1,
		eff_size = 2,
		seed = 123
	)
	sim_data <- coefs_obj |>
		simulateData(N = 30, sig_eps_sq = 5, sig_eps_c_sq = 5)
	result <- sim_data |> fetwfeWithSimulatedData()
	expect_type(result, "list")
	expect_true("att_hat" %in% names(result))
})

# ------------------------------------------------------------------------------
# Test 5: Pipe test - Verify that the output of simulateData() can be directly piped into
# fetwfeWithSimulatedData() even when d = 0.
# ------------------------------------------------------------------------------
test_that("Output from simulateData can be piped into fetwfeWithSimulatedData", {
	coefs_obj <- genCoefs(
		R = 5,
		T = 6,
		d = 0,
		density = 0.1,
		eff_size = 2,
		seed = 123
	)
	sim_data <- coefs_obj |>
		simulateData(N = 30, sig_eps_sq = 5, sig_eps_c_sq = 5)
	result <- sim_data |> fetwfeWithSimulatedData()
	expect_type(result, "list")
	expect_true("att_hat" %in% names(result))
})


# ------------------------------------------------------------------------------
# Test 6: When distribution is "uniform", covariates are bounded between -sqrt(3) and sqrt(3)
# ------------------------------------------------------------------------------
test_that("Covariates are uniformly bounded when distribution = 'uniform'", {
	N <- 30
  T_val <- 6
	R_val <- 5
  d_val <- 4
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	# For gen_ints = FALSE, expected p = R + (T-1) + d + num_treats.
	p_expected <- compute_p_int(T_val, R_val, d_val)
	coefs_obj <- genCoefs(
		R = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 123
	)

	sim_data <- coefs_obj |>
		simulateData(
			N = N,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
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
# Test 7: Covariates are constant over time for each unit.
# ------------------------------------------------------------------------------
test_that("Covariates are constant over time within each unit", {
	N <- 30
	T_val <- 4
	R_val <- 3
	d_val <- 5
	sig_eps_sq <- 2
	sig_eps_c_sq <- 2
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)

	coefs_obj <- genCoefs(
		R = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 456
	)

	sim_data <- coefs_obj |>
		simulateData(
			N = N,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
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
# Test 8: Error when R >= T
# ------------------------------------------------------------------------------
test_that("simulateData errors when R >= T", {
	# For example, set T = 5 and R = 5 (since R should be <= T-1).
	N <- 30
	T_val <- 5
	R_val <- 5
	d_val <- 2
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	obj <- list(beta = beta, R = R_val, T = T_val, d = d_val, seed = 123)
	class(obj) <- "FETWFE_coefs"

	expect_error(
		simulateData(coefs_obj = obj, N = N, sig_eps_sq = 5, sig_eps_c_sq = 2),
		regexp = "R <= T - 1" # Expect an error message indicating R must be <= T-1
	)
})

# ------------------------------------------------------------------------------
# Test 9: Error when T < 3
# ------------------------------------------------------------------------------
test_that("simulateData errors when T < 3", {
	# For instance, T = 2 should trigger an error (we require at least T = 3).
	N <- 30
	T_val <- 2
	R_val <- 1
	d_val <- 2
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	obj <- list(beta = beta, R = R_val, T = T_val, d = d_val, seed = 123)
	class(obj) <- "FETWFE_coefs"

	expect_error(
		simulateData(coefs_obj = obj, N = N, sig_eps_sq = 5, sig_eps_c_sq = 2),
		regexp = "T >= 3" # Expect an error message indicating T must be at least 3
	)
})

# ------------------------------------------------------------------------------
# Test 10: Error when N < R
# ------------------------------------------------------------------------------
test_that("simulateData errors when N < R", {
	# For example, if there are 3 units and 4 treated cohorts, that's impossible.
	N <- 3
	T_val <- 5
	R_val <- 4
	d_val <- 2
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	p_expected <- compute_p_int(T_val, R_val, d_val)
	beta <- rnorm(p_expected)

	obj <- list(beta = beta, R = R_val, T = T_val, d = d_val, seed = 123)
	class(obj) <- "FETWFE_coefs"

	expect_error(
		simulateData(coefs_obj = obj, N = N, sig_eps_sq = 5, sig_eps_c_sq = 2),
		regexp = "N >= R" # Expect an error message indicating N must be at least R
	)
})
