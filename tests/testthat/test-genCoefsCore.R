library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Helper functions: expected dimensions for genCoefsCore() output.
# ------------------------------------------------------------------------------
compute_num_treats <- function(T, R) {
	T * R - (R * (R + 1)) / 2
}

compute_p <- function(T, R, d) {
	num_treats <- compute_num_treats(T, R)
	R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d
}

# ------------------------------------------------------------------------------
# Test 1: Check that genCoefsCore returns a list with components beta and theta,
#         and that both are numeric vectors.
# ------------------------------------------------------------------------------
test_that("genCoefsCore returns expected output structure", {
	res <- genCoefsCore(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2)

	expect_type(res, "list")
	expect_true("beta" %in% names(res))
	expect_true("theta" %in% names(res))

	expect_true(is.numeric(res$beta))
	expect_true(is.numeric(res$theta))
})

# ------------------------------------------------------------------------------
# Test 2: Check that the length of beta and theta equals the expected number.
# ------------------------------------------------------------------------------
test_that("genCoefsCore returns beta and theta of correct length", {
	R <- 5
	T <- 30
	d <- 12
	expected_length <- compute_p(T, R, d)

	res <- genCoefsCore(R = R, T = T, d = d, density = 0.1, eff_size = 2)
	expect_equal(length(res$beta), expected_length)
	expect_equal(length(res$theta), expected_length)
})

# ------------------------------------------------------------------------------
# Test 3: Reproducibility test: With the same seed, the output should be identical.
# ------------------------------------------------------------------------------
test_that("genCoefsCore is reproducible with the same seed", {
	res1 <- genCoefsCore(
		R = 5,
		T = 30,
		density = 0.1,
		eff_size = 2,
		d = 12,
		seed = 123
	)

	res2 <- genCoefsCore(
		R = 5,
		T = 30,
		density = 0.1,
		eff_size = 2,
		d = 12,
		seed = 123
	)

	expect_equal(res1$beta, res2$beta)
	expect_equal(res1$theta, res2$theta)
})

# ------------------------------------------------------------------------------
# Test 4: Verify that the beta vector returned by genCoefsCore() is a valid input
#         for simulateDataCore().
# ------------------------------------------------------------------------------
test_that("beta from genCoefsCore is a valid input for simulateDataCore", {
	R <- 5
	T <- 30
	d <- 12
	res_coefs <- genCoefsCore(R = R, T = T, d = d, density = 0.1, eff_size = 2)

	N <- 120
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5

	# Use res_coefs$beta as the beta argument for simulateDataCore().
	sim_data <- simulateDataCore(
		N = N,
		T = T,
		R = R,
		d = d,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		beta = res_coefs$beta,
		seed = 123,
		gen_ints = TRUE
	)

	expect_type(sim_data, "list")
	# Check that the returned list contains key elements (e.g., pdata, X, y)
	expected_names <- c("pdata", "X", "y", "coefs", "N", "T", "R", "d", "p")
	for (nm in expected_names) {
		expect_true(nm %in% names(sim_data))
	}

	# Check that the reported p equals the expected dimension.
	expect_equal(sim_data$p, compute_p(T, R, d))
})

# ------------------------------------------------------------------------------
# Test 5: Check that the returned theta vector has approximate sparsity equal to density.
# ------------------------------------------------------------------------------
test_that("genCoefsCore produces theta with approximately correct sparsity", {
	R <- 5
	T <- 30
	d <- 12
	density <- 0.1
	res <- genCoefsCore(R = R, T = T, d = d, density = density, eff_size = 2)
	theta <- res$theta
	total <- length(theta)
	nonzero <- sum(theta != 0)
	prop_nonzero <- nonzero / total

	# Allow some tolerance given randomness; here we allow a deviation of 0.03.
	expect_true(abs(prop_nonzero - density) < 0.03)
})


# ------------------------------------------------------------------------------
# Test 6: d = 0 case.
# ------------------------------------------------------------------------------
test_that("simulateDataCore and fetwfe work when d = 0", {
	# Set parameters for simulation
	N <- 120
	T_val <- 30
	R_val <- 5
	d_val <- 0
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5

	# For d = 0, the expected number of columns in the design matrix is:
	# p = R + (T - 1) + num_treats, where num_treats = T * R - (R*(R+1))/2.
	num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
	expected_p <- R_val + (T_val - 1) + num_treats

	# Generate coefficients with no covariates using genCoefsCore.
	coefs <- genCoefsCore(
		R = R_val,
		T = T_val,
		density = 0.1,
		eff_size = 2,
		d = d_val
	)

	# Use simulateDataCore() with gen_ints = FALSE.
	sim_data <- simulateDataCore(
		N = N,
		T = T_val,
		R = R_val,
		d = d_val,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		beta = coefs$beta,
		seed = 123,
		gen_ints = FALSE
	)

	# Check that cov_names is empty since d = 0.
	expect_equal(length(sim_data$cov_names), 0)

	# Check that the design matrix X has the expected dimensions.
	expect_equal(dim(sim_data$X), c(N * T_val, expected_p))

	# Now run fetwfe() on the simulated panel data.
	result <- fetwfe(
		pdata = sim_data$pdata,
		time_var = sim_data$time_var,
		unit_var = sim_data$unit_var,
		treatment = sim_data$treatment,
		covs = sim_data$cov_names, # should be empty
		response = sim_data$response,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		verbose = FALSE
	)

	# Verify that the returned d (number of covariates) is 0.
	expect_equal(result$d, 0)

	# Also, overall ATT (att_hat) should be numeric and non-missing.
	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
})


test_that("genCoefsCore errors when T < 3", {
	expect_error(
		genCoefsCore(R = 5, T = 2, density = 0.1, eff_size = 2, d = 12),
		regexp = "T must be a numeric value greater than or equal to 3"
	)
})

test_that("genCoefsCore errors when R < 2", {
	expect_error(
		genCoefsCore(R = 1, T = 30, density = 0.1, eff_size = 2, d = 12),
		regexp = "R must be a numeric value greater than or equal to 2"
	)
})

test_that("genCoefsCore errors when R > T - 1", {
	# For example, if T = 30 then R must be <= 29.
	expect_error(
		genCoefsCore(R = 30, T = 30, density = 0.1, eff_size = 2, d = 12),
		regexp = "R must be less than or equal to T - 1"
	)
})

test_that("genCoefsCore errors when d is negative", {
	expect_error(
		genCoefsCore(R = 5, T = 30, density = 0.1, eff_size = 2, d = -1),
		regexp = "d must be a non-negative numeric value"
	)
})

test_that("genCoefsCore errors when density is not between 0 and 1", {
	expect_error(
		genCoefsCore(R = 5, T = 30, density = -0.1, eff_size = 2, d = 12),
		regexp = "density must be"
	)
	expect_error(
		genCoefsCore(R = 5, T = 30, density = 1.1, eff_size = 2, d = 12),
		regexp = "density must be"
	)
})

test_that("genCoefsCore errors when eff_size is not numeric", {
	expect_error(
		genCoefsCore(R = 5, T = 30, density = 0.1, eff_size = "2", d = 12),
		regexp = "eff_size must be a numeric value"
	)
})
