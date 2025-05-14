library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Helper functions: expected dimensions for genCoefs() output.
# ------------------------------------------------------------------------------
compute_num_treats <- function(T, R) {
	T * R - (R * (R + 1)) / 2
}

compute_p <- function(T, R, d) {
	num_treats <- compute_num_treats(T, R)
	R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d
}

# ------------------------------------------------------------------------------
# Test 1: Check that genCoefs returns expeted output
# ------------------------------------------------------------------------------
test_that("genCoefs returns expected output structure", {
	res <- genCoefs(
		R = 5,
		T = 30,
		d = 12,
		density = 0.1,
		eff_size = 2,
		seed = 222
	)

	expect_type(res, "list")
	expect_true("beta" %in% names(res))
	expect_true("R" %in% names(res))
	expect_true("T" %in% names(res))
	expect_true("d" %in% names(res))
	expect_true("seed" %in% names(res))

	expect_true(is.numeric(res$beta))
	expect_true(is.numeric(res$R))
	expect_true(is.numeric(res$T))
	expect_true(is.numeric(res$d))
	expect_true(is.numeric(res$seed))

	expect_equal(res$R, 5)
	expect_equal(res$T, 30)
	expect_equal(res$d, 12)
	expect_equal(res$seed, 222)

	expect_equal(class(res), "FETWFE_coefs")
})

# ------------------------------------------------------------------------------
# Test 2: Check that the length of beta equals the expected number.
# ------------------------------------------------------------------------------
test_that("genCoefs returns beta of correct length", {
	R <- 5
	T <- 30
	d <- 12
	expected_length <- compute_p(T, R, d)

	res <- genCoefs(R = R, T = T, d = d, density = 0.1, eff_size = 2)
	expect_equal(length(res$beta), expected_length)
})

# ------------------------------------------------------------------------------
# Test 3: Reproducibility test: With the same seed, the output should be identical.
# ------------------------------------------------------------------------------
test_that("genCoefs is reproducible with the same seed", {
	res1 <- genCoefs(
		R = 5,
		T = 30,
		density = 0.1,
		eff_size = 2,
		d = 12,
		seed = 123
	)

	res2 <- genCoefs(
		R = 5,
		T = 30,
		density = 0.1,
		eff_size = 2,
		d = 12,
		seed = 123
	)

	expect_equal(res1$beta, res2$beta)
})

# ------------------------------------------------------------------------------
# Test 4: Verify that the beta vector returned by genCoefs() is a valid input
#         for simulateData().
# ------------------------------------------------------------------------------
test_that("beta from genCoefs is a valid input for simulateData", {
	R <- 5
	T <- 30
	d <- 12
	res_coefs <- genCoefs(R = R, T = T, d = d, density = 0.1, eff_size = 2)

	N <- 120
	sig_eps_sq <- 5
	sig_eps_c_sq <- 5

	# Use res_coefs$beta as the beta argument for simulateData().
	sim_data <- res_coefs |>
		simulateData(
			N = N,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq
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
# Test 5: d = 0 case.
# ------------------------------------------------------------------------------
test_that("simulateData and fetwfe work when d = 0", {
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

	# Generate coefficients with no covariates using genCoefs.
	coefs <- genCoefs(
		R = R_val,
		T = T_val,
		density = 0.1,
		eff_size = 2,
		d = d_val,
		seed = 123
	)

	# Use simulateData() with gen_ints = FALSE.
	sim_data <- coefs |>
		simulateData(
			N = N,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
		)

	# Check that cov_names is empty since d = 0.
	expect_equal(length(sim_data$covs), 0)

	# Check that the design matrix X has the expected dimensions.
	expect_equal(dim(sim_data$X), c(N * T_val, expected_p))

	# Now run fetwfe() on the simulated panel data.
	result <- sim_data |>
		fetwfeWithSimulatedData(
			verbose = FALSE
		)

	# Verify that the returned d (number of covariates) is 0.
	expect_equal(result$d, 0)

	# Also, overall ATT (att_hat) should be numeric and non-missing.
	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
})


test_that("genCoefs errors when T < 3", {
	expect_error(
		genCoefs(R = 5, T = 2, density = 0.1, eff_size = 2, d = 12),
		regexp = "T must be a numeric value greater than or equal to 3"
	)
})

test_that("genCoefs errors when R < 2", {
	expect_error(
		genCoefs(R = 1, T = 30, density = 0.1, eff_size = 2, d = 12),
		regexp = "R must be a numeric value greater than or equal to 2"
	)
})

test_that("genCoefs errors when R > T - 1", {
	# For example, if T = 30 then R must be <= 29.
	expect_error(
		genCoefs(R = 30, T = 30, density = 0.1, eff_size = 2, d = 12),
		regexp = "R must be less than or equal to T - 1"
	)
})

test_that("genCoefs errors when d is negative", {
	expect_error(
		genCoefs(R = 5, T = 30, density = 0.1, eff_size = 2, d = -1),
		regexp = "d must be a non-negative numeric value"
	)
})

test_that("genCoefs errors when density is not between 0 and 1", {
	expect_error(
		genCoefs(R = 5, T = 30, density = -0.1, eff_size = 2, d = 12),
		regexp = "density must be"
	)
	expect_error(
		genCoefs(R = 5, T = 30, density = 1.1, eff_size = 2, d = 12),
		regexp = "density must be"
	)
})

test_that("genCoefs errors when eff_size is not numeric", {
	expect_error(
		genCoefs(R = 5, T = 30, density = 0.1, eff_size = "2", d = 12),
		regexp = "eff_size must be a numeric value"
	)
})
