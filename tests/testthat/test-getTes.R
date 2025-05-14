library(testthat)
library(fetwfe)

# Helper function to compute the expected number of treatment effects.
compute_num_treats <- function(T, R) {
	T * R - (R * (R + 1)) / 2
}

# ------------------------------------------------------------------------------
# Test 1: Check that getTes returns a list with the expected output components.
# ------------------------------------------------------------------------------
test_that("getTes returns a list with att_true and actual_cohort_tes", {
	R_val <- 5
	T_val <- 30
	d_val <- 12

	# Call getTes() with the coefs object
	res <- genCoefs(
		R = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 234
	) |>
		getTes()

	expect_type(res, "list")
	expect_true("att_true" %in% names(res))
	expect_true("actual_cohort_tes" %in% names(res))
})

# ------------------------------------------------------------------------------
# Test 2: Check that att_true equals the mean of actual_cohort_tes.
# ------------------------------------------------------------------------------
test_that("att_true equals mean(actual_cohort_tes)", {
	R_val <- 5
	T_val <- 30
	d_val <- 12

	res <- genCoefs(
		R = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 234
	) |>
		getTes()

	expect_equal(
		res$att_true,
		as.numeric(mean(res$actual_cohort_tes)),
		tolerance = 1e-6
	)
})

# ------------------------------------------------------------------------------
# Test 3: Check that actual_cohort_tes equals the expected cohort effects computed via getActualCohortTes.
# ------------------------------------------------------------------------------
test_that("actual_cohort_tes matches expected output", {
	R_val <- 5
	T_val <- 30
	d_val <- 12

	res_coefs <- genCoefs(
		R = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 234
	)
	beta <- res_coefs$beta

	res <- getTes(res_coefs)

	# Compute number of treatment effects and determine the indices for treatment dummy coefficients.
	num_treats <- compute_num_treats(T_val, R_val)
	base_cols <- if (d_val > 0) {
		R_val + (T_val - 1) + d_val + d_val * R_val + d_val * (T_val - 1)
	} else {
		R_val + (T_val - 1)
	}
	treat_inds <- seq(from = base_cols + 1, length.out = num_treats)

	# Use getFirstInds() to get the starting indices.
	first_inds <- getFirstInds(R = R_val, T = T_val)

	actual_cohort_tes_expected <- getActualCohortTes(
		R = R_val,
		first_inds = first_inds,
		treat_inds = treat_inds,
		coefs = beta,
		num_treats = num_treats
	)

	expect_equal(res$actual_cohort_tes, actual_cohort_tes_expected)
	expect_equal(
		res$att_true,
		as.numeric(mean(actual_cohort_tes_expected)),
		tolerance = 1e-6
	)
})

# ------------------------------------------------------------------------------
# Test 4: Check that getTes errors when beta has an incorrect length.
# ------------------------------------------------------------------------------
test_that("getTes errors with invalid beta length", {
	R_val <- 5
	T_val <- 30
	d_val <- 12

	res_coefs <- genCoefs(
		R = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 234
	)

	# Create a faulty coefs object with incorrect beta length.
	faulty_coefs <- res_coefs
	faulty_coefs$beta <- faulty_coefs$beta[-1]

	expect_error(
		getTes(faulty_coefs),
		regexp = "length(beta) == p is not TRUE",
		fixed = TRUE
	)
})
