library(testthat)
library(fetwfe)

# Tests for the two micro-helpers added by issue #83 to R/utility.R:
#   - .cohort_block_inds(g, R, first_inds, num_treats) replaces 6 inline
#     `first_ind_g <- first_inds[g]; last_ind_g <- if (g < R) ... else
#     num_treats` patterns across core_funcs / fetwfe_core / sim_helpers
#     / variance_machinery.
#   - .multinomial_cov(probs) replaces 3 inline `Sigma_pi_hat <-
#     -outer(probs, probs); diag(...) <- probs * (1 - probs)` blocks
#     across fetwfe_core / variance_machinery / event_study.

# ------------------------------------------------------------------------------
# .cohort_block_inds()
# ------------------------------------------------------------------------------

test_that(".cohort_block_inds returns the correct range for an interior cohort", {
	# R = 4, first_inds = c(1, 4, 6, 9), num_treats = 10.
	# g = 2 → block is first_inds[2]:first_inds[3]-1 = 4:5.
	expect_equal(
		fetwfe:::.cohort_block_inds(
			g = 2L,
			R = 4L,
			first_inds = c(1L, 4L, 6L, 9L),
			num_treats = 10L
		),
		4L:5L
	)
})

test_that(".cohort_block_inds returns the correct range for the last cohort", {
	# R = 4, first_inds = c(1, 4, 6, 9), num_treats = 10.
	# g = R = 4 → block is first_inds[4]:num_treats = 9:10.
	expect_equal(
		fetwfe:::.cohort_block_inds(
			g = 4L,
			R = 4L,
			first_inds = c(1L, 4L, 6L, 9L),
			num_treats = 10L
		),
		9L:10L
	)
})

test_that(".cohort_block_inds returns the correct range for the singleton R = 1 case", {
	# R = 1, first_inds = c(1), num_treats = 7.
	# g = 1 = R → block is 1:7.
	expect_equal(
		fetwfe:::.cohort_block_inds(
			g = 1L,
			R = 1L,
			first_inds = 1L,
			num_treats = 7L
		),
		1L:7L
	)
})

test_that(".cohort_block_inds matches the prior inline construction for all cohorts in a panel", {
	# Reference: hand-build the 6-line inline pattern that 6 call sites
	# used pre-#83 and confirm bit-identical output across all R cohorts.
	first_inds <- c(1L, 3L, 6L, 10L)
	R <- 4L
	num_treats <- 12L
	for (g in seq_len(R)) {
		first_ind_g <- first_inds[g]
		last_ind_g <- if (g < R) first_inds[g + 1] - 1L else num_treats
		expected <- first_ind_g:last_ind_g
		expect_identical(
			fetwfe:::.cohort_block_inds(g, R, first_inds, num_treats),
			expected
		)
	}
})

# ------------------------------------------------------------------------------
# .multinomial_cov()
# ------------------------------------------------------------------------------

test_that(".multinomial_cov returns a symmetric matrix with the documented form", {
	probs <- c(0.3, 0.5)
	# Off-diagonal: -outer(probs, probs)[1, 2] = -0.3 * 0.5 = -0.15.
	# Diagonal:    probs * (1 - probs) = c(0.21, 0.25).
	expected <- matrix(c(0.21, -0.15, -0.15, 0.25), 2L, 2L)
	res <- fetwfe:::.multinomial_cov(probs)
	expect_equal(res, expected, tolerance = 1e-14)
	expect_equal(t(res), res, tolerance = 1e-14)
})

test_that(".multinomial_cov diagonal is probs * (1 - probs); off-diagonal is -outer", {
	probs <- c(0.1, 0.2, 0.4)
	res <- fetwfe:::.multinomial_cov(probs)
	expect_equal(diag(res), probs * (1 - probs), tolerance = 1e-14)
	# Off-diagonal entries.
	expect_equal(res[1L, 2L], -0.1 * 0.2, tolerance = 1e-14)
	expect_equal(res[1L, 3L], -0.1 * 0.4, tolerance = 1e-14)
	expect_equal(res[2L, 3L], -0.2 * 0.4, tolerance = 1e-14)
	# Symmetry.
	expect_equal(t(res), res, tolerance = 1e-14)
	# Dimensions.
	expect_equal(dim(res), c(3L, 3L))
})

test_that(".multinomial_cov matches the prior inline -outer/diag construction", {
	# Reference: replicate the pre-#83 inline 2-line block byte-identically.
	probs <- c(0.2, 0.3, 0.1)
	S_ref <- -outer(probs, probs)
	diag(S_ref) <- probs * (1 - probs)
	expect_equal(
		fetwfe:::.multinomial_cov(probs),
		S_ref,
		tolerance = 1e-14
	)
})

test_that(".multinomial_cov works for length-1 probs (degenerate single-cohort case)", {
	# Single-cohort: 1x1 matrix with diag = p * (1 - p), off-diag absent.
	probs <- 0.4
	res <- fetwfe:::.multinomial_cov(probs)
	expect_equal(dim(res), c(1L, 1L))
	expect_equal(res[1L, 1L], 0.24, tolerance = 1e-14)
})
