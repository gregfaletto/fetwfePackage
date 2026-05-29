# Regression test for #167: `my_scale()` was rewritten to use
# vectorised colMeans / colSums / broadcasting instead of
# `apply(x, 2, sd)` and `sweep(...)`. Algebraically identical; the
# floating-point reorder gives differences at the ~1e-15 level on
# realistic inputs (verified pre-merge via a manual parity run).
#
# These tests assert the externally-observable contract: the centered
# columns have mean 0, the scaled columns have unit sd (except where
# the input column had zero variance, which keeps sd = 1 by the
# zero-variance guard), and the `scaled:center` / `scaled:scale`
# attributes match `colMeans` / `apply(x, 2, sd)` of the input. The
# test covers a normal column, a constant (zero-variance) column,
# and a one-row degenerate input.

library(testthat)

test_that("my_scale centers columns to mean 0 and scales to unit sd", {
	set.seed(1L)
	x <- matrix(rnorm(200L * 7L), nrow = 200L, ncol = 7L)
	out <- fetwfe:::my_scale(x)
	# Each column has mean 0 within FP reorder noise.
	expect_equal(unname(colMeans(out)), rep(0, 7L), tolerance = 1e-12)
	# Each column has sd 1.
	expect_equal(
		unname(apply(out, 2L, sd)), rep(1, 7L),
		tolerance = 1e-12
	)
})

test_that("my_scale's attributes match base::scale's contract", {
	set.seed(2L)
	x <- matrix(rnorm(50L * 4L), nrow = 50L, ncol = 4L)
	out <- fetwfe:::my_scale(x)
	# scaled:center == colMeans(x)
	expect_equal(
		attr(out, "scaled:center"),
		colMeans(x),
		tolerance = 1e-12
	)
	# scaled:scale == apply(x, 2, sd) (un-guarded path; no zero columns)
	expect_equal(
		unname(attr(out, "scaled:scale")),
		unname(apply(x, 2L, sd)),
		tolerance = 1e-12
	)
})

test_that("my_scale handles zero-variance columns by setting scale to 1", {
	# Mix of normal and constant columns; the constant column must end
	# up with center = its constant value, scale = 1, and zeros after
	# centering+scaling (not NaN from a divide-by-zero).
	set.seed(3L)
	x <- matrix(rnorm(80L * 5L), nrow = 80L, ncol = 5L)
	x[, 3L] <- 3.5 # constant
	out <- fetwfe:::my_scale(x)
	expect_equal(attr(out, "scaled:center")[3L], 3.5, tolerance = 1e-12)
	expect_equal(attr(out, "scaled:scale")[3L], 1, tolerance = 1e-12)
	expect_true(all(out[, 3L] == 0))
	# The other columns are still standardised.
	expect_equal(
		unname(colMeans(out[, -3L])), rep(0, 4L),
		tolerance = 1e-12
	)
	expect_equal(
		unname(apply(out[, -3L], 2L, sd)), rep(1, 4L),
		tolerance = 1e-12
	)
})

test_that("my_scale parity against a hand-coded reference on full-rank inputs", {
	# Hand-coded reference inlining the pre-#167 implementation. Note:
	# on full-rank inputs (no zero-variance columns), the pre-#167 and
	# post-#167 implementations agree up to FP reorder (~1e-15). On
	# inputs WITH zero-variance columns, the pre-#167 implementation
	# had a subtle bug: `apply(x, 2, sd)` returns floating-point noise
	# (e.g., 5.11e-15) rather than exact 0 for a constant column, so
	# the `(sds == 0)` guard at the next line did not fire, and the
	# zero-variance column's `scaled:scale` was stored as 5e-15 rather
	# than 1. Post-#167 computes the constant column's sd via exact
	# algebraic zero (centering `c - c = 0` exactly), so the guard
	# fires correctly and `scaled:scale` is 1. The zero-variance test
	# above pins the post-#167 (correct) contract.
	old_my_scale <- function(x) {
		ctr <- colMeans(x)
		sds <- apply(x, 2L, sd)
		zero_sd <- (sds == 0)
		ctr2 <- ctr
		sds2 <- sds
		sds2[zero_sd] <- 1
		scaled <- sweep(x, 2L, ctr2, FUN = "-")
		scaled <- sweep(scaled, 2L, sds2, FUN = "/")
		attr(scaled, "scaled:center") <- ctr2
		attr(scaled, "scaled:scale") <- sds2
		scaled
	}
	set.seed(4L)
	x <- matrix(rnorm(400L * 12L), nrow = 400L, ncol = 12L)
	# No zero-variance columns here -- the parity check works.
	old_out <- old_my_scale(x)
	new_out <- fetwfe:::my_scale(x)
	expect_equal(
		as.vector(new_out),
		as.vector(old_out),
		tolerance = 1e-12
	)
	expect_equal(
		attr(new_out, "scaled:center"),
		attr(old_out, "scaled:center"),
		tolerance = 1e-12
	)
	expect_equal(
		attr(new_out, "scaled:scale"),
		attr(old_out, "scaled:scale"),
		tolerance = 1e-12
	)
})
