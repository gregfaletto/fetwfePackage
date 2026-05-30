# Regression tests for #165 / hardened in #178: the GLS whitening step in
# `.estimate_variance_and_gls()` (R/core_funcs.R) uses the block-apply
# identity (I_N kron A) %*% vec_{T,N}(M) = vec_{T,N}(A %*% M) instead of
# materialising the (N*T) x (N*T) Kronecker matrix.
#
# Block 1 (sanity prologue): the algebraic identity `(I_N kron A) %*% vec(M)
# = vec(A %*% M)` holds for a small numeric fixture. Documents the math the
# optimisation relies on; this assertion is necessary but not sufficient.
#
# Block 2 (load-bearing): directly call `fetwfe:::.estimate_variance_and_gls()`
# and assert its `y_gls` / `X_gls` are byte-equal (`tolerance = 0`) to an
# explicit-Kronecker reference computed inline. A future refactor that
# silently breaks the production reshape — e.g. uses the wrong variance
# component in `A`'s scale, drops the `matrix(..., nrow = T)` reshape, or
# any other change that makes the transform not equal to (I_N kron A) %*% y
# — fails this assertion. The pre-#178 test re-derived both sides inline
# and so passed regardless of what the production code did.
#
# Block 3 (sig_eps_c_sq = 0 edge case): same production-call shape as
# Block 2, exercising the degenerate-covariance branch where the closed
# form collapses to `(1/sqrt(sig_eps_sq)) * I_T`.

library(testthat)

test_that("block-apply identity `(I_N kron A) %*% vec(M) = vec(A %*% M)` holds (sanity prologue)", {
	# The mathematical identity that justifies the block-apply form. Holds
	# for any conformable `A` (T x T), `y` (N*T-vector in (unit, time)
	# order), `X` (N*T x p in the same order). This is a sanity check on
	# the math, NOT a check on the production code (the production code's
	# correctness is locked by Block 2 below).
	set.seed(1L)
	N <- 7L
	T <- 4L
	p <- 3L
	A <- matrix(rnorm(T * T), nrow = T, ncol = T)
	y <- rnorm(N * T)
	X <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)

	# Explicit Kronecker form.
	y_kron <- as.vector(kronecker(diag(N), A) %*% y)
	X_kron <- kronecker(diag(N), A) %*% X

	# Block-apply form.
	y_block <- as.vector(A %*% matrix(y, nrow = T))
	X_block <- matrix(A %*% matrix(X, nrow = T), nrow = N * T, ncol = p)

	expect_equal(y_block, y_kron, tolerance = 0)
	expect_equal(
		X_block,
		matrix(X_kron, nrow = N * T, ncol = p),
		tolerance = 0
	)
})

test_that(".estimate_variance_and_gls() matches explicit Kronecker GLS form (load-bearing)", {
	# Tiny fixture. Variance components supplied so the REML path in
	# `estOmegaSqrtInv()` doesn't fire; we're testing the GLS transform,
	# not the variance-component estimator.
	set.seed(42L)
	N <- 12L
	T <- 5L
	p <- 7L
	y <- rnorm(N * T)
	X_mod <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)
	# X_ints is the input to the REML estimator and only consulted when
	# variance components are NA; any conformable matrix works here.
	X_ints <- X_mod
	sig_eps_sq <- 2.3
	sig_eps_c_sq <- 1.1

	res <- fetwfe:::.estimate_variance_and_gls(
		y = y,
		X_ints = X_ints,
		X_mod = X_mod,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		N = N,
		T = T,
		p = p,
		verbose = FALSE
	)

	# Explicit-Kronecker reference. Mirrors `R/core_funcs.R:437-440`'s
	# closed form for `Omega_sqrt_inv` and `R/core_funcs.R:463`'s scaling
	# `A <- sqrt(sig_eps_sq) * Omega_sqrt_inv`.
	J_over_T <- matrix(1 / T, nrow = T, ncol = T)
	Omega_sqrt_inv <- (1 / sqrt(sig_eps_sq)) *
		(diag(T) - J_over_T) +
		(1 / sqrt(sig_eps_sq + T * sig_eps_c_sq)) * J_over_T
	A <- sqrt(sig_eps_sq) * Omega_sqrt_inv
	big_kron <- kronecker(diag(N), A)
	y_gls_ref <- big_kron %*% y # (N*T) x 1 matrix
	X_gls_ref <- big_kron %*% X_mod # (N*T) x p matrix

	# Production return shapes per `R/core_funcs.R:464-475`: `y_gls` is
	# (N*T) x 1 matrix (NOT a length-N*T vector); `X_gls` is (N*T) x p
	# matrix. `ignore_attr = TRUE` because `kronecker(...)` returns
	# matrices without dimnames whereas the production return may carry
	# names.
	expect_equal(res$y_gls, y_gls_ref, tolerance = 0, ignore_attr = TRUE)
	expect_equal(res$X_gls, X_gls_ref, tolerance = 0, ignore_attr = TRUE)
	expect_identical(res$sig_eps_sq, sig_eps_sq)
	expect_identical(res$sig_eps_c_sq, sig_eps_c_sq)
})

test_that(".estimate_variance_and_gls() handles sig_eps_c_sq = 0 (edge case)", {
	# Degenerate-covariance branch: when `sig_eps_c_sq = 0`, the closed
	# form for `Omega_sqrt_inv` reduces to `(1/sqrt(sig_eps_sq)) * I_T`
	# (per the comment at R/core_funcs.R:434-436). Verify the production
	# function still matches the explicit Kronecker form on this branch.
	set.seed(7L)
	N <- 30L
	T <- 4L
	p <- 5L
	y <- rnorm(N * T)
	X_mod <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)
	X_ints <- X_mod
	sig_eps_sq <- 1.5
	sig_eps_c_sq <- 0

	res <- fetwfe:::.estimate_variance_and_gls(
		y = y,
		X_ints = X_ints,
		X_mod = X_mod,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		N = N,
		T = T,
		p = p,
		verbose = FALSE
	)

	J_over_T <- matrix(1 / T, nrow = T, ncol = T)
	Omega_sqrt_inv <- (1 / sqrt(sig_eps_sq)) *
		(diag(T) - J_over_T) +
		(1 / sqrt(sig_eps_sq + T * sig_eps_c_sq)) * J_over_T
	A <- sqrt(sig_eps_sq) * Omega_sqrt_inv
	big_kron <- kronecker(diag(N), A)

	expect_equal(res$y_gls, big_kron %*% y, tolerance = 0, ignore_attr = TRUE)
	expect_equal(
		res$X_gls,
		big_kron %*% X_mod,
		tolerance = 0,
		ignore_attr = TRUE
	)
})
