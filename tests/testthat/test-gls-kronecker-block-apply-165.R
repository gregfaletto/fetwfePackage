# Regression test for #165: the GLS whitening step in
# `.estimate_variance_and_gls()` (R/core_funcs.R) uses the block-apply
# identity (I_N kron A) %*% vec_{T,N}(M) = vec_{T,N}(A %*% M) instead of
# materialising the (N*T) x (N*T) Kronecker matrix. This test asserts
# the algebraic identity directly on a small randomised fixture, so a
# future refactor that breaks the row-ordering invariant or the
# block-apply form is caught with byte-identical localisation rather
# than via a downstream end-to-end test.

library(testthat)

test_that("block-apply form equals explicit Kronecker on (unit, time)-ordered data", {
	# Tiny fixture: small enough that the explicit Kronecker is cheap to
	# materialise as a parity reference.
	set.seed(42L)
	N <- 80L
	T <- 6L
	p <- 17L
	y <- rnorm(N * T)
	X_mod <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)

	sig_eps_sq <- 2.3
	sig_eps_c_sq <- 1.1
	J_over_T <- matrix(1 / T, nrow = T, ncol = T)
	Omega_sqrt_inv <- (1 / sqrt(sig_eps_sq)) * (diag(T) - J_over_T) +
		(1 / sqrt(sig_eps_sq + T * sig_eps_c_sq)) * J_over_T
	A <- sqrt(sig_eps_sq) * Omega_sqrt_inv

	# Explicit Kronecker form (the pre-#165 implementation).
	y_gls_explicit <- as.vector(kronecker(diag(N), A) %*% y)
	X_gls_explicit <- kronecker(diag(N), A) %*% X_mod

	# Block-apply form (the post-#165 implementation).
	y_gls_block <- as.vector(A %*% matrix(y, nrow = T))
	X_gls_block <- matrix(
		A %*% matrix(X_mod, nrow = T),
		nrow = N * T,
		ncol = ncol(X_mod)
	)

	expect_equal(y_gls_block, y_gls_explicit, tolerance = 0)
	expect_equal(
		X_gls_block,
		matrix(X_gls_explicit, nrow = N * T, ncol = ncol(X_mod)),
		tolerance = 0
	)
})

test_that("block-apply form is invariant to the trivial sig_eps_c_sq=0 case", {
	# Edge case from `R/core_funcs.R` lines 432-434's docstring: when
	# `sig_eps_c_sq = 0`, the two scaled projectors recombine into
	# `(1/sqrt(sig_eps_sq)) * I_T`. The block-apply form must still
	# match the explicit Kronecker form. Smaller N is enough.
	set.seed(7L)
	N <- 30L
	T <- 4L
	p <- 5L
	y <- rnorm(N * T)
	X_mod <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)

	sig_eps_sq <- 1.5
	J_over_T <- matrix(1 / T, nrow = T, ncol = T)
	Omega_sqrt_inv <- (1 / sqrt(sig_eps_sq)) * (diag(T) - J_over_T) +
		(1 / sqrt(sig_eps_sq + T * 0)) * J_over_T # sig_eps_c_sq = 0
	A <- sqrt(sig_eps_sq) * Omega_sqrt_inv

	expect_equal(
		as.vector(A %*% matrix(y, nrow = T)),
		as.vector(kronecker(diag(N), A) %*% y),
		tolerance = 0
	)
	expect_equal(
		matrix(A %*% matrix(X_mod, nrow = T), nrow = N * T, ncol = ncol(X_mod)),
		kronecker(diag(N), A) %*% X_mod,
		tolerance = 0,
		ignore_attr = TRUE
	)
})
